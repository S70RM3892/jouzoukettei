'use strict';
// 京大有機の出題要素の予測と、その「当てる練習」（時系列の後ろ向き検証）。
// ある年 y を当てるときは、y より前の年のデータだけを使う。予測の仕方（パラメータ）も
// y より前の年での成績だけで選び直す（入れ子の検証）ので、的中率は答えを見ていない値になる。
// 使い方: node scripts/forecast.js [予測したい年]
const fs = require('fs');
const path = require('path');

const DATA = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'data', 'kakomon.json'), 'utf8'));
const TAGS = Object.keys(DATA.tags);
const YEARS = Object.keys(DATA.years).map(Number).sort((a, b) => a - b);
const tagsOf = (y) => new Set(DATA.years[y].tags);

// 予測: 近い年ほど重い出現率（減衰 lambda）＋しばらく出ていない要素を上げる補正（gap）
function predict(target, params) {
  const past = YEARS.filter((y) => y < target);
  const p = {};
  for (const t of TAGS) {
    let num = params.alpha, den = 2 * params.alpha, last = null;
    for (const y of past) {
      const w = params.lambda ** (target - y - 1);
      den += w;
      if (tagsOf(y).has(t)) { num += w; last = y; }
    }
    let q = num / den;
    const gap = last === null ? 5 : Math.min(5, target - last - 1);
    q *= 1 + params.gap * (gap / 5 - 0.3);
    p[t] = Math.max(0.01, Math.min(0.99, q));
  }
  return p;
}

// 評価: その年に実際に出た要素の数 k だけ上位を取り、何個当たったか（R 精度）と、確率の二乗誤差（ブライア）
function score(p, actual) {
  const k = actual.size;
  const top = Object.entries(p).sort((a, b) => b[1] - a[1]).slice(0, k).map(([t]) => t);
  const hit = top.filter((t) => actual.has(t)).length;
  const brier = TAGS.reduce((s, t) => s + (p[t] - (actual.has(t) ? 1 : 0)) ** 2, 0) / TAGS.length;
  return { rprec: hit / k, brier, hit, k, top };
}

const GRID = [];
for (const lambda of [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]) {
  for (const gap of [-0.4, -0.2, 0, 0.2, 0.4, 0.6]) {
    for (const alpha of [0.25, 0.5, 1]) GRID.push({ lambda, gap, alpha });
  }
}

const MIN_HISTORY = 4;
const evalYears = YEARS.filter((y) => YEARS.filter((z) => z < y).length >= MIN_HISTORY);

function meanScore(params, years) {
  const s = years.map((y) => score(predict(y, params), tagsOf(y)));
  return { rprec: s.reduce((a, b) => a + b.rprec, 0) / s.length, brier: s.reduce((a, b) => a + b.brier, 0) / s.length };
}

function tune(years) {
  let best = null;
  for (const g of GRID) {
    const m = meanScore(g, years);
    const obj = m.rprec - m.brier; // 当たりが多く、確率の外れが小さいもの
    if (!best || obj > best.obj) best = { params: g, obj, ...m };
  }
  return best;
}

function baselines(y) {
  const actual = tagsOf(y);
  const prev = YEARS.filter((z) => z < y).pop();
  const copy = {};
  TAGS.forEach((t) => { copy[t] = tagsOf(prev).has(t) ? 0.9 : 0.1; });
  return {
    freq: score(predict(y, { lambda: 1, gap: 0, alpha: 0.5 }), actual),
    lastYear: score(copy, actual),
  };
}

function main() {
  const rows = [];
  for (const y of evalYears) {
    const tuneYears = evalYears.filter((z) => z < y);
    const params = tuneYears.length >= 2 ? tune(tuneYears).params : { lambda: 1, gap: 0, alpha: 0.5 };
    const s = score(predict(y, params), tagsOf(y));
    const b = baselines(y);
    rows.push({ y, params, s, b });
  }
  const avg = (f) => rows.reduce((a, r) => a + f(r), 0) / rows.length;
  console.log('年   | 学習後モデル | 出現頻度のみ | 前年と同じ | 使ったパラメータ');
  for (const r of rows) {
    console.log(`${r.y} | ${r.s.hit}/${r.s.k} (${(100 * r.s.rprec).toFixed(0)}%) | ${r.b.freq.hit}/${r.b.freq.k} | ${r.b.lastYear.hit}/${r.b.lastYear.k} | λ=${r.params.lambda} gap=${r.params.gap} α=${r.params.alpha}`);
  }
  console.log(`平均 R精度: 学習後 ${(100 * avg((r) => r.s.rprec)).toFixed(1)}% / 出現頻度 ${(100 * avg((r) => r.b.freq.rprec)).toFixed(1)}% / 前年コピー ${(100 * avg((r) => r.b.lastYear.rprec)).toFixed(1)}%`);
  console.log(`平均ブライア: 学習後 ${avg((r) => r.s.brier).toFixed(3)} / 出現頻度 ${avg((r) => r.b.freq.brier).toFixed(3)}`);
  const random = avg((r) => r.s.k / TAGS.length);
  console.log(`当てずっぽう（出る数だけ無作為に選ぶ）の期待値: ${(100 * random).toFixed(1)}%`);

  // 入れ子の検証で勝った方を採用する（調整が効かないなら単純な出現頻度）
  const tunedWins = avg((r) => r.s.rprec) > avg((r) => r.b.freq.rprec);
  const target = Number(process.argv[2]) || YEARS[YEARS.length - 1] + 1;
  const final = tunedWins ? tune(evalYears) : { params: { lambda: 1, gap: 0, alpha: 0.5 } };
  console.log(`採用: ${tunedWins ? '調整したモデル' : '出現頻度モデル（調整は当てすぎで逆効果だった）'}`);
  const p = predict(target, final.params);
  const ranked = Object.entries(p).sort((a, b) => b[1] - a[1]);
  console.log(`\n${target}年度の予測（λ=${final.params.lambda} gap=${final.params.gap} α=${final.params.alpha}）`);
  ranked.slice(0, 15).forEach(([t, q]) => console.log(`  ${(100 * q).toFixed(0).padStart(3)}%  ${DATA.tags[t]}`));

  const out = {
    target,
    params: final.params,
    backtest: {
      years: rows.map((r) => ({ year: r.y, hit: r.s.hit, k: r.s.k, freqHit: r.b.freq.hit, lastYearHit: r.b.lastYear.hit })),
      rprec: avg((r) => r.s.rprec), rprecFreq: avg((r) => r.b.freq.rprec), rprecLastYear: avg((r) => r.b.lastYear.rprec),
      brier: avg((r) => r.s.brier), brierFreq: avg((r) => r.b.freq.brier), rprecRandom: random, chosen: tunedWins ? 'tuned' : 'frequency',
    },
    probs: Object.fromEntries(ranked),
  };
  fs.mkdirSync(path.join(__dirname, '..', 'generated'), { recursive: true });
  fs.writeFileSync(path.join(__dirname, '..', 'generated', 'forecast.json'), JSON.stringify(out, null, 2) + '\n');
}

if (require.main === module) main();
module.exports = { predict, TAGS };
