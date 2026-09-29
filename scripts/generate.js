'use strict';
// 京大レベルの大問を自動で作り、problems/generated.json に書き出す。
//  1. 京大の過去問の答え X を同じ手順で問題化し、難易度の物差し（採用する範囲）を決める
//  2. 型の選び方を「前の年までの京大だけで次の年を当てる」形で検証する
//  3. 日付を種にして N 問作り、検証を通ったものだけを残す
// 使い方: node scripts/generate.js [種(YYYYMMDD)] [問題数]
const fs = require('fs');
const path = require('path');
const crypto = require('crypto');
const chem = require('../src/chem');
const G = require('../src/generator');
const { predict } = require('./forecast');
const { loadRDKit, checkBig } = require('./validate');

const ROOT = path.join(__dirname, '..');
const DATA = JSON.parse(fs.readFileSync(path.join(ROOT, 'data', 'kakomon.json'), 'utf8'));

function todayJST() {
  const d = new Date(Date.now() + 9 * 3600 * 1000);
  return d.toISOString().slice(0, 10).replace(/-/g, '');
}

// 断片ライブラリは異性体の全列挙で重いので、ソースが変わらない限り使い回す
function loadLibrary(RDKit) {
  const h = crypto.createHash('sha1');
  for (const f of ['src/generator.js', 'src/enumerate.js', 'src/chem.js']) h.update(fs.readFileSync(path.join(ROOT, f)));
  const key = h.digest('hex');
  const file = path.join(ROOT, 'generated', 'library.json');
  try {
    const c = JSON.parse(fs.readFileSync(file, 'utf8'));
    if (c.key === key) return c.lib;
  } catch (e) { /* 作り直す */ }
  const t = Date.now();
  const lib = G.buildLibrary(RDKit);
  fs.mkdirSync(path.dirname(file), { recursive: true });
  fs.writeFileSync(file, JSON.stringify({ key, lib }));
  console.log(`断片ライブラリを作った（${Object.keys(lib).length} 分類・${((Date.now() - t) / 1000).toFixed(0)} 秒）`);
  return lib;
}

const yearNum = (y) => parseInt(y, 10);
const quantile = (xs, q) => {
  const s = [...xs].sort((a, b) => a - b);
  const i = (s.length - 1) * q;
  return s[Math.floor(i)] + (s[Math.ceil(i)] - s[Math.floor(i)]) * (i - Math.floor(i));
};

function kyotoSpec(RDKit, lib, X) {
  const prods = chem.hydrolyze(RDKit, chem.graphFromSmiles(RDKit, X));
  const classes = prods.map((p) => G.classify(RDKit, lib, p));
  if (classes.some((c) => !c)) return null;
  return { template: G.templateFor(classes), frags: prods.map((smiles, i) => ({ cls: classes[i], smiles })), answerX: X };
}

async function main() {
  const seedStr = process.argv[2] || todayJST();
  const N = Number(process.argv[3]) || 12;
  const RDKit = await loadRDKit();
  const lib = loadLibrary(RDKit);
  const forecast = JSON.parse(fs.readFileSync(path.join(ROOT, 'generated', 'forecast.json'), 'utf8'));
  const weights = G.cardWeights(forecast);

  // ---- 1. 京大の X を問題化（到達範囲と難易度の物差し）----
  const kyoto = [];
  for (const [y, X] of Object.entries(DATA.structures)) {
    const spec = kyotoSpec(RDKit, lib, X);
    if (!spec || !spec.template) { console.log(`  ${y}: 型の範囲外`); continue; }
    // 計算段階の乱数で難易度が揺れるので、種を変えて平均する
    const ds = [];
    let first = null;
    for (let k = 1; k <= 5; k++) {
      const p = G.buildProblem(RDKit, G.rng(k * 7919), lib, weights, { ...spec, id: `kyoto-${y}` });
      if (p && !checkBig(RDKit, p).errors.length) { ds.push(p.meta.difficulty); if (!first) first = p; }
    }
    if (!first) { console.log(`  ${y}: 問題化できない`); continue; }
    const d = ds.reduce((a, b) => a + b, 0) / ds.length;
    first.title = `京大${yearNum(y)}年度型の再現${/B$/.test(y) ? '（その2）' : ''}`;
    first.meta.kyoto = y;
    kyoto.push({ y, template: spec.template, d, p: first });
  }
  const ks = kyoto.map((k) => k.d);
  const lo = Math.round(quantile(ks, 0.25) * 10) / 10;
  const hi = Math.round(Math.max(...ks) * 1.25 * 10) / 10;
  const kstats = { min: Math.min(...ks), p25: quantile(ks, 0.25), median: quantile(ks, 0.5), p75: quantile(ks, 0.75), max: Math.max(...ks) };
  console.log(`難易度の6段階の境目: 基礎 < ${(kstats.p25 * 0.7).toFixed(1)} ≦ 標準 < ${kstats.p25.toFixed(1)} ≦ 京大下位 < ${kstats.median.toFixed(1)} ≦ 京大平均 < ${kstats.p75.toFixed(1)} ≦ 京大上位 < ${kstats.max.toFixed(1)} ≦ 京大超え`);
  console.log(`京大の答え ${Object.keys(DATA.structures).length} 個のうち ${kyoto.length} 個を生成器の型で再現できた`);
  kyoto.forEach((k) => console.log(`  ${k.y.padEnd(5)} 型=${k.template.padEnd(15)} 難易度 ${k.d.toFixed(1)}`));
  console.log(`採用する難易度: ${lo} 〜 ${hi}（京大の下位25% 〜 最大の1.25倍）`);

  // ---- 2. 型の選び方の当てる練習（前の年までの京大の型だけで次の年の型を当てる）----
  // 候補: 一様 / 過去の型の割合 / 出題予測の要素確率 / 両者の半々。実際の型に付けた確率の平均が最も高いものを使う
  const tplKeys = Object.keys(G.TEMPLATES);
  const priorFrom = (list) => {
    const pr = {};
    tplKeys.forEach((t) => { pr[t] = (list.filter((z) => z.template === t).length + 0.5) / (list.length + 0.5 * tplKeys.length); });
    return pr;
  };
  const METHODS = {
    uniform: () => tplKeys.map((t) => [t, 1]),
    prior: (fc, pr) => tplKeys.map((t) => [t, pr[t]]),
    forecast: (fc) => G.templateWeights(fc, null),
    mix: (fc, pr) => G.templateWeights(fc, pr),
  };
  const sorted = [...kyoto].sort((a, b) => yearNum(a.y) - yearNum(b.y) || a.y.localeCompare(b.y));
  const btScore = {};
  let nBt = 0;
  for (const k of sorted) {
    const past = sorted.filter((z) => yearNum(z.y) < yearNum(k.y));
    if (past.length < 2) continue;
    nBt++;
    const fc = { probs: predict(yearNum(k.y), { lambda: 1, gap: 0, alpha: 0.5 }) };
    for (const [m, f] of Object.entries(METHODS)) {
      const w = f(fc, priorFrom(past));
      const tot = w.reduce((a, [, v]) => a + v, 0);
      btScore[m] = (btScore[m] || 0) + w.find(([t]) => t === k.template)[1] / tot;
    }
  }
  let method = 'uniform';
  for (const m of Object.keys(METHODS)) if (btScore[m] > btScore[method] + 1e-9) method = m;
  console.log(`型の当てる練習（${nBt} 題、実際の型に付けた確率の平均）: ${Object.keys(METHODS).map((m) => `${m} ${(100 * btScore[m] / nBt).toFixed(0)}%`).join(' / ')} → 採用 ${method}`);
  const finalW = Object.fromEntries(METHODS[method](forecast, priorFrom(kyoto)));

  // ---- 3. 生成 ----
  const seed = parseInt(seedStr, 10) || [...seedStr].reduce((a, c) => a * 31 + c.charCodeAt(0), 7);
  const r = G.rng(seed);
  const out = [];
  const perTpl = {};
  const cap = Math.ceil(N / 4);
  const seenX = new Set(kyoto.map((k) => k.p.answer));
  const stats = { tried: 0, invalid: 0, band: 0, dup: 0 };
  const perGrade = {};
  const seenSig = new Set();
  for (let i = 0; i < N * 200 && out.length < N; i++) {
    let spec = G.sampleSpec(r, lib, finalW);
    // 7 割の問題で、断片どうしに関係が生まれるように寄せる（京大型の「A を酸化すると C」など）
    if (spec && r() < 0.7) spec = G.biasRelations(RDKit, r, lib, spec);
    if (!spec || (perTpl[spec.template] || 0) >= cap) continue;
    stats.tried++;
    spec.id = `g${seedStr}-${String(out.length + 1).padStart(2, '0')}`;
    let p;
    try { p = G.buildProblem(RDKit, r, lib, weights, spec); } catch (e) { p = null; }
    if (!p) { stats.invalid++; continue; }
    if (seenX.has(p.answer)) { stats.dup++; continue; }
    // 京大型: 化合物どうしの関係を含む問題を 6 割以上にする
    if (!(p.relations && p.relations.length) && out.length >= 2 && out.filter((x) => x.relations && x.relations.length).length < Math.ceil(out.length * 0.6)) { stats.norel = (stats.norel || 0) + 1; continue; }
    // 多様性: 同じ型で同じ分類の組み合わせ（例: フタル酸＋アルコール2つ）は1日1問まで
    const sig = spec.template + ':' + p.fragments.map((f) => (f.given ? 'g' : f.kind)).sort().join('+');
    if (seenSig.has(sig)) { stats.dup++; continue; }
    // 6段階それぞれから同じくらいの数を出す（段階ごとの上限 = N/6 を切り上げ）
    const g = G.gradeOf(p.meta.difficulty, kstats);
    if (p.meta.difficulty > hi * 1.3 || (perGrade[g] || 0) >= Math.ceil(N / 6)) { stats.band++; continue; }
    perGrade[g] = (perGrade[g] || 0) + 1;
    p.meta.grade = g;
    const v = checkBig(RDKit, p);
    if (v.errors.length) { stats.invalid++; console.log('  検証で落とした: ' + v.errors[0]); continue; }
    seenX.add(p.answer);
    seenSig.add(sig);
    perTpl[spec.template] = (perTpl[spec.template] || 0) + 1;
    out.push(p);
  }
  console.log(`生成: ${out.length} 問（試行 ${stats.tried}、作れない ${stats.invalid}、難易度が範囲外 ${stats.band}、重複 ${stats.dup}、関係なしで見送り ${stats.norel || 0}）。化合物どうしの関係を含む問題 ${out.filter((x) => x.relations && x.relations.length).length} 問、誘導体を含む問題 ${out.filter((x) => x.derived && x.derived.length).length} 問`);
  out.forEach((p) => console.log(`  ${p.id} ${p.meta.template.padEnd(15)} ${p.formula.padEnd(10)} 難易度 ${p.meta.difficulty}  ${p.answer}`));

  // 予測の上位の出題要素のうち、今日の問題で練習できるもの
  const tagsIn = new Set();
  for (const p of out) {
    tagsIn.add(G.TEMPLATES[p.meta.template].tag);
    for (const f of p.fragments) (f.clues || []).forEach((c) => tagsIn.add(G.CARD_TAG[c.card]));
    if (p.assemble) p.assemble.clues.forEach((c) => tagsIn.add(G.CARD_TAG[c.card]));
    (p.calcs || []).forEach((c) => { tagsIn.add('mole_calc'); tagsIn.add(c.key === 'h2' ? 'hydrogenation' : 'combustion_analysis'); });
  }
  const top = Object.keys(forecast.probs).slice(0, 12);
  console.log(`${forecast.target}年度の予測上位12要素のうち、今日の問題で練習できるもの ${top.filter((t) => tagsIn.has(t)).length}/12: ${top.filter((t) => tagsIn.has(t)).map((t) => DATA.tags[t]).join('、')}`);

  kyoto.forEach((k) => { k.p.meta.grade = G.gradeOf(k.p.meta.difficulty, kstats); });
  console.log('段階ごとの問題数: ' + G.GRADES.map((n, i) => `${n} ${out.filter((p) => p.meta.grade === i + 1).length}`).join('・'));
  const problems = [...kyoto.map((k) => k.p), ...out].map((p) => ({ ...p, meta: { ...p.meta, seed: seedStr, band: [lo, hi], kyotoStats: kstats } }));
  fs.writeFileSync(path.join(ROOT, 'problems', 'generated.json'), JSON.stringify(problems, null, 1) + '\n');
  console.log(`problems/generated.json に ${problems.length} 問（京大の再現 ${kyoto.length}・自動生成 ${out.length}）を書いた`);
}

main();
