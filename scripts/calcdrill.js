'use strict';
// 入試の定番の計算を、数値を変えて問題にし、problems/calc.json に書き出す（知識確認と同じ4択）。
// 誤答は、よくある間違い（×3 を忘れる・縮合で外れる水を引き忘れる・C=C を数え違えるなど）から作る。
// 原子量は入試の値（H 1.0, C 12, N 14, O 16, K 39, I 127）
// 使い方: node scripts/calcdrill.js [種]
const fs = require('fs');
const path = require('path');
const G = require('../src/generator');

const r = G.rng(parseInt(process.argv[2] || '20260930', 10));
const pick = (a) => a[Math.floor(r() * a.length)];
const sig3 = (x) => Number(x.toPrecision(3));
const f3 = (x) => { const s = x.toPrecision(3); return s.includes('e') ? String(Number(s)) : s; };
const items = [];
function add(topic, q, a, wrong, explain) {
  const w = [...new Set(wrong.filter((x) => x !== a))].slice(0, 3);
  if (w.length < 3) return; // 誤答がそろわない数値の組は使わない
  items.push({ topic, q, a, wrong: w, explain, calc: true });
}

// ---- 油脂: けん化価・ヨウ素価・水素付加 ----
const ACIDS = [
  { name: 'パルミチン酸', f: 'C₁₅H₃₁COOH', M: 256, cc: 0 },
  { name: 'ステアリン酸', f: 'C₁₇H₃₅COOH', M: 284, cc: 0 },
  { name: 'オレイン酸', f: 'C₁₇H₃₃COOH', M: 282, cc: 1 },
  { name: 'リノール酸', f: 'C₁₇H₃₁COOH', M: 280, cc: 2 },
  { name: 'リノレン酸', f: 'C₁₇H₂₉COOH', M: 278, cc: 3 },
];
const fatOf = (as) => ({ M: 92 + as.reduce((s, a) => s + a.M, 0) - 3 * 18, cc: as.reduce((s, a) => s + a.cc, 0) });
const fatName = (as) => (as.every((a) => a === as[0]) ? `${as[0].name}だけからなる油脂` : `${as.map((a) => a.name).join('・')}が1分子ずつ結合した油脂`);
for (let i = 0; i < 6; i++) {
  const as = i < 3 ? Array(3).fill(ACIDS[[1, 2, 3][i]]) : [pick(ACIDS), pick(ACIDS), pick(ACIDS)];
  const fat = fatOf(as);
  const sv = (3 * 56 * 1000) / fat.M;
  add('fat', `${fatName(as)}（分子量 ${fat.M}）のけん化価（油脂 1 g のけん化に必要な KOH の mg 数）は？`, f3(sv),
    [f3((56 * 1000) / fat.M), f3((3 * 40 * 1000) / fat.M), f3((3 * 56 * 1000) / (fat.M + 54))],
    `油脂 1 mol のけん化に KOH は 3 mol。1 g は 1/${fat.M} mol なので、KOH は 3 × 56 × 1000 / ${fat.M} = ${f3(sv)} mg。`);
  if (fat.cc) {
    const iv = (fat.cc * 254 * 100) / fat.M;
    add('fat', `${fatName(as)}（分子量 ${fat.M}、1分子に C=C が ${fat.cc} 個）のヨウ素価（油脂 100 g に付加する I₂ の g 数）は？`, f3(iv),
      [f3((fat.cc * 127 * 100) / fat.M), f3((254 * 100) / fat.M), f3(((fat.cc + 1) * 254 * 100) / fat.M)],
      `C=C 1 個に I₂ 1 分子が付加する。100 g は 100/${fat.M} mol なので、I₂ は ${fat.cc} × 254 × 100 / ${fat.M} = ${f3(iv)} g。`);
    const g = pick([10, 50, 100, 200]);
    const L = (g / fat.M) * fat.cc * 22.4;
    add('fat', `${fatName(as)}（分子量 ${fat.M}）${g} g を完全に硬化油にするのに必要な H₂ は標準状態で何 L か`, f3(L),
      [f3((g / fat.M) * 22.4), f3((g / fat.M) * fat.cc * 3 * 22.4), f3((g / fat.M) * fat.cc * 11.2)],
      `C=C は1分子に ${fat.cc} 個。${g} g は ${f3(g / fat.M)} mol なので、H₂ は ${f3(g / fat.M)} × ${fat.cc} × 22.4 = ${f3(L)} L。`);
  }
}
// けん化価から分子量と C=C の数（逆算）
for (let i = 0; i < 2; i++) {
  const as = Array(3).fill(pick(ACIDS.slice(2)));
  const fat = fatOf(as);
  const sv = sig3((3 * 56 * 1000) / fat.M);
  add('fat', `ある油脂（1種類の脂肪酸だけからなる）のけん化価は ${sv} だった。この油脂の分子量はおよそいくらか`, String(Math.round((3 * 56 * 1000) / sv)),
    [String(Math.round((56 * 1000) / sv)), String(Math.round((3 * 40 * 1000) / sv)), String(Math.round((3 * 56 * 1000) / sv) - 54)],
    `分子量 M の油脂 1 g に KOH は 3 × 56 × 1000 / M mg。${sv} = 168000 / M より M ≈ ${Math.round(168000 / sv)}。`);
}

// ---- アミノ酸: 等電点・電気泳動・イオンの割合 ----
const AA = [
  { name: 'グリシン', pk: [2.34, 9.60], kind: '中性' },
  { name: 'アラニン', pk: [2.34, 9.69], kind: '中性' },
  { name: 'フェニルアラニン', pk: [1.83, 9.13], kind: '中性' },
  { name: 'グルタミン酸', pk: [2.19, 4.25, 9.67], kind: '酸性' },
  { name: 'アスパラギン酸', pk: [1.88, 3.65, 9.60], kind: '酸性' },
  { name: 'リシン', pk: [2.18, 8.95, 10.53], kind: '塩基性' },
];
const pI = (a) => (a.kind === '塩基性' ? (a.pk[1] + a.pk[2]) / 2 : (a.pk[0] + a.pk[1]) / 2);
for (const a of AA) {
  const p = pI(a);
  const avg = a.pk.reduce((s, x) => s + x, 0) / a.pk.length;
  const wrongPair = a.kind === '塩基性' ? (a.pk[0] + a.pk[1]) / 2 : a.pk.length === 3 ? (a.pk[1] + a.pk[2]) / 2 : (a.pk[0] + 7) / 2;
  add('amino_acid', `${a.name}（${a.kind}アミノ酸）の電離に関する pKa は ${a.pk.join('、')} である。等電点はいくらか`, p.toFixed(2),
    [avg.toFixed(2), wrongPair.toFixed(2), (a.pk[a.pk.length - 1] - a.pk[0]).toFixed(2)],
    a.kind === '中性' ? `–COOH と –NH₃⁺ の pKa の平均。(${a.pk[0]} + ${a.pk[1]}) / 2 = ${p.toFixed(2)}。`
      : a.kind === '酸性' ? `電荷 0 の双性イオンをはさむ2つの pKa（α-COOH と側鎖の COOH）の平均。(${a.pk[0]} + ${a.pk[1]}) / 2 = ${p.toFixed(2)}。`
      : `電荷 0 の双性イオンをはさむ2つの pKa（α-NH₃⁺ と側鎖の NH₃⁺）の平均。(${a.pk[1]} + ${a.pk[2]}) / 2 = ${p.toFixed(2)}。`);
}
for (const pH of [3, 6, 11]) {
  const a = pick(AA);
  const p = pI(a);
  const ans = pH < p - 0.3 ? '陰極側' : pH > p + 0.3 ? '陽極側' : 'ほとんど移動しない';
  add('amino_acid', `${a.name}（等電点 ${p.toFixed(2)}）を pH ${pH} の緩衝液中で電気泳動すると？`, ans,
    ['陰極側', '陽極側', 'ほとんど移動しない', '両方の極に分かれる'],
    `等電点より pH が低いと陽イオンが多くなり陰極へ、高いと陰イオンが多くなり陽極へ動く。pH ${pH} と等電点 ${p.toFixed(2)} を比べる。`);
}
for (const pH of [1.34, 2.34, 3.34]) {
  const ratio = 10 ** (2.34 - pH);
  add('amino_acid', `グリシンの –COOH の電離定数 K₁ = 10^−2.34 mol/L とする。pH ${pH} の水溶液中で、陽イオン H₃N⁺CH₂COOH と双性イオン H₃N⁺CH₂COO⁻ の濃度の比 [陽イオン]/[双性イオン] は？`,
    f3(ratio), [f3(1 / ratio), f3(ratio * 10), f3(2.34 / pH)],
    `K₁ = [双性イオン][H⁺]/[陽イオン] より、[陽イオン]/[双性イオン] = [H⁺]/K₁ = 10^(2.34 − ${pH}) = ${f3(ratio)}。`);
}

// ---- ペプチドの分子量 ----
const RES = { グリシン: 75, アラニン: 89, セリン: 105, システイン: 121, フェニルアラニン: 165, チロシン: 181, リシン: 146, グルタミン酸: 147 };
for (let i = 0; i < 5; i++) {
  const n = 2 + (i % 3);
  const seq = Array.from({ length: n }, () => pick(Object.keys(RES)));
  const sum = seq.reduce((s, x) => s + RES[x], 0);
  const M = sum - 18 * (n - 1);
  add('protein', `${seq.join('・')}（分子量 ${seq.map((x) => RES[x]).join('・')}）が この順にペプチド結合でつながった${['ジ', 'トリ', 'テトラ'][n - 2]}ペプチドの分子量は？`, String(M),
    [String(sum), String(sum - 18 * n), String(sum + 18 * (n - 1))],
    `ペプチド結合が1つできるごとに H₂O（18）が外れる。${n} 個のアミノ酸では ${n - 1} 個。${sum} − 18 × ${n - 1} = ${M}。`);
}
// タンパク質中の窒素の割合から含有量（ケルダール法の考え方）
for (const pct of [16, 15]) {
  const g = pick([0.20, 0.35, 0.48]);
  add('protein', `ある食品のタンパク質は質量の ${pct}% が窒素である。食品 10 g から窒素が ${g.toFixed(2)} g 検出された（窒素はすべてタンパク質由来）。食品中のタンパク質の質量の割合は何 %？`,
    f3((g / (pct / 100)) / 10 * 100), [f3(g * pct / 10), f3(g / 10 * 100), f3((g / (pct / 100)))],
    `タンパク質の質量 = ${g.toFixed(2)} / ${pct / 100} = ${f3(g / (pct / 100))} g。10 g あたりなので ${f3((g / (pct / 100)) / 10 * 100)}%。`);
}

// ---- ビニロン: アセタール化 ----
for (let i = 0; i < 4; i++) {
  const mass = pick([44, 88, 100, 132]);
  const x = pick([30, 35, 40, 50]);
  const nOH = mass / 44;
  const gain = (nOH * x / 100 / 2) * 12;
  add('vinylon', `ポリビニルアルコール ${mass} g の –OH の ${x}% をホルムアルデヒドでアセタール化した。得られたビニロンの質量は？`, f3(mass + gain),
    [f3(mass + nOH * x / 100 * 12), f3(mass + (nOH * x / 100 / 2) * 30), f3(mass - (nOH * x / 100 / 2) * 18)],
    `–OH 2 個と HCHO 1 分子から H₂O 1 分子が外れて –O–CH₂–O– になる。1か所あたり +30 − 18 = +12。–OH は ${f3(nOH)} mol、アセタール化は ${f3(nOH * x / 100 / 2)} mol なので +${f3(gain)} g。`);
}
{
  const mass = 88, gain = 3.6;
  const x = (gain / 12) * 2 / (mass / 44) * 100;
  add('vinylon', `ポリビニルアルコール ${mass} g をホルムアルデヒドで処理すると、質量が ${gain} g 増えた。–OH の何 % がアセタール化されたか`, f3(x),
    [f3(x / 2), f3((gain / 30) * 2 / (mass / 44) * 100), f3((gain / 12) / (mass / 44) * 100 * 4)],
    `アセタール化1か所で +12、–OH 2 個を使う。${gain} / 12 = ${f3(gain / 12)} mol のアセタール → –OH ${f3(gain / 6)} mol。全 –OH ${mass / 44} mol に対して ${f3(x)}%。`);
}

// ---- 元素分析 ----
const COMP = [['CH₂O', 'C₂H₄O₂', 60, 1, 2, 1], ['C₂H₆O', 'C₂H₆O', 46, 2, 6, 1], ['C₃H₆O', 'C₃H₆O', 58, 3, 6, 1], ['C₂H₄O', 'C₄H₈O₂', 88, 2, 4, 1], ['C₇H₆O₂', 'C₇H₆O₂', 122, 7, 6, 2]];
for (const [emp, mol, M, c, h, o] of COMP) {
  const mm = 12 * c + 1 * h + 16 * o;
  const k = M / mm;
  const w = pick([2.0, 3.0, 4.4]) * 10;
  const nC = (w / M) * c * k, nH = (w / M) * h * k;
  const co2 = sig3(nC * 44), h2o = sig3(nH / 2 * 18);
  add('elemental_analysis', `C・H・O からなる化合物 ${w.toFixed(1)} mg を完全燃焼させると、CO₂ ${co2} mg と H₂O ${h2o} mg が生じた。分子量が ${M} のとき分子式は？`, mol,
    [emp === mol ? `C${c * 2}H${h * 2}O${o * 2}`.replace(/(\d+)/g, (d) => [...d].map((x) => '₀₁₂₃₄₅₆₇₈₉'[x]).join('')) : emp, `C${c + 1}H${h}O${o}`.replace(/(\d+)/g, (d) => [...d].map((x) => '₀₁₂₃₄₅₆₇₈₉'[x]).join('')), `C${c}H${h + 2}O${o}`.replace(/(\d+)/g, (d) => [...d].map((x) => '₀₁₂₃₄₅₆₇₈₉'[x]).join(''))],
    `C = ${co2} × 12/44、H = ${h2o} × 2/18、O は残り。物質量の比から組成式 ${emp}（式量 ${mm}）。分子量 ${M} ÷ ${mm} = ${k} なので ${mol}。`);
}

// ---- 高分子: 重合度 ----
for (const [name, unit, Mu] of [['ポリエチレン', '–CH₂–CH₂–', 28], ['ポリ塩化ビニル', '–CH₂–CHCl–', 62.5], ['ポリスチレン', '–CH₂–CH(C₆H₅)–', 104]]) {
  const n = pick([500, 1000, 2000, 5000]);
  const M = Mu * n;
  add('polymer_basics', `平均分子量 ${M.toLocaleString('en')} の${name}（繰り返し単位 ${unit}、式量 ${Mu}）の平均重合度は？`, String(n),
    [String(n * 2), String(Math.round(n / 2)), String(Math.round(M / 14))],
    `平均重合度 = 平均分子量 / 繰り返し単位の式量 = ${M} / ${Mu} = ${n}。`);
}

const out = items.map((x, i) => ({ id: `c-${String(i + 1).padStart(3, '0')}`, ...x }));
fs.writeFileSync(path.join(__dirname, '..', 'problems', 'calc.json'), JSON.stringify(out, null, 1) + '\n');
console.log(`problems/calc.json に ${out.length} 問: ${[...new Set(out.map((x) => x.topic))].map((t) => `${t} ${out.filter((x) => x.topic === t).length}`).join('・')}`);
