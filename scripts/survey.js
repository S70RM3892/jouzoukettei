'use strict';
// 異性体俯瞰型の問題を作り、problems/survey.json に書き出す。
// 分子式の構造異性体を環・三重結合まで全部候補にし、性質と反応の手がかりを2枚以上組み合わせて1つに決める絞り込み問題。
//  1. data/survey.json（京大以外の入試3年分の出題要素）のタグの出た年数から、カードの重みを決める
//  2. data/survey.json の答えの化合物を答えにした問題を先に作り、残りは分子式と答えを変えて作る
// 3年分では京大のような難易度の物差しや出題予測は作れないので、候補の数と手がかりの数で3段階に分ける
// 使い方: node scripts/survey.js [種] [問題数]
const fs = require('fs');
const path = require('path');
const chem = require('../src/chem');
const G = require('../src/generator');
const { enumerate } = require('../src/enumerate');
const { loadRDKit, checkNarrow } = require('./validate');

const ROOT = path.join(__dirname, '..');
const D = JSON.parse(fs.readFileSync(path.join(ROOT, 'data', 'survey.json'), 'utf8'));

// 使うカードと、そのカードで練習できる出題要素（京大で「与えられた規則」として出る反応のうち、教科書の範囲のものだけ）
const CARD_TAGS = {
  acetylide: ['alkyne'], alkyne_hydration: ['alkyne', 'addition_selectivity'], markovnikov: ['addition_selectivity'],
  h2_uptake: ['hydrogenation', 'ring_compound'], hydrogenation: ['hydrogenation', 'isomer_enumeration'], bromine_addition: ['stereo_3d', 'isomer_enumeration'],
  stereo_count: ['stereo_3d', 'optical_resolution'], cis_trans: ['stereo_3d'], chiral: ['optical_resolution'], bromine: ['hydrogenation'],
  silver_mirror: ['sugar'], iodoform: [], sodium: [], mild_oxidation: ['industrial_process'], kmno4: [], dehydration: [], ozonolysis: [], cl_sub: ['isomer_enumeration'],
};
const years = Object.values(D.years);
const freq = (t) => years.filter((y) => y.tags.includes(t)).length / years.length;
const WEIGHTS = Object.fromEntries(Object.entries(CARD_TAGS).map(([c, ts]) => [c, 0.3 + ts.reduce((a, t) => a + freq(t), 0)]));
const CARDS = Object.keys(CARD_TAGS);
// 候補の数が 4〜24 になる分子式（紙に全部書き出せる大きさ）
const FORMULAS = ['C4H6', 'C5H8', 'C4H8', 'C5H10', 'C3H4', 'C3H6O', 'C4H8O', 'C4H10O', 'C3H8O', 'C4H6O'];
const sub = (f) => f.replace(/\d+/g, (d) => [...d].map((x) => '₀₁₂₃₄₅₆₇₈₉'[x]).join(''));

const levelOf = (n, k) => { const s = Math.log2(n) + k; return s < 4.5 ? 1 : s < 6 ? 2 : 3; };

function makeProblem(RDKit, r, pool, answer, id, allowOne = false) {
  const table = G.valueTable(RDKit, pool, CARDS);
  const a = pool.indexOf(answer);
  if (a < 0 || !Object.keys(table).length) return null;
  // 1枚で決まるカードは使わない（性質と反応を組み合わせて考えさせる）。入試の答えを答えにした問題だけは、組がなければ1枚も許す
  const set = G.chooseClues(r, table, a, WEIGHTS, { lo: 2, hi: 3 }) || (allowOne ? G.chooseClues(r, table, a, WEIGHTS, { lo: 1, hi: 3 }) : null);
  if (!set) return null;
  const formula = chem.formula(chem.graphFromSmiles(RDKit, answer));
  const p = {
    id, mode: 'narrow', kind: 'survey', title: `${sub(formula)} の異性体`, level: levelOf(pool.length, set.length), formula, answer,
    candidates: pool, clues: set.map((card) => ({ card, result: chem.evaluate(RDKit, card, answer) })),
  };
  return checkNarrow(RDKit, p).errors.length ? null : p;
}

async function main() {
  const seedStr = process.argv[2] || '20260930';
  const N = Number(process.argv[3]) || 28;
  const RDKit = await loadRDKit();
  const r = G.rng(parseInt(seedStr, 10) || 11);
  const pools = {};
  for (const f of FORMULAS) {
    const pool = enumerate(RDKit, f, { rings: 1, maxRing: 6 }).map((s) => chem.canonical(RDKit, s));
    if (pool.length >= 4 && pool.length <= 24) pools[f] = pool;
  }
  const out = [];
  const used = new Set();
  const nextId = () => `iso-${String(out.length + 1).padStart(3, '0')}`;
  // 1. 入試の答えの化合物を答えにする
  for (const s of Object.values(D.structures)) {
    const X = chem.canonical(RDKit, s);
    const f = chem.formula(chem.graphFromSmiles(RDKit, X));
    if (!pools[f] || used.has(X)) continue;
    const p = makeProblem(RDKit, r, pools[f], X, nextId(), true);
    if (p) { out.push(p); used.add(X); }
  }
  const nKnown = out.length;
  // 2. 分子式と答えを変える（分子式ごとに同じくらいの数）
  const fs_ = Object.keys(pools);
  for (let i = 0; out.length < N && i < N * 40; i++) {
    const f = fs_[i % fs_.length];
    const cand = pools[f].filter((s) => !used.has(s));
    if (!cand.length) continue;
    const X = cand[Math.floor(r() * cand.length)];
    const p = makeProblem(RDKit, r, pools[f], X, nextId());
    if (!p) continue;
    used.add(X);
    out.push(p);
  }
  fs.writeFileSync(path.join(ROOT, 'problems', 'survey.json'), JSON.stringify(out, null, 1) + '\n');
  console.log(`カードの重み: ${Object.entries(WEIGHTS).sort((a, b) => b[1] - a[1]).map(([c, w]) => `${c} ${w.toFixed(2)}`).join(' / ')}`);
  console.log(`problems/survey.json に ${out.length} 問（入試の答えを答えにしたもの ${nKnown}）。分子式: ${fs_.map((f) => `${f}(${pools[f].length})`).join(' ')}`);
  out.forEach((p) => console.log(`  ${p.id} ${p.formula.padEnd(7)} ★${p.level} 候補${p.candidates.length} ${p.answer.padEnd(12)} ${p.clues.map((c) => c.card).join(', ')}`));
}

main();
