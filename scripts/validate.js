'use strict';
// 問題データの自動検証。仕様書「問題データ形式と自動検証」の項目＋データ整合性を確認する。
// 使い方: node scripts/validate.js            （problems/*.json をすべて検証）
const fs = require('fs');
const path = require('path');
const chem = require('../src/chem');
const calc = require('../src/calc');

const PROBLEM_FILES = ['narrow.json', 'big.json', 'count.json', 'polymer.json', 'generated.json'];

async function loadRDKit() {
  const init = require('@rdkit/rdkit');
  return init();
}

function loadProblems() {
  const dir = path.join(__dirname, '..', 'problems');
  return PROBLEM_FILES.flatMap((f) => {
    const file = path.join(dir, f);
    return fs.existsSync(file) ? JSON.parse(fs.readFileSync(file, 'utf8')) : [];
  });
}

// 候補すべてに全カードを当てた結果表を作る
function applyClues(RDKit, cands, clues, err) {
  const table = [];
  const expected = [];
  clues.forEach((c, ci) => {
    if (!chem.CARDS[c.card]) {
      err(`unknown card ${c.card}`);
      return;
    }
    try {
      expected[ci] = chem.normalizeResult(RDKit, c.card, c.result);
    } catch (e) {
      err(`clue ${c.card}: bad result (${e.message})`);
      return;
    }
    table[ci] = cands.map((s) => {
      try {
        return chem.evaluate(RDKit, c.card, s);
      } catch (e) {
        err(`clue ${c.card} cannot be applied to ${s} (${e.message})`);
        return null;
      }
    });
  });
  return { table, expected };
}

function canonList(RDKit, list, err) {
  try {
    return list.map((s) => chem.canonical(RDKit, s));
  } catch (e) {
    err(e.message);
    return null;
  }
}

// 絞り込み型1問（大問の各段階もこれで検証する）
function checkNarrow(RDKit, p) {
  const errors = [];
  const err = (m) => errors.push(`${p.id}: ${m}`);
  if (!Array.isArray(p.candidates) || p.candidates.length < 2) err('candidates must have 2+ entries');
  if (!Array.isArray(p.clues) || p.clues.length < 1) err('clues must have 1+ entries');
  if (errors.length) return { errors };

  const cands = canonList(RDKit, p.candidates, err);
  if (!cands) return { errors };
  const answer = chem.canonical(RDKit, p.answer);

  // 1. 全候補の分子式が formula と一致する
  cands.forEach((s) => {
    const f = chem.formula(chem.graphFromSmiles(RDKit, s));
    if (f !== p.formula) err(`candidate ${s} is ${f}, expected ${p.formula}`);
  });
  // 2. 候補に重複がない
  if (new Set(cands).size !== cands.length) err('duplicate candidates (canonical SMILES)');
  if (!cands.includes(answer)) err(`answer ${answer} is not among candidates`);

  const { table, expected } = applyClues(RDKit, cands, p.clues, err);
  if (errors.length) return { errors };

  const ok = (ci, k) => chem.consistent(p.clues[ci].card, table[ci][k], expected[ci]);
  const ai = cands.indexOf(answer);
  p.clues.forEach((c, ci) => {
    if (!ok(ci, ai)) err(`clue ${c.card}: answer gives ${JSON.stringify(table[ci][ai])}, data says ${JSON.stringify(expected[ci])}`);
  });
  // 3. 全手がかりを適用すると answer だけが残る
  const survivors = cands.filter((_, k) => p.clues.every((c, ci) => ok(ci, k)));
  if (survivors.length !== 1 || survivors[0] !== answer) {
    err(`after all clues, remaining = ${JSON.stringify(survivors)}`);
  }
  // 4. どの手がかりも単独で候補を1つ以上消す
  p.clues.forEach((c, ci) => {
    if (cands.every((_, k) => ok(ci, k))) err(`clue ${c.card} eliminates nothing`);
  });
  return { errors, cands, answer, table, expected, ok };
}

function atomCounts(RDKit, smiles) {
  const g = chem.graphFromSmiles(RDKit, smiles);
  const c = {};
  g.atoms.forEach((a, i) => {
    c[a.el] = (c[a.el] || 0) + 1;
    c.H = (c.H || 0) + chem.hCount(g, i);
  });
  return c;
}

// 大問型: 加水分解で分けて、断片ごとに絞り込み、最後に組み立てる
function checkBig(RDKit, P) {
  const errors = [];
  const err = (m) => errors.push(`${P.id}: ${m}`);
  let X;
  try {
    X = chem.canonical(RDKit, P.answer);
  } catch (e) {
    err(e.message);
    return { errors };
  }
  const fx = chem.formula(chem.graphFromSmiles(RDKit, X));
  if (fx !== P.formula) err(`X is ${fx}, expected ${P.formula}`);

  // 分解段階: 加水分解の生成物と断片の対応
  const products = chem.hydrolyze(RDKit, chem.graphFromSmiles(RDKit, X));
  if (products.length < 2) err('X must hydrolyze into 2+ molecules');
  const bonds = products.length - 1;
  const counts = {};
  products.forEach((s) => { counts[s] = (counts[s] || 0) + 1; });
  // 原子の収支: X + n H2O = 生成物の合計
  const lhs = atomCounts(RDKit, X);
  lhs.H = (lhs.H || 0) + 2 * bonds;
  lhs.O = (lhs.O || 0) + bonds;
  const rhs = {};
  products.forEach((s) => Object.entries(atomCounts(RDKit, s)).forEach(([el, n]) => { rhs[el] = (rhs[el] || 0) + n; }));
  if (JSON.stringify(Object.entries(lhs).sort()) !== JSON.stringify(Object.entries(rhs).sort())) {
    err(`atom balance mismatch: X+${bonds}H2O=${JSON.stringify(lhs)} vs products ${JSON.stringify(rhs)}`);
  }

  const frags = [];
  const seen = new Set();
  for (const f of P.fragments || []) {
    let a;
    try {
      a = chem.canonical(RDKit, f.answer);
    } catch (e) {
      err(e.message);
      continue;
    }
    if (!counts[a]) err(`fragment ${f.label} (${a}) is not a hydrolysis product of X`);
    if (seen.has(a)) err(`fragment ${f.label} duplicates another fragment`);
    seen.add(a);
    const formula = chem.formula(chem.graphFromSmiles(RDKit, a));
    let r = null;
    if (!f.given) {
      r = checkNarrow(RDKit, { id: `${P.id}/${f.label}`, formula, answer: a, candidates: f.candidates, clues: f.clues });
      errors.push(...r.errors);
    }
    frags.push({ label: f.label, answer: a, formula, count: counts[a] || 0, given: !!f.given, note: f.note || '', r, clues: f.clues || [] });
  }
  Object.keys(counts).forEach((s) => { if (!seen.has(s)) err(`hydrolysis product ${s} has no fragment entry`); });

  let assemble = null;
  if (P.assemble) {
    const r = checkNarrow(RDKit, { id: `${P.id}/X`, formula: P.formula, answer: X, candidates: P.assemble.candidates, clues: P.assemble.clues });
    errors.push(...r.errors);
    // 組み立ての候補は、どれも同じ断片に分かれるものに限る（つなぎ方だけが違う）
    (r.cands || []).forEach((s) => {
      const ps = chem.hydrolyze(RDKit, chem.graphFromSmiles(RDKit, s));
      if (JSON.stringify(ps) !== JSON.stringify(products)) err(`assemble candidate ${s} does not hydrolyze into the same fragments`);
    });
    assemble = { r, clues: P.assemble.clues };
  }
  // 計算段階: 問題文の数値から答えを計算し直して一致を確かめる
  for (const c of P.calcs || []) {
    if (!c.choices || !c.choices.includes(c.answer) || new Set(c.choices).size !== c.choices.length) err(`calc ${c.key}: choices must contain the answer once`);
    if (c.key === 'h2') {
      const M = calc.mass(P.formula);
      const n = Math.round((c.data.V / 22.4) / (c.data.m / M));
      if (n !== c.answer || n !== chem.h2Uptake(chem.graphFromSmiles(RDKit, X))) err(`calc h2: data gives ${n}, answer ${c.answer}`);
    } else if (c.key === 'combustion' || c.key === 'combustion_x') {
      const f = calc.molecularFromEmpirical(calc.empiricalFromCombustion(c.data.sample, c.data.co2, c.data.h2o), { mw: c.data.M });
      const target = c.key === 'combustion_x' ? f === P.formula : frags.some((fr) => fr.formula === f);
      if (f !== c.answer || !target) err(`calc ${c.key}: data gives ${f}, answer ${c.answer}`);
      // 選択肢の中で、データから一意に決まること
      if (c.choices.filter((x) => x === f).length !== 1) err(`calc ${c.key}: answer not uniquely among choices`);
    } else if (c.key === 'n_x') {
      if (!P.assemble || c.answer !== P.assemble.candidates.length) err(`calc n_x: answer ${c.answer} but ${P.assemble ? P.assemble.candidates.length : 0} assembly candidates`);
    } else err(`unknown calc ${c.key}`);
  }
  return { errors, X, bonds, frags, assemble, calcs: P.calcs || [] };
}

// 数え上げ型: 母集団（その分子式のすべての異性体）から条件に合うものを数える
function checkCount(RDKit, P) {
  const errors = [];
  const err = (m) => errors.push(`${P.id}: ${m}`);
  const pool = canonList(RDKit, P.pool || [], err);
  if (!pool) return { errors };
  pool.forEach((s) => {
    const f = chem.formula(chem.graphFromSmiles(RDKit, s));
    if (f !== P.formula) err(`pool member ${s} is ${f}, expected ${P.formula}`);
  });
  if (new Set(pool).size !== pool.length) err('duplicate pool members');
  if (pool.length !== P.poolSize) err(`pool has ${pool.length} members, scope says ${P.poolSize}`);
  const { table, expected } = applyClues(RDKit, pool, P.clues || [], err);
  if (errors.length) return { errors };
  const ok = (ci, k) => chem.consistent(P.clues[ci].card, table[ci][k], expected[ci]);
  const answers = pool.map((_, k) => k).filter((k) => P.clues.every((c, ci) => ok(ci, k)));
  if (answers.length === 0 || answers.length === pool.length) err(`condition selects ${answers.length} of ${pool.length}`);
  P.clues.forEach((c, ci) => {
    if (pool.every((_, k) => ok(ci, k))) err(`clue ${c.card} eliminates nothing`);
  });
  const stereo = pool.map((s, k) => {
    if (!answers.includes(k)) return null;
    try {
      return chem.stereoCount(chem.graphFromSmiles(RDKit, s));
    } catch (e) {
      err(`stereo count for ${s}: ${e.message}`);
      return null;
    }
  });
  const nStereo = stereo.reduce((a, b) => a + (b || 0), 0);
  if (P.expect) {
    if (P.expect[0] !== answers.length) err(`expected ${P.expect[0]} structural isomers, got ${answers.length}`);
    if (P.expect[1] !== nStereo) err(`expected ${P.expect[1]} incl. stereoisomers, got ${nStereo}`);
  }
  return { errors, pool, table, expected, answers, stereo, nStereo };
}

function formulaCounts(f) {
  const c = {};
  for (const [, el, n] of f.matchAll(/([A-Z][a-z]?)(\d*)/g)) c[el] = (c[el] || 0) + (n ? +n : 1);
  return c;
}

// 高分子型: 繰り返し単位と単量体の原子の収支、重合度の計算
function checkPolymer(RDKit, P) {
  const errors = [];
  const err = (m) => errors.push(`${P.id}: ${m}`);
  let unitFormula;
  try {
    const g = chem.graphFromSmiles(RDKit, P.unit);
    if (g.atoms.filter((a) => a.el === 'R').length !== 2) err('unit must have exactly two * ends');
    unitFormula = chem.formula(g);
  } catch (e) {
    err(e.message);
    return { errors };
  }
  const unitMass = chem.formulaMass(unitFormula);
  if (!Number.isInteger(P.n) || P.n < 10) err('n must be an integer >= 10');
  const mw = Math.round(unitMass * P.n * 10) / 10;
  // 単量体の合計 = 繰り返し単位 + (縮合で除かれた H2O)
  const lhs = formulaCounts(unitFormula);
  lhs.H = (lhs.H || 0) + 2 * P.bondsPerUnit;
  lhs.O = (lhs.O || 0) + P.bondsPerUnit;
  const rhs = {};
  const mons = [];
  for (const m of P.monomers) {
    let a;
    try {
      a = chem.canonical(RDKit, m.answer);
    } catch (e) {
      err(e.message);
      continue;
    }
    const formula = chem.formula(chem.graphFromSmiles(RDKit, a));
    Object.entries(formulaCounts(formula)).forEach(([el, n]) => { rhs[el] = (rhs[el] || 0) + n; });
    let r = null;
    if (!m.given) {
      r = checkNarrow(RDKit, { id: `${P.id}/${m.label}`, formula, answer: a, candidates: m.candidates, clues: m.clues });
      errors.push(...r.errors);
    }
    mons.push({ label: m.label, answer: a, formula, given: !!m.given, note: m.note || '', r, clues: m.clues || [] });
  }
  const norm = (o) => JSON.stringify(Object.entries(o).filter(([, n]) => n).sort());
  if (norm(lhs) !== norm(rhs)) err(`atom balance: unit ${unitFormula} + ${P.bondsPerUnit} H2O = ${norm(lhs)} but monomers give ${norm(rhs)}`);
  return { errors, unitFormula, unitMass, mw, mons };
}

function checkAny(RDKit, p) {
  if (p.mode === 'polymer') return checkPolymer(RDKit, p);
  if (p.mode === 'narrow') return checkNarrow(RDKit, p);
  if (p.mode === 'big') return checkBig(RDKit, p);
  if (p.mode === 'count') return checkCount(RDKit, p);
  return { errors: [`${p.id}: unknown mode ${p.mode}`] };
}

function checkAll(RDKit, problems) {
  const errors = [];
  const results = [];
  const ids = new Set();
  for (const p of problems) {
    if (ids.has(p.id)) errors.push(`${p.id}: duplicate id`);
    ids.add(p.id);
    const r = checkAny(RDKit, p);
    errors.push(...r.errors);
    results.push(r);
  }
  return { errors, results };
}

async function main() {
  const problems = loadProblems();
  const RDKit = await loadRDKit();
  const { errors } = checkAll(RDKit, problems);
  if (errors.length) {
    console.error(`NG: ${errors.length} error(s)`);
    errors.forEach((e) => console.error('  ' + e));
    process.exit(1);
  }
  const by = {};
  problems.forEach((p) => { by[p.mode] = (by[p.mode] || 0) + 1; });
  console.log(`OK: ${problems.length} problems passed (${Object.entries(by).map(([m, n]) => `${m} ${n}`).join(', ')})`);
}

if (require.main === module) main();

module.exports = { loadRDKit, loadProblems, checkNarrow, checkBig, checkCount, checkPolymer, checkAny, checkAll };
