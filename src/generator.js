'use strict';
// 京大レベルの大問（構造決定）を自動で作る。
//  1. 断片ライブラリ: 分類ごとに分子式の異性体を全部列挙する（候補集合に漏れがない）
//  2. 型を選んで断片をつなぎ、化合物 X を作る（型の重みは過去問からの出題予測）
//  3. 断片ごとに「答えがちょうど1つに決まる最小の手がかりの組」を探す。1枚で決まる組は捨てる
//  4. つなぎ方が複数あれば、組み立て段階を部分加水分解などで決めさせる
//  5. 水素付加量・燃焼分析の計算段階を付ける
//  6. 難易度を測り、京大の実物を同じ物差しで測った範囲に入るものだけ採用する
const chem = require('./chem');
const calc = require('./calc');
const { enumerate, hasMatch } = require('./enumerate');

// ---------- 乱数（種つき） ----------
function rng(seed) {
  let s = seed >>> 0 || 1;
  return () => {
    s ^= s << 13; s >>>= 0;
    s ^= s >> 17;
    s ^= s << 5; s >>>= 0;
    return s / 4294967296;
  };
}
const pick = (r, arr) => arr[Math.floor(r() * arr.length)];
function weighted(r, items) {
  const tot = items.reduce((a, [, w]) => a + w, 0);
  let x = r() * tot;
  for (const [v, w] of items) { x -= w; if (x <= 0) return v; }
  return items[items.length - 1][0];
}

// ---------- 1. 断片ライブラリ ----------
const BENZ = { seed: 'c1ccccc1', rings: 0 };
const ACYC = { rings: 0 };
const LIB_SPECS = [
  { cls: 'alcohol', label: 'アルコール', formulas: ['CH4O', 'C2H6O', 'C3H8O', 'C4H10O', 'C5H12O', 'C6H14O', 'C4H8O', 'C5H10O', 'C6H12O'], opts: ACYC, only: ['[CX4][OX2H1]'], none: ['C(=O)', 'c', '[OX2]([#6])[#6]'], nOH: 1 },
  { cls: 'aralcohol', label: 'アルコール', formulas: ['C7H8O', 'C8H10O'], opts: BENZ, only: ['[CX4][OX2H1]'], none: ['c[OX2H1]', '[OX2]([#6])[#6]'], nOH: 1 },
  { cls: 'phenol', label: 'フェノール類', formulas: ['C6H6O', 'C7H8O', 'C8H10O'], opts: BENZ, only: ['c[OX2H1]'], none: ['[CX4][OX2H1]', '[OX2]([#6])[#6]'], nOH: 1 },
  { cls: 'diol', label: '二価アルコール', formulas: ['C2H6O2', 'C3H8O2', 'C4H10O2', 'C5H12O2'], opts: ACYC, only: ['[CX4][OX2H1]'], none: ['[OX2]([#6])[#6]'], nOH: 2 },
  { cls: 'acid', label: 'カルボン酸', formulas: ['C2H4O2', 'C3H6O2', 'C4H8O2', 'C5H10O2', 'C3H4O2', 'C4H6O2', 'C5H8O2'], opts: ACYC, only: ['[CX3](=O)[OX2H1]'], nCOOH: 1, extraO: 0 },
  { cls: 'aracid', label: '芳香族カルボン酸', formulas: ['C7H6O2', 'C8H8O2', 'C9H8O2', 'C9H10O2'], opts: BENZ, only: ['[CX3](=O)[OX2H1]'], nCOOH: 1, extraO: 0 },
  { cls: 'naphacid', label: 'ナフタレンカルボン酸', formulas: ['C11H8O2'], opts: { seed: 'c1ccc2ccccc2c1', rings: 0 }, only: ['[CX3](=O)[OX2H1]'], nCOOH: 1, extraO: 0 },
  { cls: 'hydroxyaracid', label: '芳香族ヒドロキシ酸', formulas: ['C7H6O3'], opts: BENZ, only: ['[CX3](=O)[OX2H1]', 'c[OX2H1]'], nCOOH: 1, extraO: 1 },
  { cls: 'diacid', label: '二価カルボン酸', formulas: ['C4H4O4', 'C4H6O4', 'C5H8O4', 'C6H10O4'], opts: ACYC, only: ['[CX3](=O)[OX2H1]'], nCOOH: 2, extraO: 0 },
  { cls: 'ardiacid', label: '芳香族二価カルボン酸', formulas: ['C8H6O4'], opts: BENZ, only: ['[CX3](=O)[OX2H1]'], nCOOH: 2, extraO: 0 },
  { cls: 'amine', label: 'アミン', formulas: ['C2H7N', 'C3H9N', 'C4H11N'], opts: ACYC, only: ['[NX3;H1,H2]'] },
  { cls: 'aniline', label: '芳香族アミン', formulas: ['C6H7N', 'C7H9N'], opts: BENZ, only: ['c[NX3H2]'] },
  { cls: 'aminoaracid', label: '芳香族アミノ酸', formulas: ['C7H7NO2'], opts: BENZ, only: ['[CX3](=O)[OX2H1]', 'c[NX3H2]'], nCOOH: 1 },
  { cls: 'carbonyl', label: 'カルボニル化合物', formulas: ['C3H6O', 'C4H8O', 'C5H10O'], opts: ACYC, only: ['[CX3;!$(C(=O)O)]=O'] },
];

function countMatches(RDKit, s, smarts) {
  const m = RDKit.get_mol(s);
  const q = RDKit.get_qmol(smarts);
  try { const r = JSON.parse(m.get_substruct_matches(q)); return Array.isArray(r) ? r.length : 0; } finally { m.delete(); q.delete(); }
}

function buildLibrary(RDKit) {
  const lib = {};
  for (const spec of LIB_SPECS) {
    for (const f of spec.formulas) {
      let list = enumerate(RDKit, f, spec.opts);
      list = list.filter((s) => spec.only.every((sm) => hasMatch(RDKit, s, sm)));
      if (spec.none) list = list.filter((s) => spec.none.every((sm) => !hasMatch(RDKit, s, sm)));
      if (spec.nOH) list = list.filter((s) => countMatches(RDKit, s, '[OX2H1][#6;!$(C=O)]') === spec.nOH);
      if (spec.nCOOH) list = list.filter((s) => countMatches(RDKit, s, '[CX3](=O)[OX2H1]') === spec.nCOOH);
      if (spec.extraO === 0) list = list.filter((s) => countMatches(RDKit, s, '[OX2H1][#6;!$(C=O)]') === 0 && !hasMatch(RDKit, s, '[OX2]([#6])[#6;!$(C=O)]') && countMatches(RDKit, s, '[#6]=O') === spec.nCOOH);
      if (!list.length) continue;
      const pool = [...new Set(list.map((x) => chem.canonical(RDKit, x)))].sort();
      // 手がかりの組で1つに決められる答えだけを出題に使う（決められないものは問題にならない）
      let solvable = pool.map((_, i) => i);
      if (pool.length >= 3) {
        const table = valueTable(RDKit, pool, FRAG_CARDS);
        const r = rng(7);
        solvable = solvable.filter((a) => fragClues(r, table, a, pool.length, {}));
      }
      lib[`${spec.cls}:${f}`] = { cls: spec.cls, label: spec.label, formula: f, pool, solvable };
    }
  }
  return lib;
}

// ---------- 2. 断片をつなぐ ----------
function toGraph(RDKit, s) { return chem.graphFromSmiles(RDKit, s); }
const nbOf = (g, i) => chem.neighbors(g, i);

function acidSites(g) {
  const out = [];
  g.atoms.forEach((a, c) => {
    if (a.el !== 'C') return;
    const nb = nbOf(g, c);
    const dO = nb.find((x) => x.order === 2 && g.atoms[x.atom].el === 'O');
    const oh = nb.find((x) => x.order === 1 && g.atoms[x.atom].el === 'O' && chem.hCount(g, x.atom) === 1);
    if (dO && oh) out.push({ kind: 'acid', c, o: oh.atom });
  });
  return out;
}
function nucSites(g, cls) {
  const out = [];
  g.atoms.forEach((a, i) => {
    if (a.el === 'O' && chem.hCount(g, i) === 1) {
      const c = nbOf(g, i)[0].atom;
      if (nbOf(g, c).some((x) => x.order === 2 && g.atoms[x.atom].el === 'O')) return; // COOH の OH は除く
      out.push({ kind: 'OH', n: i });
    }
    if (a.el === 'N' && chem.hCount(g, i) >= 1) out.push({ kind: 'NH', n: i });
  });
  if (cls === 'carbonyl') {
    // エノール形の O（ビニルエステルになる）
    g.atoms.forEach((a, c) => {
      const dO = nbOf(g, c).find((x) => x.order === 2 && g.atoms[x.atom].el === 'O');
      if (a.el !== 'C' || !dO) return;
      const alpha = nbOf(g, c).find((x) => g.atoms[x.atom].el === 'C' && chem.hCount(g, x.atom) >= 1);
      if (alpha) out.push({ kind: 'enol', n: dO.atom, c, alpha: alpha.atom });
    });
  }
  return out;
}

function merge(parts) {
  const g = { atoms: [], bonds: [] };
  const offs = [];
  for (const p of parts) {
    offs.push(g.atoms.length);
    p.atoms.forEach((a) => g.atoms.push({ ...a }));
    p.bonds.forEach((b) => g.bonds.push({ a: b.a + offs[offs.length - 1], b: b.b + offs[offs.length - 1], order: b.order }));
  }
  return { g, offs };
}
function removeAtoms(g, dead) {
  const map = new Map();
  const atoms = [];
  g.atoms.forEach((a, i) => { if (!dead.has(i)) { map.set(i, atoms.length); atoms.push(a); } });
  const bonds = g.bonds.filter((b) => !dead.has(b.a) && !dead.has(b.b)).map((b) => ({ a: map.get(b.a), b: map.get(b.b), order: b.order }));
  return { atoms, bonds };
}

// 断片（SMILES と分類）の並びと、酸の位置→求核の位置の対応から X を作る
function assemble(RDKit, frags, pairs) {
  const graphs = frags.map((f) => toGraph(RDKit, f.smiles));
  const { g, offs } = merge(graphs);
  const dead = new Set();
  for (const { ai, as, ni, ns } of pairs) {
    const A = frags[ai].acid[as];
    const N = frags[ni].nuc[ns];
    const ac = A.c + offs[ai];
    dead.add(A.o + offs[ai]);
    const n = N.n + offs[ni];
    if (N.kind === 'enol') {
      const c = N.c + offs[ni];
      const al = N.alpha + offs[ni];
      g.bonds.find((b) => (b.a === c && b.b === n) || (b.a === n && b.b === c)).order = 1;
      g.bonds.find((b) => (b.a === c && b.b === al) || (b.a === al && b.b === c)).order = 2;
    }
    g.bonds.push({ a: ac, b: n, order: 1 });
  }
  const h = removeAtoms(g, dead);
  const smi = chem.canonical(RDKit, chem.toMolblock(h));
  return smi.includes('.') ? null : smi;
}

// すべての完全対応（酸の位置と求核の位置）で X を作り、重複を除く
function allAssemblies(RDKit, frags) {
  const acids = [];
  const nucs = [];
  frags.forEach((f, i) => {
    f.acid.forEach((_, k) => acids.push({ i, k }));
    f.useNuc.forEach((k) => nucs.push({ i, k }));
  });
  if (acids.length !== nucs.length) return [];
  const out = new Set();
  const used = new Array(nucs.length).fill(false);
  const cur = [];
  const rec = (d) => {
    if (d === acids.length) {
      const pairs = cur.map((ni, di) => ({ ai: acids[di].i, as: acids[di].k, ni: nucs[ni].i, ns: nucs[ni].k }));
      if (pairs.some((p) => p.ai === p.ni)) return;
      try { const s = assemble(RDKit, frags, pairs); if (s) out.add(s); } catch (e) { /* つなげない組み合わせ */ }
      return;
    }
    for (let j = 0; j < nucs.length; j++) {
      if (used[j]) continue;
      used[j] = true; cur.push(j);
      rec(d + 1);
      used[j] = false; cur.pop();
    }
  };
  rec(0);
  return [...out];
}

// ---------- 3. 手がかりの自動選択 ----------
const CARD_TAG = {
  silver_mirror: 'silver_mirror', iodoform: 'iodoform', fecl3: 'fecl3', sodium: 'mole_calc', nahco3: 'extraction', naoh: 'extraction', hcl: 'extraction',
  kmno4: 'side_chain_oxidation', mild_oxidation: 'alcohol_oxidation', ozonolysis: 'ozonolysis', kmno4_cleave: 'kmno4_cleavage',
  dehydration: 'dehydration', dehydration_count: 'dehydration', dehydration_ozonolysis: 'dehydration', h2_uptake: 'hydrogenation', hydrogenation: 'hydrogenation',
  chiral: 'chiral', cis_trans: 'cis_trans', carbon_env: 'symmetry_carbons', cl_sub: 'symmetry_carbons', ring_cl: 'symmetry_carbons', anhydride: 'anhydride',
  periodate: 'novel_rule', bromine: 'hydrogenation', ninhydrin: 'amino_acid', alpha_amino: 'amino_acid', partial_hydrolysis: 'partial_hydrolysis',
};
const FRAG_CARDS = ['silver_mirror', 'iodoform', 'fecl3', 'kmno4', 'mild_oxidation', 'ozonolysis', 'kmno4_cleave', 'dehydration', 'dehydration_count',
  'dehydration_ozonolysis', 'h2_uptake', 'hydrogenation', 'chiral', 'cis_trans', 'carbon_env', 'cl_sub', 'ring_cl', 'anhydride', 'periodate', 'bromine', 'naoh', 'hcl'];

function valueTable(RDKit, pool, cards) {
  const table = {};
  for (const card of cards) {
    try {
      const vals = pool.map((s) => JSON.stringify(chem.evaluate(RDKit, card, s)));
      if (new Set(vals).size > 1) table[card] = vals;
    } catch (e) { /* この候補集合には使えないカード */ }
  }
  return table;
}

// 答え a を1つに決める最小の組（大きさ lo..hi）をすべて探し、重みで1つ選ぶ
function chooseClues(r, table, a, weights, { lo = 2, hi = 3 } = {}) {
  const cards = Object.keys(table);
  const n = table[cards[0]] ? table[cards[0]].length : 0;
  const elim = {};
  cards.forEach((c) => { elim[c] = new Set([...Array(n).keys()].filter((k) => table[c][k] !== table[c][a])); });
  const others = [...Array(n).keys()].filter((k) => k !== a);
  const covers = (set) => others.every((k) => set.some((c) => elim[c].has(k)));
  // 1枚で決まるカードは「京大レベル」にならないので使わない（組み合わせて考えさせる）
  const usable = lo >= 2 ? cards.filter((c) => !covers([c])) : cards;
  const found = [];
  const rec = (start, chosen) => {
    if (chosen.length >= lo && covers(chosen)) {
      // 最小性: どれか1枚を抜くと決まらない
      if (chosen.every((c) => !covers(chosen.filter((x) => x !== c)))) found.push([...chosen]);
      return;
    }
    if (chosen.length === hi) return;
    for (let i = start; i < usable.length; i++) {
      if (elim[usable[i]].size === 0) continue;
      chosen.push(usable[i]);
      rec(i + 1, chosen);
      chosen.pop();
    }
  };
  rec(0, []);
  if (!found.length) return null;
  const score = (set) => set.reduce((s, c) => s + (weights[c] || 0.3), 0) + new Set(set.map((c) => CARD_TAG[c])).size * 0.4 + r() * 0.8;
  found.sort((x, y) => score(y) - score(x));
  return found[0];
}

// 断片の手がかり: まず「1枚では決まらない組」を探し、無理なら1枚で決まるカードも許す
function fragClues(r, table, a, n, weights) {
  if (!Object.keys(table).length) return null;
  const lo = n >= 4 ? 2 : 1;
  return chooseClues(r, table, a, weights, { lo, hi: 3 })
    || chooseClues(r, table, a, weights, { lo, hi: 4 })
    || (lo > 1 ? chooseClues(r, table, a, weights, { lo: 1, hi: 3 }) : null);
}

// ---------- 5. 計算段階 ----------
function sig3(x) { return Number(x.toPrecision(3)); }
function h2Stage(r, X, formula, n) {
  const m = pick(r, [5.0, 10.0, 12.0, 15.0]);
  const M = calc.mass(formula);
  const V = sig3((m / M) * n * 22.4);
  return {
    key: 'h2', data: { m, V }, prompt: `X ${m.toFixed(1)} g に白金触媒で水素を付加させると、標準状態で ${V} L の H₂ が消費された。X 1分子がもつ C=C の数はいくつか（ベンゼン環には付加しない）`,
    answer: n, choices: [1, 2, 3, 4], unit: '個',
    explain: `X の分子量は ${M}。${m.toFixed(1)} ÷ ${M} = ${(m / M).toFixed(4)} mol、${V} ÷ 22.4 = ${(V / 22.4).toFixed(4)} mol なので、1 mol あたり ${n} mol の H₂。`,
  };
}
function combustionStage(r, label, formula) {
  const sample = pick(r, [10.0, 15.0, 20.0, 25.0]);
  const c = calc.counts(formula);
  const M = calc.mass(formula);
  const co2 = sig3((sample / M) * c.C * 44);
  const h2o = sig3((sample / M) * (c.H / 2) * 18);
  const emp = calc.empiricalFromCombustion(sample, co2, h2o);
  const mol = calc.molecularFromEmpirical(emp, { mw: M });
  if (mol !== formula) return null; // 丸めた数値で一意に戻らないなら使わない
  const alts = new Set([formula]);
  const tweak = [{ C: 1, H: 2 }, { C: -1, H: -2 }, { O: 1 }, { H: 2 }, { H: -2 }, { O: -1 }];
  while (alts.size < 4) {
    const t = pick(r, tweak);
    const d = { ...c };
    Object.entries(t).forEach(([e, k]) => { d[e] = (d[e] || 0) + k; });
    if (Object.values(d).some((v) => v < 0) || !d.C) continue;
    const order = ['C', 'H', ...Object.keys(d).filter((e) => e !== 'C' && e !== 'H').sort()];
    alts.add(order.filter((e) => d[e]).map((e) => e + (d[e] > 1 ? d[e] : '')).join(''));
  }
  return {
    key: 'combustion', data: { sample, co2, h2o, M: Math.round(M), formula }, prompt: `加水分解で得た${label} ${sample.toFixed(1)} mg を完全燃焼させると CO₂ ${co2} mg と H₂O ${h2o} mg が得られた（分子量は ${Math.round(M)}）。${label} の分子式はどれか`,
    answer: formula, choices: [...alts].sort(), unit: '',
    explain: (() => {
      const mC = co2 * 12 / 44, mH = h2o * 2 / 18, mO = sample - mC - mH;
      const nC = mC / 12, nH = mH / 1.0, nO = mO / 16;
      const base = Math.min(nC, nH, ...(nO > 0.01 ? [nO] : []));
      const ratio = [nC, nH, nO].map((x) => (x / base).toFixed(2));
      return `C = ${co2} × 12/44 = ${mC.toFixed(2)} mg、H = ${h2o} × 2/18 = ${mH.toFixed(2)} mg、O = ${sample.toFixed(1)} − ${mC.toFixed(2)} − ${mH.toFixed(2)} = ${mO.toFixed(2)} mg。`
        + `C : H : O = ${mC.toFixed(2)}/12 : ${mH.toFixed(2)}/1.0 : ${mO.toFixed(2)}/16 = ${ratio.join(' : ')} より組成式 ${subf(emp)}。分子量 ${Math.round(M)} から ${subf(formula)}。`;
    })(),
  };
}

// ---------- 6. 難易度 ----------
function difficulty(p) {
  const fr = p.fragments.filter((f) => !f.given);
  const poolBits = fr.reduce((s, f) => s + Math.log2(f.candidates.length), 0);
  const clues = fr.reduce((s, f) => s + f.clues.length, 0) + (p.assemble ? p.assemble.clues.length : 0);
  const cardTypes = new Set([...fr.flatMap((f) => f.clues.map((c) => c.card)), ...(p.assemble ? p.assemble.clues.map((c) => c.card) : [])]);
  const stereo = [...cardTypes].filter((c) => c === 'chiral' || c === 'cis_trans').length;
  const novel = [...cardTypes].filter((c) => CARD_TAG[c] === 'novel_rule').length;
  const nC = calc.counts(p.formula).C;
  return Math.round((poolBits + 1.2 * clues + 0.6 * cardTypes.size + 1.0 * stereo + 2.0 * novel + (p.assemble ? 2 : 0) + (p.calcs || []).length + 0.15 * nC) * 10) / 10;
}

// ---------- 組み立て全体 ----------
const TEMPLATES = {
  diester_diacid: { title: 'ジカルボン酸のジエステル', tag: 'ester_hydrolysis', parts: [['diacid', 'ardiacid'], ['alcohol', 'aralcohol', 'phenol'], ['alcohol', 'aralcohol', 'phenol']] },
  diester_diol: { title: 'ジオールのジエステル', tag: 'ester_hydrolysis', parts: [['diol'], ['acid', 'aracid', 'naphacid'], ['acid', 'aracid', 'hydroxyaracid', 'naphacid']] },
  ester_amide: { title: 'エステルとアミド', tag: 'amide_hydrolysis', parts: [['aminoaracid'], ['acid'], ['alcohol', 'aralcohol']] },
  diamide: { title: '2つのアミド結合', tag: 'amide_hydrolysis', parts: [['aminoaracid'], ['acid'], ['amine']] },
  triester: { title: 'グリセリンのトリエステル', tag: 'glycerol_ester', parts: [['glycerol'], ['acid', 'aracid', 'hydroxyaracid'], ['acid', 'aracid', 'hydroxyaracid'], ['acid', 'aracid', 'hydroxyaracid']] },
  vinyl: { title: 'ジカルボン酸のジエステル', tag: 'enol_tautomer', parts: [['ardiacid', 'diacid'], ['alcohol'], ['carbonyl']] },
};

function fragmentInfo(RDKit, smiles, cls) {
  const g = toGraph(RDKit, smiles);
  const acid = acidSites(g);
  const nuc = nucSites(g, cls);
  let useNuc;
  if (cls === 'carbonyl') useNuc = [nuc.findIndex((x) => x.kind === 'enol')];
  else if (cls === 'aminoaracid' || cls === 'hydroxyaracid') useNuc = nuc.map((x, i) => (x.kind === 'NH' ? i : -1)).filter((i) => i >= 0).slice(0, 1);
  else if (/acid/.test(cls)) useNuc = [];
  else useNuc = nuc.map((_, i) => i).filter((i) => nuc[i].kind !== 'enol');
  if (cls === 'hydroxyaracid') useNuc = []; // フェノール性 OH は残す（FeCl₃ で見分ける材料）
  return { smiles, cls, acid, nuc, useNuc };
}

const LABELS = ['A', 'B', 'C', 'D', 'E'];
const subf = (f) => f.replace(/\d/g, (d) => '₀₁₂₃₄₅₆₇₈₉'[d]);

function fail(spec, why) { spec.why = why; return null; }

// 型と断片から大問を1つ作る（answerX を渡すと、その X に固定して問題化する＝過去問の再現用）
function buildProblem(RDKit, r, lib, weights, spec) {
  const { frags, id, answerX } = spec;
  const infos = frags.map((f) => fragmentInfo(RDKit, f.smiles, f.cls));
  const alts = allAssemblies(RDKit, infos);
  if (!alts.length) return fail(spec, 'つなげない');
  const X = answerX ? chem.canonical(RDKit, answerX) : pick(r, alts);
  if (!alts.includes(X)) return fail(spec, 'Xが組み立て候補にない');
  const gx = toGraph(RDKit, X);
  const formula = chem.formula(gx);
  const prods = chem.hydrolyze(RDKit, gx);

  // 断片の段階（同じ化合物は1段階）
  const seen = new Set();
  const fragments = [];
  let li = 0;
  for (const f of frags) {
    const can = chem.canonical(RDKit, f.smiles);
    if (!prods.includes(can)) {
      // ビニルエステルはカルボニル化合物として出てくる
      if (!(f.cls === 'carbonyl' && prods.includes(can))) { /* 下で確認 */ }
    }
    if (seen.has(can)) continue;
    seen.add(can);
    const label = LABELS[li++];
    const entry = f.cls === 'glycerol' ? null : lib[`${f.cls}:${chem.formula(chem.graphFromSmiles(RDKit, can))}`];
    if (!entry || entry.pool.length < 3) {
      fragments.push({ label, answer: can, given: true, note: f.cls === 'glycerol' ? 'グリセリン（1,2,3-プロパントリオール）' : `${entry ? entry.label : ''}（この分子式では構造がほぼ1つに決まるので与える）` });
      continue;
    }
    const pool = entry.pool;
    const a = pool.indexOf(can);
    if (a < 0) return fail(spec, `断片${label}が候補集合にない`);
    const table = valueTable(RDKit, pool, FRAG_CARDS);
    const set = fragClues(r, table, a, pool.length, weights);
    if (!set) {
      // カードで決められない断片は、京大と同じく名称を与える（自動出題では選ばれないが、過去問の再現では使う）
      if (!answerX) return fail(spec, `断片${label}の手がかりが作れない(${pool.length})`);
      fragments.push({ label, answer: can, given: true, note: `${entry.label}（名称が与えられる）` });
      continue;
    }
    // 必須ではないが効くカードを1枚まぜる（最少枚数を見極める練習）
    const extra = Object.keys(table).filter((c) => !set.includes(c) && table[c].some((v, k) => k !== a && v !== table[c][a]));
    const cards = extra.length && r() < 0.7 ? [...set, pick(r, extra)] : set;
    const clues = cards.map((c) => ({ card: c, result: chem.evaluate(RDKit, c, can) }))
      .sort(() => r() - 0.5);
    fragments.push({ label, answer: can, candidates: pool, clues, kind: entry.label });
  }
  // 加水分解生成物すべてに断片の段階がある（ビニルエステル由来のカルボニル化合物も含む）
  for (const p of prods) {
    if (!fragments.some((f) => f.answer === p)) return fail(spec, `生成物${p}に段階がない`);
  }

  // 組み立て
  let assemble = null;
  if (alts.length > 1) {
    const pool = alts;
    const a = pool.indexOf(X);
    const cardsX = {};
    try {
      const partials = chem.partialProducts(RDKit, gx);
      for (const pp of partials) {
        const vals = pool.map((s) => {
          try { return chem.partialProducts(RDKit, toGraph(RDKit, s)).includes(pp) ? 'y' : 'n'; } catch (e) { return 'x'; }
        });
        if (new Set(vals).size > 1 && !vals.includes('x')) cardsX[`partial_hydrolysis|${pp}`] = vals;
      }
    } catch (e) { /* 部分加水分解が使えない */ }
    for (const c of ['chiral', 'cis_trans', 'carbon_env']) {
      try {
        const vals = pool.map((s) => JSON.stringify(chem.evaluate(RDKit, c, s)));
        if (new Set(vals).size > 1) cardsX[c] = vals;
      } catch (e) { /* 使えない */ }
    }
    if (!Object.keys(cardsX).length) return fail(spec, '組み立てを決めるカードがない');
    const w2 = { ...weights };
    Object.keys(cardsX).forEach((k) => { if (k.startsWith('partial')) w2[k] = weights.partial_hydrolysis || 0.5; });
    const set = chooseClues(r, cardsX, a, w2, { lo: 1, hi: 2 });
    if (!set) return fail(spec, '組み立ての手がかりが作れない');
    assemble = {
      candidates: pool,
      clues: set.map((k) => (k.startsWith('partial_hydrolysis|') ? { card: 'partial_hydrolysis', result: k.split('|')[1] } : { card: k, result: chem.evaluate(RDKit, k, X) })),
    };
  }

  // 計算段階
  const calcs = [];
  const nH2 = chem.h2Uptake(gx);
  if (nH2 >= 1 && nH2 <= 4) calcs.push(h2Stage(r, X, formula, nH2));
  const small = fragments.filter((f) => !f.given && !/N/.test(f.answer)).sort((p, q) => p.answer.length - q.answer.length)[0];
  if (small && r() < 0.8) {
    const st = combustionStage(r, `化合物 ${small.label}`, chem.formula(toGraph(RDKit, small.answer)));
    if (st) calcs.push(st);
  }

  const bonds = prods.length - 1;
  const hasAmide = /N/.test(formula);
  const names = fragments.map((f) => f.label).join('・');
  const story = `化合物 X（分子式 ${subf(formula)}）は中性の化合物で、${hasAmide ? 'エステル結合とアミド結合' : 'エステル結合'}をもつ。${nH2 ? 'X は臭素水を脱色する。' : ''}X を完全に加水分解すると、化合物 ${names} が得られた。`;
  const problem = {
    id, mode: 'big', level: 3, generated: true,
    title: `自動生成・${(TEMPLATES[spec.template] || {}).title || '大問'}`,
    story, formula, answer: X, fragments, assemble, calcs,
    meta: { bonds, template: spec.template || 'kakomon', frags: frags.map((f) => f.cls) },
  };
  problem.meta.difficulty = difficulty(problem);
  return problem;
}

// 型の重み: 京大で実際に出た型の割合（prior）と、出題予測でのその型の要素の確率を半々で混ぜる
function templateWeights(forecast, prior) {
  const tagP = (t) => (forecast && forecast.probs && forecast.probs[t]) || 0.4;
  const keys = Object.keys(TEMPLATES);
  const tagSum = keys.reduce((a, k) => a + tagP(TEMPLATES[k].tag), 0);
  return keys.map((k) => [k, 0.5 * ((prior && prior[k]) || 0) + 0.5 * (tagP(TEMPLATES[k].tag) / tagSum)]);
}

// tplWeights: { 型: 重み }（generate.js で当てる練習の成績がいちばん良かった選び方）
function sampleSpec(r, lib, tplWeights) {
  const tpl = weighted(r, Object.keys(TEMPLATES).map((k) => [k, (tplWeights && tplWeights[k]) || 1]));
  const frags = [];
  for (const choices of TEMPLATES[tpl].parts) {
    const cls = pick(r, choices);
    if (cls === 'glycerol') { frags.push({ cls, smiles: 'OCC(O)CO' }); continue; }
    const keys = Object.keys(lib).filter((k) => k.startsWith(cls + ':') && lib[k].solvable.length);
    if (!keys.length) return null;
    const entry = lib[pick(r, keys)];
    frags.push({ cls, smiles: entry.pool[pick(r, entry.solvable)] });
  }
  return { template: tpl, frags };
}

function cardWeights(forecast) {
  const w = {};
  for (const [card, tag] of Object.entries(CARD_TAG)) w[card] = 0.3 + ((forecast && forecast.probs && forecast.probs[tag]) || 0.3);
  return w;
}

// 過去問の X を分解して断片の分類を推定（再現と到達範囲の確認用）
function classify(RDKit, lib, smiles) {
  const can = chem.canonical(RDKit, smiles);
  if (can === chem.canonical(RDKit, 'OCC(O)CO')) return 'glycerol';
  const f = chem.formula(chem.graphFromSmiles(RDKit, can));
  for (const [k, e] of Object.entries(lib)) {
    if (e.formula === f && e.pool.some((s) => chem.canonical(RDKit, s) === can)) return e.cls;
  }
  return null;
}

// 過去問の断片の分類の組が、どれかの型でつくれるか（生成できる範囲に入っているか）
function templateFor(classes) {
  for (const [name, t] of Object.entries(TEMPLATES)) {
    if (t.parts.length !== classes.length) continue;
    const used = new Array(classes.length).fill(false);
    const rec = (i) => {
      if (i === t.parts.length) return true;
      for (let j = 0; j < classes.length; j++) {
        if (used[j] || !t.parts[i].includes(classes[j])) continue;
        used[j] = true;
        if (rec(i + 1)) return true;
        used[j] = false;
      }
      return false;
    };
    if (rec(0)) return name;
  }
  return null;
}

module.exports = { CARD_TAG, templateWeights, templateFor, chooseClues, valueTable, FRAG_CARDS, rng, buildLibrary, buildProblem, sampleSpec, cardWeights, difficulty, classify, TEMPLATES };
