'use strict';
// 京大レベルの大問（構造決定）を自動で作る。
//  1. 断片ライブラリ: 分類ごとに分子式の異性体を全部列挙する（候補集合に漏れがない）
//  2. 型を選んで断片をつなぎ、化合物 X を作る（型の重みは過去問からの出題予測）
//  3. 断片ごとに「答えがちょうど1つに決まる最小の手がかりの組」を探す。1枚で決まる組は捨てる
//  4. つなぎ方が複数あれば、組み立て段階を部分加水分解などで決めさせる
//  5. 水素付加量・燃焼分析の計算段階を付ける
//  6. 難易度を測り、京大の実物を同じ物差しで測った範囲に入るものだけ採用する
//  ブレ: 数値（試料の質量）、問題文の言い回し、X の分子式を伏せるか、計算段階の組み合わせ、
//        つなぎ方の候補の数を問うか、を毎回変える。同じ構造でも同じ問題にならないようにする
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
  { cls: 'alcohol', label: 'アルコール', formulas: ['CH4O', 'C2H6O', 'C3H8O', 'C4H10O', 'C5H12O', 'C6H14O'], opts: ACYC, only: ['[CX4][OX2H1]'], none: ['C(=O)', 'c', '[OX2]([#6])[#6]'], nOH: 1 },
  // 不飽和度 1 のアルコール: C=C をもつ鎖状のものと環状のものを同じ候補集合に入れる（京大2024・2025型。H₂ の付加量や臭素水で見分ける）
  { cls: 'alcohol', label: 'アルコール（不飽和または環状）', formulas: ['C4H8O', 'C5H10O', 'C6H12O'], opts: { rings: 1, maxRing: 6 }, only: ['[CX4][OX2H1]'], none: ['C(=O)', 'c', '[OX2]([#6])[#6]', '[r3]', '[r4]'], nOH: 1 },
  { cls: 'hydroxyacid', label: 'ヒドロキシ酸', formulas: ['C2H4O3', 'C3H6O3', 'C4H8O3'], opts: ACYC, only: ['[CX3](=O)[OX2H1]', '[CX4][OX2H1]'], none: ['[CX3;!$(C(=O)O)]=O', '[OX2]([#6])[#6;!$(C=O)]'], nCOOH: 1, nOH: 1 },
  { cls: 'aralcohol', label: 'アルコール', formulas: ['C7H8O', 'C8H10O'], opts: BENZ, only: ['[CX4][OX2H1]'], none: ['c[OX2H1]', '[OX2]([#6])[#6]'], nOH: 1 },
  { cls: 'phenol', label: 'フェノール類', formulas: ['C6H6O', 'C7H8O', 'C8H10O'], opts: BENZ, only: ['c[OX2H1]'], none: ['[CX4][OX2H1]', '[OX2]([#6])[#6]'], nOH: 1 },
  { cls: 'diol', label: '二価アルコール', formulas: ['C2H6O2', 'C3H8O2', 'C4H10O2', 'C5H12O2'], opts: ACYC, only: ['[CX4][OX2H1]'], none: ['[OX2]([#6])[#6]'], nOH: 2 },
  { cls: 'acid', label: 'カルボン酸', formulas: ['C2H4O2', 'C3H6O2', 'C4H8O2', 'C5H10O2', 'C3H4O2', 'C4H6O2', 'C5H8O2'], opts: ACYC, only: ['[CX3](=O)[OX2H1]'], nCOOH: 1, extraO: 0 },
  { cls: 'aracid', label: '芳香族カルボン酸', formulas: ['C7H6O2', 'C8H8O2', 'C9H8O2', 'C9H10O2'], opts: BENZ, only: ['[CX3](=O)[OX2H1]'], nCOOH: 1, extraO: 0 },
  { cls: 'pyacid', label: 'ピリジンカルボン酸', formulas: ['C6H5NO2', 'C7H7NO2'], opts: { seed: 'c1ccncc1', rings: 0 }, only: ['[CX3](=O)[OX2H1]', 'n'], nCOOH: 1, extraO: 0 },
  { cls: 'artriacid', label: '芳香族三価カルボン酸', formulas: ['C9H6O6'], opts: BENZ, only: ['[CX3](=O)[OX2H1]'], nCOOH: 3, extraO: 0 },
  { cls: 'naphacid', label: 'ナフタレンカルボン酸', formulas: ['C11H8O2'], opts: { seed: 'c1ccc2ccccc2c1', rings: 0 }, only: ['[CX3](=O)[OX2H1]'], nCOOH: 1, extraO: 0 },
  { cls: 'hydroxyaracid', label: '芳香族ヒドロキシ酸', formulas: ['C7H6O3'], opts: BENZ, only: ['[CX3](=O)[OX2H1]', 'c[OX2H1]'], nCOOH: 1, extraO: 1 },
  { cls: 'diacid', label: '二価カルボン酸', formulas: ['C4H4O4', 'C4H6O4', 'C5H8O4', 'C6H10O4'], opts: ACYC, only: ['[CX3](=O)[OX2H1]'], nCOOH: 2, extraO: 0 },
  { cls: 'ardiacid', label: '芳香族二価カルボン酸', formulas: ['C8H6O4'], opts: BENZ, only: ['[CX3](=O)[OX2H1]'], nCOOH: 2, extraO: 0 },
  { cls: 'amine', label: 'アミン', formulas: ['C2H7N', 'C3H9N', 'C4H11N'], opts: ACYC, only: ['[NX3;H1,H2]'] },
  { cls: 'aniline', label: '芳香族アミン', formulas: ['C6H7N', 'C7H9N', 'C8H11N', 'C9H13N'], opts: BENZ, only: ['c[NX3H2]'] },
  { cls: 'aminoaracid', label: '芳香族アミノ酸', formulas: ['C7H7NO2'], opts: BENZ, only: ['[CX3](=O)[OX2H1]', 'c[NX3H2]'], nCOOH: 1 },
  { cls: 'carbonyl', label: 'カルボニル化合物', formulas: ['C3H6O', 'C4H8O', 'C5H10O'], opts: ACYC, only: ['[CX3;!$(C(=O)O)]=O'] },
  // 環状のケトン・アルデヒド（京大2003: 六員環のエノールエステル）
  { cls: 'cyclocarbonyl', label: '環状のカルボニル化合物', formulas: ['C6H10O', 'C7H12O'], opts: { rings: 1, maxRing: 6 }, only: ['[CX3;!$(C(=O)O)]=O', '[R]'], none: ['C=C', 'C#C', '[r3]', '[r4]'] }, // 高校で扱う五員環・六員環だけ
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
  if (cls === 'carbonyl' || cls === 'cyclocarbonyl') {
    // エノール形の O（ビニルエステルになる）。α炭素ごとに別のエノールになる（京大2003: 2-メチルシクロヘキサノンの2つのエノールエステル）
    g.atoms.forEach((a, c) => {
      const dO = nbOf(g, c).find((x) => x.order === 2 && g.atoms[x.atom].el === 'O');
      if (a.el !== 'C' || !dO) return;
      nbOf(g, c).filter((x) => g.atoms[x.atom].el === 'C' && chem.hCount(g, x.atom) >= 1)
        .forEach((alpha) => out.push({ kind: 'enol', n: dO.atom, c, alpha: alpha.atom }));
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
  periodate: 'novel_rule', markovnikov: 'addition_selectivity', bromine: 'hydrogenation',
  nitration: 'nitration_reduction', bromine_water: 'symmetry_carbons', acetylation_primary: 'novel_rule', acetonide: 'novel_rule',
  bromine_addition: 'addition_selectivity', stereo_count: 'stereo_3d', ninhydrin: 'amino_acid', alpha_amino: 'amino_acid', partial_hydrolysis: 'partial_hydrolysis',
};
const FRAG_CARDS = ['silver_mirror', 'iodoform', 'fecl3', 'kmno4', 'mild_oxidation', 'ozonolysis', 'kmno4_cleave', 'dehydration', 'dehydration_count',
  'dehydration_ozonolysis', 'h2_uptake', 'hydrogenation', 'markovnikov', 'chiral', 'cis_trans', 'carbon_env', 'cl_sub', 'ring_cl', 'anhydride', 'periodate', 'bromine', 'naoh', 'hcl',
  'nitration', 'bromine_water', 'acetylation_primary', 'acetonide', 'bromine_addition', 'stereo_count'];

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
const f3 = (x) => x.toPrecision(3); // 表示は有効数字3けた（7.70 を 7.7 と書かない）
function h2Stage(r, X, formula, n) {
  // 質量は毎回ばらす（3〜30 g、有効数字3けた）。覚えた数値では解けないようにする
  const m = sig3(3 + r() * 27);
  const M = calc.mass(formula);
  const V = sig3((m / M) * n * 22.4);
  return {
    key: 'h2', data: { m, V }, prompt: `X ${f3(m)} g に白金触媒で水素を付加させると、標準状態で ${f3(V)} L の H₂ が消費された。X 1分子がもつ C=C の数はいくつか（ベンゼン環には付加しない）`,
    answer: n, choices: [1, 2, 3, 4], unit: '個',
    explain: `X の分子量は ${M}。${m} ÷ ${M} = ${(m / M).toFixed(4)} mol、${V} ÷ 22.4 = ${(V / 22.4).toFixed(4)} mol なので、1 mol あたり ${n} mol の H₂。`,
  };
}
function combustionStage(r, label, formula, key = 'combustion') {
  const sample = sig3(5 + r() * 35);
  const c = calc.counts(formula);
  const M = calc.mass(formula);
  const co2 = sig3((sample / M) * c.C * 44);
  const h2o = sig3((sample / M) * (c.H / 2) * 18);
  const emp = calc.empiricalFromCombustion(sample, co2, h2o);
  if (!emp) return null;
  const mol = calc.molecularFromEmpirical(emp, { mw: M });
  if (mol !== formula) return null; // 丸めた数値で一意に戻らないなら使わない
  // 選択肢: 分子量が同じになる組み替え（CH₄ ⇄ O）を優先し、分子量だけでは選べないようにする
  const alts = new Set([formula]);
  const fmt = (d) => ['C', 'H', ...Object.keys(d).filter((e) => e !== 'C' && e !== 'H').sort()].filter((e) => d[e]).map((e) => e + (d[e] > 1 ? d[e] : '')).join('');
  const ok = (d) => Object.values(d).every((v) => v >= 0) && d.C && d.H % 2 === 0 && d.H <= 2 * d.C + 2;
  const apply = (t) => { const d = { ...c }; Object.entries(t).forEach(([e, k]) => { d[e] = (d[e] || 0) + k; }); return d; };
  for (const t of [{ C: 1, H: 4, O: -1 }, { C: -1, H: -4, O: 1 }, { C: 2, H: 8, O: -2 }, { C: -2, H: -8, O: 2 }]) {
    const d = apply(t);
    if (ok(d) && alts.size < 4) alts.add(fmt(d));
  }
  const tweak = [{ C: 1, H: 2 }, { C: -1, H: -2 }, { O: 1 }, { H: 2 }, { H: -2 }, { O: -1 }];
  for (let guard = 0; alts.size < 4 && guard < 50; guard++) {
    const d = apply(pick(r, tweak));
    if (ok(d)) alts.add(fmt(d));
  }
  return {
    key, data: { sample, co2, h2o, M: Math.round(M), formula }, prompt: `${label} ${f3(sample)} mg を完全燃焼させると CO₂ ${f3(co2)} mg と H₂O ${f3(h2o)} mg が得られた（分子量は ${Math.round(M)}）。${label.replace(/^加水分解で得た/, '')} の分子式はどれか`,
    answer: formula, choices: [...alts].sort(), unit: '',
    explain: (() => {
      const mC = co2 * 12 / 44, mH = h2o * 2 / 18, mO = sample - mC - mH;
      const nC = mC / 12, nH = mH / 1.0, nO = mO / 16;
      const base = Math.min(nC, nH, ...(nO > 0.01 ? [nO] : []));
      const ratio = [nC, nH, nO].map((x) => (x / base).toFixed(2));
      return `C = ${co2} × 12/44 = ${mC.toFixed(2)} mg、H = ${h2o} × 2/18 = ${mH.toFixed(2)} mg、O = ${sample} − ${mC.toFixed(2)} − ${mH.toFixed(2)} = ${mO.toFixed(2)} mg。`
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
  // 化合物どうしの関係と誘導体は、全体を同時に考えさせるので重くする
  const ders = (p.derived || []).length + (p.derived || []).reduce((a, d) => a + d.clues.length, 0) * 0.8;
  const rels = (p.relations || []).length * 2.5;
  return Math.round((poolBits + 1.2 * clues + 0.6 * cardTypes.size + 1.0 * stereo + 2.0 * novel + (p.assemble ? 2 : 0) + (p.calcs || []).length + 0.15 * nC + 1.5 * ders + rels) * 10) / 10;
}

// ---------- 組み立て全体 ----------
const TEMPLATES = {
  diester_diacid: { title: 'ジカルボン酸のジエステル', tag: 'ester_hydrolysis', parts: [['diacid', 'ardiacid'], ['alcohol', 'aralcohol', 'phenol'], ['alcohol', 'aralcohol', 'phenol']] },
  diester_diol: { title: 'ジオールのジエステル', tag: 'ester_hydrolysis', parts: [['diol'], ['acid', 'aracid', 'naphacid', 'pyacid'], ['acid', 'aracid', 'hydroxyaracid', 'naphacid', 'pyacid']] },
  ester_amide: { title: 'エステルとアミド', tag: 'amide_hydrolysis', parts: [['aminoaracid'], ['acid'], ['alcohol', 'aralcohol']] },
  diamide: { title: '2つのアミド結合', tag: 'amide_hydrolysis', parts: [['aminoaracid'], ['acid'], ['amine']] },
  triester: { title: 'グリセリンのトリエステル', tag: 'glycerol_ester', parts: [['glycerol'], ['acid', 'aracid', 'hydroxyaracid', 'pyacid'], ['acid', 'aracid', 'hydroxyaracid'], ['acid', 'aracid', 'hydroxyaracid']] },
  vinyl: { title: 'ジカルボン酸のジエステル', tag: 'enol_tautomer', parts: [['ardiacid', 'diacid'], ['alcohol'], ['carbonyl', 'cyclocarbonyl']] },
  // 京大2020: ヒドロキシ酸がジカルボン酸とアルコールの間をつなぐ。部分加水分解の生成物でつなぎ方を決める
  linker: { title: 'ヒドロキシ酸でつないだエステル', tag: 'partial_hydrolysis', parts: [['diacid', 'ardiacid'], ['hydroxyacid'], ['alcohol', 'aralcohol'], ['alcohol', 'aralcohol', 'phenol']] },
  // 京大2009: 芳香族三価カルボン酸にフェノール・アルコール・アミンがつく
  triacid: { title: '三価カルボン酸のエステルとアミド', tag: 'amide_hydrolysis', parts: [['artriacid'], ['phenol', 'alcohol', 'aralcohol'], ['phenol', 'alcohol', 'aralcohol'], ['aniline', 'amine', 'phenol']] },
  // 京大の定番: アルコールと、それを酸化したカルボン酸（または同じ骨格の別の官能基）からできたエステル
  ester_pair: { title: 'エステル（酸とアルコールの関係）', tag: 'alcohol_oxidation', parts: [['hydroxyacid'], ['acid'], ['alcohol']] },
  // 京大2003: 環状ケトンのエノールエステル。どちら側のα炭素でエノールになったかを不斉炭素などで決める
  enol_ring: { title: '環状ケトンのエノールエステル', tag: 'enol_tautomer', parts: [['acid', 'aracid'], ['cyclocarbonyl']] },
};

function fragmentInfo(RDKit, smiles, cls) {
  const g = toGraph(RDKit, smiles);
  const acid = acidSites(g);
  const nuc = nucSites(g, cls);
  let useNuc;
  // エノールになる位置が複数あれば、そのどれか1つを使う（variants で全部試す）
  let variants = null;
  if (cls === 'carbonyl' || cls === 'cyclocarbonyl') {
    variants = nuc.map((x, i) => (x.kind === 'enol' ? i : -1)).filter((i) => i >= 0).map((i) => [i]);
    useNuc = variants[0] || [];
  } else if (cls === 'hydroxyacid') useNuc = nuc.map((x, i) => (x.kind === 'OH' ? i : -1)).filter((i) => i >= 0);
  else if (cls === 'aminoaracid' || cls === 'hydroxyaracid') useNuc = nuc.map((x, i) => (x.kind === 'NH' ? i : -1)).filter((i) => i >= 0).slice(0, 1);
  else if (/acid/.test(cls)) useNuc = [];
  else useNuc = nuc.map((_, i) => i).filter((i) => nuc[i].kind !== 'enol');
  if (cls === 'hydroxyaracid') useNuc = []; // フェノール性 OH は残す（FeCl₃ で見分ける材料）
  return { smiles, cls, acid, nuc, useNuc, variants };
}

// エノールの位置の選び方をすべて組み合わせて X の候補を集める
function allAssembliesWithVariants(RDKit, infos) {
  let combos = [infos];
  infos.forEach((f, i) => {
    if (!f.variants || f.variants.length < 2) return;
    combos = combos.flatMap((c) => f.variants.map((v) => c.map((g, j) => (j === i ? { ...g, useNuc: v } : g))));
  });
  const out = new Set();
  for (const c of combos) allAssemblies(RDKit, c).forEach((x) => out.add(x));
  return [...out];
}

const LABELS = ['A', 'B', 'C', 'D', 'E'];
const subf = (f) => f.replace(/\d/g, (d) => '₀₁₂₃₄₅₆₇₈₉'[d]);

function fail(spec, why) { spec.why = why; return null; }

// 型と断片から大問を1つ作る（answerX を渡すと、その X に固定して問題化する＝過去問の再現用）
function buildProblem(RDKit, r, lib, weights, spec) {
  const { frags, id, answerX } = spec;
  const infos = frags.map((f) => fragmentInfo(RDKit, f.smiles, f.cls));
  const alts = allAssembliesWithVariants(RDKit, infos);
  if (!alts.length) return fail(spec, 'つなげない');
  // 立体の表記（/ \\ @）は外す。組み立ての候補は構造異性体として作っている
  const X = answerX ? chem.canonical(RDKit, answerX.replace(/[/\\@]/g, '')) : pick(r, alts);
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

  // 計算段階（ブレ: どれを出すか・数値は毎回変える）
  const calcs = [];
  // 京大2020・2024型: X の分子式を伏せ、分子量と燃焼分析から求めさせる（N を含む X は燃焼で N が決まらないので除く）
  let hideFormula = false;
  if (!/N/.test(formula) && r() < 0.35) {
    const st = combustionStage(r, 'X', formula, 'combustion_x');
    if (st) { calcs.push(st); hideFormula = true; }
  }
  const nH2 = chem.h2Uptake(gx);
  if (nH2 >= 1 && nH2 <= 4 && r() < 0.85) calcs.push(h2Stage(r, X, formula, nH2));
  const small = fragments.filter((f) => !f.given && !/N/.test(f.answer)).sort((p, q) => p.answer.length - q.answer.length)[0];
  if (small && !hideFormula && r() < 0.7) {
    const st = combustionStage(r, `加水分解で得た化合物 ${small.label}`, chem.formula(toGraph(RDKit, small.answer)));
    if (st) calcs.push(st);
  }
  // 京大2001型: 断片がそろった段階で、つなぎ方の候補が何通りあるかを数えさせる
  if (assemble && alts.length >= 3 && r() < 0.6) {
    const n = alts.length;
    const choices = [...new Set([n - 2, n - 1, n, n + 1, n + 2].filter((x) => x >= 1))].slice(0, 4);
    if (!choices.includes(n)) choices[choices.length - 1] = n;
    calcs.push({
      key: 'n_x', at: 'assemble', data: { n },
      prompt: `断片 ${fragments.map((f) => f.label).join('・')} の構造がすべて決まった。つなぎ方だけが違う X の候補は何種類あるか（立体異性体は区別しない）`,
      answer: n, choices: choices.sort((a, b) => a - b), unit: '種類',
      explain: `酸の側（カルボキシ基）と、アルコール・フェノール・アミン・エノールの側を1対1に対応させる組み合わせを、同じ化合物になるものを除いて数えると ${n} 種類。`,
    });
  }

  const bonds = prods.length - 1;
  const hasAmide = /N/.test(formula);
  const names = fragments.map((f) => f.label).join('・');
  const M = Math.round(calc.mass(formula));
  const head = hideFormula ? `化合物 X（分子量 ${M}）` : pick(r, [`化合物 X（分子式 ${subf(formula)}）`, `分子式 ${subf(formula)} の化合物 X`]);
  const neutral = pick(r, ['は中性の化合物で、', 'は NaHCO₃ 水溶液にも希塩酸にも溶けない化合物で、', 'は水に溶けにくい中性の化合物で、']);
  const bondsTxt = hasAmide ? 'エステル結合とアミド結合' : pick(r, ['エステル結合', '複数のエステル結合']);
  const unsat = nH2 ? pick(r, ['X は臭素水を脱色する。', 'X に臭素水を加えると、臭素の色が消えた。', '']) : pick(r, ['X は臭素水を脱色しない。', '']);
  const story = `${head}${neutral}${bondsTxt}をもつ。${unsat}X を完全に加水分解すると、化合物 ${names} が得られた。`;
  const problem = {
    id, mode: 'big', level: 3, generated: true,
    title: `自動生成・${(TEMPLATES[spec.template] || {}).title || '大問'}`,
    story, formula, answer: X, fragments, assemble, calcs, hideFormula,
    meta: { bonds, template: spec.template || 'kakomon', frags: frags.map((f) => f.cls) },
  };
  // 京大型: 化合物どうしの関係・誘導体を混ぜて、全体を同時に考えないと決まらない手がかりの組に作り直す
  if (spec.chain !== false) {
    const Ch = require('./chain');
    let d = null;
    for (let t = 0; t < 4 && !d; t++) d = Ch.design(RDKit, r, problem);
    if (d && (d.derived.length || d.relations.length)) {
      problem.fragments = d.fragments;
      problem.derived = d.derived;
      problem.relations = d.relations;
      problem.chain = true;
    }
  }
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

// 化合物どうしの関係（A を酸化すると C、A と B を水素付加すると同じ化合物）が生まれるように、断片の1つを選び直す
const REL = ['kmno4', 'hydrogenation', 'dehydration', 'markovnikov', 'kmno4_cleave', 'mild_oxidation'];
function biasRelations(RDKit, r, lib, spec) {
  const out = { ...spec, frags: spec.frags.map((f) => ({ ...f })) };
  const opOf = (o, s) => { try { const v = chem.evaluate(RDKit, o, s); return Array.isArray(v) && v.length === 1 ? chem.canonical(RDKit, v[0]) : null; } catch (e) { return null; } };
  // 候補が3つ以上ある（与えられない）断片の分類の、出題できる構造
  const poolOf = (cls) => Object.values(lib).filter((e) => e.cls === cls && e.pool.length >= 3).flatMap((e) => e.solvable.map((k) => e.pool[k]));
  const tries = [];
  out.frags.forEach((fi, i) => out.frags.forEach((fj, j) => { if (i !== j && fi.cls !== 'glycerol' && fj.cls !== 'glycerol') REL.forEach((o) => tries.push([i, j, o])); }));
  tries.sort(() => r() - 0.5);
  for (const [i, j, o] of tries) {
    const self = chem.canonical(RDKit, out.frags[i].smiles);
    const other = chem.canonical(RDKit, out.frags[j].smiles);
    // 向き1: B を置き換える。A を o すると B（yields）か、A と B を o すると同じ化合物（same）
    const pi = opOf(o, self);
    if (pi) {
      const hits = poolOf(out.frags[j].cls).filter((m) => m !== self && (m === pi || (o !== 'mild_oxidation' && opOf(o, m) === pi)));
      if (hits.length) { out.frags[j].smiles = pick(r, hits); return out; }
    }
    // 向き2: A を置き換える。B（与えられた酸などでもよい）が A から o で得られるように
    const src = poolOf(out.frags[i].cls).filter((m) => m !== other && opOf(o, m) === other);
    if (src.length) { out.frags[i].smiles = pick(r, src); return out; }
  }
  return spec;
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

module.exports = { biasRelations, fragmentInfo, allAssembliesWithVariants, CARD_TAG, templateWeights, templateFor, chooseClues, valueTable, FRAG_CARDS, rng, buildLibrary, buildProblem, sampleSpec, cardWeights, difficulty, classify, TEMPLATES };
