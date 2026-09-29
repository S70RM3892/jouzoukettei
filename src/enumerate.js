'use strict';
// 構造異性体の列挙（標準的な正規化による重複除去つきの原子追加法）。
// 「すべて書け」「何種類あるか」の母集団と、自動生成する問題の候補集合をここから作る。
// 立体異性体は区別しない（構造異性体だけ）。
const chem = require('./chem');

const VAL = { C: 4, N: 3, O: 2, S: 2, Cl: 1, Br: 1 };

function parseFormula(f) {
  const c = {};
  for (const [, el, n] of f.matchAll(/([A-Z][a-z]?)(\d*)/g)) c[el] = (c[el] || 0) + (n ? +n : 1);
  return c;
}

function freeVal(g, i) {
  return VAL[g.atoms[i].el] - g.bonds.reduce((s, b) => s + (b.a === i || b.b === i ? b.order : 0), 0);
}

function key(RDKit, g) {
  return chem.canonical(RDKit, chem.toMolblock(g));
}

// 高校化学で扱わない・不安定な構造を除く
function acceptable(g, opts) {
  const nb = (i) => g.bonds.filter((b) => b.a === i || b.b === i).map((b) => ({ atom: b.a === i ? b.b : b.a, order: b.order }));
  for (let i = 0; i < g.atoms.length; i++) {
    const el = g.atoms[i].el;
    const ns = nb(i);
    if (el === 'O' || el === 'S') {
      if (ns.some((x) => g.atoms[x.atom].el === 'O' || g.atoms[x.atom].el === 'N' || g.atoms[x.atom].el === 'S')) return false; // 過酸化物など
      // エノール（C=C-OH）は除く
      if (ns.length === 1 && ns[0].order === 1 && freeVal(g, i) === 1 && !g.atoms[ns[0].atom].arom) {
        const c = ns[0].atom;
        if (nb(c).some((y) => y.order === 2 && g.atoms[y.atom].el === 'C')) return false;
      }
    }
    if (el === 'N') {
      if (ns.some((x) => g.atoms[x.atom].el === 'N' || g.atoms[x.atom].el === 'O')) return false;
      // イミン・ニトリルは扱わない（ピリジン環の N は芳香環の一部なので許す）
      if (!g.atoms[i].arom && ns.some((x) => x.order > 1) && !opts.allowImine) return false;
    }
    if (el === 'C') {
      const dbl = ns.filter((x) => x.order === 2).length;
      if (dbl >= 2) return false; // アレン・ケテン
      // gem-ジオール、ヘミアセタール（OH とエーテル O を同じ C に）、C(OH)=... は除く
      const oh = ns.filter((x) => g.atoms[x.atom].el === 'O' && x.order === 1 && freeVal(g, x.atom) === 1).length;
      const or = ns.filter((x) => g.atoms[x.atom].el === 'O' && x.order === 1).length;
      const co = ns.filter((x) => g.atoms[x.atom].el === 'O' && x.order === 2).length;
      if (!opts.allowAcetal && !co && or >= 2) return false; // アセタール・ヘミアセタール・gem-ジオール
      if (oh && ns.some((x) => g.atoms[x.atom].el === 'N')) return false; // ヘミアミナール
      if (ns.some((x) => x.order === 3) && ns.some((x) => g.atoms[x.atom].el !== 'C')) return false;
    }
  }
  return true;
}

/**
 * 分子式 formula の構造異性体を列挙する。
 * opts.rings: 許す環の数（既定 1）、opts.maxRing: 環の最大員数（既定 6）
 * opts.seed: 出発骨格の SMILES（例 'c1ccccc1'）。指定するとその骨格を含むものだけを作る
 * opts.filter: 最後に通す関数 (smiles) => boolean
 */
function enumerate(RDKit, formula, opts = {}) {
  const o = { rings: 1, maxRing: 6, allowImine: false, allowAcetal: false, ...opts };
  const target = parseFormula(formula);
  const hTarget = target.H || 0;
  const heavy = [];
  for (const [el, n] of Object.entries(target)) if (el !== 'H') for (let k = 0; k < n; k++) heavy.push(el);

  let start;
  const remaining = { ...target };
  delete remaining.H;
  if (o.seed) {
    start = chem.graphFromSmiles(RDKit, o.seed);
    start.atoms.forEach((a) => { remaining[a.el] -= 1; });
    if (Object.values(remaining).some((n) => n < 0)) return [];
    start.fixed = start.bonds.length; // 骨格の結合は変えない
  } else {
    const first = remaining.C ? 'C' : Object.keys(remaining)[0];
    remaining[first] -= 1;
    start = { atoms: [{ el: first }], bonds: [], fixed: 0 };
  }
  const seedRings = start.bonds.length - start.atoms.length + 1;
  const unsat = (target.C || 0) - hTarget / 2 + (target.N || 0) / 2 + 1 - ((target.Cl || 0) + (target.Br || 0)) / 2;
  if (unsat < 0 || !Number.isInteger(unsat)) return [];

  const rem = (g) => {
    const r = { ...target };
    delete r.H;
    g.atoms.forEach((a) => { r[a.el] -= 1; });
    return r;
  };
  const usedUnsat = (g) => {
    const rings = g.bonds.length - g.atoms.length + 1;
    return rings + g.bonds.reduce((s, b) => s + (b.order - 1), 0);
  };
  const ringCount = (g) => g.bonds.length - g.atoms.length + 1;

  let level = new Map([[key(RDKit, start), start]]);
  const done = new Map();
  const total = heavy.length;
  const clone = (g) => ({ atoms: g.atoms.map((a) => ({ ...a })), bonds: g.bonds.map((b) => ({ ...b })), fixed: g.fixed });

  // 1 段階ずつ: 原子を足すか、結合次数を上げるか、環を閉じる
  for (let guard = 0; guard < total * 4 + 10 && level.size; guard++) {
    const next = new Map();
    for (const g of level.values()) {
      const r = rem(g);
      const u = usedUnsat(g);
      if (g.atoms.length === total && u === unsat) {
        if (acceptable(g, o)) done.set(key(RDKit, g), g);
        continue;
      }
      const push = (h) => {
        if (usedUnsat(h) > unsat) return;
        if (!acceptable(h, { ...o, partial: true })) { /* 途中段階では許す */ }
        let k;
        try { k = key(RDKit, h); } catch (e) { return; }
        if (!next.has(k)) next.set(k, h);
      };
      // 原子を足す
      for (const [el, n] of Object.entries(r)) {
        if (n <= 0) continue;
        for (let i = 0; i < g.atoms.length; i++) {
          if (freeVal(g, i) < 1) continue;
          const h = clone(g);
          h.atoms.push({ el });
          h.bonds.push({ a: i, b: h.atoms.length - 1, order: 1 });
          push(h);
        }
      }
      if (u < unsat) {
        // 結合次数を上げる（骨格の結合は除く）
        g.bonds.forEach((b, bi) => {
          if (bi < g.fixed || b.order >= 3) return;
          if (freeVal(g, b.a) < 1 || freeVal(g, b.b) < 1) return;
          const h = clone(g);
          h.bonds[bi].order += 1;
          push(h);
        });
        // 環を閉じる
        if (ringCount(g) - seedRings < o.rings) {
          for (let i = 0; i < g.atoms.length; i++) {
            for (let j = i + 1; j < g.atoms.length; j++) {
              if (freeVal(g, i) < 1 || freeVal(g, j) < 1) continue;
              if (g.bonds.some((b) => (b.a === i && b.b === j) || (b.a === j && b.b === i))) continue;
              const d = pathLen(g, i, j);
              if (d < 2 || d + 1 > o.maxRing) continue;
              const h = clone(g);
              h.bonds.push({ a: i, b: j, order: 1 });
              push(h);
            }
          }
        }
      }
    }
    level = next;
  }
  let out = [...done.keys()];
  if (o.seed) {
    const q = RDKit.get_qmol(o.seed);
    out = out.filter((s) => {
      const m = RDKit.get_mol(s);
      try { return JSON.parse(m.get_substruct_match(q)).atoms !== undefined; } finally { m.delete(); }
    });
    q.delete();
  }
  if (o.filter) out = out.filter(o.filter);
  return out.filter((s) => chem.formula(chem.graphFromSmiles(RDKit, s)) === formulaString(target)).sort();
}

function pathLen(g, s, t) {
  const dist = new Map([[s, 0]]);
  const q = [s];
  while (q.length) {
    const x = q.shift();
    if (x === t) return dist.get(x);
    for (const b of g.bonds) {
      const y = b.a === x ? b.b : b.b === x ? b.a : -1;
      if (y >= 0 && !dist.has(y)) { dist.set(y, dist.get(x) + 1); q.push(y); }
    }
  }
  return Infinity;
}

function formulaString(c) {
  const order = ['C', 'H', ...Object.keys(c).filter((e) => e !== 'C' && e !== 'H').sort()];
  return order.filter((e) => c[e]).map((e) => e + (c[e] > 1 ? c[e] : '')).join('');
}

// よく使う分類
function hasMatch(RDKit, smiles, smarts) {
  const m = RDKit.get_mol(smiles);
  const q = RDKit.get_qmol(smarts);
  try { return JSON.parse(m.get_substruct_match(q)).atoms !== undefined; } finally { m.delete(); q.delete(); }
}

module.exports = { enumerate, parseFormula, hasMatch };
