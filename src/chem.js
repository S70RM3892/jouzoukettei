'use strict';
// 手がかりカードの判定ロジック（Node の検証・ビルド専用）。
// 各カードは「候補1つ → 結果」の関数。結果が問題データの値と一致しない候補が除外される。
// 対象外の構造（例: C=C を持つ分子に KMnO4）では Unsupported を投げ、その問題は検証で落ちる。

class Unsupported extends Error {}

const VALENCE = { C: 4, O: 2, N: 3, S: 2, F: 1, Cl: 1, Br: 1, I: 1, H: 1 };

// ---------- 分子グラフ ----------

function parseMolblock(mb) {
  const lines = mb.split('\n');
  const counts = lines[3];
  const na = parseInt(counts.slice(0, 3), 10);
  const nb = parseInt(counts.slice(3, 6), 10);
  const atoms = [];
  for (let i = 0; i < na; i++) {
    const l = lines[4 + i];
    atoms.push({ el: l.slice(31, 34).trim() });
  }
  const bonds = [];
  for (let i = 0; i < nb; i++) {
    const l = lines[4 + na + i];
    const a = parseInt(l.slice(0, 3), 10) - 1;
    const b = parseInt(l.slice(3, 6), 10) - 1;
    const order = parseInt(l.slice(6, 9), 10);
    if (order < 1 || order > 3) throw new Unsupported('aromatic/query bond in molblock');
    bonds.push({ a, b, order });
  }
  for (const l of lines) {
    if (l.startsWith('M  CHG')) throw new Unsupported('charged species');
  }
  return { atoms, bonds };
}

function graphFromSmiles(RDKit, smiles) {
  const mol = RDKit.get_mol(smiles);
  if (!mol || !mol.is_valid()) throw new Error(`invalid SMILES: ${smiles}`);
  try {
    return parseMolblock(mol.get_molblock());
  } finally {
    mol.delete();
  }
}

function cloneGraph(g) {
  return { atoms: g.atoms.map((a) => ({ ...a })), bonds: g.bonds.map((b) => ({ ...b })) };
}

function neighbors(g, i) {
  const out = [];
  for (const b of g.bonds) {
    if (b.a === i) out.push({ atom: b.b, order: b.order, bond: b });
    else if (b.b === i) out.push({ atom: b.a, order: b.order, bond: b });
  }
  return out;
}

function hCount(g, i) {
  const v = VALENCE[g.atoms[i].el];
  if (v === undefined) throw new Unsupported(`element ${g.atoms[i].el}`);
  const used = neighbors(g, i).reduce((s, n) => s + n.order, 0);
  const h = v - used;
  if (h < 0) throw new Unsupported('hypervalent atom');
  return h;
}

function formula(g) {
  const c = {};
  const add = (el, n) => { c[el] = (c[el] || 0) + n; };
  g.atoms.forEach((a, i) => { add(a.el, 1); add('H', hCount(g, i)); });
  const order = ['C', 'H', ...Object.keys(c).filter((e) => e !== 'C' && e !== 'H').sort()];
  return order.filter((e) => c[e]).map((e) => e + (c[e] > 1 ? c[e] : '')).join('');
}

function toMolblock(g) {
  const pad = (s, n) => String(s).padStart(n);
  const lines = ['', '  jouzoukettei', ''];
  lines.push(`${pad(g.atoms.length, 3)}${pad(g.bonds.length, 3)}  0  0  0  0  0  0  0  0999 V2000`);
  for (const a of g.atoms) {
    lines.push(`    0.0000    0.0000    0.0000 ${a.el.padEnd(3)} 0  0  0  0  0  0  0  0  0  0  0  0`);
  }
  for (const b of g.bonds) lines.push(`${pad(b.a + 1, 3)}${pad(b.b + 1, 3)}${pad(b.order, 3)}  0`);
  lines.push('M  END');
  return lines.join('\n');
}

function components(g) {
  const seen = new Array(g.atoms.length).fill(-1);
  let k = 0;
  for (let s = 0; s < g.atoms.length; s++) {
    if (seen[s] !== -1) continue;
    const stack = [s];
    seen[s] = k;
    while (stack.length) {
      const x = stack.pop();
      for (const n of neighbors(g, x)) if (seen[n.atom] === -1) { seen[n.atom] = k; stack.push(n.atom); }
    }
    k++;
  }
  const parts = [];
  for (let c = 0; c < k; c++) {
    const idx = [];
    seen.forEach((v, i) => { if (v === c) idx.push(i); });
    const map = new Map(idx.map((old, nw) => [old, nw]));
    parts.push({
      atoms: idx.map((i) => ({ ...g.atoms[i] })),
      bonds: g.bonds.filter((b) => map.has(b.a)).map((b) => ({ a: map.get(b.a), b: map.get(b.b), order: b.order })),
    });
  }
  return parts;
}

function canonical(RDKit, smilesOrMolblock) {
  const mol = RDKit.get_mol(smilesOrMolblock);
  if (!mol || !mol.is_valid()) throw new Error(`invalid structure: ${smilesOrMolblock}`);
  try {
    return mol.get_smiles();
  } finally {
    mol.delete();
  }
}

function graphToSmilesList(RDKit, g) {
  return components(g).map((p) => canonical(RDKit, toMolblock(p))).sort();
}

// ---------- 対称性（色の精密化） ----------
// 水素を明示した分子グラフで Weisfeiler-Lehman 精密化を行い、同じ色 = 同じ置換基とみなす。
// fixed に入れた原子は固有色にして、その原子を固定した上での同値性を見る。

function refine(g, fixed) {
  const nodes = g.atoms.map((a) => a.el);
  const adj = g.atoms.map(() => []);
  g.bonds.forEach((b) => { adj[b.a].push([b.b, b.order]); adj[b.b].push([b.a, b.order]); });
  const heavy = g.atoms.length;
  for (let i = 0; i < heavy; i++) {
    const h = hCount(g, i);
    for (let k = 0; k < h; k++) {
      const hi = nodes.length;
      nodes.push('H');
      adj.push([[i, 1]]);
      adj[i].push([hi, 1]);
    }
  }
  let color = nodes.map((el, i) => (fixed.includes(i) ? `*${fixed.indexOf(i)}` : el));
  for (let it = 0; it < nodes.length; it++) {
    const sig = color.map((c, i) => c + '|' + adj[i].map(([j, o]) => o + color[j]).sort().join(','));
    const keys = [...new Set(sig)].sort();
    const next = sig.map((s) => String(keys.indexOf(s)));
    const stable = new Set(next).size === new Set(color).size;
    color = next;
    if (stable) break;
  }
  return { color, adj };
}

function chiralCount(g) {
  let n = 0;
  g.atoms.forEach((a, i) => {
    if (a.el !== 'C') return;
    const nb = neighbors(g, i);
    if (nb.some((x) => x.order !== 1)) return;
    if (nb.length + hCount(g, i) !== 4 || hCount(g, i) > 1) return;
    const { color, adj } = refine(g, [i]);
    const cs = adj[i].map(([j]) => color[j]);
    if (new Set(cs).size === 4) n++;
  });
  return n;
}

function bondInSmallRing(g, bond, maxSize) {
  // bond を除いて a→b の最短経路を探す
  const dist = new Map([[bond.a, 0]]);
  const q = [bond.a];
  while (q.length) {
    const x = q.shift();
    for (const n of neighbors(g, x)) {
      if (n.bond === bond || dist.has(n.atom)) continue;
      dist.set(n.atom, dist.get(x) + 1);
      q.push(n.atom);
    }
  }
  return dist.has(bond.b) && dist.get(bond.b) + 1 < maxSize;
}

function hasCisTrans(g) {
  for (const bond of g.bonds) {
    if (bond.order !== 2) continue;
    const { a, b } = bond;
    if (g.atoms[a].el !== 'C' || g.atoms[b].el !== 'C') continue;
    if (bondInSmallRing(g, bond, 8)) continue;
    const { color, adj } = refine(g, [a, b]);
    const ok = [a, b].every((x) => {
      const other = adj[x].filter(([j]) => j !== a && j !== b);
      if (other.length !== 2) return false; // アレン等は除外
      return color[other[0][0]] !== color[other[1][0]];
    });
    if (ok) return true;
  }
  return false;
}

// ---------- 変換反応 ----------

function addAtom(g, el) {
  g.atoms.push({ el });
  return g.atoms.length - 1;
}

function hasCC(g, order) {
  return g.bonds.some((b) => b.order === order && g.atoms[b.a].el === 'C' && g.atoms[b.b].el === 'C');
}

function isCarbonyl(g, c) {
  return neighbors(g, c).some((n) => n.order === 2 && g.atoms[n.atom].el === 'O');
}

// 硫酸酸性 KMnO4: 第一級アルコール・アルデヒド → カルボン酸、第二級アルコール → ケトン、第三級は反応しない
function oxidize(RDKit, g0) {
  const g = cloneGraph(g0);
  if (hasCC(g, 2) || hasCC(g, 3)) throw new Unsupported('KMnO4 with C=C / C#C');
  let changed = false;
  const n0 = g.atoms.length;
  for (let c = 0; c < n0; c++) {
    if (g.atoms[c].el !== 'C') continue;
    const nb = neighbors(g, c);
    const carbonNb = nb.filter((x) => g.atoms[x.atom].el === 'C').length;
    const h = hCount(g, c);
    if (isCarbonyl(g, c)) {
      // アルデヒド基 → カルボキシ基
      const singleO = nb.filter((x) => x.order === 1 && g.atoms[x.atom].el === 'O');
      if (h === 1 && singleO.length === 0) {
        if (carbonNb === 0) throw new Unsupported('formaldehyde oxidation');
        const o = addAtom(g, 'O');
        g.bonds.push({ a: c, b: o, order: 1 });
        changed = true;
      } else if (h === 1 && singleO.length === 1) {
        throw new Unsupported('formic acid / formate oxidation');
      }
      continue;
    }
    const oh = nb.filter((x) => x.order === 1 && g.atoms[x.atom].el === 'O' && hCount(g, x.atom) === 1);
    if (oh.length === 0) continue;
    if (oh.length > 1) throw new Unsupported('gem-diol');
    if (nb.some((x) => x.order !== 1)) continue;
    if (carbonNb === 0) throw new Unsupported('methanol oxidation');
    if (h >= 2) {
      // 第一級 → カルボン酸
      const o = addAtom(g, 'O');
      g.bonds.push({ a: c, b: o, order: 2 });
      changed = true;
    } else if (h === 1) {
      oh[0].bond.order = 2; // 第二級 → ケトン
      changed = true;
    }
  }
  return changed ? graphToSmilesList(RDKit, g) : [];
}

// エステルの加水分解: R-CO-O-R' → R-COOH + R'-OH
function hydrolyze(RDKit, g0) {
  const g = cloneGraph(g0);
  let changed = false;
  const n0 = g.atoms.length;
  for (let c = 0; c < n0; c++) {
    if (g.atoms[c].el !== 'C' || !isCarbonyl(g, c)) continue;
    for (const x of neighbors(g, c)) {
      if (x.order !== 1 || g.atoms[x.atom].el !== 'O') continue;
      const oe = x.atom;
      const alkyl = neighbors(g, oe).filter((y) => y.atom !== c);
      if (alkyl.length !== 1 || g.atoms[alkyl[0].atom].el !== 'C') continue;
      const r = alkyl[0].atom;
      if (isCarbonyl(g, r)) throw new Unsupported('anhydride');
      if (neighbors(g, r).some((y) => y.order !== 1)) throw new Unsupported('enol ester');
      g.bonds.splice(g.bonds.indexOf(x.bond), 1);
      const o = addAtom(g, 'O');
      g.bonds.push({ a: c, b: o, order: 1 });
      changed = true;
      break;
    }
  }
  return changed ? graphToSmilesList(RDKit, g) : [];
}

// オゾン分解: C=C → C=O + O=C
function ozonolyze(RDKit, g0) {
  const g = cloneGraph(g0);
  if (hasCC(g, 3)) throw new Unsupported('alkyne ozonolysis');
  const dbl = g.bonds.filter((b) => b.order === 2 && g.atoms[b.a].el === 'C' && g.atoms[b.b].el === 'C');
  if (dbl.length === 0) return [];
  for (const b of dbl) {
    const i = g.bonds.indexOf(b);
    g.bonds.splice(i, 1);
    for (const c of [b.a, b.b]) {
      const o = addAtom(g, 'O');
      g.bonds.push({ a: c, b: o, order: 2 });
    }
  }
  return graphToSmilesList(RDKit, g);
}

// ---------- 官能基判定（SMARTS） ----------

function matchesAny(RDKit, smiles, smartsList) {
  const mol = RDKit.get_mol(smiles);
  try {
    return smartsList.some((s) => {
      const q = RDKit.get_qmol(s);
      try {
        const m = JSON.parse(mol.get_substruct_match(q));
        return Array.isArray(m.atoms) && m.atoms.length > 0;
      } finally {
        q.delete();
      }
    });
  } finally {
    mol.delete();
  }
}

const ALDEHYDE = ['[CX3H1](=O)', '[CH2]=O'];

// kind: bool = 陽性/陰性, count = 個数, products = 生成物の構造
const CARDS = {
  silver_mirror: {
    name: '銀鏡反応',
    action: 'アンモニア性硝酸銀水溶液を加えて温める',
    kind: 'bool', yes: '銀が析出した', no: '変化なし',
    smarts: ALDEHYDE,
  },
  fehling: {
    name: 'フェーリング液',
    action: 'フェーリング液を加えて加熱する',
    kind: 'bool', yes: '赤色沈殿が生じた', no: '変化なし',
    smarts: ALDEHYDE,
  },
  iodoform: {
    name: 'ヨードホルム反応',
    action: 'ヨウ素と NaOH 水溶液を加えて温める',
    kind: 'bool', yes: '黄色沈殿が生じた', no: '変化なし',
    smarts: ['[CH3][CX3](=O)[#6]', '[CH3][CX3H1]=O', '[CH3][CX4H1]([OX2H1])[#6]', '[CH3][CH2][OX2H1]'],
  },
  sodium: {
    name: '金属Na',
    action: '金属ナトリウムを加える',
    kind: 'bool', yes: '気体が発生した', no: '変化なし',
    smarts: ['[OX2H1]'],
  },
  nahco3: {
    name: 'NaHCO₃',
    action: '炭酸水素ナトリウム水溶液を加える',
    kind: 'bool', yes: '気体が発生した', no: '変化なし',
    smarts: ['[CX3](=O)[OX2H1]'],
  },
  fecl3: {
    name: 'FeCl₃',
    action: '塩化鉄(III) 水溶液を加える',
    kind: 'bool', yes: '呈色した', no: '呈色しない',
    smarts: ['c[OX2H1]'],
  },
  bromine: {
    name: '臭素水',
    action: '臭素水を加える',
    kind: 'bool', yes: '臭素の色が消えた', no: '色は消えない',
    smarts: ['C=C', 'C#C'],
  },
  kmno4: {
    name: 'KMnO₄ 酸化',
    action: '硫酸酸性の KMnO₄ 水溶液で十分に酸化する',
    kind: 'products', none: '酸化されなかった',
    transform: oxidize,
  },
  hydrolysis: {
    name: '加水分解',
    action: '希硫酸を加えて加熱し、加水分解する',
    kind: 'products', none: '加水分解されなかった',
    transform: hydrolyze,
  },
  ozonolysis: {
    name: 'オゾン分解',
    action: 'オゾン分解する',
    kind: 'products', none: '反応しない',
    transform: ozonolyze,
  },
  chiral: {
    name: '不斉炭素',
    action: '不斉炭素原子の数を調べる',
    kind: 'count',
    compute: (RDKit, g) => chiralCount(g),
  },
  cis_trans: {
    name: 'シス-トランス異性',
    action: 'シス-トランス異性体が存在するか調べる',
    kind: 'bool', yes: '存在する', no: '存在しない',
    compute: (RDKit, g) => hasCisTrans(g),
  },
};

function evaluate(RDKit, card, smiles) {
  const def = CARDS[card];
  if (!def) throw new Error(`unknown card: ${card}`);
  if (def.smarts) return matchesAny(RDKit, smiles, def.smarts);
  const g = graphFromSmiles(RDKit, smiles);
  if (def.transform) return def.transform(RDKit, g);
  return def.compute(RDKit, g);
}

function normalizeResult(RDKit, card, result) {
  if (CARDS[card].kind === 'products') return result.map((s) => canonical(RDKit, s)).sort();
  return result;
}

function sameResult(a, b) {
  return JSON.stringify(a) === JSON.stringify(b);
}

module.exports = {
  CARDS, Unsupported, evaluate, normalizeResult, sameResult,
  canonical, graphFromSmiles, formula, hCount, chiralCount, hasCisTrans,
};
