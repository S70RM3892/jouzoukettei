'use strict';
// 手がかりカードの判定ロジック（Node の検証・ビルド専用）。
// 各カードは「候補1つ → 結果」の関数。結果が問題データの値と一致しない候補が除外される。
// 対象外の構造（例: C=C を持つ分子に KMnO4）では Unsupported を投げ、その問題は検証で落ちる。

class Unsupported extends Error {}

const VALENCE = { C: 4, O: 2, N: 3, S: 2, F: 1, Cl: 1, Br: 1, I: 1, H: 1, R: 1 }; // R = 高分子の繰り返し単位の端（*）
const MASS = { H: 1.0, C: 12, N: 14, O: 16, S: 32, Cl: 35.5 }; // 入試で与えられる原子量

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

// ペプチドの略記 "pep:Gly-Ala-Phe"（N末端 → C末端）を SMILES に展開する
const RESIDUES = {
  Gly: '', Ala: 'C', Ser: 'CO', Cys: 'CS', Phe: 'Cc1ccccc1', Tyr: 'Cc1ccc(O)cc1',
  Lys: 'CCCCN', Glu: 'CCC(=O)O', Asp: 'CC(=O)O', Val: 'C(C)C', Leu: 'CC(C)C',
};

function expand(s) {
  if (!s.startsWith('pep:')) return s;
  const seq = s.slice(4).split('-');
  return 'N' + seq.map((r, i) => {
    if (!(r in RESIDUES)) throw new Error(`unknown residue ${r}`);
    const side = RESIDUES[r] ? `(${RESIDUES[r]})` : '';
    return `C${side}C(=O)` + (i === seq.length - 1 ? 'O' : 'N');
  }).join('');
}

function graphFromSmiles(RDKit, smiles) {
  const mol = RDKit.get_mol(expand(smiles));
  if (!mol || !mol.is_valid()) throw new Error(`invalid SMILES: ${smiles}`);
  try {
    const g = parseMolblock(mol.get_molblock());
    const q = RDKit.get_qmol('a');
    try {
      const m = JSON.parse(mol.get_substruct_matches(q));
      (Array.isArray(m) ? m : []).forEach((x) => { g.atoms[x.atoms[0]].arom = true; });
    } finally {
      q.delete();
    }
    return g;
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
  g.atoms.forEach((a, i) => { if (a.el !== 'R') { add(a.el, 1); add('H', hCount(g, i)); } });
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
  const mol = RDKit.get_mol(expand(smilesOrMolblock));
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
  g.bonds.forEach((b) => {
    const o = g.atoms[b.a].arom && g.atoms[b.b].arom ? 'a' : b.order;
    adj[b.a].push([b.b, o]);
    adj[b.b].push([b.a, o]);
  });
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

// シス-トランス異性を生じる C=C の数
function stereoBondCount(g) {
  let n = 0;
  for (const bond of g.bonds) {
    if (bond.order !== 2) continue;
    const { a, b } = bond;
    if (g.atoms[a].el !== 'C' || g.atoms[b].el !== 'C' || g.atoms[a].arom) continue;
    if (bondInSmallRing(g, bond, 8)) continue;
    const { color, adj } = refine(g, [a, b]);
    const ok = [a, b].every((x) => {
      const other = adj[x].filter(([j]) => j !== a && j !== b);
      if (other.length !== 2) return false; // アレン等は除外
      return color[other[0][0]] !== color[other[1][0]];
    });
    if (ok) n++;
  }
  return n;
}

function hasCisTrans(g) {
  return stereoBondCount(g) > 0;
}

// ---------- 変換反応 ----------

function addAtom(g, el) {
  g.atoms.push({ el });
  return g.atoms.length - 1;
}

// ベンゼン環の中の結合は数えない
function isAromBond(g, b) {
  return g.atoms[b.a].arom && g.atoms[b.b].arom;
}

function hasCC(g, order) {
  return g.bonds.some((b) => b.order === order && !isAromBond(g, b)
    && g.atoms[b.a].el === 'C' && g.atoms[b.b].el === 'C');
}

function isCarbonyl(g, c) {
  return neighbors(g, c).some((n) => n.order === 2 && g.atoms[n.atom].el === 'O');
}

function isPhenolO(g, o) {
  return g.atoms[o].el === 'O' && hCount(g, o) === 1 && neighbors(g, o).some((n) => g.atoms[n.atom].arom);
}

// 生成物のうち、除外する原子を含む成分（CO2 や H2O になった断片）を取り除いて SMILES にする
function productsExcluding(RDKit, g, drop) {
  // components は原子をコピーするので、元の添字で判定するために印を付けておく
  g.atoms.forEach((a, i) => { a._i = i; });
  const keep = components(g).filter((p) => !p.atoms.some((a) => drop.has(a._i)));
  g.atoms.forEach((a) => { delete a._i; });
  return keep.map((p) => canonical(RDKit, toMolblock(p))).sort();
}

// 硫酸酸性 KMnO4
// - ベンゼン環の側鎖: 環に直結した炭素に H があれば、側鎖全体が -COOH になる
// - 第一級アルコール・アルデヒド → カルボン酸、第二級アルコール → ケトン、第三級は反応しない
function oxidize(RDKit, g0) {
  const g = cloneGraph(g0);
  if (hasCC(g, 2) || hasCC(g, 3)) throw new Unsupported('KMnO4 with C=C / C#C');
  if (g.atoms.some((a, i) => isPhenolO(g, i))) throw new Unsupported('KMnO4 with phenol');
  let changed = false;
  const drop = new Set();
  const n0 = g.atoms.length;

  for (let c = 0; c < n0; c++) {
    const a = g.atoms[c];
    if (a.el !== 'C' || a.arom) continue;
    const nb = neighbors(g, c);
    const ring = nb.filter((x) => g.atoms[x.atom].arom);
    if (ring.length !== 1) continue;
    if (isCarbonyl(g, c)) {
      const singleHet = nb.some((x) => x.order === 1 && g.atoms[x.atom].el !== 'C');
      if (!singleHet && hCount(g, c) === 0) throw new Unsupported('aryl ketone oxidation');
      continue; // アルデヒドは下で、カルボン酸・エステル・アミドはそのまま
    }
    if (nb.some((x) => x.order !== 1)) continue;
    if (hCount(g, c) === 0) continue; // tert-ブチル基などは酸化されない
    // 側鎖を切り離して -COOH にする
    for (const x of nb) {
      if (x.atom === ring[0].atom) continue;
      g.bonds.splice(g.bonds.indexOf(x.bond), 1);
      markFragment(g, x.atom, drop);
    }
    const o1 = addAtom(g, 'O');
    const o2 = addAtom(g, 'O');
    g.bonds.push({ a: c, b: o1, order: 2 }, { a: c, b: o2, order: 1 });
    changed = true;
  }

  for (let c = 0; c < n0; c++) {
    if (g.atoms[c].el !== 'C' || g.atoms[c].arom || drop.has(c)) continue;
    const nb = neighbors(g, c);
    const carbonNb = nb.filter((x) => g.atoms[x.atom].el === 'C').length;
    const h = hCount(g, c);
    if (isCarbonyl(g, c)) {
      // アルデヒド基 → カルボキシ基
      const singleHet = nb.filter((x) => x.order === 1 && g.atoms[x.atom].el !== 'C');
      if (h === 1 && singleHet.length === 0) {
        if (carbonNb === 0) throw new Unsupported('formaldehyde oxidation');
        const o = addAtom(g, 'O');
        g.bonds.push({ a: c, b: o, order: 1 });
        changed = true;
      } else if (h >= 1 && singleHet.length >= 1) {
        throw new Unsupported('formic acid / formate / formamide oxidation');
      } else if (h === 2) {
        throw new Unsupported('formaldehyde oxidation');
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
  return changed ? productsExcluding(RDKit, g, drop) : [];
}

// start から結合をたどれる原子すべてに印を付ける（切り離した側鎖）
function markFragment(g, start, set) {
  const stack = [start];
  set.add(start);
  while (stack.length) {
    const x = stack.pop();
    for (const n of neighbors(g, x)) if (!set.has(n.atom)) { set.add(n.atom); stack.push(n.atom); }
  }
}

// 加水分解で切れる結合: エステル R-CO-O-R' と アミド R-CO-NH-R'
function hydrolyzable(g) {
  const out = [];
  g.atoms.forEach((a, c) => {
    if (a.el !== 'C' || !isCarbonyl(g, c)) return;
    for (const x of neighbors(g, c)) {
      if (x.order !== 1) continue;
      const het = g.atoms[x.atom].el;
      if (het !== 'O' && het !== 'N') continue;
      const other = neighbors(g, x.atom).filter((y) => y.atom !== c);
      if (het === 'O') {
        if (other.length !== 1 || g.atoms[other[0].atom].el !== 'C') continue; // カルボン酸
        const r = other[0].atom;
        if (isCarbonyl(g, r)) throw new Unsupported('anhydride');
        const dbl = g.atoms[r].arom ? [] : neighbors(g, r).filter((y) => y.order !== 1);
        if (dbl.length > 1 || (dbl.length === 1 && (dbl[0].order !== 2 || g.atoms[dbl[0].atom].el !== 'C'))) throw new Unsupported('unusual enol ester');
        // ビニルエステル: 生じるビニルアルコール型はすぐにカルボニル化合物へ変わる
        if (dbl.length === 1) { out.push({ c, bond: x.bond, enol: { r, o: x.atom, s: dbl[0].atom } }); continue; }
      } else {
        if (other.some((y) => isCarbonyl(g, y.atom))) throw new Unsupported('imide');
        if (other.length === 0) continue; // 第一級アミド R-CONH2 は扱わない
      }
      out.push({ c, bond: x.bond });
    }
  });
  return out;
}

function breakBonds(RDKit, g0, list) {
  const g = cloneGraph(g0);
  const bonds = list.map((h) => g.bonds.find((b) => b.a === h.bond.a && b.b === h.bond.b));
  for (let i = 0; i < list.length; i++) {
    g.bonds.splice(g.bonds.indexOf(bonds[i]), 1);
    const o = addAtom(g, 'O');
    g.bonds.push({ a: list[i].c, b: o, order: 1 });
    const e = list[i].enol;
    if (e) {
      const find = (p, q) => g.bonds.find((b) => (b.a === p && b.b === q) || (b.a === q && b.b === p));
      find(e.r, e.o).order = 2;
      find(e.r, e.s).order = 1;
    }
  }
  return graphToSmilesList(RDKit, g);
}

// 完全な加水分解
function hydrolyze(RDKit, g) {
  const hs = hydrolyzable(g);
  return hs.length ? breakBonds(RDKit, g, hs) : [];
}

// 部分的な加水分解で得られうる化合物すべて（全部は切れない切り方）
function partialProducts(RDKit, g) {
  const hs = hydrolyzable(g);
  if (hs.length < 2) throw new Unsupported('partial hydrolysis needs 2+ hydrolyzable bonds');
  if (hs.length > 6) throw new Unsupported('too many hydrolyzable bonds');
  const set = new Set();
  for (let mask = 1; mask < (1 << hs.length) - 1; mask++) {
    const pick = hs.filter((_, i) => mask & (1 << i));
    breakBonds(RDKit, g, pick).forEach((s) => set.add(s));
  }
  return [...set].sort();
}

// 分子内脱水（濃硫酸・加熱）で生じるアルケン（構造異性体のみ、重複なし）
function dehydrate(RDKit, g0) {
  const alcohols = [];
  g0.atoms.forEach((a, c) => {
    if (a.el !== 'C' || a.arom || isCarbonyl(g0, c)) return;
    const nb = neighbors(g0, c);
    if (nb.some((x) => x.order !== 1)) return;
    nb.forEach((x) => {
      if (g0.atoms[x.atom].el === 'O' && hCount(g0, x.atom) === 1) alcohols.push({ c, o: x.atom, bond: x.bond });
    });
  });
  if (alcohols.length === 0) return [];
  if (alcohols.length > 1) throw new Unsupported('dehydration of polyol');
  const { c, o } = alcohols[0];
  const set = new Set();
  for (const x of neighbors(g0, c)) {
    const b = x.atom;
    if (g0.atoms[b].el !== 'C' || g0.atoms[b].arom || isCarbonyl(g0, b)) continue;
    if (neighbors(g0, b).some((y) => y.order !== 1)) continue;
    if (hCount(g0, b) === 0) continue;
    const g = cloneGraph(g0);
    const ob = g.bonds.find((y) => (y.a === c && y.b === o) || (y.a === o && y.b === c));
    g.bonds.splice(g.bonds.indexOf(ob), 1);
    const cb = g.bonds.find((y) => (y.a === c && y.b === b) || (y.a === b && y.b === c));
    cb.order = 2;
    productsExcluding(RDKit, g, new Set([o])).forEach((s) => set.add(s));
  }
  if (set.size === 0) throw new Unsupported('alcohol without beta hydrogen');
  return [...set].sort();
}

// 脱水で生じるアルケンの種類（シス-トランス異性体も別に数える）
function dehydrationCount(RDKit, g) {
  return dehydrate(RDKit, g).reduce((n, s) => n + 2 ** stereoBondCount(graphFromSmiles(RDKit, s)), 0);
}

// 脱水で得たアルケンをオゾン分解して得られる化合物（混合物全体）
function dehydrateOzonolyze(RDKit, g) {
  const set = new Set();
  for (const s of dehydrate(RDKit, g)) ozonolyze(RDKit, graphFromSmiles(RDKit, s)).forEach((p) => set.add(p));
  return [...set].sort();
}

// オゾン分解: C=C → C=O + O=C
function ozonolyze(RDKit, g0) {
  const g = cloneGraph(g0);
  if (hasCC(g, 3)) throw new Unsupported('alkyne ozonolysis');
  const dbl = g.bonds.filter((b) => b.order === 2 && !isAromBond(g, b)
    && g.atoms[b.a].el === 'C' && g.atoms[b.b].el === 'C');
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

// 炭素上の H を1つ Cl に置き換えた化合物の種類（構造異性体）。onlyRing ならベンゼン環の H だけ
function chloroCount(RDKit, g, onlyRing) {
  const set = new Set();
  g.atoms.forEach((a, i) => {
    if (a.el !== 'C' || (onlyRing && !a.arom) || hCount(g, i) === 0) return;
    const h = cloneGraph(g);
    const cl = addAtom(h, 'Cl');
    h.bonds.push({ a: i, b: cl, order: 1 });
    set.add(canonical(RDKit, toMolblock(h)));
  });
  if (onlyRing && set.size === 0) throw new Unsupported('no aromatic C-H');
  return set.size;
}

// 加熱で分子内脱水して環状の酸無水物（5員環・6員環）になるか
function formsAnhydride(g) {
  const acids = g.atoms.map((a, i) => i).filter((c) => g.atoms[c].el === 'C' && isCarbonyl(g, c)
    && neighbors(g, c).some((x) => x.order === 1 && g.atoms[x.atom].el === 'O' && hCount(g, x.atom) === 1));
  for (let i = 0; i < acids.length; i++) {
    for (let j = i + 1; j < acids.length; j++) {
      const path = shortestPath(g, acids[i], acids[j]);
      const size = path.length + 1; // 経路上の原子 + 橋かけの O
      if (size !== 5 && size !== 6) continue;
      for (let k = 0; k + 1 < path.length; k++) {
        const b = g.bonds.find((y) => (y.a === path[k] && y.b === path[k + 1]) || (y.b === path[k] && y.a === path[k + 1]));
        if (b.order === 2 && !isAromBond(g, b)) throw new Unsupported('anhydride across C=C needs cis/trans');
      }
      return true;
    }
  }
  return false;
}

function shortestPath(g, s, t) {
  const prev = new Map([[s, -1]]);
  const q = [s];
  while (q.length) {
    const x = q.shift();
    if (x === t) break;
    for (const n of neighbors(g, x)) if (!prev.has(n.atom)) { prev.set(n.atom, x); q.push(n.atom); }
  }
  const path = [];
  for (let x = t; x !== -1; x = prev.get(x)) path.unshift(x);
  return path;
}

// 立体異性体を含めた数（不斉炭素2個以上はメソ体の判定が要るので扱わない）
function stereoCount(g) {
  const n = chiralCount(g);
  if (n > 1) throw new Unsupported('2+ chiral centers (meso check needed)');
  return 2 ** (n + stereoBondCount(g));
}

// ---------- 官能基判定（SMARTS） ----------

function matchesAny(RDKit, smiles, smartsList) {
  const mol = RDKit.get_mol(expand(smiles));
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
  // ---- Phase 2: 芳香族・脱水・部分加水分解 ----
  naoh: {
    name: 'NaOH 水溶液',
    action: 'NaOH 水溶液に加えて振り混ぜる',
    kind: 'bool', yes: '溶けた', no: '溶けない',
    smarts: ['[CX3](=O)[OX2H1]', 'c[OX2H1]'],
  },
  hcl: {
    name: '希塩酸',
    action: '希塩酸に加えて振り混ぜる',
    kind: 'bool', yes: '溶けた', no: '溶けない',
    smarts: ['[NX3;!$(N-C=O);!$(N=*);!$(N#*)]'],
  },
  dehydration: {
    name: '脱水',
    action: '濃硫酸を加えて加熱し、分子内で脱水する',
    kind: 'products', none: '脱水されなかった',
    transform: dehydrate,
  },
  dehydration_count: {
    name: '脱水生成物の数',
    action: '分子内脱水で生じるアルケンの種類を数える（シス-トランス異性体は別に数える）',
    kind: 'count',
    compute: dehydrationCount,
  },
  dehydration_ozonolysis: {
    name: '脱水 → オゾン分解',
    action: '分子内脱水で得たアルケンをすべてオゾン分解する',
    kind: 'products', none: '反応しない',
    transform: dehydrateOzonolyze,
  },
  ring_cl: {
    name: '環の塩素置換体',
    action: 'ベンゼン環の H 原子1個を Cl 原子に置き換えた化合物の種類を数える',
    kind: 'count',
    compute: (RDKit, g) => chloroCount(RDKit, g, true),
  },
  cl_sub: {
    name: '塩素置換体',
    action: '炭素原子に結合した H 原子1個を Cl 原子に置き換えた化合物の種類を数える（構造異性体のみ）',
    kind: 'count',
    compute: (RDKit, g) => chloroCount(RDKit, g, false),
  },
  anhydride: {
    name: '加熱脱水',
    action: '加熱して、分子内で脱水した環状の酸無水物になるか調べる',
    kind: 'bool', yes: '酸無水物になった', no: 'ならなかった',
    compute: (RDKit, g) => formsAnhydride(g),
  },
  partial_hydrolysis: {
    name: '部分加水分解',
    action: '穏やかに加水分解し、途中の段階の生成物を調べる',
    kind: 'contains', yes: 'が得られた', no: 'は得られない',
    transform: partialProducts,
  },
  // ---- Phase 3: アミノ酸・ペプチド ----
  ninhydrin: {
    name: 'ニンヒドリン反応',
    action: 'ニンヒドリン水溶液を加えて温める',
    kind: 'bool', yes: '赤紫色になった', no: '変化なし',
    smarts: ['[NX3H2][CX4]'],
  },
  alpha_amino: {
    name: 'α-アミノ酸か',
    action: '同じ炭素原子にアミノ基とカルボキシ基が結合しているか調べる',
    kind: 'bool', yes: 'α-アミノ酸である', no: 'α-アミノ酸ではない',
    smarts: ['[NX3;!$(NC=O)][CX4][CX3](=O)[OX2H1]'],
  },
  xanthoprotein: {
    name: 'キサントプロテイン反応',
    action: '濃硝酸を加えて加熱し、冷やしてアンモニア水を加える',
    kind: 'bool', yes: '黄色 → 橙黄色になった', no: '変化なし',
    smarts: ['c'],
  },
  sulfur: {
    name: '硫黄の検出',
    action: 'NaOH を加えて加熱し、酢酸鉛(II) 水溶液を加える',
    kind: 'bool', yes: '黒色沈殿が生じた', no: '変化なし',
    smarts: ['[SX2H1]'],
  },
  biuret: {
    name: 'ビウレット反応',
    action: 'NaOH 水溶液と少量の CuSO₄ 水溶液を加える',
    kind: 'bool', yes: '赤紫色になった', no: '青色のまま',
    compute: (RDKit, g) => hydrolyzable(g).filter((h) => g.atoms[h.bond.a === h.c ? h.bond.b : h.bond.a].el === 'N').length >= 2,
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
  const kind = CARDS[card].kind;
  if (kind === 'products') return result.map((s) => canonical(RDKit, s)).sort();
  if (kind === 'contains') return canonical(RDKit, result);
  return result;
}

function sameResult(a, b) {
  return JSON.stringify(a) === JSON.stringify(b);
}

// 候補にカードを当てた値 value が、観察された結果 expected と矛盾しないか
function consistent(card, value, expected) {
  if (CARDS[card].kind === 'contains') return Array.isArray(value) && value.includes(expected);
  return sameResult(value, expected);
}

// 分子式の式量（入試の原子量で）
function formulaMass(f) {
  let m = 0;
  for (const [, el, n] of f.matchAll(/([A-Z][a-z]?)(\d*)/g)) m += MASS[el] * (n ? +n : 1);
  return Math.round(m * 10) / 10;
}

module.exports = {
  formulaMass,
  CARDS, Unsupported, evaluate, normalizeResult, sameResult, consistent, expand,
  canonical, graphFromSmiles, formula, hCount, chiralCount, hasCisTrans, stereoCount, hydrolyze,
};
