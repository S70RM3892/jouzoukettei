'use strict';
// 立体を区別する反応と判定。RDKit に座標つきの構造を渡して、立体配置を読ませる。
//  - アンチ付加（臭素）とシン付加（白金触媒の水素）: 京大2024 III(b)
//  - フィッシャー投影からの構造、硝酸酸化（両端を COOH）、炭素を1つ減らす反応、メソ体の判定: 京大2019 IV・2020 IV
// 立体を含む SMILES（@, @@, /, \）をそのまま受け取り、立体つきの正規 SMILES を返す。

function isoCanonical(RDKit, s) {
  const m = RDKit.get_mol(s);
  if (!m || !m.is_valid()) throw new Error(`invalid SMILES: ${s}`);
  try { return m.get_smiles(); } finally { m.delete(); }
}

// 鏡像: @ と @@ を入れ替える（/ \ は鏡に映しても変わらない）
function mirror(RDKit, s) {
  return isoCanonical(RDKit, s.replace(/@@|@/g, (x) => (x === '@@' ? '@' : '@@')));
}
function hasStereocenter(s) { return /@/.test(s); }
// 光学活性か（鏡像と重ならない）。メソ体や不斉炭素のない化合物は false
function opticallyActive(RDKit, s) {
  const c = isoCanonical(RDKit, s);
  return hasStereocenter(c) && mirror(RDKit, c) !== c;
}

// ---------- フィッシャー投影 ----------
// top: 'CHO' | 'COOH' | 'CH2OH'、centers: 上から順に OH が 'R'（右）か 'L'（左）、bottom: 'CH2OH' | 'COOH' | 'CHO'
// 上から下へ SMILES を書くと、各中心の隣は（上, H, OH, 下）の順。OH が右なら @、左なら @@
// （D-グリセルアルデヒド = (R) になることを RDKit の CIP で確かめてある）
const TOP = { CHO: 'O=C', COOH: 'OC(=O)', CH2OH: 'OC' };
const BOTTOM = { CH2OH: 'CO', COOH: 'C(=O)O', CHO: 'C=O' };
function fischer(RDKit, top, centers, bottom) {
  const s = TOP[top] + centers.map((c) => (c === 'R' ? '[C@H](O)' : '[C@@H](O)')).join('') + BOTTOM[bottom];
  return isoCanonical(RDKit, s);
}

// ---------- 座標つき mol ファイルの読み書き ----------
function readMol(mb) {
  const lines = mb.split('\n');
  const na = parseInt(lines[3].slice(0, 3), 10);
  const nb = parseInt(lines[3].slice(3, 6), 10);
  const atoms = [];
  for (let i = 0; i < na; i++) {
    const l = lines[4 + i];
    atoms.push({ x: +l.slice(0, 10), y: +l.slice(10, 20), z: +l.slice(20, 30), el: l.slice(31, 34).trim(), chg: 0 });
  }
  const bonds = [];
  for (let i = 0; i < nb; i++) {
    const l = lines[4 + na + i];
    bonds.push({ a: +l.slice(0, 3) - 1, b: +l.slice(3, 6) - 1, order: +l.slice(6, 9), stereo: +(l.slice(9, 12) || 0) });
  }
  for (const l of lines) {
    if (!l.startsWith('M  CHG')) continue;
    const n = +l.slice(6, 9);
    for (let k = 0; k < n; k++) atoms[+l.slice(10 + 8 * k, 13 + 8 * k) - 1].chg = +l.slice(14 + 8 * k, 17 + 8 * k);
  }
  return { atoms, bonds };
}
function writeMol(m, dim = '2D') {
  const pad = (s, n) => String(s).padStart(n);
  // 2 行目の 21〜22 桁目が次元（2D/3D）。ここがずれると RDKit がくさびや 3D 座標を読まない
  const out = ['', `  jouzou  0101262000${dim}`, ''];
  out.push(`${pad(m.atoms.length, 3)}${pad(m.bonds.length, 3)}  0  0  0  0  0  0  0  0999 V2000`);
  for (const a of m.atoms) out.push(`${pad(a.x.toFixed(4), 10)}${pad(a.y.toFixed(4), 10)}${pad(a.z.toFixed(4), 10)} ${a.el.padEnd(3)} 0  0  0  0  0  0  0  0  0  0  0  0`);
  for (const b of m.bonds) out.push(`${pad(b.a + 1, 3)}${pad(b.b + 1, 3)}${pad(b.order, 3)}${pad(b.stereo || 0, 3)}`);
  const ch = m.atoms.map((a, i) => [i, a.chg]).filter(([, c]) => c);
  if (ch.length) out.push(`M  CHG${pad(ch.length, 3)}${ch.map(([i, c]) => ` ${pad(i + 1, 3)} ${pad(c, 3)}`).join('')}`);
  out.push('M  END');
  return out.join('\n');
}
function molOf(RDKit, s) {
  const m = RDKit.get_mol(s);
  try { return readMol(m.get_molblock()); } finally { m.delete(); }
}
function toSmiles(RDKit, m, dim) {
  const mol = RDKit.get_mol(writeMol(m, dim));
  if (!mol || !mol.is_valid()) throw new Error('invalid edited structure');
  try { return mol.get_smiles().split('.').map((x) => isoCanonical(RDKit, x)).sort(); } finally { mol.delete(); }
}
const nbrs = (m, i) => m.bonds.filter((b) => b.a === i || b.b === i).map((b) => ({ atom: b.a === i ? b.b : b.a, order: b.order, bond: b }));
const VAL = { C: 4, O: 2, N: 3, Br: 1, Cl: 1 };
const hOf = (m, i) => (VAL[m.atoms[i].el] || 0) - nbrs(m, i).reduce((s, x) => s + x.order, 0);
function addAtomNear(m, parent, el, dx = 0.8, dy = 0.8, dz = 0) {
  const p = m.atoms[parent];
  m.atoms.push({ x: p.x + dx, y: p.y + dy, z: p.z + dz, el, chg: 0 });
  return m.atoms.length - 1;
}
function dropAtoms(m, dead) {
  const map = new Map();
  const atoms = [];
  m.atoms.forEach((a, i) => { if (!dead.has(i)) { map.set(i, atoms.length); atoms.push(a); } });
  return { atoms, bonds: m.bonds.filter((b) => !dead.has(b.a) && !dead.has(b.b)).map((b) => ({ ...b, a: map.get(b.a), b: map.get(b.b) })) };
}

// ---------- 付加の立体（京大2024 III(b)） ----------
// C=C を含む平面の上下から付加させる。anti: 一方は上・他方は下、syn: 両方とも同じ側。
// 平面の表と裏のどちらから近づくかは同じ確率なので、両方の生成物（鏡像の組）を返す
function addition(RDKit, s, x, mode) {
  const base = molOf(RDKit, s);
  const dbl = base.bonds.filter((b) => b.order === 2 && base.atoms[b.a].el === 'C' && base.atoms[b.b].el === 'C');
  if (dbl.length !== 1) throw new Error('one C=C is required');
  const out = new Set();
  for (const face of [1, -1]) {
    const m = { atoms: base.atoms.map((a) => ({ ...a, z: 0 })), bonds: base.bonds.map((b) => ({ ...b, stereo: 0 })) };
    const b = m.bonds[base.bonds.indexOf(dbl[0])];
    b.order = 1;
    addAtomNear(m, b.a, x, 0, 0, 1.4 * face);
    addAtomNear(m, b.b, x, 0, 0, (mode === 'anti' ? -1.4 : 1.4) * face);
    m.bonds.push({ a: b.a, b: m.atoms.length - 2, order: 1, stereo: 0 });
    m.bonds.push({ a: b.b, b: m.atoms.length - 1, order: 1, stereo: 0 });
    toSmiles(RDKit, m, '3D').forEach((p) => out.add(p));
  }
  return [...out].sort();
}
const antiAddition = (RDKit, s, x = 'Br') => addition(RDKit, s, x, 'anti');
const synAddition = (RDKit, s, x = 'H') => addition(RDKit, s, x, 'syn');

// ---------- 糖の反応（京大2019 IV） ----------
// 反応1: 硝酸で酸化。アルデヒド基と末端の CH₂OH がどちらも -COOH になる（不斉炭素の配置は変わらない）
function nitricOxidation(RDKit, s) {
  const m = molOf(RDKit, s);
  const n0 = m.atoms.length;
  for (let c = 0; c < n0; c++) {
    if (m.atoms[c].el !== 'C') continue;
    const nb = nbrs(m, c);
    const dO = nb.find((x) => x.order === 2 && m.atoms[x.atom].el === 'O');
    const oh = nb.filter((x) => x.order === 1 && m.atoms[x.atom].el === 'O' && hOf(m, x.atom) === 1);
    const cs = nb.filter((x) => m.atoms[x.atom].el === 'C');
    if (dO && hOf(m, c) === 1) { m.bonds.push({ a: c, b: addAtomNear(m, c, 'O'), order: 1, stereo: 0 }); continue; }
    if (!dO && oh.length === 1 && cs.length === 1 && hOf(m, c) === 2) m.bonds.push({ a: c, b: addAtomNear(m, c, 'O', -0.8, 0.8), order: 2, stereo: 0 });
  }
  return toSmiles(RDKit, m, '2D');
}
// 反応2: アルデヒド基の炭素が外れ、隣の炭素がアルデヒド基になる（炭素が1つ減る）
function degrade(RDKit, s) {
  const m = molOf(RDKit, s);
  const c1 = m.atoms.findIndex((a, c) => a.el === 'C' && hOf(m, c) === 1 && nbrs(m, c).some((x) => x.order === 2 && m.atoms[x.atom].el === 'O'));
  if (c1 < 0) throw new Error('no aldehyde');
  const o1 = nbrs(m, c1).find((x) => x.order === 2).atom;
  const c2 = nbrs(m, c1).find((x) => m.atoms[x.atom].el === 'C').atom;
  const oh2 = nbrs(m, c2).find((x) => m.atoms[x.atom].el === 'O' && x.order === 1);
  if (!oh2) throw new Error('C2 has no OH');
  oh2.bond.order = 2;
  m.bonds.forEach((b) => { if (b.a === c2 || b.b === c2) b.stereo = 0; });
  return toSmiles(RDKit, dropAtoms(m, new Set([c1, o1])), '2D');
}

// 立体異性体の数（鏡像異性体・メソ体・シス-トランスを区別して数える）。
// 不斉炭素ごとにくさび（上・下）を付けた構造をすべて作り、RDKit で同じものをまとめる。C=C のシス-トランスは 2 倍ずつ掛ける
function countStereoisomers(RDKit, s) {
  const chem = require('./chem');
  // 原子の番号をそろえるため、グラフと座標は同じ（立体を外した）SMILES から作る
  const flat = isoCanonical(RDKit, s.replace(/[@/\\]/g, ''));
  const g = chem.graphFromSmiles(RDKit, flat);
  const centers = chem.chiralCenters(g);
  const ez = chem.stereoBondCount(g);
  if (centers.length > 12) throw new Error('too many centers');
  const base = molOf(RDKit, flat);
  // くさびを付ける結合: 中心から出る単結合（環の外を優先）
  const pickBond = (c) => {
    const bs = base.bonds.filter((b) => (b.a === c || b.b === c) && b.order === 1);
    return bs.find((b) => !centers.includes(b.a === c ? b.b : b.a)) || bs[0];
  };
  const seen = new Set();
  for (let mask = 0; mask < 1 << centers.length; mask++) {
    const m = { atoms: base.atoms.map((a) => ({ ...a })), bonds: base.bonds.map((b) => ({ ...b, stereo: 0 })) };
    centers.forEach((c, k) => {
      const b = m.bonds[base.bonds.indexOf(pickBond(c))];
      if (b.a !== c) { const t = b.a; b.a = b.b; b.b = t; }
      b.stereo = mask & (1 << k) ? 1 : 6;
    });
    seen.add(toSmiles(RDKit, m, '2D').join('.'));
  }
  return seen.size * 2 ** ez;
}

module.exports = {
  countStereoisomers, isoCanonical, mirror, opticallyActive, fischer, antiAddition, synAddition, nitricOxidation, degrade };
