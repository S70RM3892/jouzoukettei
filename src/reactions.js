'use strict';
// 京大で「問題文で与えられた規則」として出た反応と、芳香族の置換反応（配向性）。
// どれも分子グラフ（立体なし）の変換。立体を扱う反応は stereo.js。
// 出典: 2001（ジアゾカップリング）、2010（第一級 OH の選択的アセチル化）、2016（ジアリールエーテルの水素化分解）、
//       2018（配向性と合成経路）、2019（イミドの穏やかな加水分解）、2020・2026（アセタール・アセトニド）、
//       2025（小員環の水素化による開環）、2026 IV（メチル化分析）
const chem = require('./chem');

const { neighbors, hCount, cloneGraph, Unsupported } = chem;

function addAtom(g, el, chg) {
  g.atoms.push(chg ? { el, chg } : { el });
  return g.atoms.length - 1;
}
function addBond(g, a, b, order = 1) { g.bonds.push({ a, b, order }); }
function removeBond(g, a, b) {
  const i = g.bonds.findIndex((x) => (x.a === a && x.b === b) || (x.a === b && x.b === a));
  if (i >= 0) g.bonds.splice(i, 1);
}
function bondBetween(g, a, b) { return g.bonds.find((x) => (x.a === a && x.b === b) || (x.a === b && x.b === a)); }
function removeAtoms(g, dead) {
  const map = new Map();
  const atoms = [];
  g.atoms.forEach((a, i) => { if (!dead.has(i)) { map.set(i, atoms.length); atoms.push(a); } });
  return { atoms, bonds: g.bonds.filter((b) => !dead.has(b.a) && !dead.has(b.b)).map((b) => ({ a: map.get(b.a), b: map.get(b.b), order: b.order })) };
}
const el = (g, i) => g.atoms[i].el;
const isCarbonylC = (g, c) => el(g, c) === 'C' && neighbors(g, c).some((x) => x.order === 2 && el(g, x.atom) === 'O');
const smilesOf = (RDKit, g) => chem.canonical(RDKit, chem.toMolblock(g));
function productsOf(RDKit, g, dropSmiles = []) {
  const out = chem.graphToSmilesList(RDKit, g);
  const drop = new Set(dropSmiles.map((s) => chem.canonical(RDKit, s)));
  return out.filter((s) => !drop.has(s));
}
const uniqSorted = (xs) => [...new Set(xs)].sort();

// ---------- 芳香族の置換反応と配向性（2018・2009・2016・2004・2010） ----------
// 置換基の強さ: 正 = オルト・パラ配向（大きいほど強く活性化）、負 = メタ配向
function directorStrength(g, ring, r, x) {
  const e = el(g, x);
  if (e === 'N') {
    if (g.atoms[x].chg === 1) return -1; // ニトロ基
    if (neighbors(g, x).some((y) => y.atom !== r && isCarbonylC(g, y.atom))) return 2; // -NHCOCH₃（アセトアミド）
    return 3; // -NH₂, -NHR
  }
  if (e === 'O') {
    if (hCount(g, x) === 1 || g.atoms[x].chg === -1) return 3; // -OH, -O⁻
    if (neighbors(g, x).some((y) => y.atom !== r && isCarbonylC(g, y.atom))) return 1.5; // -OCOR
    return 2; // -OR
  }
  if (e === 'Cl' || e === 'Br' || e === 'I' || e === 'F') return 0.5; // ハロゲン: 不活性化だがオルト・パラ配向
  if (e === 'S') return -1; // -SO₃H
  if (e === 'C') {
    if (ring.has(x)) return 1;
    if (g.atoms[x].arom) return 1; // フェニル基
    const nb = neighbors(g, x);
    if (nb.some((y) => y.order === 2 && el(g, y.atom) === 'O')) return -1; // -COOH, -CHO, -COR, -COOR, -CONH₂
    if (nb.some((y) => y.order === 3 && el(g, y.atom) === 'N')) return -1; // -CN
    return 1; // アルキル基・ビニル基
  }
  throw new Unsupported(`substituent ${e}`);
}

function benzeneRing(g) {
  // 6 個の芳香族炭素からなる環が1つだけの分子に限る（ナフタレン・ピリジン・複数のベンゼン環は対象外）
  const arom = g.atoms.map((a, i) => (a.arom ? i : -1)).filter((i) => i >= 0);
  if (!arom.length) throw new Unsupported('no benzene ring');
  if (arom.length !== 6 || arom.some((i) => el(g, i) !== 'C')) throw new Unsupported('only one benzene ring is supported');
  return new Set(arom);
}

function ringDistance(g, ring, a, b) {
  const dist = new Map([[a, 0]]);
  const q = [a];
  while (q.length) {
    const x = q.shift();
    for (const n of neighbors(g, x)) {
      if (!ring.has(n.atom) || dist.has(n.atom)) continue;
      dist.set(n.atom, dist.get(x) + 1);
      q.push(n.atom);
    }
  }
  return dist.get(b);
}

// 置換が起こる位置（環の炭素の添字）を返す。exhaustive = 臭素水のように強い活性化基のオルト・パラ全部
function easPositions(g, { exhaustive = false } = {}) {
  const ring = benzeneRing(g);
  const subs = [];
  for (const r of ring) {
    for (const n of neighbors(g, r)) if (!ring.has(n.atom)) subs.push({ r, s: directorStrength(g, ring, r, n.atom) });
  }
  const free = [...ring].filter((r) => hCount(g, r) > 0);
  if (!free.length) return [];
  const op = (p, r) => { const d = ringDistance(g, ring, p, r); return d === 1 || d === 3; };
  if (exhaustive) {
    const top = subs.filter((x) => x.s >= 3);
    if (!top.length) return null; // 臭素水では反応しない
    return [free.filter((p) => top.some((x) => op(p, x.r)))];
  }
  let targets;
  const act = subs.filter((x) => x.s > 0);
  if (act.length) {
    const max = Math.max(...act.map((x) => x.s));
    const top = act.filter((x) => x.s === max);
    targets = free.filter((p) => top.some((x) => op(p, x.r)));
    // 立体障害: 2つの置換基に挟まれた位置は、ほかに候補があれば使わない
    const between = (p) => subs.filter((x) => ringDistance(g, ring, p, x.r) === 1).length >= 2;
    if (targets.some((p) => !between(p))) targets = targets.filter((p) => !between(p));
  } else if (subs.length) {
    targets = free.filter((p) => subs.every((x) => ringDistance(g, ring, p, x.r) === 2));
  } else targets = free;
  return targets.map((p) => [p]);
}

// 置換の生成物（オルトとパラが両方できるなら両方＝混合物）
function eas(group) {
  return (RDKit, g0, opts = {}) => {
    const sets = easPositions(g0, opts);
    if (sets === null) return [];
    const out = [];
    for (const ps of sets) {
      const g = cloneGraph(g0);
      for (const p of ps) group(g, p);
      out.push(smilesOf(RDKit, g));
    }
    return uniqSorted(out);
  };
}
const attachNitro = (g, p) => { const n = addAtom(g, 'N', 1); addBond(g, p, n); addBond(g, n, addAtom(g, 'O'), 2); addBond(g, n, addAtom(g, 'O', -1)); };
const attachHalogen = (x) => (g, p) => addBond(g, p, addAtom(g, x));
const attachSulfo = (g, p) => { const s = addAtom(g, 'S'); addBond(g, p, s); addBond(g, s, addAtom(g, 'O'), 2); addBond(g, s, addAtom(g, 'O'), 2); addBond(g, s, addAtom(g, 'O')); };
const nitration = eas(attachNitro);
const bromination = eas(attachHalogen('Br'));
const chlorination = eas(attachHalogen('Cl'));
const sulfonation = eas(attachSulfo);
function bromineWater(RDKit, g0) { return eas(attachHalogen('Br'))(RDKit, g0, { exhaustive: true }); }

// ニトロ基の還元（スズと塩酸、または H₂/Pd）→ アミノ基
function nitroReduction(RDKit, g0) {
  const g = cloneGraph(g0);
  const dead = new Set();
  g.atoms.forEach((a, i) => {
    if (a.el !== 'N' || a.chg !== 1) return;
    const os = neighbors(g, i).filter((x) => el(g, x.atom) === 'O');
    if (os.length !== 2) return;
    os.forEach((x) => dead.add(x.atom));
    delete a.chg;
  });
  if (!dead.size) return [];
  return chem.graphToSmilesList(RDKit, removeAtoms(g, dead));
}

// ---------- アセチル化（2015・2018・2010） ----------
function acylate(g, o) {
  const c = addAtom(g, 'C');
  addBond(g, o, c);
  addBond(g, c, addAtom(g, 'O'), 2);
  addBond(g, c, addAtom(g, 'C'));
}
function acetylSites(g) {
  const out = [];
  g.atoms.forEach((a, i) => {
    if (a.el === 'O' && hCount(g, i) === 1) {
      const c = neighbors(g, i)[0].atom;
      if (!isCarbonylC(g, c)) out.push({ i, c, kind: g.atoms[c].arom ? 'phenol' : 'alcohol' });
    }
    if (a.el === 'N' && !a.chg && hCount(g, i) >= 1 && !neighbors(g, i).some((x) => isCarbonylC(g, x.atom))) out.push({ i, kind: 'amine' });
  });
  return out;
}
// 無水酢酸で OH と NH₂ をすべてアセチル化
function acetylation(RDKit, g0) {
  const g = cloneGraph(g0);
  const sites = acetylSites(g);
  if (!sites.length) return [];
  sites.forEach((s) => acylate(g, s.i));
  return chem.graphToSmilesList(RDKit, g);
}
// 京大2010: 同じ物質量の無水酢酸では第一級アルコールの OH だけがアセチル化される
function acetylationPrimary(RDKit, g0) {
  const g = cloneGraph(g0);
  const sites = acetylSites(g).filter((s) => s.kind === 'alcohol' && hCount(g, s.c) >= 2);
  if (!sites.length) return [];
  if (sites.length > 1) throw new Unsupported('several primary OH');
  acylate(g, sites[0].i);
  return chem.graphToSmilesList(RDKit, g);
}

// ---------- ジアゾ化（2001・2017・2018） ----------
function arylAmines(g) {
  return g.atoms.map((a, i) => i).filter((i) => el(g, i) === 'N' && !g.atoms[i].chg && hCount(g, i) === 2 && neighbors(g, i).length === 1 && g.atoms[neighbors(g, i)[0].atom].arom);
}
// ジアゾ化して H₃PO₂ で還元: Ar–NH₂ → Ar–H
function deamination(RDKit, g0) {
  const ns = arylAmines(g0);
  if (!ns.length) return [];
  return chem.graphToSmilesList(RDKit, removeAtoms(g0, new Set(ns)));
}
// ジアゾ化して水と温める: Ar–NH₂ → Ar–OH
function diazoHydrolysis(RDKit, g0) {
  const ns = arylAmines(g0);
  if (!ns.length) return [];
  const g = cloneGraph(g0);
  ns.forEach((n) => { g.atoms[n] = { el: 'O' }; });
  return chem.graphToSmilesList(RDKit, g);
}
// カップリングの位置: フェノール（またはアニリン）の OH のパラ位、ふさがっていればオルト位
function couplingPosition(g) {
  const ring = benzeneRing(g);
  let best = null;
  for (const r of ring) {
    const oh = neighbors(g, r).find((x) => el(g, x.atom) === 'O' && hCount(g, x.atom) === 1);
    if (!oh) continue;
    const free = [...ring].filter((p) => hCount(g, p) > 0);
    const para = free.filter((p) => ringDistance(g, ring, p, r) === 3);
    const ortho = free.filter((p) => ringDistance(g, ring, p, r) === 1);
    best = para.length ? para : ortho;
    break;
  }
  if (!best || !best.length) throw new Unsupported('no coupling position');
  return best;
}
// 芳香族アミンをジアゾ化し、フェノール類とカップリングさせる
function azoCoupling(RDKit, amine, partner = 'Oc1ccccc1') {
  const ga = chem.graphFromSmiles(RDKit, amine);
  const gp = chem.graphFromSmiles(RDKit, partner);
  const ns = arylAmines(ga);
  if (ns.length !== 1) throw new Unsupported('need one aromatic NH2');
  const out = [];
  for (const p of couplingPosition(gp)) {
    const g = cloneGraph(ga);
    const off = g.atoms.length;
    gp.atoms.forEach((a) => g.atoms.push({ ...a }));
    gp.bonds.forEach((b) => g.bonds.push({ a: b.a + off, b: b.b + off, order: b.order }));
    const n2 = addAtom(g, 'N');
    addBond(g, ns[0], n2, 2);
    addBond(g, n2, p + off);
    out.push(smilesOf(RDKit, g));
  }
  return uniqSorted(out);
}

// ---------- 京大2019: イミドの穏やかな加水分解 ----------
// C(=O)–N–C(=O) の C–N 結合が1つだけ切れてアミドとカルボン酸になる。どちら側が切れるかは決まらないので両方の生成物（混合物）
function imideHydrolysis(RDKit, g0) {
  const out = [];
  let found = false;
  g0.atoms.forEach((a, n) => {
    if (a.el !== 'N') return;
    const acyl = neighbors(g0, n).filter((x) => x.order === 1 && isCarbonylC(g0, x.atom));
    if (acyl.length !== 2) return;
    found = true;
    for (const c of acyl) {
      const g = cloneGraph(g0);
      removeBond(g, n, c.atom);
      addBond(g, c.atom, addAtom(g, 'O'));
      out.push(...chem.graphToSmilesList(RDKit, g));
    }
  });
  if (!found) return [];
  return uniqSorted(out);
}

// ---------- 京大2016: ジアリールエーテルの水素化分解 ----------
// Ar–O–Ar' + H₂ → Ar–OH + Ar'–H。どちらの C–O が切れるかで2通り（両方の生成物を返す）
function etherHydrogenolysis(RDKit, g0) {
  const out = [];
  g0.atoms.forEach((a, o) => {
    if (a.el !== 'O') return;
    const nb = neighbors(g0, o);
    if (nb.length !== 2 || !nb.every((x) => g0.atoms[x.atom].arom)) return;
    for (const x of nb) {
      const g = cloneGraph(g0);
      removeBond(g, o, x.atom);
      out.push(...chem.graphToSmilesList(RDKit, g));
    }
  });
  return uniqSorted(out);
}

// ---------- 京大2025: 小員環のひずみと水素化による開環 ----------
// 三員環・四員環の C–C 結合が1本切れて H が2つ付く。切れる結合は、生成物に残る小員環がいちばん少なくなるもの
function smallRingBonds(g) {
  return g.bonds.filter((b) => b.order === 1 && el(g, b.a) === 'C' && el(g, b.b) === 'C' && !g.atoms[b.a].arom && ringSize(g, b) <= 4);
}
function ringSize(g, bond) {
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
  return dist.has(bond.b) ? dist.get(bond.b) + 1 : Infinity;
}
function ringHydrogenolysis(RDKit, g0) {
  const cands = smallRingBonds(g0);
  if (!cands.length) return [];
  const results = cands.map((b) => {
    const g = cloneGraph(g0);
    removeBond(g, b.a, b.b);
    return { g, strain: smallRingBonds(g).length };
  });
  const min = Math.min(...results.map((r) => r.strain));
  return uniqSorted(results.filter((r) => r.strain === min).map((r) => smilesOf(RDKit, r.g)));
}

// ---------- アセタール（京大2026 III・2020 IV・2013） ----------
// アセタール炭素: C=O をもたず、単結合の O が2つ（OR または OH）
function acetalCarbons(g) {
  return g.atoms.map((a, c) => c).filter((c) => el(g, c) === 'C' && !isCarbonylC(g, c)
    && neighbors(g, c).filter((x) => x.order === 1 && el(g, x.atom) === 'O').length === 2
    && neighbors(g, c).every((x) => x.order === 1));
}
// グリコシド結合: アセタール炭素の2つの O のうち、1つが環の中、もう1つが環の外で炭素とつながっているもの。
// 加水分解では環の外側の C–O が切れる（還元性を示す炭素に OH が戻る）。スクロースのように両側がアセタール炭素のときは片方だけ切る
function glycosideBonds(g) {
  const out = [];
  const seenO = new Set();
  for (const c of acetalCarbons(g)) {
    const os = neighbors(g, c).filter((x) => el(g, x.atom) === 'O');
    const ring = os.filter((x) => ringSize(g, x.bond) !== Infinity);
    const exo = os.filter((x) => ringSize(g, x.bond) === Infinity && neighbors(g, x.atom).filter((y) => el(g, y.atom) === 'C').length === 2);
    if (ring.length !== 1 || exo.length !== 1 || seenO.has(exo[0].atom)) continue;
    seenO.add(exo[0].atom);
    out.push({ c, bond: exo[0].bond });
  }
  return out;
}

// すべてのアセタール・ヘミアセタールをカルボニル化合物とアルコールに戻す
function openAcetals(g0) {
  const g = cloneGraph(g0);
  for (const c of acetalCarbons(g0)) {
    const os = neighbors(g, c).filter((x) => el(g, x.atom) === 'O').map((x) => x.atom);
    const oh = os.find((o) => hCount(g, o) === 1);
    if (oh !== undefined) {
      os.filter((o) => o !== oh).forEach((o) => removeBond(g, c, o));
      bondBetween(g, c, oh).order = 2;
    } else {
      os.forEach((o) => removeBond(g, c, o));
      addBond(g, c, addAtom(g, 'O'), 2);
    }
  }
  return g;
}
// 酸と水でアセタールを加水分解（アセトニドの脱保護にも使う）
function acetalHydrolysis(RDKit, g0) {
  if (!acetalCarbons(g0).length) return [];
  return productsOf(RDKit, openAcetals(g0), ['O']);
}
// 京大2026: 酸と大過剰のアルコール R'OH の中での平衡。
// アルデヒドは、同じ分子の OH と五員環・六員環をつくれるなら環状アセタール（もう一方は OR'）、つくれなければ R'O–CH–OR'。
// 放出される単純なアルコールと溶媒は生成物に数えない
function acetalExchange(solvent) {
  return (RDKit, g0) => {
    if (!acetalCarbons(g0).length && !g0.atoms.some((a, c) => isCarbonylC(g0, c) && hCount(g0, c) === 1)) return [];
    const opened = openAcetals(g0);
    const out = [];
    for (const part of chem.graphToSmilesList(RDKit, opened)) {
      const g = chem.graphFromSmiles(RDKit, part);
      const ald = g.atoms.map((a, c) => c).filter((c) => isCarbonylC(g, c) && hCount(g, c) === 1 && neighbors(g, c).length === 2);
      if (!ald.length) {
        // 単純なアルコール（O が1つだけの化合物）は放出されたものとして除く
        const nO = g.atoms.filter((a) => a.el === 'O').length;
        if (!(nO === 1 && g.atoms.some((a, i) => a.el === 'O' && hCount(g, i) === 1))) out.push(part);
        continue;
      }
      if (ald.length > 1) throw new Unsupported('several aldehydes');
      const c = ald[0];
      const oCarbonyl = neighbors(g, c).find((x) => x.order === 2).atom;
      const ohs = g.atoms.map((a, i) => i).filter((i) => el(g, i) === 'O' && hCount(g, i) === 1 && !isCarbonylC(g, neighbors(g, i)[0].atom));
      const dist = pathLengths(g, c);
      const ringOH = ohs.filter((o) => dist.get(o) === 4 || dist.get(o) === 5); // 環の員数 = 経路の原子数
      const variants = ringOH.length ? ringOH.map((o) => ({ o })) : [{ o: null }];
      for (const v of variants) {
        const h = cloneGraph(g);
        bondBetween(h, c, oCarbonyl).order = 1;
        // カルボニルの O を溶媒のアルキル基につなぐ（OR'）
        attachAlkyl(h, oCarbonyl, solvent);
        if (v.o !== null) addBond(h, c, v.o);
        else { const o2 = addAtom(h, 'O'); addBond(h, c, o2); attachAlkyl(h, o2, solvent); }
        out.push(smilesOf(RDKit, h));
      }
    }
    return uniqSorted(out);
  };
}
function attachAlkyl(g, o, solvent) {
  // solvent: 'C'（メタノール）か 'CC'（エタノール）
  let prev = o;
  for (let k = 0; k < solvent.length; k++) { const c = addAtom(g, 'C'); addBond(g, prev, c); prev = c; }
}
function pathLengths(g, s) {
  const dist = new Map([[s, 0]]);
  const q = [s];
  while (q.length) {
    const x = q.shift();
    for (const n of neighbors(g, x)) if (!dist.has(n.atom)) { dist.set(n.atom, dist.get(x) + 1); q.push(n.atom); }
  }
  return dist;
}

// 京大2020 IV: アセトンによるヒドロキシ基の保護（アセトニド）。
// 隣り合う炭素の OH（1,2 → 五員環）か1つおいた炭素の OH（1,3 → 六員環）。OH が3つ以上なら、より近い組（1,2）が優先
function acetonide(RDKit, g0) {
  const ohs = g0.atoms.map((a, i) => i).filter((i) => el(g0, i) === 'O' && hCount(g0, i) === 1 && !g0.atoms[neighbors(g0, i)[0].atom].arom && !isCarbonylC(g0, neighbors(g0, i)[0].atom));
  const pairs = [];
  for (let i = 0; i < ohs.length; i++) {
    for (let j = i + 1; j < ohs.length; j++) {
      const ci = neighbors(g0, ohs[i])[0].atom, cj = neighbors(g0, ohs[j])[0].atom;
      const d = pathLengths(g0, ci).get(cj);
      if (d === 1 || d === 2) pairs.push({ a: ohs[i], b: ohs[j], d });
    }
  }
  if (!pairs.length) return [];
  const best = Math.min(...pairs.map((p) => p.d));
  const out = pairs.filter((p) => p.d === best).map((p) => {
    const g = cloneGraph(g0);
    const k = addAtom(g, 'C');
    addBond(g, k, p.a); addBond(g, k, p.b);
    addBond(g, k, addAtom(g, 'C')); addBond(g, k, addAtom(g, 'C'));
    return smilesOf(RDKit, g);
  });
  return uniqSorted(out);
}

// ---------- 京大2026 IV・糖の構造決定: メチル化分析 ----------
// すべての OH をメチル化（-OCH₃）してから、グリコシド結合（環のアセタール）だけを加水分解する。
// 結合に使われていた位置が OH として残るので、つながり方がわかる
function methylationAnalysis(RDKit, g0) {
  const g = cloneGraph(g0);
  const ohs = g.atoms.map((a, i) => i).filter((i) => el(g, i) === 'O' && hCount(g, i) === 1 && !isCarbonylC(g, neighbors(g, i)[0].atom));
  if (!ohs.length) return [];
  ohs.forEach((o) => addBond(g, o, addAtom(g, 'C')));
  // グリコシド: アセタール炭素のうち、一方の O が同じ環の中にあるもの。環の外側の C–O を切って OH にする
  for (const c of acetalCarbons(g)) {
    const os = neighbors(g, c).filter((x) => el(g, x.atom) === 'O');
    const exo = os.filter((x) => ringSize(g, x.bond) === Infinity);
    if (exo.length !== 1) continue;
    removeBond(g, c, exo[0].atom);
    addBond(g, c, addAtom(g, 'O'));
  }
  return productsOf(RDKit, g, ['CO']); // メタノールは除く
}

// 臭素の付加（生成物。立体は stereo.js）
function bromineAddition(RDKit, g0) {
  const g = cloneGraph(g0);
  const dbl = g.bonds.filter((b) => b.order === 2 && el(g, b.a) === 'C' && el(g, b.b) === 'C' && !(g.atoms[b.a].arom && g.atoms[b.b].arom));
  if (!dbl.length) return [];
  for (const b of dbl) {
    b.order = 1;
    addBond(g, b.a, addAtom(g, 'Br'));
    addBond(g, b.b, addAtom(g, 'Br'));
  }
  return chem.graphToSmilesList(RDKit, g);
}

module.exports = {
  easPositions, nitration, bromination, chlorination, sulfonation, bromineWater, nitroReduction,
  acetylation, acetylationPrimary, deamination, diazoHydrolysis, azoCoupling, imideHydrolysis, etherHydrogenolysis,
  ringHydrogenolysis, acetalHydrolysis, acetalExchange, acetonide, methylationAnalysis, bromineAddition, glycosideBonds,
};
