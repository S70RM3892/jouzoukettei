'use strict';
// 答えのリーク（問題文・図に出す情報が、設問の答えそのもの・答えを読み取れるものになっていないか）を判定する。
//  same:        図の化合物が答えと同じ（立体を除いて）
//  substructure: 答えの重原子と結合がそのまま図の中にある（例: D のニトロ化生成物の図に D が丸ごと入っている）
//  skeleton:    答えの炭素骨格がそのまま図の中にあり、炭素数も同じ（例: B に臭素を付加した化合物の図から B が読める）
// 実験の生成物がリークになるときは、構造式の代わりに「分子式と種類の数」だけを示す（shownValue）。
const chem = require('./chem');

const cache = new Map();
function graph(RDKit, s) {
  const k = s.replace(/[/\\]/g, '').replace(/@+/g, '');
  if (cache.has(k)) return cache.get(k);
  const g = chem.graphFromSmiles(RDKit, k);
  if (cache.size > 20000) cache.clear();
  cache.set(k, g);
  return g;
}

function adjacency(g) {
  const adj = g.atoms.map(() => new Map());
  for (const b of g.bonds) {
    const o = g.atoms[b.a].arom && g.atoms[b.b].arom ? 'a' : b.order;
    adj[b.a].set(b.b, o);
    adj[b.b].set(b.a, o);
  }
  return adj;
}

// small の原子（keep で選んだもの）が big の原子に単射で写り、small の結合が big の結合になるか
function embeds(small, big, keepS, keepB, atomEq, bondEq) {
  const sIdx = small.atoms.map((_, i) => i).filter(keepS);
  const bIdx = big.atoms.map((_, i) => i).filter(keepB);
  if (sIdx.length > bIdx.length) return false;
  if (!sIdx.length) return true;
  const sa = adjacency(small);
  const ba = adjacency(big);
  // 次数の大きい原子から、つながった順に並べる（探索を早く打ち切る）
  const order = [];
  const seen = new Set();
  const deg = (i) => [...sa[i].keys()].filter(keepS).length;
  while (order.length < sIdx.length) {
    const start = sIdx.filter((i) => !seen.has(i)).sort((x, y) => deg(y) - deg(x))[0];
    const q = [start];
    seen.add(start);
    while (q.length) {
      const x = q.shift();
      order.push(x);
      for (const n of sa[x].keys()) if (keepS(n) && !seen.has(n)) { seen.add(n); q.push(n); }
    }
  }
  const map = new Map();
  const used = new Set();
  const rec = (d) => {
    if (d === order.length) return true;
    const s = order[d];
    for (const b of bIdx) {
      if (used.has(b) || !atomEq(small.atoms[s], big.atoms[b])) continue;
      if ([...sa[s].keys()].filter(keepS).length > [...ba[b].keys()].filter(keepB).length) continue;
      let ok = true;
      for (const [n, o] of sa[s]) {
        if (!keepS(n) || !map.has(n)) continue;
        const bo = ba[b].get(map.get(n));
        if (bo === undefined || !bondEq(o, bo)) { ok = false; break; }
      }
      if (!ok) continue;
      map.set(s, b); used.add(b);
      if (rec(d + 1)) return true;
      map.delete(s); used.delete(b);
    }
    return false;
  };
  return rec(0);
}

const nC = (g) => g.atoms.filter((a) => a.el === 'C').length;

// 答え ans が図 fig から読み取れるか。読み取れるなら理由（same / substructure / skeleton）
function leakKind(RDKit, ans, fig) {
  let a, f;
  try { a = graph(RDKit, ans); f = graph(RDKit, fig); } catch (e) { return null; }
  const strip = (s) => chem.canonical(RDKit, s.replace(/[/\\]/g, '').replace(/@+/g, ''));
  try { if (strip(ans) === strip(fig)) return 'same'; } catch (e) { /* 比べられない */ }
  const all = () => true;
  // 小さな答え（アセトアルデヒドなど）はどんな図にも部分として入るので、答えが図の半分以上を占めるときだけ数える
  if (a.atoms.length >= 5 && a.atoms.length * 2 >= f.atoms.length && embeds(a, f, all, all, (x, y) => x.el === y.el && !!x.arom === !!y.arom, (x, y) => x === y)) return 'substructure';
  const isC = (g) => (i) => g.atoms[i].el === 'C';
  if (nC(a) >= 3 && nC(a) === nC(f) && embeds(a, f, isC(a), isC(f), () => true, () => true)) return 'skeleton';
  return null;
}

// 生成物の構造式を見せると、反応させた化合物 subject が読めてしまうか
function productsLeak(RDKit, subject, products) {
  return products.some((p) => leakKind(RDKit, subject, p));
}

const DIPEP = '[NX3;H2][CX4][CX3](=O)[NX3;H1][CX4][CX3](=O)[OX2H1]';
function isDipeptide(RDKit, s) {
  const m = RDKit.get_mol(chem.expand(s));
  const q = RDKit.get_qmol(DIPEP);
  try { const r = JSON.parse(m.get_substruct_matches(q)); return Array.isArray(r) && r.length > 0; } finally { m.delete(); q.delete(); }
}

function formulaOf(RDKit, s) { return chem.formula(graph(RDKit, s)); }

// 問題文に出す形の結果。生成物の構造式がリークになるときは { n: 種類の数, formulas: 分子式 } に弱める
function shownValue(RDKit, card, subject, value) {
  const kind = chem.CARDS[card] && chem.CARDS[card].kind;
  if (kind === 'products' && Array.isArray(value) && value.length && productsLeak(RDKit, subject, value)) {
    return { n: value.length, formulas: value.map((p) => formulaOf(RDKit, p)).sort() };
  }
  if (kind === 'contains' && typeof value === 'string') {
    // ペプチドの部分加水分解（ジペプチドの配列）は、アミノ酸の記号（化合物の記号）の並びで書く（構造式は出さない。exam.js）
    if (isDipeptide(RDKit, value)) return value;
    // エステルの部分加水分解の生成物は断片をそのまま含むので、分子式・不斉炭素の数・ヨードホルム反応だけを示す
    let iodo = false;
    try { iodo = !!chem.evaluate(RDKit, 'iodoform', value); } catch (e) { /* 判定できない */ }
    return { formula: formulaOf(RDKit, value), chiral: chem.chiralCount(graph(RDKit, value)), iodoform: iodo };
  }
  return value;
}
const isWeak = (v) => v && typeof v === 'object' && !Array.isArray(v);

module.exports = { leakKind, productsLeak, shownValue, isWeak, embeds, formulaOf };
