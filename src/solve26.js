'use strict';
// 京大2026型の大問を「問題用紙に出ている情報だけ」から解く、生成器とは別の実装（別インスタンス解答）。
//  - 答え・解説・生成器の内部情報は publicView で取り除いてから渡す
//  - 規則の適用は reactions.js を使わず、規則の文（R1〜R4）どおりに一から組み直す
//  - 逆算は、前提に書かれた範囲の構造異性体を全部作り、書かれた事実で絞る（残りの数も返す）
//  - 計算は、問題文の数値と反応式の下の分子量だけから計算する
const chem = require('./chem');
const { enumerate } = require('./enumerate');

function publicView(P) {
  const strip = (q) => {
    const { answer, explain, refs, unique, model, rubric, ...rest } = q;
    const out = { ...rest };
    if (rest.calc) { const { answer: a, ...c } = rest.calc; out.calc = c; }
    return out;
  };
  return {
    id: P.id, parts: P.parts, rules: P.rules,
    experiments: P.experiments.map((e) => ({ label: e.label, part: e.part, group: e.group, text: e.text, figs: e.figs, fact: e.fact })),
    questions: P.questions.map(strip),
  };
}

// ---------- 小さなグラフ操作（この実装だけで使う） ----------
function toG(RDKit, s) {
  const g = chem.graphFromSmiles(RDKit, s);
  return { atoms: g.atoms.map((a) => ({ el: a.el, arom: !!a.arom, chg: a.chg || 0 })), bonds: g.bonds.map((b) => ({ ...b })) };
}
const nbs = (g, i) => g.bonds.filter((b) => b.a === i || b.b === i).map((b) => ({ j: b.a === i ? b.b : b.a, o: b.order, b }));
const hOf = (g, i) => chem.hCount(g, i);
function cut(g, i, j) { g.bonds = g.bonds.filter((b) => !((b.a === i && b.b === j) || (b.a === j && b.b === i))); }
function link(g, i, j, o = 1) { g.bonds.push({ a: i, b: j, order: o }); }
function atom(g, el) { g.atoms.push({ el, arom: false, chg: 0 }); return g.atoms.length - 1; }
function pieces(RDKit, g) {
  // 孤立した原子（外した O）を除いて、連結成分ごとの SMILES
  const alive = g.atoms.map((a, i) => i).filter((i) => !g.atoms[i].dead);
  const seen = new Set();
  const out = [];
  for (const s of alive) {
    if (seen.has(s)) continue;
    const comp = [];
    const st = [s];
    seen.add(s);
    while (st.length) { const x = st.pop(); comp.push(x); for (const n of nbs(g, x)) if (!seen.has(n.j) && !g.atoms[n.j].dead) { seen.add(n.j); st.push(n.j); } }
    const map = new Map(comp.map((x, k) => [x, k]));
    const sub = { atoms: comp.map((x) => g.atoms[x]), bonds: g.bonds.filter((b) => map.has(b.a) && map.has(b.b)).map((b) => ({ a: map.get(b.a), b: map.get(b.b), order: b.order })) };
    out.push(chem.canonical(RDKit, chem.toMolblock(sub)));
  }
  return out;
}

// アルデヒドの炭素（C=O と H）か、アセタール・ヘミアセタールの炭素（単結合の O が2つと H）を探す
function carbonylCenters(g) {
  return g.atoms.map((a, i) => i).filter((i) => {
    if (g.atoms[i].el !== 'C' || g.atoms[i].arom) return false;
    const os = nbs(g, i).filter((n) => g.atoms[n.j].el === 'O');
    if (hOf(g, i) < 1) return false; // ケトン由来のものは規則の対象外
    if (os.length === 1 && os[0].o === 2) return true;
    return os.length === 2 && os.every((n) => n.o === 1);
  });
}
// R1・R2 の逆: 中心の炭素の O をすべて外してアルデヒドに戻す。OR の O は相手の炭素側に OH として残す
function openAll(g) {
  const centers = carbonylCenters(g);
  for (const c of centers) {
    for (const n of nbs(g, c).filter((x) => g.atoms[x.j].el === 'O')) {
      cut(g, c, n.j);
      if (!nbs(g, n.j).length) g.atoms[n.j].dead = true; // OH・=O だった O は水として外れる
    }
    link(g, c, atom(g, 'O'), 2);
  }
  return g;
}
function dist(g, s) {
  const d = new Map([[s, 0]]);
  const q = [s];
  while (q.length) { const x = q.shift(); for (const n of nbs(g, x)) if (!d.has(n.j) && !g.atoms[n.j].dead) { d.set(n.j, d.get(x) + 1); q.push(n.j); } }
  return d;
}

// 規則 R1〜R4 を大過剰のアルコール（炭素数 nC の直鎖）と酸触媒の条件で当てはめ、主生成物を返す
function solveApply(RDKit, smiles, nC) {
  const g = openAll(toG(RDKit, smiles));
  const ald = carbonylCenters(g).filter((c) => nbs(g, c).some((n) => n.o === 2 && g.atoms[n.j].el === 'O'));
  for (const c of ald) {
    const oDbl = nbs(g, c).find((n) => n.o === 2).j;
    cut(g, c, oDbl);
    g.atoms[oDbl].dead = true;
    const d = dist(g, c);
    // R4: 同じ分子の OH で五員環・六員環（環の原子数 = C から O までの結合数 + 1）
    const ohs = g.atoms.map((a, i) => i).filter((i) => g.atoms[i].el === 'O' && !g.atoms[i].dead && nbs(g, i).length === 1
      && !nbs(g, nbs(g, i)[0].j).some((n) => n.o === 2) && (d.get(i) === 4 || d.get(i) === 5));
    if (ohs.length > 1) throw new Error('環をつくれる OH が2つ以上ある（規則で決まらない）');
    const alk = () => { const o = atom(g, 'O'); link(g, c, o); let prev = o; for (let k = 0; k < nC; k++) { const x = atom(g, 'C'); link(g, prev, x); prev = x; } };
    if (ohs.length === 1) { link(g, c, ohs[0]); alk(); } else { alk(); alk(); }
  }
  // 放出されたアルコール（O を1つだけもち、それが OH）と水は書かない
  return pieces(RDKit, g).filter((s) => {
    const h = toG(RDKit, s);
    const os = h.atoms.map((a, i) => i).filter((i) => h.atoms[i].el === 'O');
    return !(os.length === 1 && hOf(h, os[0]) === 1) && s !== 'O';
  }).sort();
}

// 加水分解（R2 の逆）: 生成物 → 物質量
function hydrolysis(RDKit, s) {
  const g = toG(RDKit, s);
  if (!carbonylCenters(g).some((c) => nbs(g, c).filter((n) => g.atoms[n.j].el === 'O' && n.o === 1).length === 2)) return null;
  const cnt = {};
  pieces(RDKit, openAll(g)).forEach((p) => { cnt[p] = (cnt[p] || 0) + 1; });
  return cnt;
}

// 性質（カードの判定は chem.evaluate を使う）
const AW = { H: 1.0, C: 12, O: 16, N: 14 };
const mw = (f) => Math.round([...f.matchAll(/([A-Z][a-z]?)(\d*)/g)].reduce((m, [, el, n]) => m + AW[el] * (n ? +n : 1), 0) * 10) / 10;
function prop(RDKit, card, s) {
  try {
    if (card === 'mw') return mw(chem.formula(chem.graphFromSmiles(RDKit, s)));
    if (card === 'oxid') {
      const p = chem.evaluate(RDKit, 'mild_oxidation', s);
      if (!p.length) return 'none';
      return p.some((x) => chem.evaluate(RDKit, 'silver_mirror', x)) ? 'ald' : 'ket';
    }
    if (card === 'kmno4_formula') { const p = chem.evaluate(RDKit, 'kmno4', s); return p.length ? p.map((x) => chem.formula(chem.graphFromSmiles(RDKit, x))).sort().join(',') : 'none'; }
    return chem.evaluate(RDKit, card, s);
  } catch (e) { return 'x'; }
}
const same = (a, b) => JSON.stringify(a) === JSON.stringify(b);
function universe(RDKit, premise, allowAcetal) {
  const opts = premise.arom ? { seed: 'c1ccccc1', rings: 0, allowAcetal } : { rings: 0, allowAcetal };
  return [...new Set(enumerate(RDKit, premise.formula, opts).map((s) => chem.canonical(RDKit, s)))];
}
const perms = (a) => (a.length <= 1 ? [a] : a.flatMap((x, i) => perms([...a.slice(0, i), ...a.slice(i + 1)]).map((p) => [x, ...p])));

// 逆算: D（アセタール）。事実: 前提（分子式・環や多重結合）、Na と反応しない、加水分解の物質量、生成物の性質
function solveD(RDKit, V) {
  const exps = V.experiments.filter((e) => e.part === '(a)' && e.fact);
  const pre = exps.find((e) => e.fact.premise);
  const hyd = exps.find((e) => e.fact.card === 'hydrolysis_pattern');
  const counts = hyd.fact.value;
  const labels = Object.keys(counts);
  const props = exps.filter((e) => labels.includes(e.fact.subject));
  const sat = [];
  for (const s of universe(RDKit, pre.fact.premise, true)) {
    if (prop(RDKit, 'sodium', s) !== pre.fact.value) continue;
    const cnt = hydrolysis(RDKit, s);
    if (!cnt || Object.keys(cnt).length !== labels.length) continue;
    // 記号と生成物の対応で、事実と矛盾しないものをすべて集める（記号が一意に決まらなければ、その設問は答えが1つにならない）
    const ok = perms(Object.keys(cnt)).filter((ps) => labels.every((L, i) => cnt[ps[i]] === counts[L]
      && props.filter((e) => e.fact.subject === L).every((e) => same(prop(RDKit, e.fact.card, ps[i]), e.fact.value))));
    if (ok.length) sat.push({ s, byLabel: Object.fromEntries(labels.map((L, i) => [L, [...new Set(ok.map((ps) => ps[i]))]])) });
  }
  return sat;
}
function solveG(RDKit, V) {
  const partB = V.parts.find((p) => p.label === '(b)');
  const facts = V.experiments.filter((e) => e.part === '(b)' && e.fact && e.fact.subject === 'G');
  return universe(RDKit, partB.premise, false).filter((s) => prop(RDKit, 'chiral', s) === partB.premise.chiral
    && facts.every((e) => same(prop(RDKit, e.fact.card, s), e.fact.value)));
}

function solve(RDKit, V) {
  const out = {};
  const nC = V.parts[0].solventC;
  const D = solveD(RDKit, V);
  const G = solveG(RDKit, V);
  const eq = V.parts.find((p) => p.label === '(b)').equation.masses;
  for (const q of V.questions) {
    if (q.kind === 'apply') {
      const fig = V.parts[0].figs.find((f) => q.prompt.includes(`化合物 ${f.label} `));
      out[q.id] = { smiles: solveApply(RDKit, fig.smiles, nC) };
    } else if (q.id === 'q2D') out[q.id] = { smiles: D.map((x) => x.s), candidates: D.length };
    else if (q.id === 'q2E') out[q.id] = { smiles: [...new Set(D.flatMap((x) => x.byLabel.E))], candidates: D.length };
    else if (q.id === 'q3G') out[q.id] = { smiles: G, candidates: G.length };
    else if (q.type === 'number' && q.calc) {
      const c = q.calc;
      const r = (x) => Number(x.toPrecision(q.sig_figs));
      if (c.kind === 'resolution_X') out[q.id] = { value: r((c.w / eq.G / 2) * eq.H) };
      else if (c.kind === 'resolution_Y_remaining') out[q.id] = { value: r((c.w / eq.G - c.b / eq.H) * eq.G) };
      else if (c.kind === 'resolution_Y_fraction') {
        const half = c.w / eq.G / 2, nH = c.b / eq.H;
        const plus = half - nH * c.p / (c.p + c.q), minus = half - nH * c.q / (c.p + c.q);
        out[q.id] = { value: r((100 * minus) / (plus + minus)) };
      } else if (c.kind === 'hydrolysis_mass') {
        // D の分子式と、F（D の加水分解生成物）の構造は解いた D から
        const lab = q.prompt.match(/化合物 ([A-Z])/)[1];
        const hyd = V.experiments.find((e) => e.fact && e.fact.card === 'hydrolysis_pattern').fact.value;
        const cand = D.length === 1 ? D[0].byLabel[lab] : [];
        if (cand.length === 1) out[q.id] = { value: r((c.w / mw(chem.formula(chem.graphFromSmiles(RDKit, D[0].s)))) * hyd[lab] * mw(chem.formula(chem.graphFromSmiles(RDKit, cand[0])))) };
      }
    }
  }
  return out;
}

// 想定解と一致するか
function agree(RDKit, P, sol) {
  const errs = [];
  const can = (s) => chem.canonical(RDKit, s);
  for (const q of P.questions) {
    if (q.type === 'essay') continue;
    const s = sol[q.id];
    if (!s) { errs.push(`${q.id}: 別解答が答えを出せない`); continue; }
    if (q.type === 'draw') {
      const a = q.answer.smiles.map(can).sort();
      const b = s.smiles.map(can).sort();
      if (!same(a, b)) errs.push(`${q.id}: 想定解 ${a} と別解答 ${b} が違う`);
      if (s.candidates !== undefined && s.candidates !== 1) errs.push(`${q.id}: 事実を満たす構造が ${s.candidates} 個（一意でない）`);
    } else if (q.type === 'number') {
      if (s.value !== q.answer.value) errs.push(`${q.id}: 想定解 ${q.answer.value} と別解答 ${s.value} が違う`);
    }
  }
  return errs;
}

module.exports = { publicView, solve, agree, solveApply, hydrolysis };
