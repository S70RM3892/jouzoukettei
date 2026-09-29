'use strict';
// 化合物どうしの関係を使う大問（京大型）の手がかり設計と検証。
// 断片 A, B, C… のそれぞれについて、
//  - 直接の手がかり（A は銀鏡反応を示す…）
//  - 誘導体の手がかり（A を酸化すると D が得られ、D は銀鏡反応を示す…）
//  - 断片どうしの関係（A と B に水素を付加させると同じ化合物、A を酸化すると C が得られる…）
// を混ぜ、全部を同時に考えたときに初めて答えが1通りに決まる組を選ぶ。
const chem = require('./chem');

// 誘導体をつくる操作（1つの生成物になるもの）
const OPS = ['mild_oxidation', 'kmno4', 'dehydration', 'hydrogenation', 'markovnikov', 'acetylation_primary', 'bromine_addition', 'kmno4_cleave', 'ozonolysis', 'bromine_water', 'nitration'];
// 化合物どうしの関係に使う操作
const REL_OPS = ['hydrogenation', 'dehydration', 'mild_oxidation', 'kmno4', 'markovnikov', 'kmno4_cleave'];
// 手がかりにするカード（誘導体にも使う）
const CLUE_CARDS = ['silver_mirror', 'iodoform', 'fecl3', 'sodium', 'nahco3', 'bromine', 'chiral', 'cis_trans', 'carbon_env', 'cl_sub', 'ring_cl',
  'anhydride', 'mild_oxidation', 'kmno4', 'dehydration', 'h2_uptake', 'hydrogenation', 'ozonolysis', 'kmno4_cleave', 'bromine_water', 'nitration', 'stereo_count'];

function makeTools(RDKit) {
  const cache = new Map();
  const ev = (card, s) => {
    const key = card + '|' + s;
    if (cache.has(key)) return cache.get(key);
    let v;
    try { v = JSON.stringify(chem.normalizeResult(RDKit, card, chem.evaluate(RDKit, card, s))); } catch (e) { v = 'x'; }
    cache.set(key, v);
    return v;
  };
  // 操作の生成物（ちょうど1つのときだけ）
  const op = (o, s) => {
    const v = ev(o, s);
    if (v === 'x') return null;
    const a = JSON.parse(v);
    return Array.isArray(a) && a.length === 1 ? a[0] : null;
  };
  return { ev, op, can: (s) => chem.canonical(RDKit, s) };
}

// 手がかりの集合から、各断片の候補（添字）と、関係まで含めた解の数を求める
function solve(frs, feats) {
  const dom = frs.map((f, i) => f.pool.map((_, k) => k).filter((k) => feats.every((x) => x.f !== i || x.test(k))));
  const rels = feats.filter((x) => x.pair);
  const involved = [...new Set(rels.flatMap((x) => x.pair))];
  let count = 1;
  const solutions = [];
  frs.forEach((_, i) => { if (!involved.includes(i)) count *= dom[i].length; });
  if (!involved.length) return { dom, count, pairs: 1 };
  // 関係のある断片だけ総当たり
  let n = 0;
  const cur = {};
  const rec = (d) => {
    if (d === involved.length) {
      if (rels.every((x) => x.test2(cur[x.pair[0]], cur[x.pair[1]]))) { n++; if (solutions.length < 3) solutions.push({ ...cur }); }
      return;
    }
    const i = involved[d];
    for (const k of dom[i]) {
      cur[i] = k;
      rec(d + 1);
      if (n > 50) return;
    }
  };
  rec(0);
  // 関係で絞られた後の各断片の候補
  const dom2 = dom.map((ds, i) => (involved.includes(i) ? ds.filter((k) => solutions.some((s) => s[i] === k) || n > 3) : ds));
  return { dom: dom2, count: count * n, pairs: n };
}

// P（buildProblem の出力）の断片の手がかりを、関係を使う組に作り直す
function design(RDKit, r, P, { labels = 'DEFGHIJ' } = {}) {
  const T = makeTools(RDKit);
  const frs = P.fragments.map((f, i) => ({ i, f, pool: f.given ? [T.can(f.answer)] : f.candidates.map(T.can) }))
    .map((x) => ({ ...x, a: x.pool.indexOf(T.can(x.f.answer)) }));
  if (frs.some((x) => x.a < 0)) return null;
  const feats = [];
  // 直接の手がかり
  frs.forEach((x, i) => {
    if (x.f.given) return;
    for (const c of CLUE_CARDS) {
      const vals = x.pool.map((s) => T.ev(c, s));
      if (vals[x.a] === 'x' || new Set(vals).size < 2) continue;
      feats.push({ kind: 'direct', f: i, card: c, value: vals[x.a], test: (k) => vals[k] === vals[x.a], w: 1 });
    }
  });
  // 誘導体（断片ごとに1つまで）
  const derived = [];
  let li = 0;
  frs.forEach((x, i) => {
    if (x.f.given || x.pool.length < 3) return;
    const ops = OPS.map((o) => ({ o, prods: x.pool.map((s) => T.op(o, s)) }))
      .filter(({ prods }) => prods[x.a] && new Set(prods.filter(Boolean)).size >= 2);
    if (!ops.length || r() < 0.2) return;
    const { o, prods } = ops[Math.floor(r() * ops.length)];
    const label = labels[li++];
    const d = { label, from: i, op: o, answer: prods[x.a], prods };
    derived.push(d);
    // 「A を〜すると D が得られた」という事実そのもの（反応が起こる）
    feats.push({ kind: 'applies', f: i, d, test: (k) => prods[k] !== null, w: 0, fixed: true });
    for (const c of CLUE_CARDS) {
      if (c === o) continue;
      const vals = prods.map((p) => (p ? T.ev(c, p) : 'x'));
      if (vals[x.a] === 'x' || new Set(vals.filter((v) => v !== 'x')).size < 2) continue;
      feats.push({ kind: 'derived', f: i, d, card: c, value: vals[x.a], test: (k) => vals[k] === vals[x.a], w: 1.4 });
    }
  });
  // 断片どうしの関係
  for (let i = 0; i < frs.length; i++) {
    for (let j = 0; j < frs.length; j++) {
      if (i === j) continue;
      const A = frs[i], B = frs[j];
      if (A.f.given && B.f.given) continue;
      for (const o of REL_OPS) {
        const pa = T.op(o, A.pool[A.a]);
        if (!pa) continue;
        // A を o すると B が得られる
        if (pa === B.pool[B.a]) {
          feats.push({ kind: 'yields', pair: [i, j], op: o, test2: (p, q) => T.op(o, A.pool[p]) === B.pool[q], w: 2.2 });
        }
        // A と B に o を行うと同じ化合物（i < j で1回だけ）
        if (i < j && pa === T.op(o, B.pool[B.a]) && o !== 'mild_oxidation') {
          feats.push({ kind: 'same', pair: [i, j], op: o, product: pa, test2: (p, q) => { const x = T.op(o, A.pool[p]); return x !== null && x === T.op(o, B.pool[q]); }, w: 2.2 });
        }
      }
    }
  }
  // 選ぶ: 関係 → 誘導体の手がかり → 直接の手がかり の順に重く、解が1つになるまで足す
  const sel = feats.filter((x) => x.fixed);
  const rels = feats.filter((x) => x.pair).sort(() => r() - 0.5).slice(0, 2);
  sel.push(...rels);
  const total = () => solve(frs, sel).count;
  let guard = 0;
  while (total() > 1 && guard++ < 40) {
    const now = solve(frs, sel);
    let best = null;
    for (const x of feats) {
      if (sel.includes(x) || x.pair) continue;
      // すでに1つに決まった断片の手がかりは足さない
      if (now.dom[x.f].length <= 1) continue;
      const after = now.dom[x.f].filter((k) => x.test(k)).length;
      if (after === now.dom[x.f].length) continue;
      // 1枚で全候補から1つに決まる直接の手がかりは避ける（京大らしく組み合わせで決めさせる）
      const single = x.kind === 'direct' && frs[x.f].pool.filter((_, k) => x.test(k)).length === 1;
      const score = (Math.log(now.dom[x.f].length / Math.max(1, after)) + 0.3) * x.w * (single ? 0.25 : 1) + r() * 0.3;
      if (!best || score > best.score) best = { x, score };
    }
    if (!best) return null;
    sel.push(best.x);
  }
  if (total() !== 1) return null;
  // 余分な手がかりを外す（関係と「反応が起こる」は残す）
  // ただし候補が4つ以上の断片には、それ自身や誘導体についての事実を2つ以上残す（1つで決まる問題にしない）
  const facts = (set, i) => set.filter((y) => y.f === i && (y.kind === 'direct' || y.kind === 'derived')).length;
  for (const x of [...sel].sort(() => r() - 0.5)) {
    if (x.fixed || x.pair) continue;
    const rest = sel.filter((y) => y !== x);
    if (frs[x.f].pool.length >= 4 && facts(rest, x.f) < 2) continue;
    if (solve(frs, rest).count === 1) sel.splice(sel.indexOf(x), 1);
  }
  // 使われなかった誘導体は消す（手がかりがない誘導体は問題文に出さない）
  const usedD = derived.filter((d) => sel.some((x) => x.kind === 'derived' && x.d === d));
  const final = sel.filter((x) => x.kind !== 'applies' || usedD.includes(x.d));
  if (solve(frs, final).count !== 1) return null;
  // 書き戻す
  const fragments = P.fragments.map((f, i) => ({
    ...f,
    clues: f.given ? [] : final.filter((x) => x.kind === 'direct' && x.f === i).map((x) => ({ card: x.card, result: JSON.parse(x.value) })),
  }));
  const ders = usedD.map((d, k) => ({
    label: labels[k], from: P.fragments[d.from].label, op: d.op, answer: d.answer,
    candidates: [...new Set(d.prods.filter(Boolean))],
    clues: final.filter((x) => x.kind === 'derived' && x.d === d).map((x) => ({ card: x.card, result: JSON.parse(x.value) })),
  }));
  const relations = final.filter((x) => x.pair).map((x) => ({ type: x.kind, op: x.op, a: P.fragments[x.pair[0]].label, b: P.fragments[x.pair[1]].label, ...(x.product ? { product: x.product } : {}) }));
  return { fragments, derived: ders, relations };
}

// 検証: P（chain 付き）の手がかりを全部当てはめたとき、各断片の答えがただ1通りに決まるか
function check(RDKit, P) {
  const T = makeTools(RDKit);
  const errors = [];
  const byLabel = Object.fromEntries(P.fragments.map((f, i) => [f.label, i]));
  const frs = P.fragments.map((f) => ({ f, pool: f.given ? [T.can(f.answer)] : f.candidates.map(T.can) })).map((x) => ({ ...x, a: x.pool.indexOf(T.can(x.f.answer)) }));
  frs.forEach((x) => { if (x.a < 0) errors.push(`${x.f.label}: answer not in candidates`); });
  if (errors.length) return errors;
  const feats = [];
  frs.forEach((x, i) => (x.f.clues || []).forEach((c) => {
    const exp = JSON.stringify(chem.normalizeResult(RDKit, c.card, c.result));
    if (T.ev(c.card, x.pool[x.a]) !== exp) errors.push(`${x.f.label}: clue ${c.card} does not match the answer`);
    feats.push({ f: i, test: (k) => T.ev(c.card, x.pool[k]) === exp });
  }));
  for (const d of P.derived || []) {
    const i = byLabel[d.from];
    const x = frs[i];
    if (T.op(d.op, x.pool[x.a]) !== T.can(d.answer)) errors.push(`${d.label}: ${d.op}(${d.from}) is not the answer`);
    feats.push({ f: i, test: (k) => T.op(d.op, x.pool[k]) !== null });
    d.clues.forEach((c) => {
      const exp = JSON.stringify(chem.normalizeResult(RDKit, c.card, c.result));
      if (T.ev(c.card, T.can(d.answer)) !== exp) errors.push(`${d.label}: clue ${c.card} does not match`);
      feats.push({ f: i, test: (k) => { const p = T.op(d.op, x.pool[k]); return p !== null && T.ev(c.card, p) === exp; } });
    });
  }
  for (const rl of P.relations || []) {
    const i = byLabel[rl.a], j = byLabel[rl.b];
    const A = frs[i], B = frs[j];
    const test2 = rl.type === 'yields' ? (p, q) => T.op(rl.op, A.pool[p]) === B.pool[q]
      : (p, q) => { const x = T.op(rl.op, A.pool[p]); return x !== null && x === T.op(rl.op, B.pool[q]); };
    if (!test2(A.a, B.a)) errors.push(`relation ${rl.type} ${rl.a}-${rl.b} does not hold`);
    feats.push({ pair: [i, j], test2 });
  }
  const s = solve(frs, feats);
  if (s.count !== 1) errors.push(`joint solutions: ${s.count}`);
  return errors;
}

module.exports = { design, check, OPS };
