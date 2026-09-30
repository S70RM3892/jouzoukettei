'use strict';
// 京大2026年度 第3問の「設計思想」（本文は写さない）にそろえた大問を自動で作る。
//  (a) 規則の定義 → 適用問題（構造式を書く）→ 逆算問題（分子式と反応性から構造を書く）
//  (b) 光学分割（触媒 X は一方の鏡像異性体だけ、触媒 Y は両方を速さを変えて反応させる）→ 構造の逆算 → 量的設問（質量・割合・有効数字）
// 出題ルール:
//  - 先に答え（SMILES・数値）を決め、実験事実は答えから逆算して作る
//  - 問題文・図に出すものが設問の答えそのもの・答えを読み取れるものにならない（audit.js で検査）
//  - 分子量・質量は分子式から自動計算して埋める（手書きしない）
//  - 構造が決まって初めて答えが出る設問は、その構造の設問より後に置く（depends_on）
//  - 逆算の設問は、分子式の構造異性体を全部作り、書いた事実を満たすものがちょうど1つであることを確かめてから出す
const chem = require('./chem');
const rx = require('./reactions');
const { enumerate } = require('./enumerate');
const { sentence, infer } = require('./inference');

// 同じ分子式の列挙は何度も使うので覚えておく
const ENUM = new Map();
function enumCached(RDKit, f, opts) {
  const k = f + JSON.stringify(opts);
  if (!ENUM.has(k)) ENUM.set(k, [...new Set(enumerate(RDKit, f, opts).map((s) => chem.canonical(RDKit, s)))]);
  return ENUM.get(k);
}

// ---------- 乱数（種つき・再現できる） ----------
function rng(seed) {
  let s = seed >>> 0 || 1;
  return () => { s ^= s << 13; s >>>= 0; s ^= s >> 17; s ^= s << 5; s >>>= 0; return s / 4294967296; };
}
const pick = (r, a) => a[Math.floor(r() * a.length)];
const shuffle = (r, a) => a.map((x) => [r(), x]).sort((p, q) => p[0] - q[0]).map((x) => x[1]);

// ---------- 数値 ----------
const AW = { H: 1.0, C: 12, O: 16, N: 14 }; // 京大で与えられる原子量
function molarMass(f) {
  let m = 0;
  for (const [, el, n] of f.matchAll(/([A-Z][a-z]?)(\d*)/g)) m += AW[el] * (n ? +n : 1);
  return Math.round(m * 10) / 10;
}
const sig = (x, n) => Number(Number(x).toPrecision(n));
const fmt = (x, n) => Number(x).toPrecision(n); // 有効数字 n けたの表記（末尾の 0 を残す）
const subf = (f) => f.replace(/\d/g, (d) => '₀₁₂₃₄₅₆₇₈₉'[d]);
const mStr = (M) => (Number.isInteger(M) ? String(M) : M.toFixed(1));

// ---------- 難易度（1〜5）→ パラメータ ----------
// 規則の複雑さ（適用回数・優先順位の数）、逆算の段数、立体化学の有無、量的計算の段数
const LEVELS = {
  1: { nApply: 1, ring: false, trap: false, mixed: false, allowMW: true, calcA: false, Y: false, remain: false, essay: false },
  2: { nApply: 2, ring: true, trap: false, mixed: false, allowMW: true, calcA: false, Y: true, remain: false, essay: false },
  3: { nApply: 2, ring: true, trap: false, mixed: false, allowMW: false, calcA: true, Y: true, remain: false, essay: true },
  4: { nApply: 3, ring: true, trap: true, mixed: false, allowMW: false, calcA: true, Y: true, remain: true, essay: true },
  5: { nApply: 3, ring: true, trap: true, mixed: true, allowMW: false, calcA: true, Y: true, remain: true, essay: true },
};

// ---------- (a) アセタールの規則 ----------
const ALK = { C: 'メタノール', CC: 'エタノール', CCC: '1-プロパノール' };
const RULES_A = (solvent) => [
  { id: 'R1', title: 'ヘミアセタールの生成',
    text: 'アルデヒド R¹–CHO にアルコール R²–OH が付加すると、ヘミアセタール R¹–CH(OH)–O–R² が生じる。この反応は可逆で、ヘミアセタールは単離できず、もとのアルデヒド・アルコールとの平衡混合物として存在する。',
    scheme: ['O=CC1CCCCC1', '+', 'CO', '⇄', 'COC(O)C1CCCCC1'] },
  { id: 'R2', title: 'アセタールの生成と加水分解',
    text: '酸触媒があると、ヘミアセタールはもう1分子のアルコール R³–OH と反応して水を失い、アセタール R¹–CH(OR²)(OR³) になる。この反応も可逆で、酸触媒と多量の水があると、アセタールはもとのアルデヒドとアルコールに戻る（加水分解）。',
    scheme: ['COC(O)C1CCCCC1', '+', 'CO', '⇄', 'COC(OC)C1CCCCC1', '+', 'O'] },
  { id: 'R3', title: 'アセタール交換',
    text: `アセタールを、酸触媒とともに大過剰のアルコール R⁴–OH の中に十分な時間置くと、R1・R2 の平衡がくり返され、–OR² と –OR³ はすべて –OR⁴ に置き換わる（R⁴–OH が大過剰なので、平衡は R⁴–OH が結合した側に偏る）。`,
    scheme: ['COC(OC)C1CCCCC1', '→', 'CCOC(OCC)C1CCCCC1'], schemeNote: '大過剰のエタノール・酸触媒' },
  { id: 'R4', title: '環状構造の優先',
    text: '同じ分子の中に、アルデヒドの炭素と五員環または六員環をつくれる位置にヒドロキシ基があるときは、そのヒドロキシ基が R1 の R²–OH としてはたらき、環状のヘミアセタール・アセタールが主生成物となる。このとき R3 の条件でも環は開かず、環の外側の –OR だけが –OR⁴ に置き換わる。三員環・四員環・七員環以上の環はできない。',
    scheme: ['O=CCCC(O)c1ccccc1', '⇄', 'OC1CCC(c2ccccc2)O1'] },
];
const LIMIT_A = 'この問題では、R1〜R4 の反応のみが起こるものとする。酸触媒がなければ R2〜R4 は起こらない。ケトン・エーテル・エステルは反応しない。立体異性体は区別しない。';

// 鎖の番号 1 がアルデヒドの炭素。ohAt に OH（ohAt=0 なら OH なし）、methyls の位置にメチル基
function chainSmiles(spec, head) {
  let s = head;
  for (let i = 2; i <= spec.L; i++) {
    let a = 'C';
    if (spec.methyls.includes(i)) a += '(C)';
    if (i === spec.ohAt) a += spec.cyclic ? (i < spec.L ? '(O1)' : 'O1') : (i < spec.L ? '(O)' : 'O');
    // 分岐の後に鎖を続けるため、最後以外の OH は枝にする
    s += a;
  }
  return s;
}
const alkO = (a) => `O${a}`;
function reactantForm(RDKit, spec, form, alk1, alk2) {
  let smi;
  if (form === 'open') smi = chainSmiles({ ...spec, cyclic: false }, `C(${alkO(alk1)})(${alkO(alk2)})`);
  else if (form === 'hemi') smi = chainSmiles({ ...spec, cyclic: true }, 'C1(O)');
  else if (form === 'cyc') smi = chainSmiles({ ...spec, cyclic: true }, `C1(${alkO(alk1)})`);
  return chem.canonical(RDKit, smi);
}
const ringOf = (ohAt) => (ohAt === 4 ? 5 : ohAt === 5 ? 6 : 0);

function applyAcetalRules(RDKit, smiles, solvent) {
  return rx.acetalExchange(solvent)(RDKit, chem.graphFromSmiles(RDKit, smiles)).map((s) => chem.canonical(RDKit, s));
}

// 適用問題の反応物を作る。case: ring（R4 で環をつくる）/ exchange（R3 だけ）/ trap（OH があるが環をつくれない）/ hemi（環状ヘミアセタール）
function makeApply(RDKit, r, kase, solvent) {
  const others = Object.keys(ALK).filter((a) => a !== solvent);
  for (let guard = 0; guard < 60; guard++) {
    let spec;
    if (kase === 'exchange') spec = { L: 2 + Math.floor(r() * 3), ohAt: 0, methyls: [] };
    else if (kase === 'trap') { const oh = pick(r, [3, 7]); spec = { L: oh + (r() < 0.4 ? 1 : 0), ohAt: oh, methyls: [] }; }
    else { const oh = pick(r, [4, 5]); spec = { L: oh + (r() < 0.4 ? 1 : 0), ohAt: oh, methyls: [] }; }
    const nMe = Math.floor(r() * (kase === 'exchange' ? 2 : 2.4));
    for (let k = 0; k < nMe; k++) {
      const pos = 2 + Math.floor(r() * (spec.L - 1));
      if (!spec.methyls.includes(pos)) spec.methyls.push(pos);
    }
    const form = kase === 'hemi' ? 'hemi' : kase === 'ring' ? pick(r, ['open', 'open', 'cyc']) : 'open';
    const a1 = pick(r, others);
    const a2 = r() < 0.8 ? a1 : pick(r, Object.keys(ALK));
    let react;
    try { react = reactantForm(RDKit, spec, form, a1, a2); } catch (e) { continue; }
    const g = chem.graphFromSmiles(RDKit, react);
    const nC = chem.formula(g).match(/^C(\d*)/)[1];
    if (+nC > 11) continue;
    let prods;
    try { prods = applyAcetalRules(RDKit, react, solvent); } catch (e) { continue; }
    if (prods.length !== 1 || prods[0] === react) continue;
    return { case: kase, reactant: react, product: prods[0], ring: kase === 'ring' || kase === 'hemi' ? ringOf(spec.ohAt) : 0 };
  }
  return null;
}

// ---------- 逆算: 事実の評価 ----------
// 加水分解（R2 の逆）: アセタール炭素の C–O を全部切って、生成物と物質量（係数）を返す。アセタールでなければ null
function hydrolysisPattern(RDKit, s) {
  let prods;
  try { prods = rx.acetalHydrolysis(RDKit, chem.graphFromSmiles(RDKit, s)); } catch (e) { return null; }
  if (!prods || !prods.length) return null;
  const cnt = {};
  prods.forEach((p) => { cnt[p] = (cnt[p] || 0) + 1; });
  return cnt;
}
const PROP_CARDS = {
  silver_mirror: { f: (R, s) => chem.evaluate(R, 'silver_mirror', s) },
  iodoform: { f: (R, s) => chem.evaluate(R, 'iodoform', s) },
  chiral: { f: (R, s) => chem.evaluate(R, 'chiral', s) },
  fecl3: { f: (R, s) => chem.evaluate(R, 'fecl3', s) },
  sodium: { f: (R, s) => chem.evaluate(R, 'sodium', s) },
  oxid: { f: (R, s) => oxidClass(R, s) },
  mw: { f: (R, s) => molarMass(chem.formula(chem.graphFromSmiles(R, s))) },
  dehydration_count: { f: (R, s) => chem.evaluate(R, 'dehydration_count', s) },
  ring_cl: { f: (R, s) => chem.evaluate(R, 'ring_cl', s) },
  kmno4_formula: { f: (R, s) => { const p = chem.evaluate(R, 'kmno4', s); return p.length ? p.map((x) => chem.formula(chem.graphFromSmiles(R, x))).sort().join(',') : 'none'; } },
};
// 穏やかな酸化の生成物の種類: 'ald'（銀鏡反応を示す）/ 'ket'（示さない）/ 'none'（酸化されない）
function oxidClass(RDKit, s) {
  let p;
  try { p = chem.evaluate(RDKit, 'mild_oxidation', s); } catch (e) { return 'x'; }
  if (!p.length) return 'none';
  return p.some((x) => chem.evaluate(RDKit, 'silver_mirror', x)) ? 'ald' : 'ket';
}
function propValue(RDKit, card, s) {
  try { return PROP_CARDS[card].f(RDKit, s); } catch (e) { return 'x'; }
}
function propSentence(card, subj, v) {
  if (card === 'oxid') return v === 'none' ? `${subj} を硫酸酸性の二クロム酸カリウム水溶液で穏やかに酸化しようとしたが、酸化されなかった。`
    : `${subj} を硫酸酸性の二クロム酸カリウム水溶液で穏やかに酸化すると、銀鏡反応を${v === 'ald' ? '示す' : '示さない'}化合物が得られた。`;
  if (card === 'mw') return `${subj} の分子量は ${mStr(v)} であった。`;
  if (card === 'kmno4_formula') return v === 'none' ? `${subj} を過マンガン酸カリウムで十分に酸化しても、ベンゼン環の側鎖は変化しなかった。`
    : `${subj} を硫酸酸性の過マンガン酸カリウム水溶液で十分に酸化すると、分子式 ${v.split(',').map(subf).join('、')} の芳香族化合物が得られた。`;
  return sentence(card, subj, v, chem.CARDS);
}
function propHint(card, v) {
  if (card === 'oxid') return v === 'none' ? '第三級アルコール（または酸化される OH がない）' : v === 'ald' ? '第一級アルコール（–CH₂OH）' : '第二級アルコール（>CH–OH）';
  if (card === 'mw') return '分子量から分子式を決め、候補を絞る';
  if (card === 'kmno4_formula') return 'ベンゼン環に直結した炭素はすべて –COOH になる。COOH の数 = 側鎖の数';
  return infer(card, v, chem.CARDS);
}

// 事実のリストを満たす候補（labels: 生成物の記号 → 係数、props: [{label, card, value}]）
function hydroMatches(RDKit, s, pat) {
  const cnt = hydrolysisPattern(RDKit, s);
  if (!cnt) return false;
  const prods = Object.keys(cnt);
  const labels = Object.keys(pat.counts);
  if (prods.length !== labels.length) return false;
  // 記号と生成物の対応をすべて試す
  const perm = (arr) => (arr.length <= 1 ? [arr] : arr.flatMap((x, i) => perm([...arr.slice(0, i), ...arr.slice(i + 1)]).map((p) => [x, ...p])));
  return perm(prods).some((ps) => labels.every((L, i) => cnt[ps[i]] === pat.counts[L]
    && pat.props.filter((p) => p.label === L).every((p) => JSON.stringify(propValue(RDKit, p.card, ps[i])) === JSON.stringify(p.value))));
}

// ---------- (a) 逆算: アセタール D ----------
// アルデヒドは R–C=O の形（末尾が C=O）、アルコールのアルキル基は O につく原子から書く
const ALDS = ['CC=O', 'CCC=O', 'CCCC=O', 'CC(C)C=O', 'c1ccccc1C=O', 'c1ccccc1CC=O', 'CCC(C)C=O'];
const ALCS = ['C', 'CC', 'CCC', 'C(C)C'];
function makeInverseD(RDKit, r, knobs) {
  for (let guard = 0; guard < 40; guard++) {
    const ald = pick(r, ALDS);
    const a1 = pick(r, ALCS);
    const a2 = knobs.mixed ? pick(r, ALCS.filter((x) => x !== a1)) : a1;
    const D = chem.canonical(RDKit, acetalOf(ald, a1, a2));
    const f = chem.formula(chem.graphFromSmiles(RDKit, D));
    const nC = +(f.match(/^C(\d*)/)[1] || 1);
    if (nC > 10) continue;
    const arom = /c/.test(D);
    const opts = arom ? { seed: 'c1ccccc1', rings: 0, allowAcetal: true } : { rings: 0, allowAcetal: true };
    let pool;
    try { pool = enumCached(RDKit, f, opts); } catch (e) { continue; }
    if (!pool.includes(D)) continue;
    const cnt = hydrolysisPattern(RDKit, D);
    const prods = Object.keys(cnt);
    // 生成物の記号: カルボニル化合物を E、アルコールを F（, G）
    const carbonyl = prods.find((p) => propValue(RDKit, 'silver_mirror', p) === true);
    if (!carbonyl) continue;
    const alcs = prods.filter((p) => p !== carbonyl);
    const labelOf = { [carbonyl]: 'E' };
    alcs.forEach((p, i) => { labelOf[p] = 'FK'[i]; }); // G・H は (b) で使う
    const counts = Object.fromEntries(prods.map((p) => [labelOf[p], cnt[p]]));
    const pat = { counts, props: [] };
    // 事実を足していく（1つに決まるまで）。値は答えから計算する
    const cards = ['silver_mirror', 'iodoform', 'oxid', 'chiral', ...(arom ? ['fecl3'] : []), ...(knobs.allowMW ? ['mw'] : [])];
    const options = [];
    for (const p of prods) for (const c of cards) {
      if (c === 'silver_mirror' && p !== carbonyl) continue;
      if (c === 'oxid' && p === carbonyl) continue;
      options.push({ label: labelOf[p], card: c, value: propValue(RDKit, c, p) });
    }
    // 3種類の生成物がどれも 1 mol のときは、E がカルボニル化合物であることを事実として書く（記号の対応を1つに決める）
    if (knobs.mixed) pat.props.push({ label: 'E', card: 'silver_mirror', value: true });
    const base = pool.filter((s) => propValue(RDKit, 'sodium', s) === false && hydroMatches(RDKit, s, pat));
    let left = base;
    const chosen = [];
    // 分子量など決定的な事実は後回しにして、性質の組み合わせで絞らせる
    const order = shuffle(r, options).sort((x, y) => (x.card === 'mw') - (y.card === 'mw'));
    while (left.length > 1) {
      let best = null;
      for (const o of order) {
        if (chosen.includes(o)) continue;
        const trial = { counts, props: [...pat.props, o] };
        const l2 = left.filter((s) => hydroMatches(RDKit, s, trial));
        if (l2.length < left.length && (!best || l2.length < best.l.length || (l2.length === best.l.length && best.o.card === 'mw'))) best = { o, l: l2 };
      }
      if (!best) break;
      chosen.push(best.o);
      pat.props.push(best.o);
      left = best.l;
    }
    if (left.length !== 1 || left[0] !== D) continue;
    // 余分な事実を外す（1つ抜いても決まるなら抜く）
    for (const o of [...chosen]) {
      const trial = { counts, props: pat.props.filter((x) => x !== o) };
      if (pool.filter((s) => propValue(RDKit, 'sodium', s) === false && hydroMatches(RDKit, s, trial)).length === 1) pat.props = trial.props;
    }
    return { D, formula: f, arom, universe: opts, pool, pat, carbonyl, labelOf, nBase: base.length };
  }
  return null;
}
function acetalOf(ald, a1, a2) {
  // R–CHO → R–CH(OR1)(OR2)
  return `${ald.slice(0, -3)}C(O${a1})O${a2}`;
}

// ---------- (b) 光学分割: G ----------
const G_FORMULAS = [['C4H10O', false], ['C5H12O', false], ['C6H14O', false], ['C8H10O', true], ['C9H12O', true], ['C10H14O', true]];
function makeG(RDKit, r, knobs) {
  for (let guard = 0; guard < 30; guard++) {
    const [f, arom] = pick(r, G_FORMULAS);
    const opts = arom ? { seed: 'c1ccccc1', rings: 0 } : { rings: 0 };
    const pool = enumCached(RDKit, f, opts);
    // 不斉炭素を1個もつ第二級アルコールで、不斉炭素が OH の炭素（酸化するとケトンになり不斉炭素がなくなる）
    const Gs = pool.filter((s) => {
      try {
        if (chem.evaluate(RDKit, 'chiral', s) !== 1 || oxidClass(RDKit, s) !== 'ket' || chem.evaluate(RDKit, 'fecl3', s)) return false;
        const k = chem.evaluate(RDKit, 'mild_oxidation', s)[0];
        return chem.evaluate(RDKit, 'chiral', k) === 0;
      } catch (e) { return false; }
    });
    if (!Gs.length) continue;
    const G = pick(r, Gs);
    // 前提（問題文に書く）: 分子式、不斉炭素を1個もつ。事実はこの前提の上に足す
    const base = pool.filter((s) => propValue(RDKit, 'chiral', s) === 1);
    const cards = ['sodium', 'iodoform', 'oxid', 'dehydration_count', ...(arom ? ['fecl3', 'ring_cl', 'kmno4_formula'] : [])];
    const facts = [];
    let left = base;
    const order = shuffle(r, cards);
    while (left.length > 1) {
      let best = null;
      for (const c of order) {
        if (facts.some((x) => x.card === c)) continue;
        const v = propValue(RDKit, c, G);
        if (v === 'x') continue;
        const l2 = left.filter((s) => JSON.stringify(propValue(RDKit, c, s)) === JSON.stringify(v));
        if (l2.length < left.length && (!best || l2.length < best.l.length)) best = { c, v, l: l2 };
      }
      if (!best) break;
      facts.push({ card: best.c, value: best.v });
      left = best.l;
    }
    if (left.length !== 1 || left[0] !== G) continue;
    // 1枚で決まる問題は避ける（組み合わせて考えさせる）
    if (facts.length < (knobs === LEVELS[1] ? 1 : 2)) continue;
    const H = chem.canonical(RDKit, rx.acetylation(RDKit, chem.graphFromSmiles(RDKit, G))[0]);
    return { G, H, formula: f, arom, universe: opts, pool, facts, nBase: base.length };
  }
  return null;
}

// ---------- 大問を1つ作る ----------
function buildK26(RDKit, seed, level) {
  const r = rng(seed);
  const knobs = LEVELS[level];
  const solvent = pick(r, ['C', 'CC']);
  // ---- (a) 適用 ----
  const cases = ['exchange'];
  if (knobs.ring) cases.unshift(pick(r, ['ring', 'ring', 'hemi']));
  if (knobs.trap) cases.push('trap');
  while (cases.length < knobs.nApply) cases.push(pick(r, ['ring', 'exchange', 'hemi']));
  const applies = [];
  for (const c of cases.slice(0, knobs.nApply)) {
    const a = makeApply(RDKit, r, c, solvent);
    if (!a) return null;
    if (applies.some((x) => x.reactant === a.reactant || x.product === a.product)) return null;
    applies.push(a);
  }
  const Dq = makeInverseD(RDKit, r, knobs);
  if (!Dq) return null;
  const Gq = makeG(RDKit, r, knobs);
  if (!Gq) return null;

  const labelsP = ['P', 'Q', 'R'];
  const questions = [];
  const experiments = [];
  let ek = 0;
  const KANA = 'アイウエオカキクケコサシスセソ'.split('');
  const addExp = (part, group, text, hint, fact) => { const e = { label: `（${KANA[ek++]}）`, part, group, text, figs: [], hint, fact }; experiments.push(e); return e; };

  // (a) の本文
  const rulesA = RULES_A(solvent);
  const partA = {
    label: '(a)', theme: 'acetal_exchange', solventC: solvent.length,
    intro: 'アルデヒドとアルコールの反応について、次の規則 R1〜R4 が成り立つ。',
    rules: rulesA, limitation: LIMIT_A,
    figs: applies.map((a, i) => ({ label: labelsP[i], smiles: a.reactant })),
    tail: `図の化合物 ${labelsP.slice(0, applies.length).join('・')} を、それぞれ少量の硫酸を含む大過剰の${ALK[solvent]}の中に十分な時間置いた。`,
    expLead: '次に、構造のわからない化合物 D について、次の実験を行った。',
  };
  applies.forEach((a, i) => {
    questions.push({
      id: `q1${labelsP[i]}`, part: '(a)', kind: 'apply', type: 'draw', points: 4, depends_on: [],
      prompt: `化合物 ${labelsP[i]} から得られる主生成物の構造式を書け。`,
      answer: { smiles: [a.product] },
      notes: ['立体異性体は区別しない', `溶媒の${ALK[solvent]}・水・OH を1つだけもつアルコール（放出されたもの）は書かなくてよい`],
      explain: a.case === 'ring' || a.case === 'hemi'
        ? `ヒドロキシ基とアルデヒドの炭素で ${a.ring} 員環がつくれるので R4 が優先し、環状アセタールになる。環の外の –OR は R3 で溶媒の –O${solvent === 'C' ? 'CH₃' : 'C₂H₅'} に置き換わる。`
        : a.case === 'trap' ? 'OH はあるが、アルデヒドの炭素との間で五員環・六員環をつくれない（四員環や七員環以上になる）ので R4 は起こらない。R3 だけが起こり、OH はそのまま残る。'
          : 'OH がないので R3 だけが起こり、–OR はすべて溶媒のアルコールの –OR に置き換わる。',
    });
  });

  // (a) 逆算: D
  const Dsub = '化合物 D';
  const pat = Dq.pat;
  const lbls = Object.keys(pat.counts).sort();
  const cntTxt = lbls.map((L) => `${L} ${pat.counts[L]} mol`).join('、').replace(/、([^、]*)$/, ' と $1');
  const premD = `${Dsub} は分子式 ${subf(Dq.formula)} で、${Dq.arom ? 'ベンゼン環を1つもち、それ以外に環や多重結合をもたない' : '環構造や多重結合をもたない'}化合物である。`;
  const eD0 = addExp('(a)', 'D', `${premD}${Dsub} は金属ナトリウムと反応しなかった。`, 'OH も COOH ももたない。O が2つあるのにエーテルやアセタール', { subject: 'D', card: 'sodium', value: false, premise: { formula: Dq.formula, arom: Dq.arom } });
  const eD1 = addExp('(a)', 'D', `${Dsub} を希硫酸と加熱すると R2 の逆反応が完全に進み、D 1 mol から ${cntTxt} が得られた（ほかの生成物は水だけ）。`, 'アセタールの加水分解。生成物の物質量の比から、2つの –OR が同じか違うかがわかる', { subject: 'D', card: 'hydrolysis_pattern', value: pat.counts });
  const propExps = pat.props.map((p) => addExp('(a)', 'D', propSentence(p.card, `化合物 ${p.label}`, p.value), propHint(p.card, p.value), { subject: p.label, card: p.card, value: p.value }));
  questions.push({
    id: 'q2D', part: '(a)', kind: 'identify', type: 'draw', points: 6, depends_on: [],
    prompt: '化合物 D の構造式を書け。', answer: { smiles: [Dq.D] }, notes: ['立体異性体は区別しない'],
    refs: [eD0, eD1, ...propExps].map((e) => e.label),
    explain: `分子式 ${subf(Dq.formula)} の構造異性体（${Dq.pool.length} 種類）のうち、Na と反応せず加水分解で ${cntTxt} を与えるアセタールは ${Dq.nBase} 種類。${pat.props.map((p) => `${p.label}: ${propHint(p.card, p.value)}`).join('。')}。これをすべて満たすのは1つだけ。`,
    unique: { universe: Dq.universe, formula: Dq.formula, count: 1 },
  });
  questions.push({
    id: 'q2E', part: '(a)', kind: 'identify', type: 'draw', points: 4, depends_on: ['q2D'],
    prompt: '化合物 E の構造式を書け。', answer: { smiles: [Dq.carbonyl] }, notes: ['立体異性体は区別しない'],
    explain: 'D のアセタール炭素が –CHO に戻ったものが E。',
  });
  let calcA = null;
  if (knobs.calcA) {
    // 物質量の係数だけで記号が決まる生成物を問う（係数が同じ生成物が2つあれば E を問う）
    const Fl = lbls.find((L) => L !== 'E' && lbls.filter((M) => pat.counts[M] === pat.counts[L]).length === 1) || 'E';
    const Fs = Object.keys(Dq.labelOf).find((s) => Dq.labelOf[s] === Fl);
    const MD = molarMass(Dq.formula);
    const MF = molarMass(chem.formula(chem.graphFromSmiles(RDKit, Fs)));
    const w = sig(2 + r() * 18, 3);
    const n = pat.counts[Fl];
    const ans = sig((w / MD) * n * MF, 3);
    calcA = { id: 'q2m', part: '(a)', kind: 'calc', type: 'number', points: 4, depends_on: ['q2D'], unit: 'g', sig_figs: 3,
      prompt: `D ${fmt(w, 3)} g を完全に加水分解すると、化合物 ${Fl} は何 g 得られるか。有効数字 3 けたで答えよ。`,
      answer: { value: ans, sig_figs: 3 }, notes: ['有効数字 3 けた'],
      calc: { kind: 'hydrolysis_mass', w, formulaX: Dq.formula, formulaT: chem.formula(chem.graphFromSmiles(RDKit, Fs)), n, answer: ans, sig: 3 },
      explain: `D の分子量 ${mStr(MD)}、${Fl} の分子量 ${mStr(MF)}。D 1 mol から ${Fl} は ${n} mol。${fmt(w, 3)} ÷ ${mStr(MD)} × ${n} × ${mStr(MF)} = ${fmt(ans, 3)} g。` };
    questions.push(calcA);
  }

  // ---- (b) 光学分割 ----
  const MG = molarMass(Gq.formula);
  const fH = chem.formula(chem.graphFromSmiles(RDKit, Gq.H));
  const MH = molarMass(fH);
  const partB = {
    label: '(b)', theme: 'optical_resolution',
    premise: { formula: Gq.formula, arom: Gq.arom, chiral: 1 },
    intro: `化合物 G（分子式 ${subf(Gq.formula)}）は${Gq.arom ? 'ベンゼン環を1つもち、それ以外に環や多重結合をもたない化合物で' : '環構造や多重結合をもたない化合物で'}、不斉炭素原子を1個もち、1対の鏡像異性体 (+)-G と (−)-G がある。G は触媒の存在下で無水酢酸と反応して、酢酸エステル H になる（式1）。はじめに用いる G は (+)-G と (−)-G を同じ物質量ずつ含む混合物（ラセミ体）である。`,
    equation: { lhs: ['G', '(CH₃CO)₂O'], rhs: ['H', 'CH₃COOH'], masses: { G: MG, H: MH } },
    rules: [
      { id: 'X', title: '触媒 X', text: '触媒 X の存在下では、(+)-G だけが式1 の反応をする。(−)-G はまったく反応しない。' },
      { id: 'Y', title: '触媒 Y', text: '触媒 Y の存在下では、(+)-G と (−)-G の両方が式1 の反応をするが、(+)-G のほうが速く反応する。' },
    ],
    expLead: '化合物 G について、次のことがわかっている。',
    limitation: '触媒がなければ式1 の反応は起こらない。無水酢酸は十分な量を用いる。H の加水分解、G と H のラセミ化（鏡像異性体どうしの変化）は起こらない。',
  };
  const eG = Gq.facts.map((f) => addExp('(b)', 'G', propSentence(f.card, '化合物 G', f.value), propHint(f.card, f.value), { subject: 'G', card: f.card, value: f.value }));
  questions.push({
    id: 'q3G', part: '(b)', kind: 'identify', type: 'draw', points: 6, depends_on: [],
    prompt: '化合物 G の構造式を書け。', answer: { smiles: [Gq.G] }, notes: ['立体異性体（鏡像異性体）は区別しない'],
    refs: eG.map((e) => e.label),
    explain: `分子式 ${subf(Gq.formula)} の構造異性体（${Gq.pool.length} 種類）のうち、不斉炭素原子を1個もつものは ${Gq.nBase} 種類。${Gq.facts.map((f) => propHint(f.card, f.value)).join('。')}。すべて満たすのは1つだけ。`,
    unique: { universe: Gq.universe, formula: Gq.formula, count: 1 },
  });
  // X: 一方だけ反応 → H の質量
  const wX = sig(3 + r() * 27, 3);
  const ansX = sig((wX / MG / 2) * MH, 3);
  questions.push({
    id: 'q4X', part: '(b)', kind: 'calc', type: 'number', points: 5, depends_on: [], unit: 'g', sig_figs: 3,
    prompt: `ラセミ体の G ${fmt(wX, 3)} g を、触媒 X の存在下で十分な量の無水酢酸と反応させた。反応が終わったとき、得られた H は何 g か。有効数字 3 けたで答えよ。`,
    answer: { value: ansX, sig_figs: 3 }, notes: ['有効数字 3 けた'],
    calc: { kind: 'resolution_X', w: wX, formulaG: Gq.formula, formulaH: fH, answer: ansX, sig: 3 },
    explain: `G の物質量は ${fmt(wX, 3)} ÷ ${mStr(MG)} mol。反応するのは (+)-G（半分）だけなので H はその半分の物質量。${fmt(wX, 3)} ÷ ${mStr(MG)} ÷ 2 × ${mStr(MH)} = ${fmt(ansX, 3)} g。`,
  });
  let nCalc = 1 + (calcA ? 1 : 0);
  if (knobs.Y) {
    // Y: 反応した G のうち (+) の割合 p/(p+q)、転化率 c。表示する数値（有効数字3けた）から答えを計算し直す
    let wY, b, p, q, c;
    for (let g2 = 0; g2 < 50; g2++) {
      [p, q] = pick(r, [[9, 1], [4, 1], [7, 3], [17, 3], [3, 1], [19, 1]]);
      c = sig(0.2 + r() * 0.35, 2);
      if (c * p / (p + q) < 0.48) break;
    }
    wY = sig(3 + r() * 27, 3);
    b = sig((wY / MG) * c * MH, 3);
    const n0 = wY / MG, nH = b / MH;
    const remPlus = n0 / 2 - nH * p / (p + q);
    const remMinus = n0 / 2 - nH * q / (p + q);
    if (remPlus <= 0 || remMinus <= 0) return null;
    const frac = sig((100 * remMinus) / (remPlus + remMinus), 2);
    const yIntro = `ラセミ体の G ${fmt(wY, 3)} g を、触媒 Y の存在下で無水酢酸と反応させ、途中で反応を止めた。反応後の混合物から H ${fmt(b, 3)} g と、反応しなかった G を分離した。分離した H をすべて加水分解して得た G は、(+)-G と (−)-G を物質量の比 ${p} : ${q} で含んでいた。`;
    partB.tailY = yIntro;
    if (knobs.remain) {
      const mRem = sig((n0 - nH) * MG, 3);
      questions.push({
        id: 'q5m', part: '(b)', kind: 'calc', type: 'number', points: 4, depends_on: [], unit: 'g', sig_figs: 3,
        prompt: '触媒 Y を用いた実験で、反応しなかった G は何 g か。有効数字 3 けたで答えよ。',
        answer: { value: mRem, sig_figs: 3 }, notes: ['有効数字 3 けた'],
        calc: { kind: 'resolution_Y_remaining', w: wY, b, formulaG: Gq.formula, formulaH: fH, answer: mRem, sig: 3 },
        explain: `はじめの G: ${fmt(wY, 3)} ÷ ${mStr(MG)} mol。反応した G = 生じた H = ${fmt(b, 3)} ÷ ${mStr(MH)} mol。差に ${mStr(MG)} を掛けて ${fmt(mRem, 3)} g。`,
      });
      nCalc++;
    }
    questions.push({
      id: 'q5Y', part: '(b)', kind: 'calc', type: 'number', points: 6, depends_on: knobs.remain ? ['q5m'] : [], unit: '%', sig_figs: 2,
      prompt: '触媒 Y を用いた実験で、反応しなかった G に含まれる (−)-G の割合は、物質量で何 % か。有効数字 2 けたで答えよ。',
      answer: { value: frac, sig_figs: 2 }, notes: ['有効数字 2 けた', '物質量の割合（%）で答える'],
      calc: { kind: 'resolution_Y_fraction', w: wY, b, p, q, formulaG: Gq.formula, formulaH: fH, answer: frac, sig: 2 },
      explain: `はじめの (+)-G・(−)-G はそれぞれ ${fmt(wY, 3)} ÷ ${mStr(MG)} ÷ 2 mol。反応した G（= H）${fmt(b, 3)} ÷ ${mStr(MH)} mol のうち (+) が ${p}/${p + q}、(−) が ${q}/${p + q}。残った (−)-G ÷ 残った G 全体 × 100 = ${fmt(frac, 2)} %。`,
    });
    nCalc++;
  }
  if (knobs.essay) {
    questions.push({
      id: 'q6', part: '(b)', kind: 'essay', type: 'essay', points: 4, depends_on: [],
      prompt: '触媒 Y を用いて反応を長く続けると、残る G に含まれる (−)-G の割合は高くなるが、得られる (−)-G の量は少なくなる。その理由を 50 字程度で説明せよ。',
      model: '触媒 Y は (−)-G とも反応するので、反応を続けると (+)-G が減るだけでなく (−)-G も H に変わって減っていくから。',
      rubric: ['触媒 Y は (−)-G とも反応する（両方の鏡像異性体が反応する）', '(+)-G の方が速く減るので割合は上がるが、(−)-G も H に変わって量が減る'],
      notes: ['50 字程度'],
      explain: '採点基準の要素がすべて入っていれば満点。',
    });
  }

  // 並べる: (a) 適用 → (a) 逆算 → (b)
  questions.forEach((x, i) => { x.no = i + 1; });
  const params = {
    seed, difficulty: level, theme: ['acetal_exchange', 'optical_resolution'], solvent,
    rule_applications: applies.length, rule_priorities: applies.filter((a) => a.ring).length + (applies.some((a) => a.case === 'trap') ? 1 : 0),
    inverse_steps: { D: pat.props.length + 2, G: Gq.facts.length + 1 }, stereo: !!knobs.Y, calc_steps: nCalc,
    apply_cases: applies.map((a) => a.case),
  };
  return {
    id: `k26-${seed}-${level}`, mode: 'k26', level, title: `京大2026型 ・ アセタールの規則と光学分割（難易度 ${level}）`,
    theme: 'acetal_exchange+optical_resolution', params, difficulty: level,
    parts: [partA, partB], experiments, questions,
    // 規則の一覧（データスキーマ: id・text・figure_smiles）
    rules: [...rulesA.map((x) => ({ ...x, part: '(a)' })), ...partB.rules.map((x) => ({ ...x, part: '(b)' }))]
      .map((x) => ({ ...x, figure_smiles: (x.scheme || []).filter((t) => !['+', '⇄', '→'].includes(t)) })),
    hidden: { D: Dq.D, G: Gq.G, H: Gq.H, applies },
  };
}

module.exports = { buildK26, molarMass, LEVELS, applyAcetalRules, hydrolysisPattern, hydroMatches, propValue, rng };
