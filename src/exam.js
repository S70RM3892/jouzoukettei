'use strict';
// 大問（断片・手がかり・組み立て・計算）を、京大の大問と同じ形に組み直す。
//  問題文: 前置き → X についての実験・測定 → 化合物 A, B, … についての実験（あ）（い）… → 部分加水分解
//  設問: 分子式（記述）、数値（記述）、「〜だけから考えられるものをすべて選べ」、構造（選択）、
//        生成物の予測、条件を満たす異性体の数、反応名、X の構造
// 答えはすべてエンジンで計算し、選択肢の構造は候補集合から作る。
const chem = require('./chem');
const { sentence, infer, NAMES } = require('./inference');

const KANA = 'あいうえおかきくけこさしすせそたちつてとなにぬねのはひふへほまみむめも'.split('');
const subf = (f) => f.replace(/\d/g, (d) => '₀₁₂₃₄₅₆₇₈₉'[d]);

function seedOf(id) {
  let h = 2166136261;
  for (const c of id) { h ^= c.charCodeAt(0); h = Math.imul(h, 16777619) >>> 0; }
  return h || 1;
}
function rng(seed) {
  let s = seed >>> 0 || 1;
  return () => { s ^= s << 13; s >>>= 0; s ^= s >> 17; s ^= s << 5; s >>>= 0; return s / 4294967296; };
}

// 異性体の数を問うときの条件（真偽で分けられるカード）
const COND = {
  iodoform: 'ヨードホルム反応を示す', silver_mirror: '銀鏡反応を示す', fecl3: '塩化鉄(III) 水溶液で呈色する',
  bromine: '臭素水を脱色する', cis_trans: 'シス-トランス異性体が存在する', chiral: '不斉炭素原子をもつ', anhydride: '加熱すると酸無水物になる',
};
// 生成物の予測に使うカード
const PREDICT = ['mild_oxidation', 'kmno4', 'dehydration', 'hydrogenation', 'bromine_water', 'nitration', 'acetylation', 'markovnikov', 'ozonolysis', 'kmno4_cleave', 'acetylation_primary', 'bromine_addition'];

// 断片の SMILES の並び（加水分解の生成物、重複あり）から、ありうる X をすべて作る
function assemblies(RDKit, list) {
  const G = require('./generator');
  const has = (s, sm) => { const m = RDKit.get_mol(s); const qm = RDKit.get_qmol(sm); try { return JSON.parse(m.get_substruct_match(qm)).atoms !== undefined; } finally { m.delete(); qm.delete(); } };
  const clsOf = (s) => {
    const cooh = has(s, '[CX3](=O)[OX2H1]');
    if (cooh) {
      if (has(s, '[CX4][OX2H1]')) return 'hydroxyacid';
      if (has(s, 'c[NX3H2]')) return 'aminoaracid';
      if (has(s, 'c[OX2H1]')) return 'hydroxyaracid';
      return 'acid';
    }
    if (has(s, '[CX3;!$(C(=O)[O,N])]=O') && !has(s, '[OX2H1,NX3;H1,H2]')) return 'carbonyl';
    return 'alcohol';
  };
  const infos = list.map((s) => G.fragmentInfo(RDKit, s, clsOf(s)));
  return G.allAssembliesWithVariants(RDKit, infos);
}

function norm(RDKit, card, v) { return JSON.stringify(chem.normalizeResult(RDKit, card, v)); }

function buildExam(RDKit, P) {
  const r = rng(seedOf(P.id));
  const cards = chem.CARDS;
  const ev = (card, s) => chem.evaluate(RDKit, card, s);
  const can = (s) => chem.canonical(RDKit, s);
  const fml = (s) => chem.formula(chem.graphFromSmiles(RDKit, s));
  const X = can(P.answer);
  const mols = new Set([X]);
  const experiments = [];
  const questions = [];
  let k = 0;
  const addExp = (group, text, figs = [], hint = '', data = {}) => {
    const e = { label: `（${KANA[k++]}）`, group, text, figs, hint, ...data };
    figs.forEach((f) => mols.add(f));
    experiments.push(e);
    return e;
  };
  const q = (x) => { questions.push({ no: questions.length + 1, ...x }); return x; };

  const calcs = P.calcs || [];
  const cx = calcs.find((c) => c.key === 'combustion_x');
  const h2 = calcs.find((c) => c.key === 'h2');
  const f3 = (x) => Number(x).toPrecision(3);

  // ---- X についての測定 ----
  if (cx) {
    const e = addExp('X', `X ${f3(cx.data.sample)} mg を完全燃焼させると、二酸化炭素 ${f3(cx.data.co2)} mg と水 ${f3(cx.data.h2o)} mg が生じた。X の分子量は ${cx.data.M} である。`, [],
      'C の質量 = CO₂ × 12/44、H の質量 = H₂O × 2/18、O は残り。物質量の比から組成式、分子量から分子式');
    q({ type: 'formula', prompt: `化合物 X の分子式を求めよ。`, answer: P.formula, points: 6, refs: [e.label], explain: cx.explain });
  }
  if (h2) {
    const e = addExp('X', `X ${f3(h2.data.m)} g に白金触媒を用いて水素を付加させると、標準状態で ${f3(h2.data.V)} L の水素が消費された。ベンゼン環には水素は付加しないものとする。`, [],
      'X の物質量（質量 ÷ 分子量）と H₂ の物質量（体積 ÷ 22.4）の比 = X 1分子の C=C の数');
    q({ type: 'number', prompt: 'X 1分子に含まれる炭素原子間の二重結合の数を答えよ。', answer: h2.answer, unit: '個', points: 5, refs: [e.label], explain: h2.explain });
  }

  // 加水分解の生成物の数（原子の収支）
  const prods = chem.hydrolyze(RDKit, chem.graphFromSmiles(RDKit, X));
  const labels = P.fragments.map((f) => f.label);
  q({ type: 'number', prompt: `X 1 mol を完全に加水分解したとき、消費される水は何 mol か。`, answer: prods.length - 1, unit: 'mol', points: 4, refs: [],
    explain: `X ＋ n H₂O → ${labels.join(' ＋ ')}（同じ化合物が複数あればその数も）。H と O の数の収支から n = ${prods.length - 1}。エステル結合・アミド結合 1 つにつき H₂O 1 分子。` });

  // ---- 断片ごとの実験 ----
  const fragQs = [];
  const combustionFrag = calcs.find((c) => c.key === 'combustion');
  for (const f of P.fragments) {
    const subj = `化合物 ${f.label}`;
    const formula = fml(f.answer);
    mols.add(can(f.answer));
    if (f.given) {
      addExp(f.label, `${subj} は${(f.note || '').replace(/（.*）$/, '') || '構造が与えられた化合物'}であった。`, [can(f.answer)]);
      continue;
    }
    const cands = f.candidates.map(can);
    cands.forEach((c) => mols.add(c));
    const isCombustion = combustionFrag && combustionFrag.answer === formula && !fragQs.some((x) => x.formula === formula);
    if (isCombustion) {
      const e = addExp(f.label, `${subj} ${f3(combustionFrag.data.sample)} mg を完全燃焼させると、二酸化炭素 ${f3(combustionFrag.data.co2)} mg と水 ${f3(combustionFrag.data.h2o)} mg が生じた。${f.label} の分子量は ${combustionFrag.data.M} である。`, [],
        'C・H・O の質量 → 物質量の比 → 組成式 → 分子量で分子式');
      q({ type: 'formula', prompt: `化合物 ${f.label} の分子式を求めよ。`, answer: formula, points: 5, refs: [e.label], explain: combustionFrag.explain });
    } else {
      addExp(f.label, `${subj} の分子式は ${subf(formula)} であった。${f.kind ? `${subj} は${f.kind.replace(/（.*）/, '')}である。` : ''}`);
    }
    const expOfClue = f.clues.map((c) => {
      const res = chem.normalizeResult(RDKit, c.card, c.result);
      const figs = Array.isArray(res) ? res : cards[c.card].kind === 'contains' ? [res] : [];
      return addExp(f.label, sentence(c.card, subj, res, cards), figs, infer(c.card, res, cards), { card: c.card });
    });
    fragQs.push({ f, formula, cands, expOfClue });
  }

  // ---- 組み立て ----
  let asmExp = [];
  if (P.assemble) {
    asmExp = P.assemble.clues.map((c) => {
      const res = chem.normalizeResult(RDKit, c.card, c.result);
      const figs = Array.isArray(res) ? res : c.card === 'partial_hydrolysis' ? [res] : [];
      return addExp('X', sentence(c.card, '化合物 X', res, cards), figs, infer(c.card, res, cards), { card: c.card });
    });
  }

  // ---- 断片の設問 ----
  let predicted = 0;
  let counted = 0;
  for (const { f, formula, cands, expOfClue } of fragQs) {
    const ans = can(f.answer);
    const table = f.clues.map((c) => cands.map((s) => { try { return norm(RDKit, c.card, ev(c.card, s)); } catch (e) { return 'x'; } }));
    const expected = f.clues.map((c) => norm(RDKit, c.card, c.result));
    const ok = (ci, i) => (cards[f.clues[ci].card].kind === 'contains' ? table[ci][i].includes(JSON.parse(expected[ci])) : table[ci][i] === expected[ci]);
    const remainAfter = (use) => cands.filter((_, i) => use.every((ci) => ok(ci, i)));
    // 「〜だけから考えられるものをすべて選べ」: 1つ抜くと 2〜6 個残る組を選ぶ
    if (f.clues.length >= 2) {
      const opts = f.clues.map((_, drop) => {
        const use = f.clues.map((__, i) => i).filter((i) => i !== drop);
        return { use, rem: remainAfter(use) };
      });
      const small = opts.filter((o) => o.rem.length >= 2 && o.rem.length <= 6);
      const large = opts.filter((o) => o.rem.length > 6 && o.rem.length < cands.length);
      if (!small.length && large.length) {
        // 候補が多いときは数で問う: ある実験がなかったら何種類残るか（残りの条件を全部当てはめて数える）
        const o = large[Math.floor(r() * large.length)];
        const dropped = f.clues.map((_, i) => i).find((i) => !o.use.includes(i));
        q({
          type: 'number', prompt: `実験${expOfClue[dropped].label}の結果がなかったとすると、化合物 ${f.label} の構造として考えられるものは何種類になるか。`,
          answer: o.rem.length, unit: '種類', points: 6, refs: o.use.map((i) => expOfClue[i].label),
          explain: o.use.map((i) => `${expOfClue[i].label} ${expOfClue[i].hint}`).join(' ／ ') + ` → この条件だけでは ${o.rem.length} 種類が残る。${expOfClue[dropped].label} がそれを1つに絞る。`,
        });
      }
      if (small.length) {
        const opts = small;
        const o = opts[Math.floor(r() * opts.length)];
        q({
          type: 'multi', prompt: `実験${o.use.map((i) => expOfClue[i].label).join('')}の結果だけから、化合物 ${f.label} の構造として考えられるものをすべて選べ。`,
          options: cands, answer: o.rem, points: 8, refs: o.use.map((i) => expOfClue[i].label),
          explain: o.use.map((i) => `${expOfClue[i].label} ${expOfClue[i].hint}`).join(' ／ ') + ` → 条件を満たすのは ${o.rem.length} 個。`,
        });
      }
    }
    q({
      type: 'single', prompt: `化合物 ${f.label} の構造を選べ。`, options: cands, answer: [ans], points: 8, refs: expOfClue.map((e) => e.label),
      explain: expOfClue.map((e) => `${e.label} ${e.hint}`).join(' ／ ') + ` → 候補 ${cands.length} 個のうち、すべてを満たすのは1つ。`,
    });
    // 生成物の予測（問題全体で1問）
    if (!predicted) {
      const pc = PREDICT.filter((c) => !f.clues.some((x) => x.card === c)).map((c) => {
        let a;
        try { a = ev(c, ans); } catch (e) { return null; }
        if (!Array.isArray(a) || a.length !== 1) return null;
        const others = new Set();
        cands.forEach((s) => { try { const v = ev(c, s); if (Array.isArray(v) && v.length === 1) others.add(can(v[0])); } catch (e) { /* 使えない */ } });
        others.delete(can(a[0]));
        return others.size >= 2 ? { c, a: can(a[0]), others: [...others] } : null;
      }).filter(Boolean);
      if (pc.length) {
        const p = pc[Math.floor(r() * pc.length)];
        const opts = [p.a, ...p.others.sort(() => r() - 0.5).slice(0, 4)].sort(() => r() - 0.5);
        opts.forEach((s) => mols.add(s));
        const verb = sentence(p.c, `化合物 ${f.label}`, [p.a], cards).replace(/、次の化合物が得られた。$/, '');
        q({ type: 'single', prompt: `${verb}、どの化合物が得られるか。`, options: opts, answer: [p.a], points: 6, refs: [],
          explain: `${f.label} は決まった構造。${infer(p.c, [p.a], cards)}` });
        predicted++;
      }
    }
    // 条件を満たす異性体の数（問題全体で1問）
    if (!counted && cands.length >= 4) {
      const cc = Object.keys(COND).filter((c) => !f.clues.some((x) => x.card === c)).map((c) => {
        let n = 0;
        try { cands.forEach((s) => { const v = ev(c, s); if (c === 'chiral' ? v > 0 : v === true) n++; }); } catch (e) { return null; }
        return n >= 1 && n < cands.length ? { c, n } : null;
      }).filter(Boolean);
      if (cc.length) {
        const x = cc[Math.floor(r() * cc.length)];
        q({ type: 'number', prompt: `化合物 ${f.label} と同じ分子式 ${subf(formula)} をもつ${f.kind ? f.kind.replace(/（.*）/, '') : '化合物'}（候補の図の ${cands.length} 種類）のうち、${COND[x.c]}ものは何種類あるか。`, answer: x.n, unit: '種類', points: 5, refs: [],
          explain: `候補の構造を1つずつ調べる。${infer(x.c, true, cards)}` });
        counted++;
      }
    }
  }

  // ---- つなぎ方 ----
  const nx = calcs.find((c) => c.key === 'n_x');
  if (nx) q({ type: 'number', prompt: nx.prompt.replace(/。$/, '') + '。', answer: nx.answer, unit: '種類', points: 5, refs: [], explain: nx.explain });
  // X の構造: つなぎ方の候補に、断片を「惜しい別の候補」に替えて組み直したものを加える
  {
    const alts = P.assemble ? P.assemble.candidates.map(can) : [X];
    const opts = new Set(alts);
    const wrong = [];
    for (const { f, cands } of fragQs) {
      const ans = can(f.answer);
      cands.filter((c) => c !== ans).slice(0, 3).forEach((w) => {
        try {
          const list = prods.map((p) => (p === ans ? w : p));
          assemblies(RDKit, list).forEach((x) => { if (!opts.has(x)) wrong.push(x); });
        } catch (e) { /* つなげない */ }
      });
    }
    wrong.sort(() => r() - 0.5).slice(0, Math.max(0, 5 - opts.size)).forEach((x) => opts.add(x));
    const options = [...opts].sort(() => r() - 0.5);
    options.forEach((s) => mols.add(s));
    q({ type: 'single', prompt: '化合物 X の構造を選べ。', options, answer: [X], points: 10, refs: asmExp.map((e) => e.label),
      explain: (asmExp.length ? asmExp.map((e) => `${e.label} ${e.hint}`).join(' ／ ') + ' → つなぎ方が1つに決まる。' : '断片 ' + labels.join('・') + ' をエステル（アミド）結合でつなぐと1通りに決まる。') + ' ほかの選択肢は、どれかの断片が違うか、つなぎ方が違う。' });
  }

  // ---- 反応名（1問） ----
  const named = experiments.filter((e) => e.card && NAMES[e.card]);
  if (named.length) {
    const e = named[Math.floor(r() * named.length)];
    const wrong = Object.values(NAMES).filter((n) => n !== NAMES[e.card]).sort(() => r() - 0.5).slice(0, 3);
    q({ type: 'choice', prompt: `実験${e.label}で起こった反応（または操作）の名称を選べ。`, options: [NAMES[e.card], ...wrong].sort(() => r() - 0.5), answer: [NAMES[e.card]], points: 2, refs: [e.label], explain: `${NAMES[e.card]}。${e.hint}` });
  }

  return {
    intro: P.story,
    experiments,
    questions,
    total: questions.reduce((a, x) => a + x.points, 0),
    mols: [...mols],
  };
}

module.exports = { buildExam };
