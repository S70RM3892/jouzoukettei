'use strict';
// 大問（断片・誘導体・化合物どうしの関係・組み立て・計算）を、京大の大問と同じ形に組み直す。
//  問題文: 規則の段落 → 前置き → X の測定 → 化合物 A, B, … の実験 → 誘導体 D, E … → 化合物どうしの関係 → 部分加水分解
//  設問: X の分子式、C=C の数、抽出、断片の分子式、「〜だけから考えられるものをすべて選べ」／何通りか、断片と誘導体の構造、
//        条件を満たす異性体の数、X の不斉炭素、質量の計算、つなぎ方の数、X の構造、理由の記述（自己採点）
// 答えはすべてエンジンで計算し、選択肢の構造は候補集合から作る。
const chem = require('./chem');
const { sentence, infer, derivation, relation, relationHint, RULE_CARDS } = require('./inference');

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
  const derived = P.derived || [];
  const relations = P.relations || [];
  const clueList = (cl) => (cl || []).map((c) => ({ ...c, result: chem.normalizeResult(RDKit, c.card, c.result) }));

  // ---- 規則の段落（教科書にない反応・京大で問題文に与えられた規則） ----
  const used = new Set([
    ...P.fragments.flatMap((f) => (f.clues || []).map((c) => c.card)),
    ...derived.flatMap((d) => [d.op, ...d.clues.map((c) => c.card)]),
    ...relations.map((rl) => rl.op),
    ...(P.assemble ? P.assemble.clues.map((c) => c.card) : []),
  ]);
  const rules = RULE_CARDS.filter((c) => used.has(c)).map((c) => ({ name: cards[c].name, text: cards[c].action }));

  // ---- X の測定 ----
  if (cx) {
    const e = addExp('X', `X ${f3(cx.data.sample)} mg を完全燃焼させると、二酸化炭素 ${f3(cx.data.co2)} mg と水 ${f3(cx.data.h2o)} mg が生じた。X の分子量は ${cx.data.M} である。`, [],
      'C の質量 = CO₂ × 12/44、H の質量 = H₂O × 2/18、O は残り。物質量の比から組成式、分子量から分子式');
    q({ type: 'formula', prompt: '化合物 X の分子式を求めよ。', answer: P.formula, points: 6, refs: [e.label], explain: cx.explain });
  }
  if (h2) {
    const e = addExp('X', `X ${f3(h2.data.m)} g に白金触媒を用いて水素を付加させると、標準状態で ${f3(h2.data.V)} L の水素が消費された。ベンゼン環には水素は付加しないものとする。`, [],
      'X の物質量（質量 ÷ 分子量）と H₂ の物質量（体積 ÷ 22.4）の比 = X 1分子の C=C の数');
    q({ type: 'number', prompt: 'X 1分子に含まれる炭素原子間の二重結合の数を答えよ。', answer: h2.answer, unit: '個', points: 4, refs: [e.label], explain: h2.explain });
  }
  const prods = chem.hydrolyze(RDKit, chem.graphFromSmiles(RDKit, X));
  const labels = P.fragments.map((f) => f.label);

  // ---- 断片の事実と直接の実験 ----
  const combustionFrag = calcs.find((c) => c.key === 'combustion');
  const fragInfo = [];
  let combustionDone = false;
  for (const f of P.fragments) {
    const subj = `化合物 ${f.label}`;
    const formula = fml(f.answer);
    mols.add(can(f.answer));
    if (f.given) {
      addExp(f.label, `${subj} は${(f.note || '').replace(/（.*）$/, '') || '構造が与えられた化合物'}であった。`, [can(f.answer)]);
      fragInfo.push({ f, formula, cands: [can(f.answer)], exps: [] });
      continue;
    }
    const cands = f.candidates.map(can);
    cands.forEach((c) => mols.add(c));
    let formulaQ = null;
    if (combustionFrag && !combustionDone && combustionFrag.answer === formula) {
      combustionDone = true;
      const e = addExp(f.label, `${subj} ${f3(combustionFrag.data.sample)} mg を完全燃焼させると、二酸化炭素 ${f3(combustionFrag.data.co2)} mg と水 ${f3(combustionFrag.data.h2o)} mg が生じた。${f.label} の分子量は ${combustionFrag.data.M} である。`, [],
        'C・H・O の質量 → 物質量の比 → 組成式 → 分子量で分子式');
      formulaQ = { type: 'formula', prompt: `化合物 ${f.label} の分子式を求めよ。`, answer: formula, points: 5, refs: [e.label], explain: combustionFrag.explain };
    } else {
      // 分類（アルコールかどうかなど）は書かない。エステルの加水分解生成物であることと実験から自分で判断させる
      addExp(f.label, `${subj} の分子式は ${subf(formula)} であった。`);
    }
    const exps = clueList(f.clues).map((c) => {
      const figs = Array.isArray(c.result) ? c.result : cards[c.card].kind === 'contains' ? [c.result] : [];
      return addExp(f.label, sentence(c.card, subj, c.result, cards), figs, infer(c.card, c.result, cards), { card: c.card });
    });
    fragInfo.push({ f, formula, cands, exps, formulaQ });
  }

  // ---- 誘導体（A を〜すると D。D の実験） ----
  const derInfo = [];
  for (const d of derived) {
    const fi = fragInfo.find((x) => x.f.label === d.from);
    const e0 = addExp(d.label, derivation(d.op, `化合物 ${d.from}`, d.label, cards), [], infer(d.op, [], cards), { card: d.op });
    fi.exps.push(e0);
    const exps = clueList(d.clues).map((c) => {
      const figs = Array.isArray(c.result) ? c.result : [];
      return addExp(d.label, sentence(c.card, `化合物 ${d.label}`, c.result, cards), figs, infer(c.card, c.result, cards), { card: c.card });
    });
    fi.exps.push(...exps);
    derInfo.push({ d, e0, exps });
    d.candidates.forEach((c) => mols.add(can(c)));
  }
  // ---- 化合物どうしの関係 ----
  const relExps = relations.map((rl) => {
    const e = addExp('rel', relation(rl, cards), [], relationHint(rl, cards), { rel: rl });
    [rl.a, rl.b].forEach((l) => { const fi = fragInfo.find((x) => x.f.label === l); if (fi) fi.exps.push(e); });
    return e;
  });
  // ---- 部分加水分解 ----
  const asmExp = P.assemble ? clueList(P.assemble.clues).map((c) => {
    const figs = Array.isArray(c.result) ? c.result : c.card === 'partial_hydrolysis' ? [c.result] : [];
    return addExp('X', sentence(c.card, '化合物 X', c.result, cards), figs, infer(c.card, c.result, cards), { card: c.card });
  }) : [];

  // ---- 抽出（京大頻出）: 加水分解後の混合物の分け方 ----
  const layerOf = (s) => {
    let a = false, b = false, c = false;
    try { a = ev('nahco3', s); b = ev('naoh', s); c = ev('hcl', s); } catch (e) { return null; }
    if (a && c) return null; // 両性（アミノ酸）は扱わない
    if (c) return 'amine';
    if (a) return 'acid';
    if (b) return 'phenol';
    return 'neutral';
  };
  const small = (s) => { const n = (fml(s).match(/^C(\d*)/) || [])[1]; return (n === '' ? 1 : +n) <= 3 && !/c/.test(s); };
  const ex = fragInfo.filter((x) => !small(x.f.answer)).map((x) => ({ label: x.f.label, layer: layerOf(can(x.f.answer)) }));
  if (ex.length >= 2 && ex.every((x) => x.layer) && new Set(ex.map((x) => x.layer)).size >= 2) {
    const hasAmine = ex.some((x) => x.layer === 'amine');
    const steps = hasAmine
      ? [['amine', '希塩酸'], ['acid', '炭酸水素ナトリウム水溶液'], ['phenol', '水酸化ナトリウム水溶液']]
      : [['acid', '炭酸水素ナトリウム水溶液'], ['phenol', '水酸化ナトリウム水溶液']];
    const names = ['Ⅰ', 'Ⅱ', 'Ⅲ', 'Ⅳ'];
    const where = (layer) => { const i = steps.findIndex((st) => st[0] === layer); return i >= 0 ? names[i] : names[steps.length]; };
    const proc = steps.map((st, i) => `${'①②③'[i]}${i ? '残ったエーテル層に' : ''}${st[1]}を加えてよく振り混ぜ、水層${names[i]}を分けた。`).join('') + `${'①②③④'[steps.length]}エーテル層${names[steps.length]}が残った。`;
    q({
      type: 'assign', prompt: `X を完全に加水分解した後の混合物をジエチルエーテルに溶かし、次の操作で分けた。${proc}化合物 ${ex.map((x) => x.label).join('・')} は、それぞれ ${names.slice(0, steps.length + 1).join('・')} のどれに含まれるか（水によく溶ける低分子の化合物はここでは考えない）。`,
      items: ex.map((x) => x.label), options: names.slice(0, steps.length + 1), answer: ex.map((x) => where(x.layer)), points: 6, refs: [],
      explain: '塩基性のアミン → 希塩酸（塩になる）。カルボン酸 → 炭酸水素ナトリウム水溶液（炭酸より強い酸）。フェノール類 → 炭酸より弱いので NaHCO₃ には溶けず、NaOH 水溶液に溶ける。中性の化合物（アルコール・エステルなど）はエーテル層に残る。',
      layers: ex,
    });
  }

  // ---- 断片ごとの設問 ----
  let counted = 0;
  const allFeatsFor = (fi) => fi.exps;
  for (const fi of fragInfo) {
    const { f, cands, formulaQ } = fi;
    if (f.given) continue;
    if (formulaQ) q(formulaQ);
    const ans = can(f.answer);
    // 直接の手がかりで「すべて選べ」（関係のない断片のとき）
    const direct = clueList(f.clues);
    const isChained = derived.some((d) => d.from === f.label) || relations.some((rl) => rl.a === f.label || rl.b === f.label);
    if (direct.length >= 1 && (isChained || direct.length >= 2)) {
      const use = isChained ? direct.map((_, i) => i) : direct.map((_, i) => i).slice(0, direct.length - 1);
      const rem = cands.filter((s) => use.every((i) => { try { return JSON.stringify(chem.normalizeResult(RDKit, direct[i].card, ev(direct[i].card, s))) === JSON.stringify(direct[i].result); } catch (e) { return false; } }));
      const refs = use.map((i) => fi.exps[i].label);
      if (rem.length >= 2 && rem.length <= 6) {
        q({ type: 'multi', prompt: `${isChained ? `化合物 ${f.label} 自身についての実験` : '実験'}${refs.join('')}の結果だけから、化合物 ${f.label} の構造として考えられるものをすべて選べ。`,
          options: cands, answer: rem, points: 8, refs, explain: use.map((i) => `${fi.exps[i].label} ${fi.exps[i].hint}`).join(' ／ ') + ` → ${rem.length} 個が残る。${isChained ? 'これをほかの化合物との関係で1つに絞る。' : ''}` });
      } else if (rem.length > 6 && rem.length < cands.length) {
        q({ type: 'number', prompt: `実験${refs.join('')}の結果だけでは、化合物 ${f.label} の構造として考えられるものは何通りあるか（構造異性体のみ）。`,
          answer: rem.length, unit: '通り', points: 6, refs, explain: use.map((i) => `${fi.exps[i].label} ${fi.exps[i].hint}`).join(' ／ ') + ` → ${rem.length} 通り。` });
      }
    }
    q({ type: 'single', prompt: `化合物 ${f.label} の構造を選べ。`, options: cands, answer: [ans], points: 8, refs: fi.exps.map((e) => e.label),
      explain: fi.exps.map((e) => `${e.label} ${e.hint}`).join(' ／ ') + ` → 候補 ${cands.length} 個のうち、すべてを満たすのは1つ。` });
    // 誘導体の構造
    for (const di of derInfo.filter((x) => x.d.from === f.label)) {
      const opts = [...new Set(di.d.candidates.map(can))];
      const pickOpts = [can(di.d.answer), ...opts.filter((o) => o !== can(di.d.answer)).sort(() => r() - 0.5).slice(0, 7)].sort(() => r() - 0.5);
      pickOpts.forEach((o) => mols.add(o));
      q({ type: 'single', prompt: `化合物 ${di.d.label} の構造を選べ。`, options: pickOpts, answer: [can(di.d.answer)], points: 6, refs: [di.e0.label, ...di.exps.map((e) => e.label)],
        explain: `${di.e0.label} ${di.e0.hint} ／ ` + di.exps.map((e) => `${e.label} ${e.hint}`).join(' ／ ') });
    }
    // 条件を満たす異性体の数（1問）
    if (!counted && cands.length >= 4) {
      const cc = Object.keys(COND).filter((c) => !(f.clues || []).some((x) => x.card === c)).map((c) => {
        let n = 0;
        try { cands.forEach((s) => { const v = ev(c, s); if (c === 'chiral' ? v > 0 : v === true) n++; }); } catch (e) { return null; }
        return n >= 1 && n < cands.length ? { c, n } : null;
      }).filter(Boolean);
      if (cc.length) {
        const x = cc[Math.floor(r() * cc.length)];
        q({ type: 'number', prompt: `化合物 ${f.label} と同じ分子式 ${subf(fi.formula)} をもつ${f.kind ? f.kind.replace(/（.*）/, '') : '化合物'}の構造異性体（${f.label} を含む）のうち、${COND[x.c]}ものは何種類あるか。`,
          answer: x.n, unit: '種類', points: 5, refs: [], explain: `すべての構造異性体を書き出して1つずつ調べる（全部で ${cands.length} 種類）。${infer(x.c, true, cards)}` });
        counted++;
      }
    }
  }

  // ---- X の立体（不斉炭素の数） ----
  const nChiral = chem.chiralCount(chem.graphFromSmiles(RDKit, X));
  q({ type: 'number', prompt: '化合物 X に含まれる不斉炭素原子の数を答えよ。', answer: nChiral, unit: '個', points: 3, refs: [],
    explain: `X の構造で、4つの異なる原子・原子団が結合した炭素を数える → ${nChiral} 個。` });

  // ---- 質量の計算 ----
  {
    const Mx = chem.formulaMass(P.formula);
    const tgt = fragInfo.filter((x) => !x.f.given)[Math.floor(r() * fragInfo.filter((x) => !x.f.given).length)] || fragInfo[0];
    const nT = prods.filter((p) => p === can(tgt.f.answer)).length;
    const Mt = chem.formulaMass(tgt.formula);
    const w = Number((5 + r() * 20).toPrecision(3));
    const ansMass = (w / Mx) * nT * Mt;
    q({ type: 'number', prompt: `X ${f3(w)} g を完全に加水分解したとき、得られる化合物 ${tgt.f.label} は何 g か。有効数字 3 けたで答えよ。`, answer: Number(ansMass.toPrecision(3)), tol: 0.006, unit: 'g', points: 5, refs: [],
      explain: `X の分子量 ${Mx}、${tgt.f.label} の分子量 ${Mt}。X 1 mol から ${tgt.f.label} は ${nT} mol。${f3(w)} ÷ ${Mx} × ${nT} × ${Mt} = ${ansMass.toPrecision(3)} g。` });
  }

  // ---- つなぎ方の数 ----
  const nx = calcs.find((c) => c.key === 'n_x');
  if (nx) q({ type: 'number', prompt: nx.prompt.replace(/。$/, '') + '。', answer: nx.answer, unit: '種類', points: 5, refs: [], explain: nx.explain });

  // ---- X の構造（断片を惜しい候補に替えた誤答つき） ----
  {
    const alts = P.assemble ? P.assemble.candidates.map(can) : [X];
    const opts = new Set(alts);
    const wrong = [];
    for (const fi of fragInfo) {
      if (fi.f.given) continue;
      const ans = can(fi.f.answer);
      fi.cands.filter((c) => c !== ans).slice(0, 3).forEach((w) => {
        try {
          const list = prods.map((p) => (p === ans ? w : p));
          assemblies(RDKit, list).forEach((x) => { if (!opts.has(x)) wrong.push(x); });
        } catch (e) { /* つなげない */ }
      });
    }
    wrong.sort(() => r() - 0.5).slice(0, Math.max(0, 6 - opts.size)).forEach((x) => opts.add(x));
    const options = [...opts].sort(() => r() - 0.5);
    options.forEach((s) => mols.add(s));
    q({ type: 'single', prompt: '化合物 X の構造を選べ。', options, answer: [X], points: 10, refs: asmExp.map((e) => e.label),
      explain: (asmExp.length ? asmExp.map((e) => `${e.label} ${e.hint}`).join(' ／ ') + ' → つなぎ方が1つに決まる。' : `断片 ${labels.join('・')} をエステル（アミド）結合でつなぐと1通りに決まる。`) + ' ほかの選択肢は、どれかの断片が違うか、つなぎ方が違う。' });
  }

  // ---- 理由の記述（自己採点。京大は字数制限つきの説明を毎年のように出す） ----
  const essays = [];
  const extr = questions.find((x) => x.type === 'assign');
  if (extr) {
    const acid = extr.layers.find((x) => x.layer === 'acid');
    const phen = extr.layers.find((x) => x.layer === 'phenol');
    if (acid && phen) essays.push({ prompt: `抽出の操作で、化合物 ${acid.label} は炭酸水素ナトリウム水溶液の層に移り、化合物 ${phen.label} は移らなかった。その理由を 50 字程度で説明せよ。`,
      model: `${acid.label} はカルボキシ基をもち炭酸より強い酸なので NaHCO₃ と反応して塩になり水に溶けるが、${phen.label} のフェノール性ヒドロキシ基は炭酸より弱い酸で反応しないから。` });
  }
  const same = relations.find((rl) => rl.type === 'same');
  if (same) essays.push({ prompt: `実験${relExps[relations.indexOf(same)].label}から、化合物 ${same.a} と ${same.b} の構造についてどのようなことがわかるか。40 字程度で説明せよ。`,
    model: `${same.a} と ${same.b} は、この反応で変化する部分以外の構造（炭素骨格と官能基の位置）が同じである。` });
  const yields = relations.find((rl) => rl.type === 'yields');
  if (yields) essays.push({ prompt: `実験${relExps[relations.indexOf(yields)].label}から、化合物 ${yields.a} と ${yields.b} の構造の関係を 40 字程度で説明せよ。`,
    model: `${yields.b} は ${yields.a} にこの反応を行った生成物なので、反応で変わらない部分（炭素骨格など）は ${yields.a} と同じである。` });
  if (fragInfo.some((x) => (x.f.clues || []).some((c) => c.card === 'anhydride' && c.result === true))) {
    const fx = fragInfo.find((x) => (x.f.clues || []).some((c) => c.card === 'anhydride' && c.result === true));
    essays.push({ prompt: `化合物 ${fx.f.label} を加熱すると酸無水物が生じるのはなぜか。30 字程度で説明せよ。`, model: '2つのカルボキシ基が近い位置にあり、分子内で脱水して環状の酸無水物をつくれるから。' });
  }
  if (essays.length) {
    const e = essays[Math.floor(r() * essays.length)];
    q({ type: 'essay', prompt: e.prompt, model: e.model, points: 5, refs: [], explain: '模範解答と見比べて自己採点する。要点が入っていれば満点、半分なら部分点。' });
  }

  return {
    rules,
    intro: P.story,
    experiments,
    questions,
    total: questions.reduce((a, x) => a + x.points, 0),
    mols: [...mols],
  };
}

module.exports = { buildExam };
