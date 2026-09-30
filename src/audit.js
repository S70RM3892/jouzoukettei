'use strict';
// 問題（画面に出す形: 規則・問題文・実験・設問）を機械的に検査する。どれか1つでも落ちた問題は収録しない。
//  1. リーク検査: 図・名前が設問の答えそのもの／答えを読み取れるものになっていないか。答えの数値がそのまま事実として書かれていないか
//  2. 依存関係検査: 設問の depends_on が実在し、循環がなく、依存先が先に置かれているか
//  3. 注記: 構造・数の設問に、採点の曖昧さを消す注記があるか
// 一意性（事実を満たす構造がちょうど1つか）は、問題の種類ごとに checkUnique* で調べる。
const chem = require('./chem');
const { leakKind } = require('./leak');

const STRUCT = new Set(['single', 'multi', 'draw']);
const strip = (s) => s.replace(/[/\\]/g, '').replace(/@+/g, '');

function answersOf(q) {
  if (!STRUCT.has(q.type)) return [];
  const a = q.type === 'draw' ? q.answer.smiles : q.answer;
  return (Array.isArray(a) ? a : [a]).filter((x) => typeof x === 'string');
}

// sheet: { rules:[{figs?}], intro, experiments:[{label,text,figs,fact?}], questions:[{id,type,kind,answer,about?,depends_on?,notes?}] }
// names: { canonical SMILES: 名前 }（図のキャプションや本文に出る名前）
function leakCheck(RDKit, sheet, names = {}) {
  const errors = [];
  const figs = [];
  (sheet.rules || []).forEach((r) => (r.figs || []).forEach((f) => figs.push({ where: `規則 ${r.id || r.name}`, s: f })));
  (sheet.experiments || []).forEach((e) => (e.figs || []).forEach((f) => figs.push({ where: `実験${e.label}`, s: f })));
  (sheet.parts || []).forEach((p) => (p.figs || []).forEach((f) => figs.push({ where: `${p.label} の本文`, s: f.smiles || f })));
  const texts = [sheet.intro || '', ...(sheet.parts || []).map((p) => p.intro || ''), ...(sheet.experiments || []).map((e) => e.text),
    ...(sheet.rules || []).map((r) => r.text || '')].join('\n');
  for (const q of sheet.questions) {
    for (const a of answersOf(q)) {
      let ca;
      try { ca = chem.canonical(RDKit, strip(a)); } catch (e) { errors.push(`${q.id}: 答え ${a} が構造として読めない`); continue; }
      for (const f of figs) {
        const k = leakKind(RDKit, ca, f.s);
        // 同じもの（答えそのもの）はどの設問でも不可。読み取れるもの（部分構造・炭素骨格）は、構造を決めさせる設問で不可
        if (k === 'same' || (k && q.kind === 'identify')) errors.push(`${q.id}: ${f.where} の図 ${f.s} から答え ${ca} が読める（${k}）`);
      }
      const nm = names[ca];
      if (nm && q.kind === 'identify' && texts.includes(nm)) errors.push(`${q.id}: 答えの名前「${nm}」が問題文にある`);
    }
    // 事実のリーク: 設問が問う性質（化合物とカード）が、そのまま実験の事実として書かれている
    if (q.about) {
      const hit = (sheet.experiments || []).find((e) => e.fact && e.fact.subject === q.about.subject && e.fact.card === q.about.card);
      if (hit) errors.push(`${q.id}: 答えが実験${hit.label}にそのまま書かれている（${q.about.subject} の ${q.about.card}）`);
    }
  }
  return errors;
}

function dependencyCheck(questions) {
  const errors = [];
  const pos = new Map();
  questions.forEach((q, i) => {
    if (!q.id) errors.push(`問${i + 1}: id がない`);
    else if (pos.has(q.id)) errors.push(`${q.id}: id が重複`);
    else pos.set(q.id, i);
  });
  // 循環の検出（深さ優先）
  const state = new Map();
  const visit = (id, path) => {
    if (state.get(id) === 2) return;
    if (state.get(id) === 1) { errors.push(`依存が循環している: ${[...path, id].join(' → ')}`); return; }
    state.set(id, 1);
    const q = questions[pos.get(id)];
    for (const d of (q && q.depends_on) || []) if (pos.has(d)) visit(d, [...path, id]);
    state.set(id, 2);
  };
  questions.forEach((q, i) => {
    for (const d of q.depends_on || []) {
      if (!pos.has(d)) errors.push(`${q.id}: 依存先 ${d} がない`);
      else if (pos.get(d) > i) errors.push(`${q.id}: 依存先 ${d} より前に置かれている`);
    }
    if (q.id) visit(q.id, []);
  });
  return errors;
}

// 採点の曖昧さを消す注記
function notesCheck(questions) {
  const errors = [];
  for (const q of questions) {
    const n = (q.notes || []).join(' ');
    if (q.kind === 'count' && !/区別/.test(n)) errors.push(`${q.id}: 数える設問に、立体異性体（シス-トランス・鏡像）の扱いの注記がない`);
    if (STRUCT.has(q.type) && q.kind !== 'classify' && !/立体異性体/.test(n)) errors.push(`${q.id}: 構造の設問に「立体異性体」の扱いの注記がない`);
    if (q.type === 'draw' && Array.isArray(q.answer.smiles) && q.answer.smiles.length > 1 && !/順序/.test(n)) errors.push(`${q.id}: 複数の構造を答える設問に「順序は問わない」の注記がない`);
    if (q.type === 'number' && q.sig_figs && !/有効数字/.test(n)) errors.push(`${q.id}: 有効数字の注記がない`);
    if (q.type === 'essay' && !(q.rubric && q.rubric.length >= 2 && q.rubric.length <= 3)) errors.push(`${q.id}: 記述の採点基準（必須要素 2〜3 個）がない`);
  }
  return errors;
}

function auditSheet(RDKit, sheet, names) {
  return [...leakCheck(RDKit, sheet, names), ...dependencyCheck(sheet.questions), ...notesCheck(sheet.questions)];
}

// 大問（京大形式）: 画面に出る形の事実だけで、各化合物の答えがちょうど1つに決まるか
function checkUniqueBig(RDKit, P, exam) {
  const { shownValue } = require('./leak');
  const errors = [];
  // 化合物の記号が重なっていない（断片 D と誘導体 D が同じ記号だと問題文が読めない）
  const labels = [...P.fragments.map((f) => f.label), ...(P.derived || []).map((d) => d.label)];
  if (new Set(labels).size !== labels.length) errors.push(`化合物の記号が重なっている（${labels.join('・')}）`);
  const shown = (card, s, v) => JSON.stringify(shownValue(RDKit, card, s, v));
  const val = (card, s) => { try { return shown(card, s, chem.normalizeResult(RDKit, card, chem.evaluate(RDKit, card, s))); } catch (e) { return 'x'; } };
  // 与えられた断片: 分子式と官能基の数で絞った全異性体から、書いた事実で1つに決まる
  for (const e of exam.experiments.filter((x) => x.given)) {
    const facts = e.given.cards;
    const left = e.given.pool.filter((s) => facts.every((c) => val(c.card, s) === shown(c.card, s, c.result)));
    if (left.length !== 1) errors.push(`${e.label}: 与えた断片が事実で1つに決まらない（${left.length}）`);
  }
  if (!P.chain) {
    for (const f of P.fragments.filter((x) => !x.given)) {
      const ans = chem.canonical(RDKit, f.answer);
      const left = f.candidates.map((s) => chem.canonical(RDKit, s)).filter((s) => (f.clues || []).every((c) => val(c.card, s) === shown(c.card, ans, chem.normalizeResult(RDKit, c.card, c.result))));
      if (left.length !== 1 || left[0] !== ans) errors.push(`${f.label}: 問題文の事実だけでは ${left.length} 通り`);
    }
  }
  // chain は chain.check が問題文の形の事実で解の数を数える（validate の checkBig で実行）
  if (P.assemble) {
    const X = chem.canonical(RDKit, P.answer);
    const sigs = (s) => { try { return chem.partialProducts(RDKit, chem.graphFromSmiles(RDKit, s)).map((pp) => shown('partial_hydrolysis', s, pp)); } catch (e) { return []; } };
    const left = P.assemble.candidates.map((s) => chem.canonical(RDKit, s)).filter((s) => P.assemble.clues.every((c) => (c.card === 'partial_hydrolysis'
      ? sigs(s).includes(shown('partial_hydrolysis', X, chem.canonical(RDKit, c.result)))
      : val(c.card, s) === shown(c.card, X, chem.normalizeResult(RDKit, c.card, c.result)))));
    if (left.length !== 1 || left[0] !== X) errors.push(`X: 問題文の事実だけではつなぎ方が ${left.length} 通り`);
  }
  return errors;
}

module.exports = { auditSheet, leakCheck, dependencyCheck, notesCheck, checkUniqueBig, answersOf };
