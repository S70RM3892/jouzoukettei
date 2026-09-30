'use strict';
// 京大2026型の大問1問を、検証パイプライン（要望書4章）に通す。
//  1. 一意性: 別解答（solve26.js）が、前提の範囲の全異性体から事実で絞り、ちょうど1つを返す
//  2. リーク: 規則の図・問題の図・名前が設問の答え（読み取れるもの）になっていない（audit.js）
//  3. 依存関係: depends_on が非循環で、依存先が先に置かれている（audit.js）
//  4. 数値検算: Python（scripts/numcheck.py）で分数計算し直す（calcs を集めて一括で）
//  5. 別インスタンス解答: 答えを取り除いた問題（publicView）だけから別の実装で解き、想定解と一致する
//  6. 採点基準: 記述には必須要素 2〜3 個（audit.js の notesCheck）
const A = require('../src/audit');
const S = require('../src/solve26');

const isMol = (x) => typeof x === 'string' && !['+', '⇄', '→'].includes(x);
function sheetOf(P) {
  return {
    rules: P.rules.map((r) => ({ id: r.id, text: r.text, figs: (r.scheme || []).filter(isMol) })),
    parts: P.parts,
    intro: P.parts.map((p) => [p.intro, p.tail || '', p.tailY || '', p.limitation || ''].join('')).join('\n'),
    experiments: P.experiments,
    questions: P.questions,
  };
}

function checkK26(RDKit, P, names = {}) {
  const errors = [];
  const sheet = sheetOf(P);
  const leak = A.leakCheck(RDKit, sheet, names);
  const deps = A.dependencyCheck(P.questions);
  const notes = A.notesCheck(P.questions);
  let agreeErr = [];
  try {
    const sol = S.solve(RDKit, S.publicView(P));
    agreeErr = S.agree(RDKit, P, sol);
  } catch (e) { agreeErr = [`別解答が失敗: ${e.message}`]; }
  const uniqueErr = agreeErr.filter((m) => /一意/.test(m));
  errors.push(...leak.map((m) => `leak ${m}`), ...deps.map((m) => `deps ${m}`), ...notes.map((m) => `notes ${m}`), ...agreeErr.map((m) => `solver ${m}`));
  // 3要素（規則定義ブロック・構造式を書く設問・量的設問）
  if (!P.rules.length) errors.push('規則定義ブロックがない');
  if (!P.questions.some((q) => q.type === 'draw')) errors.push('構造式を書く設問がない');
  if (!P.questions.some((q) => q.type === 'number' && q.calc)) errors.push('量的設問がない');
  return {
    errors: errors.map((m) => `${P.id}: ${m}`),
    validation: { unique: !uniqueErr.length && !agreeErr.length, leak_free: !leak.length, deps_ok: !deps.length, notes_ok: !notes.length, solver_agree: !agreeErr.length },
  };
}

function calcsOf(P) {
  return P.questions.filter((q) => q.calc).map((q) => ({ id: `${P.id}/${q.id}`, calc: q.calc }));
}

module.exports = { checkK26, calcsOf, sheetOf };
