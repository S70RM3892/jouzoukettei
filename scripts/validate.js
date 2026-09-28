'use strict';
// 問題データの自動検証。仕様書「問題データ形式と自動検証」の4項目＋データ整合性を確認する。
// 使い方: node scripts/validate.js [problems/narrow.json]
const fs = require('fs');
const path = require('path');
const chem = require('../src/chem');

async function loadRDKit() {
  const init = require('@rdkit/rdkit');
  return init();
}

// 1問を検証し、エラー一覧と「候補×手がかり」の結果表を返す
function checkProblem(RDKit, p) {
  const errors = [];
  const err = (m) => errors.push(`${p.id}: ${m}`);
  if (p.mode !== 'narrow') err(`unsupported mode ${p.mode}`);
  if (!Array.isArray(p.candidates) || p.candidates.length < 2) err('candidates must have 2+ entries');
  if (!Array.isArray(p.clues) || p.clues.length < 1) err('clues must have 1+ entries');
  if (errors.length) return { errors };

  let cands;
  try {
    cands = p.candidates.map((s) => chem.canonical(RDKit, s));
  } catch (e) {
    err(e.message);
    return { errors };
  }
  const answer = chem.canonical(RDKit, p.answer);

  // 1. 全候補の分子式が formula と一致する
  cands.forEach((s) => {
    const f = chem.formula(chem.graphFromSmiles(RDKit, s));
    if (f !== p.formula) err(`candidate ${s} is ${f}, expected ${p.formula}`);
  });
  // 2. 候補に重複がない
  if (new Set(cands).size !== cands.length) err('duplicate candidates (canonical SMILES)');
  if (!cands.includes(answer)) err(`answer ${answer} is not among candidates`);

  const table = []; // table[clue][cand] = 結果
  const expected = [];
  p.clues.forEach((c, ci) => {
    if (!chem.CARDS[c.card]) {
      err(`unknown card ${c.card}`);
      return;
    }
    try {
      expected[ci] = chem.normalizeResult(RDKit, c.card, c.result);
    } catch (e) {
      err(`clue ${c.card}: bad result (${e.message})`);
      return;
    }
    table[ci] = cands.map((s) => {
      try {
        return chem.evaluate(RDKit, c.card, s);
      } catch (e) {
        err(`clue ${c.card} cannot be applied to ${s} (${e.message})`);
        return null;
      }
    });
  });
  if (errors.length) return { errors };

  const ai = cands.indexOf(answer);
  p.clues.forEach((c, ci) => {
    if (!chem.sameResult(table[ci][ai], expected[ci])) {
      err(`clue ${c.card}: answer gives ${JSON.stringify(table[ci][ai])}, data says ${JSON.stringify(expected[ci])}`);
    }
  });
  // 3. 全手がかりを適用すると answer だけが残る
  const survivors = cands.filter((_, k) => p.clues.every((c, ci) => chem.sameResult(table[ci][k], expected[ci])));
  if (survivors.length !== 1 || survivors[0] !== answer) {
    err(`after all clues, remaining = ${JSON.stringify(survivors)}`);
  }
  // 4. どの手がかりも単独で候補を1つ以上消す
  p.clues.forEach((c, ci) => {
    const removed = cands.filter((_, k) => !chem.sameResult(table[ci][k], expected[ci]));
    if (removed.length === 0) err(`clue ${c.card} eliminates nothing`);
  });
  return { errors, cands, answer, table, expected };
}

function checkAll(RDKit, problems) {
  const errors = [];
  const results = [];
  const ids = new Set();
  for (const p of problems) {
    if (ids.has(p.id)) errors.push(`${p.id}: duplicate id`);
    ids.add(p.id);
    const r = checkProblem(RDKit, p);
    errors.push(...r.errors);
    results.push(r);
  }
  return { errors, results };
}

async function main() {
  const file = process.argv[2] || path.join(__dirname, '..', 'problems', 'narrow.json');
  const problems = JSON.parse(fs.readFileSync(file, 'utf8'));
  const RDKit = await loadRDKit();
  const { errors } = checkAll(RDKit, problems);
  if (errors.length) {
    console.error(`NG: ${errors.length} error(s)`);
    errors.forEach((e) => console.error('  ' + e));
    process.exit(1);
  }
  console.log(`OK: ${problems.length} problems passed`);
}

if (require.main === module) main();

module.exports = { loadRDKit, checkProblem, checkAll };
