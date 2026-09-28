'use strict';
// 判定ロジックの単体テスト。教科書に載っている典型例で各カードの答えを確認する。
const chem = require('../src/chem');
const { loadRDKit } = require('./validate');

const CASES = [
  // [card, SMILES, 期待値, 物質名]
  ['silver_mirror', 'CC=O', true, 'アセトアルデヒド'],
  ['silver_mirror', 'C=O', true, 'ホルムアルデヒド'],
  ['silver_mirror', 'O=CO', true, 'ギ酸'],
  ['silver_mirror', 'COC=O', true, 'ギ酸メチル'],
  ['silver_mirror', 'CC(C)=O', false, 'アセトン'],
  ['silver_mirror', 'CC(=O)O', false, '酢酸'],
  ['iodoform', 'CCO', true, 'エタノール'],
  ['iodoform', 'CC=O', true, 'アセトアルデヒド'],
  ['iodoform', 'CC(C)=O', true, 'アセトン'],
  ['iodoform', 'CC(C)O', true, '2-プロパノール'],
  ['iodoform', 'CC(O)C(=O)O', true, '乳酸'],
  ['iodoform', 'CO', false, 'メタノール'],
  ['iodoform', 'CCCO', false, '1-プロパノール'],
  ['iodoform', 'CC(=O)O', false, '酢酸'],
  ['iodoform', 'CCOC(C)=O', false, '酢酸エチル'],
  ['iodoform', 'CC(C)(C)O', false, '2-メチル-2-プロパノール'],
  ['iodoform', 'CCC(=O)CC', false, '3-ペンタノン'],
  ['sodium', 'CCO', true, 'エタノール'],
  ['sodium', 'CC(=O)O', true, '酢酸'],
  ['sodium', 'CCOCC', false, 'ジエチルエーテル'],
  ['sodium', 'CC(C)=O', false, 'アセトン'],
  ['nahco3', 'CC(=O)O', true, '酢酸'],
  ['nahco3', 'Oc1ccccc1', false, 'フェノール'],
  ['nahco3', 'CCO', false, 'エタノール'],
  ['fecl3', 'Oc1ccccc1', true, 'フェノール'],
  ['fecl3', 'OCc1ccccc1', false, 'ベンジルアルコール'],
  ['bromine', 'C=CC', true, 'プロペン'],
  ['bromine', 'C#C', true, 'アセチレン'],
  ['bromine', 'CC(C)=O', false, 'アセトン'],
  ['kmno4', 'CCCO', ['CCC(=O)O'], '1-プロパノール → プロピオン酸'],
  ['kmno4', 'CC(C)O', ['CC(C)=O'], '2-プロパノール → アセトン'],
  ['kmno4', 'CC(C)(C)O', [], '2-メチル-2-プロパノール → 酸化されない'],
  ['kmno4', 'CC=O', ['CC(=O)O'], 'アセトアルデヒド → 酢酸'],
  ['kmno4', 'CC(C)=O', [], 'アセトン → 酸化されない'],
  ['hydrolysis', 'CCOC(C)=O', ['CC(=O)O', 'CCO'], '酢酸エチル'],
  ['hydrolysis', 'COC=O', ['CO', 'O=CO'], 'ギ酸メチル'],
  ['hydrolysis', 'CCC(=O)O', [], 'プロピオン酸'],
  ['ozonolysis', 'CC=CC', ['CC=O', 'CC=O'], '2-ブテン'],
  ['ozonolysis', 'C=C(C)C', ['C=O', 'CC(C)=O'], '2-メチルプロペン'],
  ['ozonolysis', 'CCO', [], 'エタノール'],
  ['chiral', 'CCC(C)O', 1, '2-ブタノール'],
  ['chiral', 'CC(O)C(=O)O', 1, '乳酸'],
  ['chiral', 'OC(C(O)C(=O)O)C(=O)O', 2, '酒石酸（メソ体も数は2）'],
  ['chiral', 'CC(C)O', 0, '2-プロパノール'],
  ['chiral', 'OC1CCCC1', 0, 'シクロペンタノール'],
  ['chiral', 'CC1CCCC1O', 2, '2-メチルシクロペンタノール'],
  ['chiral', 'NC(C)C(=O)O', 1, 'アラニン'],
  ['chiral', 'NCC(=O)O', 0, 'グリシン'],
  ['cis_trans', 'CC=CC', true, '2-ブテン'],
  ['cis_trans', 'C=CCC', false, '1-ブテン'],
  ['cis_trans', 'CC(C)=CC', false, '2-メチル-2-ブテン'],
  ['cis_trans', 'OC(=O)C=CC(=O)O', true, 'マレイン酸／フマル酸'],
  ['cis_trans', 'C1=CCCC1', false, 'シクロペンテン'],
];

async function main() {
  const RDKit = await loadRDKit();
  let fail = 0;
  for (const [card, smi, want, label] of CASES) {
    let got;
    try {
      got = chem.evaluate(RDKit, card, smi);
    } catch (e) {
      got = `ERROR ${e.message}`;
    }
    const exp = chem.normalizeResult(RDKit, card, want);
    if (!chem.sameResult(got, exp)) {
      fail++;
      console.error(`NG ${card} ${label} (${smi}): got ${JSON.stringify(got)}, want ${JSON.stringify(exp)}`);
    }
  }
  // 分子式
  const f = chem.formula(chem.graphFromSmiles(RDKit, 'CCC(C)=O'));
  if (f !== 'C4H8O') { fail++; console.error(`NG formula: ${f}`); }
  if (fail) {
    console.error(`${fail} rule test(s) failed`);
    process.exit(1);
  }
  console.log(`OK: ${CASES.length + 1} rule tests passed`);
}

main();
