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
  ['acetylide', 'C#C', true, 'アセチレン'],
  ['acetylide', 'C#CCC', true, '1-ブチン'],
  ['acetylide', 'CC#CC', false, '2-ブチン'],
  ['acetylide', 'C=CC=C', false, '1,3-ブタジエン'],
  ['alkyne_hydration', 'C#C', ['CC=O'], 'アセチレン → アセトアルデヒド'],
  ['alkyne_hydration', 'C#CC', ['CC(C)=O'], 'プロピン → アセトン'],
  ['alkyne_hydration', 'C#CCC', ['CCC(C)=O'], '1-ブチン → 2-ブタノン'],
  ['alkyne_hydration', 'CC#CC', ['CCC(C)=O'], '2-ブチン → 2-ブタノン'],
  ['alkyne_hydration', 'CC#CCC', ['CCC(=O)CC', 'CCCC(C)=O'], '2-ペンチン → 2-ペンタノンと 3-ペンタノン'],
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
  // ---- 芳香族 ----
  ['kmno4', 'Cc1ccccc1', ['O=C(O)c1ccccc1'], 'トルエン → 安息香酸'],
  ['kmno4', 'CCc1ccccc1', ['O=C(O)c1ccccc1'], 'エチルベンゼン → 安息香酸'],
  ['kmno4', 'Cc1ccccc1C', ['O=C(O)c1ccccc1C(=O)O'], 'o-キシレン → フタル酸'],
  ['kmno4', 'CC(C)(C)c1ccccc1', [], 'tert-ブチルベンゼン → 酸化されない'],
  ['kmno4', 'OCc1ccccc1', ['O=C(O)c1ccccc1'], 'ベンジルアルコール → 安息香酸'],
  ['naoh', 'Oc1ccccc1', true, 'フェノール'],
  ['naoh', 'OCc1ccccc1', false, 'ベンジルアルコール'],
  ['hcl', 'Nc1ccccc1', true, 'アニリン'],
  ['hcl', 'CC(=O)Nc1ccccc1', false, 'アセトアニリド'],
  ['hydrolysis', 'CC(=O)Oc1ccccc1', ['CC(=O)O', 'Oc1ccccc1'], '酢酸フェニル'],
  ['hydrolysis', 'CC(=O)Nc1ccccc1', ['CC(=O)O', 'Nc1ccccc1'], 'アセトアニリド'],
  ['hydrolysis', 'CC(=O)Oc1ccccc1C(=O)O', ['CC(=O)O', 'O=C(O)c1ccccc1O'], 'アセチルサリチル酸'],
  ['ring_cl', 'Cc1ccc(C)cc1', 1, 'p-キシレン'],
  ['ring_cl', 'Cc1cccc(C)c1', 3, 'm-キシレン'],
  ['ring_cl', 'Cc1ccccc1C', 2, 'o-キシレン'],
  ['ring_cl', 'Cc1ccccc1', 3, 'トルエン'],
  ['cl_sub', 'CCCC', 2, 'ブタン'],
  ['cl_sub', 'CC(C)C', 2, '2-メチルプロパン'],
  ['cl_sub', 'CCCCC', 3, 'ペンタン'],
  ['anhydride', 'O=C(O)c1ccccc1C(=O)O', true, 'フタル酸'],
  ['anhydride', 'O=C(O)c1ccc(C(=O)O)cc1', false, 'テレフタル酸'],
  ['anhydride', 'O=C(O)CCC(=O)O', true, 'コハク酸'],
  ['anhydride', 'O=C(O)CC(=O)O', false, 'マロン酸'],
  // ---- 脱水 ----
  ['dehydration', 'CCO', ['C=C'], 'エタノール → エチレン'],
  ['dehydration', 'CCC(C)O', ['C=CCC', 'CC=CC'], '2-ブタノール → 1-ブテン + 2-ブテン'],
  ['dehydration', 'CC(C)(C)O', ['C=C(C)C'], '2-メチル-2-プロパノール'],
  ['dehydration', 'CCC(C)=O', [], '2-ブタノン'],
  ['dehydration_count', 'CCC(C)O', 3, '2-ブタノール（シス・トランスを含め3種）'],
  ['dehydration_count', 'CCC(O)CC', 2, '3-ペンタノール（2-ペンテンのシス・トランス）'],
  ['dehydration_ozonolysis', 'CCC(C)(C)O', ['C=O', 'CC(C)=O', 'CC=O', 'CCC(C)=O'], '2-メチル-2-ブタノール'],
  // ---- 部分加水分解・ペプチド ----
  ['partial_hydrolysis', 'CC(OC(C)=O)COC(=O)c1ccccc1', ['CC(=O)O', 'CC(CO)OC(C)=O', 'CC(O)COC(=O)c1ccccc1', 'O=C(O)c1ccccc1'], '1,2-プロパンジオールのジエステル'],
  ['ninhydrin', 'NCC(=O)O', true, 'グリシン'],
  ['ninhydrin', 'CNCC(=O)O', false, 'サルコシン（第二級アミン）'],
  ['alpha_amino', 'NC(C)C(=O)O', true, 'アラニン'],
  ['alpha_amino', 'NCCC(=O)O', false, 'β-アラニン'],
  ['xanthoprotein', 'pep:Gly-Phe', true, 'Gly-Phe'],
  ['xanthoprotein', 'pep:Gly-Ala', false, 'Gly-Ala'],
  ['sulfur', 'pep:Cys-Gly', true, 'Cys-Gly'],
  ['biuret', 'pep:Gly-Ala', false, 'ジペプチド'],
  ['biuret', 'pep:Gly-Ala-Phe', true, 'トリペプチド'],
  ['hydrolysis', 'pep:Gly-Ala', ['CC(N)C(=O)O', 'NCC(=O)O'], 'Gly-Ala'],
  ['chiral', 'pep:Gly-Ala-Phe', 2, 'Gly-Ala-Phe'],
  // ---- 高分子 ----
  ['hydrolysis', 'CC(=O)OC=C', ['CC(=O)O', 'CC=O'], '酢酸ビニル → 酢酸 + アセトアルデヒド'],
  ['hydrolysis', 'C=CC(=O)OC', ['C=CC(=O)O', 'CO'], 'アクリル酸メチル'],
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
    const exp = chem.CARDS[card].kind === 'contains'
      ? want.map((s) => chem.canonical(RDKit, s)).sort()
      : chem.normalizeResult(RDKit, card, want);
    if (!chem.sameResult(got, exp)) {
      fail++;
      console.error(`NG ${card} ${label} (${smi}): got ${JSON.stringify(got)}, want ${JSON.stringify(exp)}`);
    }
  }
  // 分子式
  const f = chem.formula(chem.graphFromSmiles(RDKit, 'CCC(C)=O'));
  if (f !== 'C4H8O') { fail++; console.error(`NG formula: ${f}`); }
  const u = chem.formula(chem.graphFromSmiles(RDKit, '*C(=O)CCCCC(=O)NCCCCCCN*'));
  if (u !== 'C12H22N2O2' || chem.formulaMass(u) !== 226) { fail++; console.error(`NG nylon66 unit: ${u} ${chem.formulaMass(u)}`); }
  if (fail) {
    console.error(`${fail} rule test(s) failed`);
    process.exit(1);
  }
  console.log(`OK: ${CASES.length + 1} rule tests passed`);
}

main();
