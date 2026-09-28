'use strict';
// 京大の過去問の設問を、このエンジン（異性体列挙＋判定カード＋計算器）だけで解かせ、
// 京大の正解と一致するかを確かめる。構造は自分で解いた答えを SMILES で書いている（問題文は写していない）。
// 解けない出題要素は最後に一覧にして、エンジンを伸ばす順番の材料にする。
const chem = require('../src/chem');
const calc = require('../src/calc');
const { enumerate, hasMatch } = require('../src/enumerate');
const { loadRDKit } = require('./validate');

async function main() {
  const R = await loadRDKit();
  const can = (s) => chem.canonical(R, s);
  const ev = (card, s) => chem.evaluate(R, card, s);
  const norm = (card, v) => chem.normalizeResult(R, card, v);
  const fits = (s, clues) => clues.every(([card, v]) => chem.consistent(card, ev(card, s), norm(card, v)));
  // 候補集合を列挙し、手がかりで絞り込んだ結果が京大の答えと一致するか
  const solve = (formula, opts, clues) => enumerate(R, formula, opts).filter((s) => fits(s, clues));
  const same = (a, b) => JSON.stringify([...a].map(can).sort()) === JSON.stringify([...b].map(can).sort());
  const benz = { seed: 'c1ccccc1', rings: 0 };

  const CASES = [
    ['2007 III(a) C₈H₁₀ 炭素の種類 3', () => same(solve('C8H10', benz, [['carbon_env', 3]]), ['Cc1ccc(C)cc1'])],
    ['2007 III(a) C₈H₁₀ 炭素の種類 5', () => same(solve('C8H10', benz, [['carbon_env', 5]]), ['Cc1cccc(C)c1'])],
    ['2007 III(a) C₈H₁₀ 炭素4種・KMnO₄ 生成物が酸無水物に', () => same(solve('C8H10', benz, [['carbon_env', 4]]), ['Cc1ccccc1C'])
      && ev('anhydride', ev('kmno4', 'Cc1ccccc1C')[0]) === true],
    ['2007 III(a) C₈H₁₀ 炭素の種類 6', () => same(solve('C8H10', benz, [['carbon_env', 6]]), ['CCc1ccccc1'])],
    ['2019 III(a) C₉H₁₂ 炭素の種類 3', () => same(solve('C9H12', benz, [['carbon_env', 3]]), ['Cc1cc(C)cc(C)c1'])],
    ['2019 III(a) C₉H₁₂ 炭素6種・KMnO₄ で安息香酸（クメン）', () => same(solve('C9H12', benz, [['carbon_env', 6], ['kmno4', ['O=C(O)c1ccccc1']]]), ['CC(C)c1ccccc1'])],
    ['2019 III(a) C₉H₁₂ 炭素7種・安息香酸', () => same(solve('C9H12', benz, [['carbon_env', 7], ['kmno4', ['O=C(O)c1ccccc1']]]), ['CCCc1ccccc1'])],
    ['2019 III(a) C₉H₁₂ 炭素7種・テレフタル酸', () => same(solve('C9H12', benz, [['carbon_env', 7], ['kmno4', ['O=C(O)c1ccc(C(=O)O)cc1']]]), ['CCc1ccc(C)cc1'])],
    ['2019 III(a) C₉H₁₂ 等価な炭素のないもの3種', () => solve('C9H12', benz, [['carbon_env', 9]]).length === 3],
    ['2008 III(a) 鎖状 C₄H₈O の A（Br₂付加・Na・不斉炭素）', () => same(solve('C4H8O', { rings: 0 }, [['h2_uptake', 1], ['sodium', true], ['chiral', 1]]), ['C=CC(C)O'])],
    ['2008 III(a) A を酸化した D に H₂ で C（臭素と反応しない）', () => {
      const d = ev('mild_oxidation', 'C=CC(C)O')[0];
      const c = ev('hydrogenation', d)[0];
      return can(d) === can('C=CC(C)=O') && can(c) === can('CCC(C)=O') && ev('bromine', c) === false;
    }],
    ['2008 III(a) エーテル B をすべて（Br₂付加・不斉炭素なし）: 構造4・シストランス込み5', () => {
      const bs = solve('C4H8O', { rings: 0 }, [['h2_uptake', 1], ['sodium', false], ['chiral', 0]]).filter((s) => hasMatch(R, s, '[OX2]([#6])[#6]'));
      const withStereo = bs.reduce((n, s) => n + (chem.hasCisTrans(chem.graphFromSmiles(R, s)) ? 2 : 1), 0);
      return bs.length === 4 && withStereo === 5;
    }],
    ['2016 III(a) 燃焼分析 18.5 mg → C₄H₁₀O', () => calc.molecularFromEmpirical(calc.empiricalFromCombustion(18.5, 43.9, 22.5), { maxMW: 100 }) === 'C4H10O'],
    ['2016 III(a) ジエステル A・B が同じ分子式、不斉炭素1つずつ', () => {
      const A = 'CC(C)OC(=O)C=CC(=O)OC(C)CC';
      const B = 'CCOC(=O)C=CC(=O)OCC(C)CC';
      return chem.formula(chem.graphFromSmiles(R, A)) === chem.formula(chem.graphFromSmiles(R, B))
        && ev('chiral', A) === 1 && ev('chiral', B) === 1
        && same(ev('hydrolysis', A), ['CC(C)O', 'CCC(C)O', 'O=C(O)C=CC(=O)O'])
        && same(ev('mild_oxidation', 'CC(C)O'), ['CC(C)=O'])
        && ev('iodoform', 'CCC(C)CO') === false;
    }],
    ['2013 III(a) C₁₆H₁₈O₄: H₂ 2 mol、加水分解でフタル酸・ブタナール（ビニルエステル）・3-ブテン-1-オール', () => {
      const A = 'C=CCCOC(=O)c1ccccc1C(=O)OC=CCC';
      return chem.formula(chem.graphFromSmiles(R, A)) === 'C16H18O4' && ev('h2_uptake', A) === 2
        && same(ev('hydrolysis', A), ['O=C(O)c1ccccc1C(=O)O', 'CCCC=O', 'C=CCCO'])
        && ev('cis_trans', 'C=CCCO') === false;
    }],
    ['2013 III(a) 水素付加した B の加水分解で 1-ブタノール、脱水で 1-ブテン、酸化でブタナール', () => {
      const B = ev('hydrogenation', 'C=CCCOC(=O)c1ccccc1C(=O)OC=CCC')[0];
      const prods = ev('hydrolysis', B);
      return same(prods, ['O=C(O)c1ccccc1C(=O)O', 'CCCCO', 'CCCCO'])
        && same(ev('dehydration', 'CCCCO'), ['C=CCC']) && same(ev('mild_oxidation', 'CCCCO'), ['CCCC=O']);
    }],
    ['2017 III(a) C₁₆H₁₈O₆: H₂ 0.732 L/10 g、加水分解でグリセリン・酢酸2・不飽和酸、不斉炭素なし', () => {
      const A = 'CC(=O)OCC(COC(C)=O)OC(=O)C=Cc1ccccc1';
      const L = calc.liters((10 / calc.mass('C16H18O6')) * ev('h2_uptake', A));
      return chem.formula(chem.graphFromSmiles(R, A)) === 'C16H18O6' && Math.abs(L - 0.732) < 0.002 && ev('chiral', A) === 0
        && same(ev('hydrolysis', A), ['OCC(O)CO', 'CC(=O)O', 'CC(=O)O', 'O=C(O)C=Cc1ccccc1'])
        && ev('ozonolysis', 'O=C(O)C=Cc1ccccc1').map(can).includes(can('O=Cc1ccccc1'))
        && ev('cis_trans', 'O=C(O)C=Cc1ccccc1') === true;
    }],
    ['2012 III(a) O₂ 8.5 mol と質量%から C₆H₁₀、KMnO₄ 開裂でアジピン酸（シクロヘキセン）', () => {
      const emp = calc.empiricalFromPercent({ C: 87.8, H: 12.2 });
      const f = ['C3H5', 'C6H10', 'C9H15'].find((x) => calc.o2ForCombustion(x) === 8.5);
      return emp === 'C3H5' && f === 'C6H10' && same(ev('kmno4_cleave', 'C1=CCCCC1'), ['O=C(O)CCCCC(=O)O']);
    }],
    ['2012 III(a) D のオゾン分解生成物を還元した二価アルコール（1-メチルシクロペンテン）', () => same(ev('ozonolysis', 'CC1=CCCC1'), ['CC(=O)CCCC=O'])],
    ['2012 III(a) C（メチル2つ）に H₂ 2分子でメチル4つのアルカン', () => ev('h2_uptake', 'C=C(C)C(C)=C') === 2 && same(ev('hydrogenation', 'C=C(C)C(C)=C'), ['CC(C)C(C)C'])],
    ['2011 III(a) C₉H₁₀ の E を酸化 → p-トルイル酸（さらに空気酸化でテレフタル酸）', () => same(ev('kmno4_cleave', 'C=Cc1ccc(C)cc1'), ['Cc1ccc(C(=O)O)cc1']) && same(ev('kmno4', 'Cc1ccc(C(=O)O)cc1'), ['O=C(O)c1ccc(C(=O)O)cc1'])],
    ['2006 III C₄H₁₁N のうちアミド D になり得ないもの（第三級アミン）', () => {
      const all = enumerate(R, 'C4H11N', { rings: 0 });
      const tert = all.filter((s) => hasMatch(R, s, '[NX3H0]'));
      return all.length === 8 && same(tert, ['CCN(C)C']);
    }],
    ['2022 III(a) 酢酸のアミドが C₅H₁₁NO になる第一級アミン B をすべて（2つ）', () => {
      const bs = enumerate(R, 'C3H9N', { rings: 0 }).filter((s) => hasMatch(R, s, '[NX3H2]'));
      return same(bs, ['CCCN', 'CC(C)N']);
    }],
    ['2018 IV 問1 不飽和酸 A の KMnO₄ 開裂でプロピオン酸1・マロン酸2', () => same(ev('kmno4_cleave', 'CCC=CCC=CCC(=O)O'), ['CCC(=O)O', 'O=C(O)CC(=O)O', 'O=C(O)CC(=O)O'])],
    ['2022 IV 過ヨウ素酸開裂（与えられた規則）: グリセルアルデヒド → ギ酸2・ホルムアルデヒド', () => same(ev('periodate', 'OCC(O)C=O'), ['O=CO', 'O=CO', 'C=O'])],
    ['2022 IV 過ヨウ素酸開裂: グルコース（鎖状）→ ギ酸5・ホルムアルデヒド1', () => same(ev('periodate', 'OCC(O)C(O)C(O)C(O)C=O'), ['O=CO', 'O=CO', 'O=CO', 'O=CO', 'O=CO', 'C=O'])],
    ['2015 III(b) ヒドロキシ酸2分子の環状エステル（ラクチド）: 分子式と、加水分解で乳酸2分子', () => chem.formula(chem.graphFromSmiles(R, 'CC1OC(=O)C(C)OC1=O')) === 'C6H8O4'
      && same(ev('hydrolysis', 'CC1OC(=O)C(C)OC1=O'), ['CC(O)C(=O)O', 'CC(O)C(=O)O'])],
  ];

  // エンジンでは扱えない出題要素（次に伸ばす候補）
  const UNSUPPORTED = [
    '2026 III アセタール交換の平衡・光学分割の速度（アセタールの規則が未実装）',
    '2026 IV メチル化分析・糖の立体（糖の立体を区別する表現が未実装）',
    '2019 IV アルドースの立体・メソ体（フィッシャー投影が未実装）',
    '2019 III(b) イミドの穏やかな加水分解（位置選択の規則が未実装）',
    '2010 III(a) 第一級 OH の選択的アセチル化（与えられた規則が未実装）',
    '2018 III 配向性・合成経路の設計（置換反応が未実装）',
    '2016 III(b) ジアリールエーテルの水素化分解（反応規則が未実装）',
  ];

  let ok = 0;
  for (const [label, f] of CASES) {
    let r;
    try { r = f(); } catch (e) { r = false; console.log('   ERROR ' + e.message); }
    console.log(`${r ? 'OK' : 'NG'} ${label}`);
    if (r) ok++;
  }
  console.log(`\n過去問の設問 ${CASES.length} 問中 ${ok} 問をエンジンだけで解けた。`);
  console.log(`まだ解けない出題要素 ${UNSUPPORTED.length} 件:`);
  UNSUPPORTED.forEach((u) => console.log('  - ' + u));
  if (ok !== CASES.length) process.exit(1);
}
main();
