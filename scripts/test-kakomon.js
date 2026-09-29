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

  // 2000〜2005・2009・2020・2023〜2025年度（あとから追加した年）
  const ring = { rings: 1, maxRing: 6 };
  const skeleton = (sm) => { // O を除き、すべて単結合にした炭素骨格（還元・脱水・水素付加のあとと同じ）
    const g = chem.graphFromSmiles(R, sm);
    const keep = g.atoms.map((a, i) => (a.el === 'C' ? i : -1)).filter((i) => i >= 0);
    const idx = new Map(keep.map((i, k) => [i, k]));
    const h = { atoms: keep.map(() => ({ el: 'C' })), bonds: g.bonds.filter((b) => idx.has(b.a) && idx.has(b.b)).map((b) => ({ a: idx.get(b.a), b: idx.get(b.b), order: 1 })) };
    return can(chem.toMolblock(h));
  };
  const monochloro = (sm) => { // H を1つ Cl に置き換えた構造異性体すべて
    const g = chem.graphFromSmiles(R, sm);
    const out = new Set();
    g.atoms.forEach((a, i) => {
      if (a.el !== 'C' || chem.hCount(g, i) === 0) return;
      const h = chem.cloneGraph(g);
      h.atoms.push({ el: 'Cl' });
      h.bonds.push({ a: i, b: h.atoms.length - 1, order: 1 });
      out.add(can(chem.toMolblock(h)));
    });
    return [...out];
  };
  CASES.push(
    ['2024 III(a) 燃焼分析 25.8 mg → CO₂ 66.0・H₂O 27.0、分子量 86 → C₅H₁₀O', () => calc.molecularFromEmpirical(calc.empiricalFromCombustion(25.8, 66.0, 27.0), { mw: 86 }) === 'C5H10O'],
    ['2024 III(a) C₅H₁₀O の5物質（鎖状・環状）を銀鏡・ヨードホルム・不斉炭素・非等価な炭素の数・Na で決める', () => {
      const all = enumerate(R, 'C5H10O', ring);
      // (い) は不斉炭素のないメチルケトン。候補は2つ残り、(う) の条件（水素付加で非等価な炭素が1種類減り、酸化すると (い)）で1つに決まる
      const Iall = all.filter((x) => ev('iodoform', x) && hasMatch(R, x, '[CX3](=O)') && ev('chiral', x) === 0 && !hasMatch(R, x, 'C=C'));
      const U = all.filter((x) => hasMatch(R, x, 'C=C') && !hasMatch(R, x, '[#8]C=C') && !hasMatch(R, x, 'C=O') && hasMatch(R, x, '[OX2H]')
        && ev('carbon_env', ev('hydrogenation', x)[0]) === ev('carbon_env', x) - 1 && Iall.some((y) => same(ev('mild_oxidation', ev('hydrogenation', x)[0]), [y])));
      const I = [...new Set(U.map((x) => can(ev('mild_oxidation', ev('hydrogenation', x)[0])[0])))];
      const A = all.filter((x) => ev('silver_mirror', x) && ev('chiral', x) === 0 && I.some((y) => skeleton(y) === skeleton(x)));
      const E = all.filter((x) => ev('carbon_env', x) === 3 && !hasMatch(R, x, '[CH3]') && !hasMatch(R, x, '*=*') && !ev('sodium', x));
      const O = all.filter((x) => ev('sodium', x) && ev('carbon_env', x) === 3 && (() => { try { return ev('dehydration', x).some((d) => ev('carbon_env', ev('hydrogenation', d)[0]) === 1); } catch (e) { return false; } })());
      return same(I, ['CC(=O)C(C)C']) && same(A, ['CC(C)CC=O']) && same(U, ['C=C(C)C(C)O']) && same(E, ['C1CCOCC1']) && same(O, ['OC1CCCC1']);
    }],
    ['2023 III 問2 ピリジン環をもつ C₇H₉N（環の N に置換基なし）は 9 種類', () => enumerate(R, 'C7H9N', { seed: 'c1ccncc1', rings: 0 }).length === 9],
    ['2025 III ラクトン C₈H₆O₂（フタリド）を加水分解・酸化すると酸無水物をつくるフタル酸。C₈H₈O₂ で NaHCO₃ に溶け、酸化でテレフタル酸になるのは p-トルイル酸だけ', () => {
      const F = ev('hydrolysis', 'O=C1OCc2ccccc21');
      const G = ev('kmno4', F[0]);
      const H = enumerate(R, 'C8H8O2', benz).filter((x) => ev('nahco3', x)).filter((x) => { try { return same(ev('kmno4', x), ['O=C(O)c1ccc(C(=O)O)cc1']); } catch (e) { return false; } });
      return same(F, ['O=C(O)c1ccccc1CO']) && same(G, ['O=C(O)c1ccccc1C(=O)O']) && ev('anhydride', G[0]) === true && same(H, ['Cc1ccc(C(=O)O)cc1']);
    }],
    ['2020 III 1.00 mmol の燃焼で H₂O 216 mg・CO₂ 836 mg、分子量 348 → C₁₉H₂₄O₆。X の部分加水分解物に C₇H₁₄O₃ と組成式 C₃H₃O の酸', () => {
      const X = 'Cc1ccc(COC(=O)C=CC(=O)OCC(=O)OC(C)C(C)C)cc1';
      const parts = chem.partialProducts(R, chem.graphFromSmiles(R, X)).map((x) => chem.formula(chem.graphFromSmiles(R, x)));
      return calc.molecularFromEmpirical(calc.empiricalFromCombustion(348, 836, 216), { mw: 348 }) === 'C19H24O6'
        && chem.formula(chem.graphFromSmiles(R, X)) === 'C19H24O6' && ev('chiral', X) === 1
        && parts.includes('C7H14O3') && parts.includes('C12H12O4') && same(ev('markovnikov', 'C=CC(C)C'), ['CC(O)C(C)C']);
    }],
    ['2009 III(a) C₃₀H₂₅NO₅ を加水分解すると フェノール2・ベンゼントリカルボン酸・トリメチルアニリン（2:1:1）', () => {
      const X = 'O=C(Oc1ccccc1)c1cc(C(=O)Oc2ccccc2)cc(C(=O)Nc2c(C)cc(C)cc2C)c1';
      return chem.formula(chem.graphFromSmiles(R, X)) === 'C30H25NO5' && same(ev('hydrolysis', X), ['Oc1ccccc1', 'Oc1ccccc1', 'O=C(O)c1cc(C(=O)O)cc(C(=O)O)c1', 'Cc1cc(C)c(N)c(C)c1'])
        && same(ev('kmno4', 'Cc1cc(C)cc(C)c1'), ['O=C(O)c1cc(C(=O)O)cc(C(=O)O)c1']) && ev('ring_cl', 'Cc1cc(C)cc(C)c1') === 1;
    }],
    ['2009 III(b) C₆H₁₂ のアルケンに H₂ で 5 種のアルカン、その一塩素置換体で不斉炭素のないものは合計 8', () => {
      const alkanes = [...new Set(enumerate(R, 'C6H12', { rings: 0 }).map((x) => ev('hydrogenation', x)[0]).map(can))];
      const n = alkanes.reduce((a, x) => a + monochloro(x).filter((y) => ev('chiral', y) === 0).length, 0);
      return alkanes.length === 5 && n === 8;
    }],
    ['2005 III(c) C₄H₈ に HCl を付加（マルコフニコフ則）して不斉炭素ができないのは 2-メチルプロペンだけ', () => {
      const n = Math.round(36.5 / (0.65 * 14));
      const ok = enumerate(R, `C${n}H${2 * n}`, { rings: 0 }).filter((x) => ev('markovnikov', x).every((y) => ev('chiral', y) === 0));
      return n === 4 && same(ok, ['C=C(C)C']);
    }],
    ['2004 III 段階的な水素付加: C₈H₆ → スチレン → エチルベンゼン → C₈H₁₆（燃焼 7.0 g → CO₂ 22.0 g・H₂O 9.0 g）', () => calc.molecularFromEmpirical(calc.empiricalFromCombustion(7.0, 22.0, 9.0), { maxMW: 120 }) === 'C8H16'
      && ev('h2_uptake', 'C#Cc1ccccc1') === 2 && same(ev('hydrogenation', 'C=Cc1ccccc1'), ['CCc1ccccc1'])],
    ['2004 III スチレンに水を付加して不斉炭素をもつ E（C₈H₁₀O）、スルホン化体のうち臭素置換体が最も少ないのは p 体', () => {
      const E = ev('markovnikov', 'C=Cc1ccccc1');
      const iso = ['CCc1ccccc1S(=O)(=O)O', 'CCc1cccc(S(=O)(=O)O)c1', 'CCc1ccc(S(=O)(=O)O)cc1'];
      const min = iso.reduce((a, x) => (ev('cl_sub', x) < ev('cl_sub', a) ? x : a));
      return same(E, ['CC(O)c1ccccc1']) && ev('chiral', E[0]) === 1 && can(min) === can('CCc1ccc(S(=O)(=O)O)cc1')
        && calc.molecularFromEmpirical(calc.empiricalFromCombustion(6.1, 17.6, 4.5), { mw: 122 }) === 'C8H10O';
    }],
    ['2003 III 六員環の酢酸エノールエステル C₉H₁₄O₂: 加水分解で 2-メチルシクロヘキサノン（不斉炭素あり）とシクロヘキサンカルボアルデヒド（フェーリング陽性）', () => {
      const C = 'CC(=O)OC1=C(C)CCCC1', D = 'CC(=O)OC1=CCCCC1C', E = 'CC(=O)OC=C1CCCCC1';
      return [C, D, E].every((x) => chem.formula(chem.graphFromSmiles(R, x)) === 'C9H14O2')
        && same(ev('hydrolysis', C), ['CC(=O)O', 'CC1CCCCC1=O']) && same(ev('hydrolysis', D), ['CC(=O)O', 'CC1CCCCC1=O'])
        && ev('chiral', C) === 0 && ev('chiral', D) === 1 && ev('chiral', 'CC1CCCCC1=O') === 1
        && same(ev('hydrolysis', E), ['CC(=O)O', 'O=CC1CCCCC1']) && ev('fehling', 'O=CC1CCCCC1') === true && ev('chiral', 'O=CC1CCCCC1') === 0;
    }],
    ['2003 III 問3 六員環と C=C をもつ C₇H₁₂ は 4 種類', () => enumerate(R, 'C7H12', ring).filter((x) => hasMatch(R, x, '[r6]') && hasMatch(R, x, 'C=C')).length === 4],
    ['2002 III C₅H₁₀ の鎖状異性体 5・環状異性体 5', () => enumerate(R, 'C5H10', { rings: 0 }).length === 5 && enumerate(R, 'C5H10', { rings: 1, maxRing: 5 }).filter((x) => hasMatch(R, x, '[R]')).length === 5],
    ['2000 III 質量百分率 C 85.94%・H 11.94% → C₅₄H₉₀O、100 g に I₂ 337 g 付加 → C=C 10 個', () => calc.empiricalFromPercent({ C: 85.94, H: 11.94 }) === 'C54H90O'
      && Math.round((337 / 253.8) / (100 / calc.mass('C54H90O'))) === 10],
  );

  // エンジンでは扱えない出題要素（次に伸ばす候補）
  const UNSUPPORTED = [
    '2026 III アセタール交換の平衡・光学分割の速度（アセタールの規則が未実装）',
    '2026 IV メチル化分析・糖の立体（糖の立体を区別する表現が未実装）',
    '2019 IV アルドースの立体・メソ体（フィッシャー投影が未実装）',
    '2019 III(b) イミドの穏やかな加水分解（位置選択の規則が未実装）',
    '2010 III(a) 第一級 OH の選択的アセチル化（与えられた規則が未実装）',
    '2025 III 小員環のひずみによる水素化での開環（規則が未実装）',
    '2024 III(b)・2002 III 臭素のアンチ付加・水素のシン付加、シクロプロパンの立体（三次元の立体配置が未実装）',
    '2020 IV・2026 IV アセトンによるジオールの保護（アセタールの規則が未実装）',
    '2001 III ジアゾカップリング（反応規則が未実装）',
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
