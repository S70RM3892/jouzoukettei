'use strict';
// ブラウザ用の SMILES 読み取り・比較（src/smiles-lite.js）が、RDKit の正規化と同じ判定をするか確かめる。
// 構造式の記述の採点に使うので、同じ構造は「同じ」、構造異性体は「違う」と判定できなければならない
const L = require('../src/smiles-lite');
const chem = require('../src/chem');
const { loadRDKit } = require('./validate');

(async () => {
  const RDKit = await loadRDKit();
  let bad = 0;
  const same = [
    ['c1ccccc1O', 'Oc1ccccc1'], ['C1=CC=CC=C1O', 'Oc1ccccc1'], ['CC(=O)O', 'OC(C)=O'], ['O=[N+]([O-])c1ccccc1', 'O=N(=O)c1ccccc1'],
    ['COC1CCCCO1', 'O1CCCCC1OC'], ['C[C@H](O)CC', 'CC(O)CC'], ['CCOC(C)OCC', 'C(C)(OCC)OCC'], ['C/C=C/C', 'CC=CC'],
  ];
  const diff = [['CCO', 'COC'], ['CC(O)CC', 'CCCCO'], ['COC1CCCCO1', 'COC1CCCO1C'], ['Cc1ccccc1O', 'Cc1ccc(O)cc1'], ['CC=CC', 'C=CCC']];
  for (const [a, b] of same) if (!L.same(L.parse(a), L.parse(b))) { bad++; console.error(`同じはず: ${a} ${b}`); }
  for (const [a, b] of diff) if (L.same(L.parse(a), L.parse(b))) { bad++; console.error(`違うはず: ${a} ${b}`); }
  // 問題に出る構造（京大2026型の答え・図）すべてで、RDKit の正規 SMILES を読んで分子式が合い、自分自身と同じと判定できる
  const fs = require('fs');
  const path = require('path');
  const file = path.join(__dirname, '..', 'problems', 'k26.json');
  const k26 = fs.existsSync(file) ? JSON.parse(fs.readFileSync(file, 'utf8')) : [];
  const all = new Set();
  k26.forEach((P) => {
    P.questions.filter((q) => q.type === 'draw').forEach((q) => q.answer.smiles.forEach((s) => all.add(s)));
    P.parts.forEach((pt) => (pt.figs || []).forEach((f) => all.add(f.smiles)));
  });
  for (const s of all) {
    const g = L.parse(s);
    const f = chem.formula(chem.graphFromSmiles(RDKit, s));
    if (L.formula(g) !== f) { bad++; console.error(`分子式が違う: ${s} ${L.formula(g)} ≠ ${f}`); }
    // RDKit で原子の順番を変えた SMILES（ランダム）と同じと判定できる
    const m = RDKit.get_mol(s);
    const alt = m.get_smiles(JSON.stringify({ doRandom: true, rootedAtAtom: Math.max(0, g.atoms.length - 1) }));
    m.delete();
    if (!L.same(g, L.parse(alt))) { bad++; console.error(`同じと判定できない: ${s} / ${alt}`); }
  }
  // 答えどうし（違う構造）は違うと判定する
  const list = [...all];
  for (let i = 0; i < list.length; i++) for (let j = i + 1; j < list.length; j++) {
    if (L.same(L.parse(list[i]), L.parse(list[j])) !== (chem.canonical(RDKit, list[i]) === chem.canonical(RDKit, list[j]))) { bad++; console.error(`判定が RDKit と違う: ${list[i]} ${list[j]}`); }
  }
  if (bad) { console.error(`NG: SMILES の判定 ${bad} 件`); process.exit(1); }
  console.log(`OK: SMILES の読み取りと比較（${same.length + diff.length} 組 + 問題の構造 ${all.size} 個）`);
})();
