'use strict';
// 異性体列挙の検算。教科書・過去問で知られている数と一致するか確かめる。
const { enumerate, hasMatch } = require('../src/enumerate');
const { loadRDKit } = require('./validate');

async function main() {
  const R = await loadRDKit();
  const acyclic = { rings: 0 };
  const benz = { seed: 'c1ccccc1', rings: 0 };
  const has = (sm) => (s) => hasMatch(R, s, sm);
  const not = (sm) => (s) => !hasMatch(R, s, sm);
  const CASES = [
    ['C4H10O', {}, 7, 'C₄H₁₀O（アルコール4＋エーテル3）'],
    ['C5H12O', {}, 14, 'C₅H₁₂O（アルコール8＋エーテル6）'],
    ['C5H12O', { filter: has('[OX2H1]') }, 8, 'C₅H₁₂O のアルコール'],
    ['C4H8', { maxRing: 6 }, 5, 'C₄H₈（環式を含む）'],
    ['C5H10', {}, 10, 'C₅H₁₀（環式を含む）'],
    ['C4H8O2', { ...acyclic, filter: has('[CX3](=O)[OX2]') }, 6, 'C₄H₈O₂ のカルボン酸・エステル'],
    ['C5H10O2', { ...acyclic, filter: has('[CX3](=O)[OX2]') }, 13, 'C₅H₁₀O₂ のカルボン酸・エステル'],
    ['C5H10O', { ...acyclic, filter: has('[CX3]=O') }, 7, 'C₅H₁₀O のカルボニル化合物'],
    ['C4H8O', acyclic, 11, '鎖状の C₄H₈O（京大2008：エノールを除く）'],
    ['C8H10', benz, 4, 'C₈H₁₀ の芳香族（京大2007）'],
    ['C9H12', benz, 8, 'C₉H₁₂ の芳香族（京大2019）'],
    ['C8H10O', benz, 19, 'C₈H₁₀O のベンゼン環をもつ化合物'],
    ['C7H8O', benz, 5, 'C₇H₈O の芳香族（クレゾール3・ベンジルアルコール・アニソール）'],
    ['C4H11N', { ...acyclic }, 8, 'C₄H₁₁N のアミン（京大2006）'],
  ];
  let fail = 0;
  for (const [f, o, want, label] of CASES) {
    const t = Date.now();
    const got = enumerate(R, f, o);
    const ok = got.length === want;
    if (!ok) fail++;
    console.log(`${ok ? 'OK' : 'NG'} ${label}: ${got.length}（期待 ${want}） ${Date.now() - t}ms${ok ? '' : '\n   ' + got.join(' ')}`);
  }
  if (fail) process.exit(1);
}
main();
