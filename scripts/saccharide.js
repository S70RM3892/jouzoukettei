'use strict';
// 糖の構造決定の問題を作り、problems/sugar.json に書き出す。
// 二糖・三糖の候補（単糖の種類・α/β・結合の位置の組み合わせ）から、還元性・加水分解の生成物・メチル化分析・酵素（基質特異性）の
// 手がかりを2枚以上組み合わせて1つに決める。糖は立体でしか区別できないので、候補は sac: の略記（src/sugar.js）で書く。
// 使い方: node scripts/saccharide.js [種]
const fs = require('fs');
const path = require('path');
const chem = require('../src/chem');
const G = require('../src/generator');
const { loadRDKit, checkNarrow } = require('./validate');

const ROOT = path.join(__dirname, '..');
const CARDS = ['fehling', 'hydrolysis', 'methylation_analysis', 'maltase', 'invertase', 'lactase', 'cellobiase'];
const WEIGHTS = { fehling: 1.0, hydrolysis: 0.9, methylation_analysis: 1.2, maltase: 0.8, invertase: 0.8, lactase: 0.8, cellobiase: 0.8 };

// 二糖: 非還元末端（Glc・Gal の α/β）× 還元末端グルコースの 4 位・6 位、と還元性を示す炭素どうしの結合（トレハロース・スクロース）、フルクトースが還元末端のもの
const DI = [];
for (const u of ['Glc', 'Gal']) for (const a of ['a', 'b']) for (const n of [4, 6]) DI.push(`sac:${u}(${a}1-${n})Glc`);
DI.push('sac:Glc(a1-1a)Glc', 'sac:Glc(a1-2b)Fru', 'sac:Gal(b1-4)Fru', 'sac:Glc(a1-6)Fru');
// 三糖: マルトトリオース・パノース・枝分かれ・ラフィノースなど
const TRI = ['sac:Glc(a1-4)Glc(a1-4)Glc', 'sac:Glc(a1-6)Glc(a1-4)Glc', 'sac:Glc(a1-4)Glc(a1-6)Glc', 'sac:Glc(a1-4)[Glc(a1-6)]Glc',
  'sac:Gal(a1-6)Glc(a1-2b)Fru', 'sac:Glc(b1-4)Glc(b1-4)Glc', 'sac:Gal(b1-4)Glc(a1-2b)Fru', 'sac:Glc(a1-4)Glc(a1-2b)Fru', 'sac:Gal(b1-4)[Glc(a1-6)]Glc'];

async function main() {
  const RDKit = await loadRDKit();
  const r = G.rng(parseInt(process.argv[2] || '20260930', 10));
  const out = [];
  for (const [label, raw] of [['二糖', DI], ['三糖', TRI]]) {
    const pool = raw.map((s) => chem.canonical(RDKit, s));
    const table = G.valueTable(RDKit, pool, CARDS);
    for (const answer of pool) {
      const a = pool.indexOf(answer);
      const set = G.chooseClues(r, table, a, WEIGHTS, { lo: 2, hi: 3 }) || G.chooseClues(r, table, a, WEIGHTS, { lo: 2, hi: 4 });
      if (!set) continue;
      const formula = chem.formula(chem.graphFromSmiles(RDKit, answer));
      const p = {
        id: `sugar-${String(out.length + 1).padStart(3, '0')}`, mode: 'narrow', kind: 'sugar', title: `${label}の構造`,
        level: label === '二糖' ? (set.length <= 2 ? 1 : 2) : 3, formula, answer, candidates: pool,
        clues: set.map((card) => ({ card, result: chem.evaluate(RDKit, card, answer) })),
      };
      const v = checkNarrow(RDKit, p);
      if (v.errors.length) { console.log('  落とした: ' + v.errors[0]); continue; }
      out.push(p);
    }
  }
  fs.writeFileSync(path.join(ROOT, 'problems', 'sugar.json'), JSON.stringify(out, null, 1) + '\n');
  console.log(`problems/sugar.json に ${out.length} 問（候補 二糖 ${DI.length}・三糖 ${TRI.length}）`);
  out.forEach((p) => console.log(`  ${p.id} ★${p.level} ${p.answer.padEnd(30)} ${p.clues.map((c) => c.card).join(', ')}`));
}

main();
