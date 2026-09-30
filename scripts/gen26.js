'use strict';
// 京大2026型の大問を作り、検証パイプラインを全部通ったものだけを problems/k26.json に書く。
// 使い方: node scripts/gen26.js [種] [問題数] [難易度（1〜5。省略すると 1〜5 を順に回す）]
//  問題ごとに パラメータ（難易度・テーマ・乱数の種）を params に記録し、同じ種と難易度で同じ問題が再現できる
const fs = require('fs');
const path = require('path');
const { execFileSync } = require('child_process');
const K = require('../src/kyoto26');
const { loadRDKit, nameMap } = require('./validate');
const { checkK26, calcsOf } = require('./k26check');

const ROOT = path.join(__dirname, '..');

function todayJST() {
  return new Date(Date.now() + 9 * 3600 * 1000).toISOString().slice(0, 10).replace(/-/g, '');
}

// Python で数値検算（calcs をまとめて渡す）
function pyCheck(calcs) {
  const tmp = path.join(ROOT, 'generated', 'calcs-k26-tmp.json');
  fs.writeFileSync(tmp, JSON.stringify(calcs));
  try {
    const out = execFileSync('python3', [path.join(__dirname, 'numcheck.py'), tmp, '--json'], { encoding: 'utf8' });
    return JSON.parse(out);
  } catch (e) {
    if (e.stdout) return JSON.parse(e.stdout);
    throw e;
  } finally { fs.rmSync(tmp, { force: true }); }
}

async function main() {
  const base = Number(process.argv[2] || todayJST());
  const N = Number(process.argv[3]) || 10;
  if (process.argv[4] && !K.LEVELS[process.argv[4]]) { console.error('難易度は 1〜5'); process.exit(1); }
  const RDKit = await loadRDKit();
  const names = nameMap(RDKit);
  const out = [];
  const log = [];
  for (let i = 0; i < N; i++) {
    const level = process.argv[4] ? Number(process.argv[4]) : (i % 5) + 1;
    let got = null;
    let tries = 0;
    for (let k = 0; k < 40 && !got; k++) {
      const seed = (base % 100000) * 100 + i * 7 + k * 1009;
      tries++;
      const P = K.buildK26(RDKit, seed, level);
      if (!P) continue;
      const r = checkK26(RDKit, P, names);
      if (r.errors.length) { log.push(...r.errors); continue; }
      const py = pyCheck(calcsOf(P));
      if (py.errors.length) { log.push(...py.errors); continue; }
      // 再現性: 記録したパラメータ（種・難易度）から作り直して同じ問題になるか
      const again = K.buildK26(RDKit, P.params.seed, P.params.difficulty);
      if (!again || JSON.stringify(again) !== JSON.stringify(P)) { log.push(`${P.id}: 同じパラメータで再現できない`); continue; }
      P.validation = { ...r.validation, numeric_checked: true, reproducible: true, rubric: P.questions.filter((q) => q.type === 'essay').every((q) => q.rubric && q.rubric.length >= 2) };
      delete P.hidden;
      got = P;
    }
    if (!got) { console.error(`難易度 ${level}: ${tries} 回試して作れなかった`); process.exit(1); }
    out.push(got);
    console.log(`  ${got.id}  難易度 ${level}  規則の適用 ${got.params.rule_applications}・優先 ${got.params.rule_priorities}  逆算 D ${got.params.inverse_steps.D}/G ${got.params.inverse_steps.G}  計算 ${got.params.calc_steps} 段${got.params.stereo ? '・立体' : ''}（${tries} 回目）`);
  }
  fs.writeFileSync(path.join(ROOT, 'problems', 'k26.json'), JSON.stringify(out, null, 1) + '\n');
  console.log(`problems/k26.json: ${out.length} 問（検証で落とした候補 ${log.length} 件）`);
  if (process.env.VERBOSE) log.forEach((m) => console.log('   ' + m));
}

if (require.main === module) main();
