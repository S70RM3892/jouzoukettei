'use strict';
// 高校有機化学の範囲（data/curriculum.json）をアプリがどこまで扱えているかを確かめる。
//  1. すべての節に知識確認の問題（problems/knowledge.json）が 10 問以上ある（計算問題 problems/calc.json は数えない）
//  2. 知識確認の問題の形が正しい（節が存在する・正解と誤答が重ならない・選択肢が4つ）
//  3. 節に書いた判定カードが実在し、判定カードはどれかの節に入っている
//  4. 構造決定などの問題で練習できない節（知識確認だけの節）を一覧にする
const fs = require('fs');
const path = require('path');
const chem = require('../src/chem');

const ROOT = path.join(__dirname, '..');
const C = JSON.parse(fs.readFileSync(path.join(ROOT, 'data', 'curriculum.json'), 'utf8'));
const K = ['knowledge.json', 'calc.json'].flatMap((f) => JSON.parse(fs.readFileSync(path.join(ROOT, 'problems', f), 'utf8')));
const MIN = 10;

const errors = [];
const topics = C.sections.flatMap((s) => s.topics.map((t) => ({ ...t, section: s.name })));
const ids = new Set(topics.map((t) => t.id));
if (ids.size !== topics.length) errors.push('節の id が重複している');

const seen = new Set();
for (const k of K) {
  if (seen.has(k.id)) errors.push(`${k.id}: id が重複`);
  seen.add(k.id);
  if (!ids.has(k.topic)) errors.push(`${k.id}: 節 ${k.topic} がない`);
  if (!k.q || !k.a || !k.explain) errors.push(`${k.id}: 問い・正解・解説のどれかが空`);
  if (!Array.isArray(k.wrong) || k.wrong.length !== 3) errors.push(`${k.id}: 誤答は3つ`);
  if (!k.calc && K.some((x) => x !== k && !x.calc && x.q === k.q)) errors.push(`${k.id}: 同じ問いがある`);
  const all = [k.a, ...(k.wrong || [])];
  if (new Set(all).size !== all.length) errors.push(`${k.id}: 選択肢が重なっている`);
}
const cardsUsed = new Set();
for (const t of topics) {
  const n = K.filter((k) => k.topic === t.id && !k.calc).length;
  if (n < MIN) errors.push(`${t.name}: 知識確認が ${n} 問（${MIN} 問以上いる）`);
  for (const c of t.engine) {
    if (!chem.CARDS[c]) errors.push(`${t.name}: 判定カード ${c} がない`);
    cardsUsed.add(c);
  }
}
const orphan = Object.keys(chem.CARDS).filter((c) => !cardsUsed.has(c));
// 与えられた規則のカード（京大の問題文で与えられる反応）は教科書の節に入れない
const RULES = new Set(require('../src/inference').RULE_CARDS.concat(['partial_hydrolysis', 'dehydration_ozonolysis', 'bromine_water']));
orphan.filter((c) => !RULES.has(c)).forEach((c) => errors.push(`判定カード ${c}（${chem.CARDS[c].name}）がどの節にも入っていない`));

const knowOnly = topics.filter((t) => !t.engine.length && !t.modes.length);
console.log(`範囲: ${C.sections.length} 章 ${topics.length} 節、知識確認 ${K.filter((k) => !k.calc).length} 問・計算 ${K.filter((k) => k.calc).length} 問`);
C.sections.forEach((s) => console.log(`  ${s.name}: ${s.topics.map((t) => `${t.name}(${K.filter((k) => k.topic === t.id).length}${t.engine.length || t.modes.length ? '' : '・知識のみ'})`).join('、')}`));
console.log(`構造決定などの問題では練習できず、知識確認だけで扱う節 ${knowOnly.length}: ${knowOnly.map((t) => t.name).join('、')}`);
if (errors.length) {
  errors.forEach((e) => console.error('NG ' + e));
  process.exit(1);
}
console.log('OK: すべての節に知識確認がある');
