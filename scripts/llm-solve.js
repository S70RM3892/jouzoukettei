'use strict';
// 別インスタンス解答（任意）: 答えを取り除いた問題（publicView）を文章にして Claude に解かせ、想定解と一致するか確かめる。
// 決まった手順の別解答（src/solve26.js）は npm test で毎回通す。こちらは API の認証情報があるときに手で走らせる。
// 使い方: node scripts/llm-solve.js [問題の id ...]   （結果は generated/llm-check.json）
// 認証: ANTHROPIC_API_KEY など SDK が読む設定。なければ何もせずに終わる
const fs = require('fs');
const path = require('path');
const Anthropic = require('@anthropic-ai/sdk');
const chem = require('../src/chem');
const S = require('../src/solve26');
const { loadRDKit, loadK26 } = require('./validate');

const ROOT = path.join(__dirname, '..');
const MODEL = 'claude-opus-5-5';

// 画面に出るのと同じ情報だけを文章にする（構造式の図は SMILES で渡す）
function problemText(V) {
  const out = [];
  for (const pt of V.parts) {
    out.push(`■ ${pt.label}`);
    out.push(pt.intro);
    if (pt.equation) out.push(`式1: ${pt.equation.lhs.join(' + ')} → ${pt.equation.rhs.join(' + ')}（分子量 ${Object.entries(pt.equation.masses).map(([k, v]) => `${k} = ${v}`).join('、')}）`);
    pt.rules.forEach((r) => out.push(`[${r.id}] ${r.title || ''}: ${r.text}${r.scheme ? `（図: ${r.scheme.join(' ')}）` : ''}`));
    if (pt.limitation) out.push(pt.limitation);
    (pt.figs || []).forEach((f) => out.push(`化合物 ${f.label} の構造（SMILES）: ${f.smiles}`));
    if (pt.tail) out.push(pt.tail);
    if (pt.expLead) out.push(pt.expLead);
    V.experiments.filter((e) => e.part === pt.label).forEach((e) => out.push(`${e.label} ${e.text}`));
    if (pt.tailY) out.push(pt.tailY);
  }
  out.push('原子量: H = 1.0, C = 12, O = 16。');
  out.push('');
  out.push('【設問】');
  V.questions.filter((q) => q.type !== 'essay').forEach((q) => out.push(`${q.id}: ${q.prompt}${q.notes && q.notes.length ? `（注: ${q.notes.join('／')}）` : ''}`));
  return out.join('\n');
}

const SCHEMA = {
  type: 'object',
  properties: {
    answers: {
      type: 'array',
      items: {
        type: 'object',
        properties: {
          id: { type: 'string' },
          smiles: { type: 'array', items: { type: 'string' } },
          value: { type: ['number', 'null'] },
        },
        required: ['id', 'smiles', 'value'],
        additionalProperties: false,
      },
    },
  },
  required: ['answers'],
  additionalProperties: false,
};

async function main() {
  const client = new Anthropic();
  const RDKit = await loadRDKit();
  const only = process.argv.slice(2);
  const results = [];
  for (const P of loadK26().filter((p) => !only.length || only.includes(p.id))) {
    const V = S.publicView(P);
    let msg;
    try {
      msg = await client.beta.messages.stream({
        model: MODEL,
        max_tokens: 64000,
        betas: ['server-side-fallback-2026-07-01'],
        fallbacks: 'default',
        output_config: { effort: 'high', format: { type: 'json_schema', schema: SCHEMA } },
        system: '高校化学（京都大学の入試）の有機化学の問題を解く。問題文の規則と事実だけを使う。構造式は SMILES（立体は書かない）で、数値は設問の有効数字で答える。構造の設問は value を null に、数値の設問は smiles を空の配列にする。',
        messages: [{ role: 'user', content: problemText(V) }],
      }).finalMessage();
    } catch (e) {
      if (e instanceof Anthropic.AuthenticationError || !(e instanceof Anthropic.APIError)) { console.log(`API を使えないので、別インスタンス解答（Claude）は飛ばす（${e.message}）`); return; }
      console.error(`${P.id}: API エラー ${e.status}: ${e.message}`);
      continue;
    }
    if (msg.stop_reason === 'refusal') { console.error(`${P.id}: 回答が拒否された`); continue; }
    const text = msg.content.filter((b) => b.type === 'text').map((b) => b.text).join('');
    let ans;
    try { ans = JSON.parse(text).answers; } catch (e) { console.error(`${P.id}: JSON を読めない`); continue; }
    const sol = {};
    for (const a of ans) {
      try { sol[a.id] = a.value === null ? { smiles: a.smiles.map((s) => chem.canonical(RDKit, s)) } : { value: a.value }; } catch (e) { sol[a.id] = { smiles: [] }; }
    }
    const errs = S.agree(RDKit, P, sol);
    results.push({ id: P.id, model: msg.model, agree: !errs.length, errors: errs });
    console.log(`${P.id}: ${errs.length ? `不一致 ${errs.length}` : '一致'}${errs.length ? '\n  ' + errs.join('\n  ') : ''}`);
  }
  if (results.length) fs.writeFileSync(path.join(ROOT, 'generated', 'llm-check.json'), JSON.stringify(results, null, 1) + '\n');
}

main();
