'use strict';
// 検証を通った問題だけを収録し、構造式SVG・候補ごとの判定結果を前計算して
// 単一ファイルの dist/index.html を作る。ブラウザ側は化学ライブラリを読み込まない。
// 使い方: node scripts/build.js
const fs = require('fs');
const path = require('path');
const chem = require('../src/chem');
const { loadRDKit, loadProblems, checkAny } = require('./validate');

const ROOT = path.join(__dirname, '..');

function carbonLabels(RDKit, smiles) {
  // 高校の構造式に寄せて、炭素に CH3 / CH2 / CH を明記する
  const g = chem.graphFromSmiles(RDKit, smiles);
  const labels = {};
  g.atoms.forEach((a, i) => {
    if (a.el !== 'C' || a.arom) return; // ベンゼン環は骨格式のまま
    const h = chem.hCount(g, i);
    labels[i] = h === 0 ? 'C' : h === 1 ? 'CH' : `CH<sub>${h}</sub>`;
  });
  return labels;
}

function drawSvg(RDKit, smiles) {
  const mol = RDKit.get_mol(chem.expand(smiles));
  try {
    mol.set_new_coords && mol.set_new_coords(true);
    const details = {
      width: 240,
      height: 170,
      bondLineWidth: 1.6,
      fixedBondLength: 38,
      minFontSize: 13,
      clearBackground: false,
      atomLabels: carbonLabels(RDKit, smiles),
    };
    let svg = mol.get_svg_with_highlights(JSON.stringify(details));
    svg = svg.replace(/<\?xml[^>]*>\s*/, '').replace(/<!-- END OF HEADER -->\s*/, '');
    return compactSvg(svg.replace(/\n/g, ''));
  } finally {
    mol.delete();
  }
}

// RDKit の SVG は線ごとに同じ style を繰り返すので、見た目を変えずに縮める（問題が増えてもファイルを軽く保つ）
// 線の fill:none と太さは svg 要素に置いて継承させ、文字の path は自分の fill を持つのでそのまま
function compactSvg(svg) {
  const short = (c) => c.replace(/^#([0-9A-F])\1([0-9A-F])\2([0-9A-F])\3$/i, '#$1$2$3');
  return svg
    .replace(/<svg[^>]*viewBox='([^']*)'[^>]*>/, "<svg xmlns='http://www.w3.org/2000/svg' viewBox='$1' fill='none' stroke-width='1.6'>")
    .replace(/<rect[^>]*\/>/g, '')
    .replace(/ class='[^']*'/g, '')
    .replace(/ style='fill:none;(?:fill-rule:evenodd;)?stroke:(#[0-9A-F]{6});stroke-width:1\.6px;[^']*'/gi, (m, c) => ` stroke='${short(c)}'`)
    .replace(/ fill='(#[0-9A-F]{6})'/gi, (m, c) => ` fill='${short(c)}'`)
    .replace(/ d='([^']*)'/g, (m, d) => ` d='${d.replace(/, /g, ' ').replace(/ ?([MLQZC]) /g, '$1').trim()}'`)
    .replace(/\s*\/>/g, '/>');
}

// 正解を1つに決めるのに最低限必要な手がかりの枚数（総当たり）
function minCards(r, n) {
  const ai = r.cands.indexOf(r.answer);
  for (let k = 1; k <= n; k++) {
    for (let mask = 1; mask < 1 << n; mask++) {
      const pick = [...Array(n).keys()].filter((i) => mask & (1 << i));
      if (pick.length !== k) continue;
      const left = r.cands.map((_, c) => c).filter((c) => pick.every((ci) => r.ok(ci, c)));
      if (left.length === 1 && left[0] === ai) return k;
    }
  }
  return n;
}

// 問題データ中の "pep:Gly-Ala" 表記を集めて、配列をそのまま名前にする
function collectPeptideNames(RDKit, node, names) {
  if (typeof node === 'string') {
    if (node.startsWith('pep:')) names[chem.canonical(RDKit, node)] = node.slice(4);
  } else if (node && typeof node === 'object') {
    Object.values(node).forEach((v) => collectPeptideNames(RDKit, v, names));
  }
}

function collectSmiles(RDKit, p, set) {
  const add = (x) => {
    if (typeof x === 'string') { try { set.add(chem.canonical(RDKit, x)); } catch (e) { /* 構造でない値 */ } }
    else if (Array.isArray(x)) x.forEach(add);
  };
  add(p.answer);
  const walk = (stage) => {
    add(stage.answer);
    add(stage.candidates || []);
    (stage.clues || []).forEach((c) => { if (typeof c.result === 'string' || Array.isArray(c.result)) add(c.result); });
  };
  (p.fragments || []).forEach(walk);
  if (p.assemble) walk(p.assemble);
  // 候補に結果表で出てくる生成物も含める
  (p.fragments || []).concat(p.assemble ? [p.assemble] : []).forEach((st) => (st.candidates || []).forEach((c) => {
    (st.clues || []).forEach((cl) => {
      try { add(chem.evaluate(RDKit, cl.card, c)); } catch (e) { /* 使えない */ }
    });
  }));
}

// 22600 → 2.26×10⁴（入試の表記）
function sci(x) {
  const e = Math.floor(Math.log10(x));
  const sup = String(e).split('').map((d) => '⁰¹²³⁴⁵⁶⁷⁸⁹'[+d]).join('');
  return `${(x / 10 ** e).toFixed(2)}×10${sup}`;
}

async function main() {
  const RDKit = await loadRDKit();
  const problems = loadProblems();
  const rawNames = JSON.parse(fs.readFileSync(path.join(ROOT, 'problems', 'names.json'), 'utf8'));
  const names = {};
  for (const [s, n] of Object.entries(rawNames)) names[chem.canonical(RDKit, s)] = n;
  collectPeptideNames(RDKit, problems, names);

  const molecules = {};
  const addMol = (s) => {
    if (!molecules[s]) {
      molecules[s] = { name: names[s] || null, svg: drawSvg(RDKit, s) };
      // 名前のない構造（自動生成の断片など）は分子式で呼ぶ
      if (!molecules[s].name) { try { molecules[s].formula = chem.formula(chem.graphFromSmiles(RDKit, chem.expand(s))); } catch (e) { /* なし */ } }
    }
    return s;
  };

  // 絞り込みの結果表を画面用に整える。contains 型は「得られる/得られない」だけ持つ
  const shapeClues = (clues, r) => clues.map((c, ci) => {
    const kind = chem.CARDS[c.card].kind;
    if (kind === 'products') {
      r.expected[ci].forEach(addMol);
      r.table[ci].forEach((list) => list.forEach(addMol));
    }
    if (kind === 'contains') {
      addMol(r.expected[ci]);
      return { card: c.card, result: r.expected[ci], byCandidate: r.table[ci].map((v) => (chem.consistent(c.card, v, r.expected[ci]) ? r.expected[ci] : false)) };
    }
    return { card: c.card, result: r.expected[ci], byCandidate: r.table[ci] };
  });
  const narrowData = (r, clues, extra) => {
    r.cands.forEach(addMol);
    return { ...extra, candidates: r.cands, answer: r.cands.indexOf(r.answer), minCards: minCards(r, clues.length), clues: shapeClues(clues, r) };
  };

  const out = [];
  const rejected = [];
  for (const p of problems) {
    const r = checkAny(RDKit, p);
    if (r.errors.length) {
      rejected.push(...r.errors);
      continue; // 検証に落ちた問題は収録しない
    }
    const base = { id: p.id, mode: p.mode, level: p.level || 1, formula: p.formula };
    if (p.mode === 'narrow') {
      out.push(narrowData(r, p.clues, base));
    } else if (p.mode === 'big') {
      addMol(r.X);
      // 自動生成の大問: 水素付加量・燃焼分析の計算段階を先に解く（分解の表で断片の分子式が見える前に）
      const stages = r.calcs.map((c) => ({ type: 'calc', key: c.key, prompt: c.prompt, answer: c.answer, choices: c.choices, unit: c.unit, explain: c.explain }));
      stages.push({
        type: 'split',
        bonds: r.bonds,
        products: r.frags.map((f) => ({ label: f.label, formula: f.formula, count: f.count })),
      });
      for (const f of r.frags) {
        addMol(f.answer);
        if (f.given) stages.push({ type: 'given', label: f.label, smiles: f.answer, note: f.note });
        else stages.push(narrowData(f.r, f.clues, { type: 'narrow', label: f.label, id: `${p.id}/${f.label}`, formula: f.formula, level: p.level }));
      }
      if (r.assemble) {
        stages.push(narrowData(r.assemble.r, r.assemble.clues, { type: 'narrow', label: 'X', id: `${p.id}/X`, formula: p.formula, level: p.level }));
      }
      if (p.generated) {
        // 京大の過去問の再現と自動生成は「京大レベル」モードにまとめる
        out.push({ ...base, mode: 'gen', title: p.title, story: p.story, answerSmiles: r.X, stages,
          difficulty: p.meta.difficulty, kyoto: p.meta.kyoto || null, band: p.meta.band, seed: p.meta.seed });
      } else out.push({ ...base, title: p.title, story: p.story, answerSmiles: r.X, stages });
    } else if (p.mode === 'polymer') {
      const unit = addMol(chem.canonical(RDKit, p.unit));
      const choices = (a) => [...new Set([Math.round(a / 2), a, a * 2, a * 4])].sort((x, y) => x - y);
      const stages = [
        { type: 'calc', key: 'dp', prompt: '平均重合度 n はおよそいくつか', answer: p.n, choices: choices(p.n),
          explain: `繰り返し単位 ${r.unitFormula} の式量は ${r.unitMass}。${r.unitMass} × n = ${r.mw} より n = ${p.n}` },
        { type: 'calc', key: 'per', prompt: p.perUnit.prompt, answer: p.perUnit.count * p.n, unit: p.perUnit.unit,
          choices: choices(p.perUnit.count * p.n),
          explain: `繰り返し単位 1 つあたり ${p.perUnit.count}、それが n = ${p.n} 個で ${p.perUnit.count * p.n}` },
      ];
      for (const m of r.mons) {
        addMol(m.answer);
        if (m.given) stages.push({ type: 'given', label: m.label, smiles: m.answer, note: m.note });
        else stages.push(narrowData(m.r, m.clues, { type: 'narrow', label: m.label, id: `${p.id}/${m.label}`, formula: m.formula, level: p.level }));
      }
      out.push({ ...base, formula: r.unitFormula, title: p.title, story: p.story.replace('{M}', sci(r.mw)), answerSmiles: unit,
        unitFormula: r.unitFormula, stages });
    } else if (p.mode === 'count') {
      r.pool.forEach(addMol);
      out.push({
        ...base,
        scope: p.scope,
        pool: r.pool,
        clues: p.clues.map((c, ci) => ({ card: c.card, result: r.expected[ci], byCandidate: r.table[ci] })),
        answers: r.answers,
        stereo: r.stereo,
        nStereo: r.nStereo,
      });
    }
  }

  // 自動生成の断片には名前のないものがある（構造式だけを見せる）。手で作った問題だけ名前の抜けを警告する
  const genMols = new Set();
  problems.filter((p) => p.generated).forEach((p) => collectSmiles(RDKit, p, genMols));
  const missing = Object.keys(molecules).filter((s) => !molecules[s].name && !genMols.has(s));
  if (missing.length) console.warn(`name missing for: ${missing.join(' ')}`);

  const cards = {};
  for (const [k, v] of Object.entries(chem.CARDS)) {
    cards[k] = { name: v.name, action: v.action, kind: v.kind, yes: v.yes, no: v.no, none: v.none };
  }
  const data = { version: 2, cards, molecules, problems: out };

  const template = fs.readFileSync(path.join(ROOT, 'src', 'index.html'), 'utf8');
  const json = JSON.stringify(data).replace(/</g, '\\u003c');
  const html = template.replace('/*__PUZZLE_DATA__*/null', () => json);
  if (html === template) throw new Error('data placeholder not found in src/index.html');
  fs.mkdirSync(path.join(ROOT, 'dist'), { recursive: true });
  // index.html: そのまま開ける完全な文書 / artifact.html: claude.ai Artifact 用（外枠は公開時に付く）
  fs.writeFileSync(
    path.join(ROOT, 'dist', 'index.html'),
    `<!doctype html>\n<html lang="ja">\n<head>\n<meta charset="utf-8">\n${html}\n</html>\n`,
  );
  fs.writeFileSync(path.join(ROOT, 'dist', 'artifact.html'), html);

  console.log(`built dist/index.html: ${out.length} problems, ${Object.keys(molecules).length} structures, ${(html.length / 1024).toFixed(0)} KB`);
  if (rejected.length) {
    console.error(`rejected ${rejected.length} error(s):`);
    rejected.forEach((e) => console.error('  ' + e));
    process.exit(1);
  }
}

main();
