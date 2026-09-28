'use strict';
// 検証を通った問題だけを収録し、構造式SVG・候補ごとの判定結果を前計算して
// 単一ファイルの dist/index.html を作る。ブラウザ側は化学ライブラリを読み込まない。
// 使い方: node scripts/build.js
const fs = require('fs');
const path = require('path');
const chem = require('../src/chem');
const { loadRDKit, checkProblem } = require('./validate');

const ROOT = path.join(__dirname, '..');

function carbonLabels(RDKit, smiles) {
  // 高校の構造式に寄せて、炭素に CH3 / CH2 / CH を明記する
  const g = chem.graphFromSmiles(RDKit, smiles);
  const labels = {};
  g.atoms.forEach((a, i) => {
    if (a.el !== 'C') return;
    const h = chem.hCount(g, i);
    labels[i] = h === 0 ? 'C' : h === 1 ? 'CH' : `CH<sub>${h}</sub>`;
  });
  return labels;
}

function drawSvg(RDKit, smiles) {
  const mol = RDKit.get_mol(smiles);
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
    return svg.replace(/\n/g, '');
  } finally {
    mol.delete();
  }
}

async function main() {
  const RDKit = await loadRDKit();
  const problems = JSON.parse(fs.readFileSync(path.join(ROOT, 'problems', 'narrow.json'), 'utf8'));
  const rawNames = JSON.parse(fs.readFileSync(path.join(ROOT, 'problems', 'names.json'), 'utf8'));
  const names = {};
  for (const [s, n] of Object.entries(rawNames)) names[chem.canonical(RDKit, s)] = n;

  const molecules = {};
  const addMol = (s) => {
    if (!molecules[s]) molecules[s] = { name: names[s] || null, svg: drawSvg(RDKit, s) };
    return s;
  };

  const out = [];
  const rejected = [];
  for (const p of problems) {
    const r = checkProblem(RDKit, p);
    if (r.errors.length) {
      rejected.push(...r.errors);
      continue; // 検証に落ちた問題は収録しない
    }
    r.cands.forEach(addMol);
    const clues = p.clues.map((c, ci) => {
      if (chem.CARDS[c.card].kind === 'products') {
        r.expected[ci].forEach(addMol);
        r.table[ci].forEach((list) => list.forEach(addMol));
      }
      return { card: c.card, result: r.expected[ci], byCandidate: r.table[ci] };
    });
    out.push({
      id: p.id,
      level: p.level || 1,
      formula: p.formula,
      candidates: r.cands,
      answer: r.cands.indexOf(r.answer),
      clues,
    });
  }

  const missing = Object.keys(molecules).filter((s) => !molecules[s].name);
  if (missing.length) console.warn(`name missing for: ${missing.join(' ')}`);

  const cards = {};
  for (const [k, v] of Object.entries(chem.CARDS)) {
    cards[k] = { name: v.name, action: v.action, kind: v.kind, yes: v.yes, no: v.no, none: v.none };
  }
  const data = { version: 1, cards, molecules, problems: out };

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
