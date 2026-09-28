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
    return svg.replace(/\n/g, '');
  } finally {
    mol.delete();
  }
}

// 問題データ中の "pep:Gly-Ala" 表記を集めて、配列をそのまま名前にする
function collectPeptideNames(RDKit, node, names) {
  if (typeof node === 'string') {
    if (node.startsWith('pep:')) names[chem.canonical(RDKit, node)] = node.slice(4);
  } else if (node && typeof node === 'object') {
    Object.values(node).forEach((v) => collectPeptideNames(RDKit, v, names));
  }
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
    if (!molecules[s]) molecules[s] = { name: names[s] || null, svg: drawSvg(RDKit, s) };
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
    return { ...extra, candidates: r.cands, answer: r.cands.indexOf(r.answer), clues: shapeClues(clues, r) };
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
      const stages = [{
        type: 'split',
        bonds: r.bonds,
        products: r.frags.map((f) => ({ label: f.label, formula: f.formula, count: f.count })),
      }];
      for (const f of r.frags) {
        addMol(f.answer);
        if (f.given) stages.push({ type: 'given', label: f.label, smiles: f.answer, note: f.note });
        else stages.push(narrowData(f.r, f.clues, { type: 'narrow', label: f.label, id: `${p.id}/${f.label}`, formula: f.formula, level: p.level }));
      }
      if (r.assemble) {
        stages.push(narrowData(r.assemble.r, r.assemble.clues, { type: 'narrow', label: 'X', id: `${p.id}/X`, formula: p.formula, level: p.level }));
      }
      out.push({ ...base, title: p.title, story: p.story, answerSmiles: r.X, stages });
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

  const missing = Object.keys(molecules).filter((s) => !molecules[s].name);
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
