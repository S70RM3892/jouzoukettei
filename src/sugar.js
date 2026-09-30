'use strict';
// 糖の略記 "sac:..." と、その判定。環状の糖は立体（α/β、グルコースとガラクトース）でしか区別できないので、
// SMILES から立体を落とす chem.js の一般の判定ではなく、ここで記号として判定する。構造式は立体つきの SMILES に展開して描く。
//
// 書き方（IUPAC の糖の略記に合わせる。左が非還元末端、右が還元末端）
//   sac:Glc                     グルコース（還元末端の α/β は水溶液中で入れかわるので書かない）
//   sac:Glc(a1-4)Glc            マルトース（α-グルコースの 1 位が、右のグルコースの 4 位の O と結合）
//   sac:Glc(a1-2b)Fru           スクロース（右の単糖の還元性を示す位置どうしの結合は、右の α/β も書く）
//   sac:Glc(a1-4)[Glc(a1-6)]Glc 枝分かれ（[ ] の中も次の単糖につながる）
//   sac:Glc{2,3,4,6}            2,3,4,6-テトラ-O-メチルグルコース（{ } はメチル化された OH の位置）
// 単糖: Glc（グルコース）・Gal（ガラクトース）・Fru（フルクトース、五員環）
//
// 立体つきの SMILES は PubChem の α/β-D-グルコピラノース・ガラクトピラノース・β-D-フルクトフラノースから作り、
// マルトース・ラクトース・セロビオース・スクロース・トレハロースで PubChem と一致することを test-rules.js で確かめる。

// <R> 環の番号、<C> 還元性を示す炭素（アノマー炭素）、<On> n 位の O
const UNITS = {
  Glc: { name: 'グルコース', anom: 1, pos: [1, 2, 3, 4, 6], a: '[C@H]', b: '[C@@H]', none: 'C', flipLink: true,
    tpl: 'C([C@@H]<R>[C@H]([C@@H]([C@H](<C>(O<R>)<O1>)<O2>)<O3>)<O4>)<O6>' },
  Gal: { name: 'ガラクトース', anom: 1, pos: [1, 2, 3, 4, 6], a: '[C@H]', b: '[C@@H]', none: 'C', flipLink: true,
    tpl: 'C([C@@H]<R>[C@@H]([C@@H]([C@H](<C>(O<R>)<O1>)<O2>)<O3>)<O4>)<O6>' },
  Fru: { name: 'フルクトース', anom: 2, pos: [1, 2, 3, 4, 6], a: '[C@@]', b: '[C@]', none: 'C', flipLink: false,
    tpl: 'C([C@@H]<R>[C@H]([C@@H](<C>(O<R>)(C<O1>)<O2>)<O3>)<O4>)<O6>' },
};
const flip = (t) => (t.includes('@@') ? t.replace('@@', '@') : t.replace('@', '@@'));

class SugarError extends Error {}

// ---- 読み書き ----
function parse(s) {
  if (!s.startsWith('sac:')) throw new SugarError(`not a sugar: ${s}`);
  let i = 4;
  const src = s;
  const unit = () => {
    const m = /^(a|b)?(Glc|Gal|Fru)(\{([\d,]*)\})?/.exec(src.slice(i));
    if (!m) throw new SugarError(`bad sugar at ${i}: ${src}`);
    i += m[0].length;
    const U = UNITS[m[2]];
    const me = m[4] ? m[4].split(',').filter(Boolean).map(Number).sort((x, y) => x - y) : [];
    me.forEach((p) => { if (!U.pos.includes(p)) throw new SugarError(`no OH at ${p} in ${m[2]}`); });
    return { type: m[2], anomer: m[1] || null, me, children: [] };
  };
  const link = () => {
    const m = /^\((a|b)(\d)-(\d)(a|b)?\)/.exec(src.slice(i));
    if (!m) return null;
    i += m[0].length;
    return { anomer: m[1], from: +m[2], to: +m[3], toAnomer: m[4] || null };
  };
  // 左から読み、つながり待ちの単糖を次の単糖にぶら下げる
  const chain = () => {
    let pending = [];
    for (;;) {
      if (src[i] === '[') {
        i++;
        const inner = chain();
        if (src[i] !== ']' || !inner.link) throw new SugarError(`bad branch: ${src}`);
        i++;
        pending.push(inner);
        continue;
      }
      const u = unit();
      pending.forEach((c) => u.children.push(c));
      pending = [];
      const l = link();
      if (!l) return { node: u, link: null };
      if (l.from !== UNITS[u.type].anom) throw new SugarError(`link must start at the anomeric carbon: ${src}`);
      u.anomer = l.anomer;
      pending.push({ node: u, link: l });
      if (i >= src.length || src[i] === ']') return pending.length === 1 ? pending[0] : { node: u, link: l };
    }
  };
  const top = chain();
  if (i !== src.length || top.link) throw new SugarError(`bad sugar: ${src}`);
  check(top.node);
  return top.node;
}

function check(n) {
  const used = new Set();
  for (const c of n.children) {
    const U = UNITS[n.type];
    if (!U.pos.includes(c.link.to) || used.has(c.link.to) || n.me.includes(c.link.to)) throw new SugarError(`bad link position ${c.link.to} on ${n.type}`);
    used.add(c.link.to);
    if (c.link.to === U.anom) {
      if (!c.link.toAnomer) throw new SugarError('anomer of the parent is needed');
      n.anomer = c.link.toAnomer;
    }
    check(c.node);
  }
}

const unitStr = (n) => `${n.type}${n.me.length ? `{${n.me.join(',')}}` : ''}`;
const linkStr = (l) => `(${l.anomer}${l.from}-${l.to}${l.toAnomer || ''})`;
function serialize(n) {
  // 位置の大きい枝を主鎖に（決まった順で書く）
  const cs = [...n.children].sort((x, y) => x.link.to - y.link.to);
  const main = cs.pop();
  const side = cs.map((c) => `[${serialize(c.node)}${linkStr(c.link)}]`).join('');
  return `${main ? serialize(main.node) + linkStr(main.link) : ''}${side}${unitStr(n)}`;
}
const canonical = (s) => 'sac:' + serialize(parse(s));

// ---- 立体つき SMILES ----
function toSmiles(s) {
  const root = parse(s);
  let ring = 0, bond = 50;
  const frags = [];
  // parentO: この単糖の還元性を示す炭素がつながる相手の O の番号（なければ null）
  const emit = (n, rcToParent) => {
    const U = UNITS[n.type];
    const r = ++ring;
    if (r > 9) throw new SugarError('too many units');
    const O = {};
    U.pos.forEach((p) => { O[p] = n.me.includes(p) ? 'OC' : 'O'; });
    for (const c of n.children) {
      const b = bond++;
      O[c.link.to] = `O%${b}`;
      emit(c.node, b);
    }
    let C;
    if (rcToParent) {
      // 還元性を示す炭素が相手の O と直接結合する。この炭素の O は消え、環の番号が枝より前に来るので、H をもつ炭素では立体の向きが入れかわる
      const t = U[n.anomer];
      C = (U.flipLink ? flip(t) : t) + `%${rcToParent}`;
      O[U.anom] = '';
    } else {
      const usedAnom = n.children.some((c) => c.link.to === U.anom);
      C = n.anomer && (usedAnom || n.anomer) ? U[n.anomer] : U.none;
    }
    let t = U.tpl.replace(/<R>/g, String(r)).replace('<C>', C);
    U.pos.forEach((p) => { t = t.replace(`<O${p}>`, O[p]); });
    frags.push(t);
  };
  emit(root, null);
  return frags.join('.');
}

// ---- 判定 ----
const links = (n, out = []) => { n.children.forEach((c) => { out.push({ parent: n, child: c.node, link: c.link }); links(c.node, out); }); return out; };
const units = (n, out = []) => { out.push(n); n.children.forEach((c) => units(c.node, out)); return out; };
// 還元性: 還元末端の単糖の、還元性を示す炭素が結合に使われていない
const reducing = (root) => !root.children.some((c) => c.link.to === UNITS[root.type].anom);
const mono = (n, me) => `sac:${n.type}${me.length ? `{${[...me].sort((x, y) => x - y).join(',')}}` : ''}`;

function hydrolysis(root) {
  const us = units(root);
  if (us.length < 2) return [];
  return [...new Set(us.map((u) => mono(u, u.me)))].sort();
}

// メチル化分析: すべての OH をメチル化してから、単糖どうしの結合（グリコシド結合）だけを加水分解する。
// 還元性を示す炭素の OH（メチルグリコシド）も加水分解で OH に戻る
function methylation(root) {
  const us = units(root);
  const used = new Map(us.map((u) => [u, new Set()]));
  links(root).forEach(({ parent, child, link }) => { used.get(parent).add(link.to); used.get(child).add(UNITS[child.type].anom); });
  return [...new Set(us.map((u) => {
    const U = UNITS[u.type];
    const me = U.pos.filter((p) => p !== U.anom && !used.get(u).has(p));
    return mono(u, me);
  }))].sort();
}

// 酵素が切る結合（高校の教科書の基質特異性）
const ENZYMES = {
  maltase: (l) => l.child.type === 'Glc' && l.link.anomer === 'a' && l.link.to !== UNITS[l.parent.type].anom, // α-グルコシド結合（マルトースなど）
  invertase: (l) => l.parent.type === 'Fru' && l.link.to === 2 && l.link.toAnomer === 'b', // スクロースの結合
  lactase: (l) => l.child.type === 'Gal' && l.link.anomer === 'b', // β-ガラクトシド結合（ラクトース）
  cellobiase: (l) => l.child.type === 'Glc' && l.link.anomer === 'b' && l.link.to !== UNITS[l.parent.type].anom, // β-グルコシド結合（セロビオース）
};

function evaluate(card, s) {
  const root = parse(s);
  switch (card) {
    case 'silver_mirror':
    case 'fehling':
      return reducing(root);
    case 'hydrolysis':
      return hydrolysis(root);
    case 'methylation_analysis':
      return methylation(root);
    case 'maltase': case 'invertase': case 'lactase': case 'cellobiase': {
      const ls = links(root);
      if (!ls.length) return false;
      return ls.some((l) => ENZYMES[card](l));
    }
    case 'sac_units':
      return units(root).length;
    default:
      return undefined; // 一般の判定（立体なしの SMILES）に任せる
  }
}

// 画面に出す名前（名前のない糖）: α-Glc(1→4)Glc、2,3,4,6-テトラ-O-メチルグルコース
const NUM = ['', 'モノ', 'ジ', 'トリ', 'テトラ', 'ペンタ'];
function displayName(s) {
  const root = parse(s);
  if (!root.children.length) {
    const U = UNITS[root.type];
    return root.me.length ? `${root.me.join(',')}-${NUM[root.me.length]}-O-メチル${U.name}` : U.name;
  }
  const one = (n) => {
    const cs = [...n.children].sort((x, y) => x.link.to - y.link.to);
    const main = cs.pop();
    const side = cs.map((c) => `[${one(c.node)}(${c.link.from}→${c.link.to})]`).join('');
    const pre = main ? `${one(main.node)}(${main.link.from}→${main.link.to})` : '';
    return `${pre}${side}${n.anomer && n !== root ? (n.anomer === 'a' ? 'α-' : 'β-') : n.anomer ? (n.anomer === 'a' ? 'α-' : 'β-') : ''}${n.type}`;
  };
  return one(root);
}

const isSugar = (s) => typeof s === 'string' && s.startsWith('sac:');
module.exports = { isSugar, parse, canonical, toSmiles, evaluate, displayName, UNITS, SugarError };
