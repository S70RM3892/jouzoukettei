'use strict';
// ブラウザで動く小さな SMILES 読み取り・比較・描画（ネットなしの HTML で構造式の記述を採点するため）。
//  parse:  SMILES → 分子グラフ（原子の元素・芳香族・電荷・水素の数、結合）。立体の記号（@ / \）は読み飛ばす
//  same:   2つの分子が同じ構造か（重原子のつながりと、各原子の元素・水素の数で比べる。ケクレ構造と芳香族の書き方の違い、立体は区別しない）
//  toSvg:  2D 配置して骨格式の SVG を描く（入力のプレビュー）
//  formula: 分子式
// ビルド時に index.html に埋め込む。Node では require して検査（scripts/test-smiles.js）
(function (root) {
  const ORGANIC = { B: [3], C: [4], N: [3, 5], O: [2], P: [3, 5], S: [2, 4, 6], F: [1], Cl: [1], Br: [1], I: [1] };
  const AROM = new Set(['b', 'c', 'n', 'o', 'p', 's']);

  function parse(src) {
    const s = String(src).trim();
    if (!s) throw new Error('空');
    const atoms = [];
    const bonds = [];
    const stack = [];
    const rings = {};
    let prev = null;
    let bondOrder = null;
    let i = 0;
    const addAtom = (a) => {
      atoms.push(a);
      const idx = atoms.length - 1;
      if (prev !== null) {
        const arom = atoms[prev].arom && a.arom;
        bonds.push({ a: prev, b: idx, order: bondOrder || (arom ? 1.5 : 1) });
      }
      bondOrder = null;
      prev = idx;
    };
    while (i < s.length) {
      const ch = s[i];
      if (ch === '(') { if (prev === null) throw new Error('( の前に原子がない'); stack.push(prev); i++; continue; }
      if (ch === ')') { if (!stack.length) throw new Error(') が多い'); prev = stack.pop(); i++; continue; }
      if (ch === '-') { bondOrder = 1; i++; continue; }
      if (ch === '=') { bondOrder = 2; i++; continue; }
      if (ch === '#') { bondOrder = 3; i++; continue; }
      if (ch === ':') { bondOrder = 1.5; i++; continue; }
      if (ch === '/' || ch === '\\') { i++; continue; }
      if (ch === '.') { prev = null; bondOrder = null; i++; continue; }
      if (/[0-9%]/.test(ch)) {
        let num;
        if (ch === '%') { num = s.slice(i + 1, i + 3); i += 3; } else { num = ch; i++; }
        if (prev === null) throw new Error('環の番号の前に原子がない');
        if (rings[num]) {
          const r = rings[num];
          const arom = atoms[r.atom].arom && atoms[prev].arom;
          bonds.push({ a: r.atom, b: prev, order: bondOrder || r.order || (arom ? 1.5 : 1) });
          delete rings[num];
        } else rings[num] = { atom: prev, order: bondOrder };
        bondOrder = null;
        continue;
      }
      if (ch === '[') {
        const j = s.indexOf(']', i);
        if (j < 0) throw new Error('] がない');
        const body = s.slice(i + 1, j);
        const m = body.match(/^(\d*)([A-Z][a-z]?|[a-z][a-z]?)(@*)(H\d*)?([+-]\d*|[+-]+)?(:\d+)?$/);
        if (!m) throw new Error(`[${body}] が読めない`);
        let el = m[2];
        const arom = /^[a-z]/.test(el);
        if (arom) el = el[0].toUpperCase() + el.slice(1);
        const h = m[4] ? (m[4].length > 1 ? +m[4].slice(1) : 1) : 0;
        let chg = 0;
        if (m[5]) { const sign = m[5][0] === '+' ? 1 : -1; const rest = m[5].slice(1); chg = sign * (rest === '' ? 1 : /^\d+$/.test(rest) ? +rest : m[5].length); }
        addAtom({ el, arom, chg, hExplicit: h, bracket: true });
        i = j + 1;
        continue;
      }
      // 有機サブセット
      let el = null;
      if (s.startsWith('Cl', i)) el = 'Cl'; else if (s.startsWith('Br', i)) el = 'Br';
      if (el) { addAtom({ el, arom: false, chg: 0 }); i += 2; continue; }
      if ('BCNOPSFI'.includes(ch)) { addAtom({ el: ch, arom: false, chg: 0 }); i++; continue; }
      if (AROM.has(ch)) { addAtom({ el: ch.toUpperCase(), arom: true, chg: 0 }); i++; continue; }
      throw new Error(`「${ch}」が読めない（${i + 1} 文字目）`);
    }
    if (stack.length) throw new Error('( が閉じていない');
    if (Object.keys(rings).length) throw new Error(`環の番号 ${Object.keys(rings).join(',')} が閉じていない`);
    // 水素の数
    atoms.forEach((a, k) => {
      if (a.bracket) { a.h = a.hExplicit; return; }
      const nb = bonds.filter((b) => b.a === k || b.b === k);
      const sum = nb.reduce((t, b) => t + (b.order === 1.5 ? 1 : b.order), 0) + (a.arom ? 1 : 0);
      const vals = ORGANIC[a.el] || [0];
      const v = vals.find((x) => x >= sum);
      a.h = v === undefined ? 0 : v - sum;
    });
    return { atoms, bonds };
  }

  function formula(g) {
    const c = {};
    g.atoms.forEach((a) => { c[a.el] = (c[a.el] || 0) + 1; c.H = (c.H || 0) + a.h; });
    const order = ['C', 'H', ...Object.keys(c).filter((e) => e !== 'C' && e !== 'H').sort()];
    return order.filter((e) => c[e]).map((e) => e + (c[e] > 1 ? c[e] : '')).join('');
  }

  // 連結成分に分ける（「順序は問わない」複数の構造の答え）
  function components(g) {
    const seen = new Set();
    const out = [];
    for (let s = 0; s < g.atoms.length; s++) {
      if (seen.has(s)) continue;
      const comp = [];
      const st = [s];
      seen.add(s);
      while (st.length) {
        const x = st.pop(); comp.push(x);
        g.bonds.forEach((b) => { const y = b.a === x ? b.b : b.b === x ? b.a : -1; if (y >= 0 && !seen.has(y)) { seen.add(y); st.push(y); } });
      }
      const map = new Map(comp.map((x, k) => [x, k]));
      out.push({ atoms: comp.map((x) => g.atoms[x]), bonds: g.bonds.filter((b) => map.has(b.a)).map((b) => ({ a: map.get(b.a), b: map.get(b.b), order: b.order })) });
    }
    return out;
  }

  // 原子の色（元素・水素の数）とつながりで同型か調べる。まず色の細分化で候補を絞り、残りは総当たり
  const label = (a) => `${a.el}${a.h}`;
  function refine(g) {
    const adj = g.atoms.map(() => []);
    g.bonds.forEach((b) => { adj[b.a].push(b.b); adj[b.b].push(b.a); });
    let col = g.atoms.map(label);
    for (let it = 0; it < g.atoms.length; it++) {
      const next = col.map((c, i) => c + '|' + adj[i].map((j) => col[j]).sort().join(','));
      const uniq = [...new Set(next)].sort();
      const nc = next.map((x) => String(uniq.indexOf(x)));
      if (new Set(nc).size === new Set(col).size) { col = nc; break; }
      col = nc;
    }
    return { adj, col: col.map((c, i) => label(g.atoms[i]) + '#' + c) };
  }
  function same(g1, g2) {
    if (g1.atoms.length !== g2.atoms.length || g1.bonds.length !== g2.bonds.length) return false;
    if (formula(g1) !== formula(g2)) return false;
    // 色は両方の分子をまとめて細分化しないと比べられないので、合わせた1つのグラフで細分化する
    const n = g1.atoms.length;
    const u = { atoms: [...g1.atoms, ...g2.atoms], bonds: [...g1.bonds, ...g2.bonds.map((b) => ({ a: b.a + n, b: b.b + n, order: b.order }))] };
    const { adj, col } = refine(u);
    const c1 = col.slice(0, n), c2 = col.slice(n);
    if ([...c1].sort().join() !== [...c2].sort().join()) return false;
    const map = new Array(n).fill(-1);
    const used = new Array(n).fill(false);
    const order = [...Array(n).keys()].sort((x, y) => c1.filter((c) => c === c1[x]).length - c1.filter((c) => c === c1[y]).length);
    const rec = (d) => {
      if (d === n) return true;
      const x = order[d];
      for (let y = 0; y < n; y++) {
        if (used[y] || c2[y] !== c1[x]) continue;
        let ok = true;
        for (const nx of adj[x]) {
          if (map[nx] < 0) continue;
          if (!adj[y + n].includes(map[nx] + n)) { ok = false; break; }
        }
        if (!ok) continue;
        map[x] = y; used[y] = true;
        if (rec(d + 1)) return true;
        map[x] = -1; used[y] = false;
      }
      return false;
    };
    return rec(0);
  }

  // 答え（SMILES の配列）と入力（SMILES。「.」・改行・読点で区切って複数）を比べる。順序は問わない
  function grade(answers, input) {
    const parts = String(input || '').split(/[\n、,，;；]+/).flatMap((x) => x.split('.')).map((x) => x.trim()).filter(Boolean);
    let gs;
    try { gs = parts.map(parse); } catch (e) { return { hits: 0, extra: parts.length, n: answers.length, error: e.message }; }
    const ans = answers.map(parse);
    const used = new Array(ans.length).fill(false);
    let hits = 0, extra = 0;
    for (const g of gs) {
      const k = ans.findIndex((a, i) => !used[i] && same(a, g));
      if (k >= 0) { used[k] = true; hits++; } else extra++;
    }
    return { hits, extra, n: ans.length };
  }

  // ---------- 2D 配置と SVG ----------
  function rings(g) {
    // 小さな環（8 員環まで）を、結合ごとに最短の閉路として集める
    const adj = g.atoms.map(() => []);
    g.bonds.forEach((b) => { adj[b.a].push(b.b); adj[b.b].push(b.a); });
    const out = [];
    const key = new Set();
    for (const b of g.bonds) {
      // b を除いて a→b の最短経路
      const prev = new Map([[b.a, -1]]);
      const q = [b.a];
      while (q.length) {
        const x = q.shift();
        if (x === b.b) break;
        for (const y of adj[x]) {
          if ((x === b.a && y === b.b) || prev.has(y)) continue;
          prev.set(y, x); q.push(y);
        }
      }
      if (!prev.has(b.b)) continue;
      const path = [];
      for (let x = b.b; x !== -1; x = prev.get(x)) path.push(x);
      if (path.length > 8) continue;
      const k = [...path].sort((p, q2) => p - q2).join(',');
      if (!key.has(k)) { key.add(k); out.push(path); }
    }
    return { adj, rings: out };
  }
  function layout(g) {
    const n = g.atoms.length;
    const { adj, rings: rs } = rings(g);
    const pos = g.atoms.map(() => null);
    // 環はまず正多角形に置き、残りは幅優先で 120° ずつ伸ばす
    const placeRing = (ring, cx, cy, start) => {
      const k = ring.length;
      const R = 1 / (2 * Math.sin(Math.PI / k));
      ring.forEach((a, i) => { if (!pos[a]) pos[a] = [cx + R * Math.cos(start + (2 * Math.PI * i) / k), cy + R * Math.sin(start + (2 * Math.PI * i) / k)]; });
    };
    const seen = new Set();
    const queue = [];
    const startAt = (a0) => {
      pos[a0] = pos[a0] || [0, 0];
      queue.push(a0); seen.add(a0);
      while (queue.length) {
        const x = queue.shift();
        // x を含む未配置の環
        for (const ring of rs) {
          if (!ring.includes(x) || ring.every((a) => pos[a])) continue;
          const placed = ring.filter((a) => pos[a]);
          if (placed.length >= 2) {
            // 既に置かれた2点の辺の外側に環を作る
            const i0 = ring.indexOf(placed[0]);
            const rot = [...ring.slice(i0), ...ring.slice(0, i0)];
            const [p, q2] = [pos[rot[0]], pos[rot[1]] || pos[rot[rot.length - 1]]];
            const mx = (p[0] + q2[0]) / 2, my = (p[1] + q2[1]) / 2;
            const cen = centroidOfPlaced(x);
            let dx = mx - cen[0], dy = my - cen[1];
            const L = Math.hypot(dx, dy) || 1; dx /= L; dy /= L;
            const k = ring.length;
            const apo = 1 / (2 * Math.tan(Math.PI / k));
            const cx = mx + dx * apo, cy = my + dy * apo;
            const a0x = Math.atan2(p[1] - cy, p[0] - cx);
            const R = 1 / (2 * Math.sin(Math.PI / k));
            const dir = pos[rot[1]] ? Math.sign(cross(cx, cy, p, pos[rot[1]])) || 1 : 1;
            rot.forEach((a, i) => { if (!pos[a]) pos[a] = [cx + R * Math.cos(a0x + dir * (2 * Math.PI * i) / k), cy + R * Math.sin(a0x + dir * (2 * Math.PI * i) / k)]; });
          } else {
            const p = pos[x];
            const par = adj[x].find((y) => pos[y] && !ring.includes(y));
            let ang = par !== undefined ? Math.atan2(p[1] - pos[par][1], p[0] - pos[par][0]) : 0;
            const k = ring.length;
            const R = 1 / (2 * Math.sin(Math.PI / k));
            const cx = p[0] + R * Math.cos(ang), cy = p[1] + R * Math.sin(ang);
            const i0 = ring.indexOf(x);
            const rot = [...ring.slice(i0), ...ring.slice(0, i0)];
            const a0x = ang + Math.PI;
            rot.forEach((a, i) => { if (!pos[a]) pos[a] = [cx + R * Math.cos(a0x + (2 * Math.PI * i) / k), cy + R * Math.sin(a0x + (2 * Math.PI * i) / k)]; });
          }
        }
        const kids = adj[x].filter((y) => !pos[y]);
        const placedNb = adj[x].filter((y) => pos[y]);
        const base = placedNb.length ? Math.atan2(pos[x][1] - avg(placedNb, 1), pos[x][0] - avg(placedNb, 0)) : 0;
        const spread = placedNb.length ? [0, 1, -1, 2, -2].map((t) => base + t * (Math.PI / 3)) : [Math.PI / 6, -Math.PI / 6 + Math.PI, Math.PI / 2 + Math.PI, -Math.PI / 2];
        // ジグザグ: 親から来た向きを交互に振る
        const flip = (x % 2 ? 1 : -1);
        let si = 0;
        for (const y of kids) {
          let ang = spread[Math.min(si, spread.length - 1)];
          if (placedNb.length === 1 && kids.length === 1) ang = base + flip * Math.PI / 6;
          else if (placedNb.length === 1 && kids.length === 2) ang = base + (si ? -1 : 1) * Math.PI / 3;
          else if (placedNb.length === 1 && kids.length === 3) ang = base + [0, 1, -1][si] * Math.PI / 2;
          pos[y] = [pos[x][0] + Math.cos(ang), pos[x][1] + Math.sin(ang)];
          si++;
        }
        for (const y of adj[x]) if (!seen.has(y)) { seen.add(y); queue.push(y); }
      }
    };
    function avg(list, k) { return list.reduce((t, y) => t + pos[y][k], 0) / list.length; }
    function centroidOfPlaced(x) { const pl = adj[x].filter((y) => pos[y]); return pl.length ? [avg(pl, 0), avg(pl, 1)] : pos[x]; }
    function cross(cx, cy, p, q2) { return (p[0] - cx) * (q2[1] - cy) - (p[1] - cy) * (q2[0] - cx); }
    let offX = 0;
    for (let a = 0; a < n; a++) {
      if (pos[a]) continue;
      pos[a] = [offX, 0];
      startAt(a);
      const xs = pos.filter(Boolean).map((p) => p[0]);
      offX = Math.max(...xs) + 2;
    }
    // 重なりを少しだけほどく（結合していない原子どうしを離す）
    for (let it = 0; it < 60; it++) {
      for (let a = 0; a < n; a++) for (let b = a + 1; b < n; b++) {
        if (adj[a].includes(b)) continue;
        const dx = pos[b][0] - pos[a][0], dy = pos[b][1] - pos[a][1];
        const d = Math.hypot(dx, dy);
        if (d < 0.75) { const f = (0.75 - d) / 2 / (d || 1); pos[a][0] -= dx * f; pos[a][1] -= dy * f; pos[b][0] += dx * f; pos[b][1] += dy * f; }
      }
    }
    return pos;
  }

  function toSvg(g, opts = {}) {
    const pos = layout(g);
    const S = opts.scale || 32;
    const xs = pos.map((p) => p[0]), ys = pos.map((p) => p[1]);
    const minX = Math.min(...xs), minY = Math.min(...ys);
    const W = (Math.max(...xs) - minX) * S + 48, H = (Math.max(...ys) - minY) * S + 40;
    const P = (i) => [(pos[i][0] - minX) * S + 24, (pos[i][1] - minY) * S + 20];
    const showLabel = (a) => a.el !== 'C' || a.chg;
    const lines = [];
    // 芳香環の結合は 1.5。環の中で交互に二重線にする（見た目だけ）
    const aromDouble = new Set();
    {
      const done = new Set();
      g.bonds.forEach((b, k) => {
        if (b.order !== 1.5) return;
        const free = (x) => !g.bonds.some((c, j) => aromDouble.has(j) && (c.a === x || c.b === x));
        if (!done.has(k) && free(b.a) && free(b.b)) aromDouble.add(k);
        done.add(k);
      });
    }
    const cen = [pos.reduce((t, p) => t + p[0], 0) / pos.length, pos.reduce((t, p) => t + p[1], 0) / pos.length];
    g.bonds.forEach((b, k) => {
      let [x1, y1] = P(b.a), [x2, y2] = P(b.b);
      const sh = (i, x, y, ox, oy) => (showLabel(g.atoms[i]) ? [x + ox * 9, y + oy * 9] : [x, y]);
      const dx = x2 - x1, dy = y2 - y1, L = Math.hypot(dx, dy) || 1;
      const ux = dx / L, uy = dy / L;
      [x1, y1] = sh(b.a, x1, y1, ux, uy);
      [x2, y2] = sh(b.b, x2, y2, -ux, -uy);
      lines.push(`<line x1="${x1.toFixed(1)}" y1="${y1.toFixed(1)}" x2="${x2.toFixed(1)}" y2="${y2.toFixed(1)}"/>`);
      const order = b.order === 1.5 ? (aromDouble.has(k) ? 2 : 1) : b.order;
      if (order >= 2) {
        // 二重線は分子の中心側に
        let nx = -uy, ny = ux;
        const mx = (pos[b.a][0] + pos[b.b][0]) / 2 - cen[0], my = (pos[b.a][1] + pos[b.b][1]) / 2 - cen[1];
        if (nx * mx + ny * my > 0) { nx = -nx; ny = -ny; }
        const o = 5;
        const t = order === 3 ? [o, -o] : [o];
        for (const tt of t) lines.push(`<line x1="${(x1 + nx * tt + ux * 4).toFixed(1)}" y1="${(y1 + ny * tt + uy * 4).toFixed(1)}" x2="${(x2 + nx * tt - ux * 4).toFixed(1)}" y2="${(y2 + ny * tt - uy * 4).toFixed(1)}"/>`);
      }
    });
    const labels = g.atoms.map((a, i) => {
      if (!showLabel(a)) return '';
      const [x, y] = P(i);
      const hh = a.h ? `H${a.h > 1 ? `<tspan dy="3" font-size="9">${a.h}</tspan>` : ''}` : '';
      const chg = a.chg ? `<tspan dy="-6" font-size="9">${a.chg > 0 ? '+' : '−'}</tspan>` : '';
      return `<text x="${x.toFixed(1)}" y="${(y + 5).toFixed(1)}" text-anchor="middle">${a.el}${hh}${chg}</text>`;
    }).join('');
    return `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${W.toFixed(0)} ${H.toFixed(0)}" width="${W.toFixed(0)}" height="${H.toFixed(0)}" style="max-width:100%;height:auto"><g stroke="currentColor" stroke-width="1.3" stroke-linecap="round">${lines.join('')}</g><g fill="currentColor" font-size="14" font-family="serif">${labels}</g></svg>`;
  }

  const api = { parse, formula, same, grade, components, toSvg, layout };
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  else root.SmilesLite = api;
}(typeof window !== 'undefined' ? window : globalThis));
