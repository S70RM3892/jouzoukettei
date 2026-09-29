'use strict';
// 京大で毎年のように出る量的計算（燃焼分析・元素分析・付加量）の計算器。
// 生成器が問題の数値を作るときと、過去問テストで答えを確かめるときに使う。
const M = { H: 1.0, C: 12.0, N: 14.0, O: 16.0, S: 32.0 };

function counts(f) {
  const c = {};
  for (const [, el, n] of f.matchAll(/([A-Z][a-z]?)(\d*)/g)) c[el] = (c[el] || 0) + (n ? +n : 1);
  return c;
}
function mass(f) {
  return Object.entries(counts(f)).reduce((s, [el, n]) => s + M[el] * n, 0);
}
function fmt(c) {
  const order = ['C', 'H', ...Object.keys(c).filter((e) => e !== 'C' && e !== 'H').sort()];
  return order.filter((e) => c[e]).map((e) => e + (c[e] > 1 ? c[e] : '')).join('');
}
function gcd(a, b) { return b ? gcd(b, a % b) : a; }

// 燃焼分析: 試料 mg と生じた CO2・H2O の mg から組成式（C,H,O のみ）
function empiricalFromCombustion(sample, co2, h2o) {
  const mC = (co2 * 12) / 44;
  const mH = (h2o * 2) / 18;
  const mO = sample - mC - mH;
  const mol = { C: mC / 12, H: mH / 1, O: Math.max(0, mO / 16) };
  return empiricalFromMoles(mol);
}
// 質量パーセントから組成式
function empiricalFromPercent(pct) {
  const mol = {};
  for (const [el, p] of Object.entries(pct)) mol[el] = p / M[el];
  if (pct.O === undefined) {
    const rest = 100 - Object.values(pct).reduce((a, b) => a + b, 0);
    if (rest > 0.5) mol.O = rest / 16;
  }
  return empiricalFromMoles(mol);
}
function empiricalFromMoles(mol) {
  const min = Math.min(...Object.values(mol).filter((x) => x > 1e-9));
  for (let k = 1; k <= 12; k++) {
    const r = {};
    let ok = true;
    for (const [el, x] of Object.entries(mol)) {
      const v = (x / min) * k;
      // 原子数が大きいと百分率の丸めの誤差も大きくなるので、許す幅を比例して広げる（京大2000: C₅₄H₉₀O）
      // k を大きくすると何でも整数に近づくので、許す幅は k によらず一定にする
      if (Math.abs(v - Math.round(v)) > 0.1 + 0.0015 * v) { ok = false; break; }
      if (Math.round(v)) r[el] = Math.round(v);
    }
    if (ok) {
      const g = Object.values(r).reduce(gcd);
      Object.keys(r).forEach((e) => { r[e] /= g; });
      return fmt(r);
    }
  }
  return null;
}
// 組成式の整数倍で分子量の上限以下のうち最大（または分子量に一致するもの）
function molecularFromEmpirical(emp, { maxMW, mw } = {}) {
  const c = counts(emp);
  const m = mass(emp);
  const cands = [];
  for (let n = 1; n <= 20; n++) {
    const f = {};
    Object.entries(c).forEach(([e, k]) => { f[e] = k * n; });
    const mm = m * n;
    const hOk = (f.H || 0) <= 2 * (f.C || 0) + 2 + (f.N || 0);
    if (hOk && (f.H || 0) % 2 === ((f.N || 0) % 2)) cands.push({ f: fmt(f), mm });
  }
  if (mw) return (cands.find((x) => Math.abs(x.mm - mw) < 0.6) || {}).f || null;
  const ok = cands.filter((x) => x.mm <= maxMW);
  return ok.length ? ok[ok.length - 1].f : null;
}
// 1 mol の完全燃焼に要る O2 の mol
function o2ForCombustion(f) {
  const c = counts(f);
  return (c.C || 0) + (c.H || 0) / 4 - (c.O || 0) / 2;
}
// 標準状態の気体体積 L
const liters = (mol) => mol * 22.4;

module.exports = { mass, counts, empiricalFromCombustion, empiricalFromPercent, molecularFromEmpirical, o2ForCombustion, liters };
