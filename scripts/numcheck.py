#!/usr/bin/env python3
"""数値検算: 計算問題を Python で独立に計算し直し、有効数字まで想定解と一致するか確かめる。

JavaScript の生成器とは別に、分数（fractions.Fraction）で厳密に計算してから有効数字に丸める。
入力: generated/calcs.json（validate.js / gen26.js が書き出す計算問題の一覧）
使い方: python3 scripts/numcheck.py [calcs.json]
"""
import json
import re
import sys
from decimal import Decimal, ROUND_HALF_EVEN
from fractions import Fraction
from pathlib import Path

# 京大で与えられる原子量
AW = {"H": Fraction(1), "C": Fraction(12), "N": Fraction(14), "O": Fraction(16), "S": Fraction(32), "Cl": Fraction(71, 2)}


def molar_mass(formula):
    m = Fraction(0)
    for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        m += AW[el] * (int(n) if n else 1)
    return m


def round_sig(x, n):
    """有効数字 n けたに丸める（JavaScript の toPrecision と同じく最近接へ）。"""
    if x == 0:
        return Decimal(0)
    d = Decimal(x.numerator) / Decimal(x.denominator)
    e = d.adjusted()
    q = Decimal(1).scaleb(e - n + 1)
    return d.quantize(q, rounding=ROUND_HALF_EVEN)


def F(x):
    return Fraction(str(x))


def counts(formula):
    return {el: int(n) if n else 1 for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", formula)}


def check(c):
    k = c["kind"]
    sig = c.get("sig", 3)
    if k == "hydrolysis_mass":
        v = F(c["w"]) / molar_mass(c["formulaX"]) * c["n"] * molar_mass(c["formulaT"])
    elif k == "resolution_X":
        v = F(c["w"]) / molar_mass(c["formulaG"]) / 2 * molar_mass(c["formulaH"])
    elif k == "resolution_Y_remaining":
        v = (F(c["w"]) / molar_mass(c["formulaG"]) - F(c["b"]) / molar_mass(c["formulaH"])) * molar_mass(c["formulaG"])
    elif k == "resolution_Y_fraction":
        half = F(c["w"]) / molar_mass(c["formulaG"]) / 2
        nh = F(c["b"]) / molar_mass(c["formulaH"])
        p, q = c["p"], c["q"]
        plus = half - nh * Fraction(p, p + q)
        minus = half - nh * Fraction(q, p + q)
        if plus <= 0 or minus <= 0:
            return f"残りの物質量が負 ({plus}, {minus})"
        v = 100 * minus / (plus + minus)
    elif k == "h2":
        # C=C の数 = H₂ の物質量 ÷ X の物質量（整数に丸めて、ずれが 3% 以内）
        ratio = (F(c["V"]) / F("22.4")) / (F(c["m"]) / molar_mass(c["formula"]))
        n = round(ratio)
        if abs(ratio - n) > Fraction(3, 100) * n or n != c["answer"]:
            return f"H2: 比 {float(ratio):.4f} → {n}、想定解 {c['answer']}"
        return None
    elif k == "combustion":
        # 燃焼分析から組成式、分子量から分子式を独立に求める
        mC = F(c["co2"]) * 12 / 44
        mH = F(c["h2o"]) * 2 / 18
        mO = F(c["sample"]) - mC - mH
        nC, nH, nO = mC / 12, mH, mO / 16
        target = counts(c["answer"])
        M = molar_mass(c["answer"])
        if abs(M - F(c["M"])) > 1:
            return f"分子量 {c['M']} と分子式 {c['answer']} の式量 {float(M)} が合わない"
        # 分子量から 1 分子の原子数を出し、想定の分子式と一致するか
        k_ = F(c["M"]) / F(c["sample"])
        got = {"C": round(nC * k_), "H": round(nH * k_), "O": round(nO * k_) if nO > Fraction(1, 100) else 0}
        want = {"C": target.get("C", 0), "H": target.get("H", 0), "O": target.get("O", 0)}
        if got != want:
            return f"燃焼分析から {got}、想定 {want}"
        return None
    else:
        return f"未対応の計算 {k}"
    got = round_sig(v, sig)
    if got != Decimal(str(c["answer"])):  # Decimal は数値として比べる（8.30 と 8.3 は等しい）
        return f"{k}: Python {got} ≠ 想定解 {c['answer']}（有効数字 {sig}）"
    return None


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    path = Path(args[0]) if args else Path(__file__).resolve().parent.parent / "generated" / "calcs.json"
    items = json.loads(path.read_text(encoding="utf-8"))
    bad = []
    for it in items:
        msg = check(it["calc"])
        if msg:
            bad.append(f"{it['id']}: {msg}")
    out = {"checked": len(items), "errors": bad}
    if "--json" in sys.argv:
        print(json.dumps(out, ensure_ascii=False))
    else:
        for b in bad:
            print("  " + b)
        print(("NG" if bad else "OK") + f": 数値検算 {len(items)} 問（Python・分数で計算し直し）")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
