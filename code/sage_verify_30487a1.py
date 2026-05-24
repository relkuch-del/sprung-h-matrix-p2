#!/usr/bin/env python3
"""
SAGE-VERIFIKATION an 30487a1 (a_2=0 ss).

Ziel:
  Verifiziere bit-genau die Brücke
    eps^-1 * (1, 2) = (-2/3, 4/3)
  in Sage, bei einem a_2=0 ss Trigger.

Mechanik:
  Analog sage_verify_17963c1_v2.py, aber bei a_p=0:
    phi    = [[0, -1/2], [1, 0]]
    eps    = (I-phi)^2 = [[1/2, 1], [-2, 1/2]],  det(eps) = 9/4
    eps^-1 = ?
    eps^-1 * (1, 2) = (-2/3, 4/3) erwartet.

Trigger:  30487a1 (im Memory als a_2=0-Anker, Lemma D.α-II.A: det=9/4 verifiziert).

Aufruf:
  conda activate sage10 && cd /mnt/c/Temp/Forschung/paper_log_alpha && \
    python3 -u sage_verify_30487a1.py 2>&1 | tee sage_verify_30487a1_log.txt
"""

import sys
from sage.rings.padics import precision_error as _pe
sys.modules['precision_error'] = _pe

from sage.all import EllipticCurve, QQ, matrix, vector


P = 2
LABEL = "30487a1"


def main():
    print("=" * 70)
    print(f"SAGE-VERIFIKATION: {LABEL}, p={P}")
    print("=" * 70)

    E = EllipticCurve(LABEL)
    a_p = int(E.ap(P))
    print(f"  Trigger:  {LABEL}")
    print(f"  Konduktor: {E.conductor()}")
    print(f"  a_2(E):   {a_p}")
    print(f"  Rang:     {E.rank()}")
    print()

    if a_p != 0:
        print(f"WARNUNG: a_2 = {a_p} ist nicht 0. Skript fuer a_p=0 geschrieben.")
        return

    p_q = QQ(P)
    ap_q = QQ(a_p)
    I2 = matrix.identity(QQ, 2)
    phi = matrix(QQ, [[0, -1 / p_q], [1, ap_q / p_q]])
    eps = (I2 - phi) ** 2
    eps_inv = eps.inverse()
    eps_inv_T = eps_inv.transpose()

    print("=== Sage-Basis bei a_p=0 ===")
    print(f"phi:    {phi.list()}")
    print(f"eps:    {eps.list()},  det = {eps.determinant()}")
    print(f"eps^-1: {eps_inv.list()}")
    print(f"eps^-T: {eps_inv_T.list()}")
    print()

    # Lemma D.α-II.A Check
    expected_det = (QQ(3) - ap_q) ** 2 / 4
    print(f"Lemma D.alpha-II.A: det(eps) = (3 - a_p)^2 / 4 = {expected_det}")
    print(f"Verifikation:       det(eps) = {eps.determinant()}")
    print(f"Match:              {eps.determinant() == expected_det}")
    print()

    # Bruecken-Test
    print("=" * 70)
    print("BRUECKEN-TEST: v_Sage -> v_Sprung via eps^-1")
    print("=" * 70)
    v_Sage = vector(QQ, [QQ(1), QQ(2)])
    v_Sprung_predicted = eps_inv * v_Sage
    v_Sprung_expected = vector(QQ, [QQ(-2) / 3, QQ(4) / 3])

    print(f"v_Sage              = {tuple(v_Sage)}")
    print(f"eps^-1 * v_Sage     = {tuple(v_Sprung_predicted)}")
    print(f"v_Sage * eps^-T     = {tuple(v_Sage * eps_inv_T)}")
    print(f"v_Sprung erwartet   = {tuple(v_Sprung_expected)} = (-2/3, 4/3)")
    print()

    match_eps_inv = v_Sprung_predicted == v_Sprung_expected
    match_eps_inv_T = (v_Sage * eps_inv_T) == v_Sprung_expected
    print(f"Match eps^-1 * v_Sage:    {match_eps_inv}")
    print(f"Match v_Sage * eps^-T:    {match_eps_inv_T}")
    print()

    # DPVAL fuer Vollstaendigkeit (rank 30487a1 muss geprueft werden)
    print("=" * 70)
    print("DPVAL-Reihen (informativ)")
    print("=" * 70)
    try:
        Lp = E.padic_lseries(P)
        sp = Lp.Dp_valued_series(n=11, prec=20)
        print(f"{'k':>3} | {'sp[0][k]':>15} | {'sp[1][k]':>15}")
        print("-" * 45)
        for k in range(11):
            print(f"{k:>3} | {str(sp[0][k]):>15} | {str(sp[1][k]):>15}")
    except Exception as e:
        print(f"DPVAL-Berechnung scheiterte: {e}")

    print()
    print("=" * 70)
    if match_eps_inv and match_eps_inv_T:
        print("THEOREM 1.2 (a_p=0) VERIFIZIERT bit-genau.")
    else:
        print("THEOREM 1.2 (a_p=0) NICHT verifiziert — Abweichung gefunden!")
    print("=" * 70)


if __name__ == "__main__":
    main()
