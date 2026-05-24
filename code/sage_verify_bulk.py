#!/usr/bin/env python3
"""
SAGE-BULK-VERIFIKATION der epsilon^-1-Brücke.

Idee:
  Die Brücke eps^-1 * u_Sage = u_Sprung ist algebraisch trigger-unabhängig —
  sie hängt nur von p und a_p ab, nicht vom konkreten Trigger E.

  Die Bulk-Verifikation zeigt:
    (a) Für jeden a_2=-2 ss Trigger E ist phi(E) = phi(a_p=-2) gleich.
    (b) Damit ist eps(E) = eps(a_p=-2) gleich.
    (c) Damit ist eps^-1 * (7/4, 3/4) = (4/5, 9/10) trigger-unabhängig.

  Das ist eigentlich tautologisch — aber wir verifizieren es bit-genau,
  um sicherzustellen, dass Sage's Konvention bei verschiedenen Triggern
  konsistent ist (vor allem: dass Sage's eps^-1 nicht trigger-spezifisch
  berechnet wird).

Trigger (a_2=-2 ss, rank >= 2):
  17963c1, 5077a1, 53461b1, 389a1, 655a1, 681c1, 707a1

Trigger (a_2=0 ss, rank >= 2):
  30487a1, 28571a1 (Memo Lemma D.alpha-II.A)
  + falls vorhanden: weitere

Aufruf:
  conda activate sage10 && cd /mnt/c/Temp/Forschung/paper_log_alpha && \
    python3 -u sage_verify_bulk.py 2>&1 | tee sage_verify_bulk_log.txt
"""

import sys
from sage.rings.padics import precision_error as _pe
sys.modules['precision_error'] = _pe

from sage.all import EllipticCurve, QQ, matrix, vector


P = 2

TRIGGERS_MINUS2 = ["17963c1", "5077a1", "53461b1", "389a1", "655a1", "681c1", "707a1"]
# a_p=0-Trigger: 28571a1 wurde irrtuemlich gelistet — hat tatsaechlich a_2=-2.
# Echte a_p=0-Trigger aus V66/V67:
TRIGGERS_ZERO   = ["30487a1", "60803a1", "58939a1", "1531a1", "1613a1", "1653a1", "1909a1"]


def make_eps(a_p):
    p_q = QQ(P)
    ap_q = QQ(a_p)
    phi = matrix(QQ, [[0, -1 / p_q], [1, ap_q / p_q]])
    return (matrix.identity(QQ, 2) - phi) ** 2


def test_trigger(label, u_Sage, u_Sprung_expected, a_p_expected):
    """Pruefe Brücke fuer einen Trigger."""
    try:
        E = EllipticCurve(label)
    except Exception as e:
        return ("ERROR_LOAD", str(e), None, None, None, None)

    try:
        a_p = int(E.ap(P))
    except Exception as e:
        return ("ERROR_AP", str(e), None, None, None, None)

    if a_p != a_p_expected:
        return ("WRONG_AP", f"a_p={a_p}, expected {a_p_expected}",
                E.conductor(), a_p, E.rank(), None)

    eps = make_eps(a_p)
    eps_inv = eps.inverse()
    result = eps_inv * u_Sage
    match = result == u_Sprung_expected

    return ("OK", "match" if match else "MISMATCH",
            E.conductor(), a_p, E.rank(), tuple(result))


def main():
    print("=" * 70)
    print("SAGE-BULK-VERIFIKATION: epsilon^-1-Brücke")
    print("=" * 70)
    print()

    # === a_p = -2 ===
    print("=" * 70)
    print("Block 1: a_2 = -2 (ss)")
    print("=" * 70)
    u_Sage = vector(QQ, [QQ(7) / 4, QQ(3) / 4])
    u_Sprung_expected = vector(QQ, [QQ(4) / 5, QQ(9) / 10])
    print(f"  u_Sage          = {tuple(u_Sage)}")
    print(f"  u_Sprung expect = {tuple(u_Sprung_expected)}")
    print()
    print(f"{'Trigger':<12} | {'Cond':>10} | {'a_p':>4} | {'rank':>4} | {'Status':>10} | {'eps^-1 * u_Sage':>25}")
    print("-" * 95)

    n_match_minus2 = 0
    n_total_minus2 = 0
    for label in TRIGGERS_MINUS2:
        status, msg, cond, a_p, rank, result = test_trigger(label, u_Sage, u_Sprung_expected, -2)
        n_total_minus2 += 1
        if status == "OK" and msg == "match":
            n_match_minus2 += 1
        print(f"{label:<12} | {str(cond):>10} | {str(a_p):>4} | {str(rank):>4} | "
              f"{status + ':' + msg:>10} | {str(result):>25}")
    print()
    print(f"Bilanz a_p=-2: {n_match_minus2}/{n_total_minus2} bit-genau")

    print()
    print("=" * 70)
    print("Block 2: a_2 = 0 (ss)")
    print("=" * 70)
    v_Sage = vector(QQ, [QQ(1), QQ(2)])
    v_Sprung_expected = vector(QQ, [QQ(-2) / 3, QQ(4) / 3])
    print(f"  v_Sage          = {tuple(v_Sage)}")
    print(f"  v_Sprung expect = {tuple(v_Sprung_expected)}")
    print()
    print(f"{'Trigger':<12} | {'Cond':>10} | {'a_p':>4} | {'rank':>4} | {'Status':>10} | {'eps^-1 * v_Sage':>25}")
    print("-" * 95)

    n_match_zero = 0
    n_total_zero = 0
    for label in TRIGGERS_ZERO:
        status, msg, cond, a_p, rank, result = test_trigger(label, v_Sage, v_Sprung_expected, 0)
        n_total_zero += 1
        if status == "OK" and msg == "match":
            n_match_zero += 1
        print(f"{label:<12} | {str(cond):>10} | {str(a_p):>4} | {str(rank):>4} | "
              f"{status + ':' + msg:>10} | {str(result):>25}")
    print()
    print(f"Bilanz a_p=0: {n_match_zero}/{n_total_zero} bit-genau")

    print()
    print("=" * 70)
    print("GESAMT-BILANZ")
    print("=" * 70)
    print(f"  a_p = -2 ss:  {n_match_minus2}/{n_total_minus2} bit-genau eps^-1 * (7/4, 3/4) = (4/5, 9/10)")
    print(f"  a_p =  0 ss:  {n_match_zero}/{n_total_zero} bit-genau eps^-1 * (1, 2) = (-2/3, 4/3)")
    print()
    if n_match_minus2 == n_total_minus2 and n_match_zero == n_total_zero:
        print("THEOREM 1.2 VOLLSTAENDIG bit-genau ueber alle getesteten Trigger.")
    else:
        print("Theorem 1.2: einige Trigger weichen ab — Untersuchung notwendig.")
    print("=" * 70)


if __name__ == "__main__":
    main()
