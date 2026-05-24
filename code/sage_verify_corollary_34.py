#!/usr/bin/env python3
"""
sage_verify_corollary_34.py — Direkter Beweis von Korollar 3.4 im Paper.

Verifiziert algebraisch (nicht über DPVAL) die Werte

  a_p = -2:  log_α^♯|_0 = (1/4)α,           log_α^♭|_0 = -1/2 + (1/4)α
  a_p =  0:  log_α^♯|_0 = -1/2,             log_α^♭|_0 = -(1/4)α
  a_p = +2:  log_α^♯|_0 = -(1/4)α,          log_α^♭|_0 = 1/2 + (1/4)α

Methode:
  1. H(0) := A^-2  (Theorem 3.1)
  2. M := [[1, 1], [-α^-1, -β^-1]]  (Def. 6.8)
  3. Halblog-Matrix = H(0) · M
  4. Erste Spalte = (log_α^♯|_0, log_α^♭|_0)
  5. Konversion in (1, α)-Basis und Vergleich mit Korollar 3.4

Aufruf:
  conda activate sage10 && cd /mnt/c/Temp/Forschung/paper_log_alpha && \
    python3 -u sage_verify_corollary_34.py 2>&1 | tee sage_verify_corollary_34_log.txt
"""

import sys
from sage.all import (matrix, Matrix, QQ, QuadraticField, sqrt, Rational,
                      NumberField, polygen)

P = 2


def setup_field(a_p):
    """Erstellt Q(α) als NumberField mit α² - a_p·α + p = 0."""
    x = polygen(QQ, 'x')
    minpoly = x ** 2 - a_p * x + P
    K = NumberField(minpoly, 'alpha')
    alpha = K.gen()
    return K, alpha


def compute_halflog_first_column(a_p):
    """Rechne H(0)·M aus und gibt erste Spalte = (log^♯|_0, log^♭|_0)."""
    K, alpha = setup_field(a_p)
    beta = a_p - alpha  # Galois-Konjugat: α + β = a_p

    # A in Q (rationale Einträge)
    A = matrix(K, [[a_p, 1], [-P, 0]])
    A_inv2 = A.inverse() ** 2  # = H(0) nach Theorem 3.1

    # M wirkt in K = Q(α)
    M = matrix(K, [[1, 1], [-1 / alpha, -1 / beta]])

    # Halblog-Matrix
    Log_mat = A_inv2 * M

    # Erste Spalte = (log_α^♯|_0, log_α^♭|_0)
    log_sharp_0 = Log_mat[0, 0]
    log_flat_0 = Log_mat[1, 0]
    return K, alpha, beta, A_inv2, M, log_sharp_0, log_flat_0


def to_alpha_coords(element, K, alpha):
    """Zerlegung x = c0 + c1·α — gibt (c0, c1)."""
    # Element von NumberField ist Polynom in alpha — Koeffizienten sind Q-rational
    coeffs = list(element.polynomial())
    c0 = QQ(coeffs[0]) if len(coeffs) > 0 else QQ(0)
    c1 = QQ(coeffs[1]) if len(coeffs) > 1 else QQ(0)
    return c0, c1


def test_case(a_p, expected_sharp, expected_flat):
    """Verifiziert einen Fall und gibt Vergleich aus."""
    K, alpha, beta, A_inv2, M, log_sharp, log_flat = compute_halflog_first_column(a_p)

    sharp_c0, sharp_c1 = to_alpha_coords(log_sharp, K, alpha)
    flat_c0, flat_c1 = to_alpha_coords(log_flat, K, alpha)

    exp_s_c0, exp_s_c1 = expected_sharp
    exp_f_c0, exp_f_c1 = expected_flat

    match_sharp = (sharp_c0 == exp_s_c0 and sharp_c1 == exp_s_c1)
    match_flat = (flat_c0 == exp_f_c0 and flat_c1 == exp_f_c1)

    print(f"  a_p = {a_p:+d}, K = Q(α) mit α² - ({a_p})·α + {P} = 0")
    print(f"  H(0) = A^-2 = {A_inv2.list()}")
    print(f"  Halblog-Matrix erste Spalte:")
    print(f"    log_α^♯|_0 = {log_sharp}")
    print(f"               = ({sharp_c0}) + ({sharp_c1})·α")
    print(f"               Erwartet: ({exp_s_c0}, {exp_s_c1})  {'✓' if match_sharp else '✗'}")
    print(f"    log_α^♭|_0 = {log_flat}")
    print(f"               = ({flat_c0}) + ({flat_c1})·α")
    print(f"               Erwartet: ({exp_f_c0}, {exp_f_c1})  {'✓' if match_flat else '✗'}")
    print()
    return match_sharp and match_flat


def main():
    print("=" * 70)
    print("Direkter algebraischer Beweis von Korollar 3.4")
    print("(unabhängig von SageMath's Dp_valued_series)")
    print("=" * 70)
    print()

    # Erwartete Werte aus Korollar 3.4
    EXPECTED = {
        # a_p: ((sharp_c0, sharp_c1), (flat_c0, flat_c1))
        -2: ((QQ(0), QQ(1)/4), (QQ(-1)/2, QQ(1)/4)),   # ♯ = α/4, ♭ = -1/2 + α/4
        0:  ((QQ(-1)/2, QQ(0)), (QQ(0), QQ(-1)/4)),    # ♯ = -1/2, ♭ = -α/4
        2:  ((QQ(0), QQ(-1)/4), (QQ(1)/2, QQ(1)/4)),   # ♯ = -α/4, ♭ = 1/2 + α/4
    }

    all_match = True
    for a_p in [-2, 0, 2]:
        print(f"--- Fall a_p = {a_p:+d} ---")
        exp_sharp, exp_flat = EXPECTED[a_p]
        ok = test_case(a_p, exp_sharp, exp_flat)
        if not ok:
            all_match = False
        print()

    print("=" * 70)
    if all_match:
        print("KOROLLAR 3.4 ALGEBRAISCH BIT-GENAU VERIFIZIERT")
        print("(alle 6 Werte stimmen mit Paper-Tabelle in Korollar 3.4 überein)")
    else:
        print("DISKREPANZEN gefunden — Korollar 3.4 oder Skript prüfen!")
    print("=" * 70)


if __name__ == "__main__":
    main()
