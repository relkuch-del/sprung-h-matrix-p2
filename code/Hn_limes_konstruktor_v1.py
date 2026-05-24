#!/usr/bin/env python3
"""
H_n-Limes-Konstruktor v1 — Saubere Implementierung der Sprung-H-Matrix bei X=0.

Mathematischer Hintergrund:
  Sprung 2012 Definition 6.8 / Sprung 2016 Theorem 2.2:

    H(X) = lim_{n→∞} C_1(X) · C_2(X) · ... · C_n(X) · A^{-(n+2)} · B

  mit
    C_k(X) = [[a_p, 1], [-Phi_{p^k}(1+X), 0]]
    A      = [[a_p, 1], [-p, 0]]
    B      = [[-1, -1], [beta, alpha]]
  und alpha, beta = Wurzeln des Char-Polynoms X² - a_p·X + p.

  Bei p = 2:
    a_p = -2: alpha = -1+i, beta = -1-i, Q(alpha) = Q(i)
    a_p =  0: alpha = i·sqrt(2), beta = -i·sqrt(2), Q(alpha) = Q(i·sqrt(2))

Ziel dieses Skripts:
  (a) H_n(0) · A^{-(n+2)} · B für n = 0, 1, ..., 15 ausrechnen.
  (b) Zeigen: bei X = 0 ist Phi_{p^k}(1) = p für k ≥ 1, also C_k(0) = A.
      Daher H_n(0) = A^n und damit H_n(0) · A^{-(n+2)} = A^{-2} unabhängig von n.
      Der Limes ist ALGEBRAISCH GESCHLOSSEN: H(0) = A^{-2} · B.
  (c) Aus H(0) die Sprung-Halblogarithmen-Konstanten in der Sprung-(1, alpha)-Basis ablesen.
  (d) Basiswechsel zur Sage-(omega, phi·omega)-Basis via epsilon = (I - phi)^2.
  (e) Vergleich mit empirischen Sage-Werten u_Sage = (7/4, 3/4) bei a_p = -2
      bzw. v_Sage = (1, 2) bei a_p = 0.

Bilanz:
  Der Limes ist TRIVIAL (n-unabhängig nach Setup). Die eigentliche Substanz steckt
  in (i) der Basiswechsel-Brücke Sage ↔ Sprung, und (ii) der Identifikation,
  welcher Eintrag von H(0) die ♯/♭-Halblogarithmen liefert.
"""

import sympy as sp
from sympy import I, Rational, Matrix, sqrt, simplify, expand, eye, Symbol, cyclotomic_poly


def make_A(a_p, p):
    """Sprung-Frobenius-Matrix A = [[a_p, 1], [-p, 0]]."""
    return Matrix([[a_p, 1], [-p, 0]])


def make_C_at_zero(a_p, p, k):
    """C_k(0) = [[a_p, 1], [-Phi_{p^k}(1), 0]]. Bei k >= 1 ist Phi_{p^k}(1) = p."""
    if k == 0:
        # Konvention: Phi_1(1+X)|_{X=0} = Phi_1(1) = 0, also C_0 müsste degeneriert sein.
        # Sprung 2012 startet das Produkt bei k=1, daher hier nur Information.
        return None
    phi_val = sp.cyclotomic_poly(p ** k, 1)  # Phi_{p^k}(1) auswerten
    return Matrix([[a_p, 1], [-phi_val, 0]])


def make_B(a_p, p):
    """B = [[-1, -1], [beta, alpha]] mit alpha, beta = Wurzeln von X² - a_p·X + p."""
    # Wir parametrisieren über alpha symbolisch und nutzen alpha + beta = a_p, alpha*beta = p
    alpha = Symbol('alpha')
    # beta = a_p - alpha
    beta = a_p - alpha
    return Matrix([[-1, -1], [beta, alpha]]), alpha


def make_B_concrete_ap_minus2():
    """B bei a_p=-2: alpha = -1+i, beta = -1-i."""
    return Matrix([[-1, -1], [-1 - I, -1 + I]])


def make_B_concrete_ap_zero():
    """B bei a_p=0, p=2: alpha = i·sqrt(2), beta = -i·sqrt(2)."""
    s = I * sqrt(2)
    return Matrix([[-1, -1], [-s, s]])


def make_phi_sage(a_p, p):
    """Sage's Dieudonne-Frobenius auf der (omega, phi·omega)-Basis.

    Sage-Konvention: phi = [[0, -1/p], [1, a_p/p]].
    """
    return Matrix([[0, Rational(-1, p)], [1, Rational(a_p, p)]])


def make_epsilon_sage(a_p, p):
    """epsilon = (I - phi)^2 in Sage-Basis."""
    I2 = eye(2)
    phi = make_phi_sage(a_p, p)
    return (I2 - phi) ** 2


def verify_n_independence(a_p, p, n_max=15):
    """Zeige: H_n(0) · A^{-(n+2)} ist n-unabhängig (algebraische Trivialitaet)."""
    A = make_A(a_p, p)
    A_inv = A.inv()
    print(f"\n=== Verify n-independence at p={p}, a_p={a_p} ===")
    print(f"A = {A.tolist()},  det(A) = {A.det()},  tr(A) = {A.trace()}")

    reference = A_inv * A_inv  # A^{-2}
    print(f"\nReferenz A^(-2) = ")
    for row in reference.tolist():
        print(f"  {row}")

    for n in range(1, n_max + 1):
        # H_n(0) = C_1(0) · C_2(0) · ... · C_n(0) = A^n  (da C_k(0) = A für k >= 1)
        H_n = eye(2)
        for k in range(1, n + 1):
            C_k = make_C_at_zero(a_p, p, k)
            H_n = H_n * C_k

        # H_n(0) · A^{-(n+2)}
        Hn_normalized = H_n * (A_inv ** (n + 2))

        # Vergleich mit Referenz A^{-2}
        diff = simplify(Hn_normalized - reference)
        is_equal = diff == Matrix([[0, 0], [0, 0]])
        marker = "OK" if is_equal else "MISMATCH"
        print(f"  n={n:2d}: H_n(0)·A^(-(n+2)) == A^(-2)  [{marker}]")
        if not is_equal:
            print(f"       diff = {diff.tolist()}")

    return reference


def compute_H_zero_full(a_p, p):
    """Vollständige H(0) = A^{-2} · B Berechnung bei p=2."""
    print(f"\n=== Compute H(0) = A^(-2) · B at p={p}, a_p={a_p} ===")

    A = make_A(a_p, p)
    A_inv_2 = A.inv() ** 2

    if a_p == -2:
        B = make_B_concrete_ap_minus2()
        print("B (a_p=-2, alpha=-1+i):")
    elif a_p == 0:
        B = make_B_concrete_ap_zero()
        print("B (a_p=0, alpha=i·sqrt(2)):")
    else:
        raise ValueError(f"a_p={a_p} not handled")

    for row in B.tolist():
        print(f"  {row}")
    print(f"det(B) = {expand(B.det())}")

    H = simplify(A_inv_2 * B)
    print(f"\nH(0) = A^(-2) · B:")
    for row in H.tolist():
        print(f"  {row}")
    print(f"det(H(0)) = {simplify(H.det())}")

    return H


def basis_change_sprung_to_sage(H_sprung, a_p, p):
    """Basiswechsel: Sprung-(1, alpha)-Basis -> Sage-(omega, phi·omega)-Basis.

    Die Brücke ist epsilon = (I - phi)^2.
    H_sage = epsilon · H_sprung   (oder Transponiert, je nach Konvention)

    Wir testen beide Konventionen und vergleichen mit empirischen Sage-Werten.
    """
    eps = make_epsilon_sage(a_p, p)
    print(f"\n=== Basiswechsel bei a_p={a_p}, p={p} ===")
    print(f"phi (Sage) = {make_phi_sage(a_p, p).tolist()}")
    print(f"epsilon = (I-phi)^2 = ")
    for row in eps.tolist():
        print(f"  {row}")
    print(f"det(epsilon) = {expand(eps.det())}")

    H_sage_v1 = simplify(eps * H_sprung)
    H_sage_v2 = simplify(H_sprung * eps)
    H_sage_v3 = simplify(eps.inv() * H_sprung)
    H_sage_v4 = simplify(H_sprung * eps.inv())

    print(f"\nKandidat 1: eps · H(0) =")
    for row in H_sage_v1.tolist():
        print(f"  {row}")

    print(f"\nKandidat 2: H(0) · eps =")
    for row in H_sage_v2.tolist():
        print(f"  {row}")

    print(f"\nKandidat 3: eps^(-1) · H(0) =")
    for row in H_sage_v3.tolist():
        print(f"  {row}")

    print(f"\nKandidat 4: H(0) · eps^(-1) =")
    for row in H_sage_v4.tolist():
        print(f"  {row}")

    return H_sage_v1, H_sage_v2, H_sage_v3, H_sage_v4


def main():
    print("=" * 70)
    print("H_n-LIMES-KONSTRUKTOR v1")
    print("=" * 70)

    # ========== Fall a_p = -2 ==========
    print("\n" + "=" * 70)
    print("FALL: a_p = -2, p = 2")
    print("=" * 70)

    ref_minus2 = verify_n_independence(a_p=-2, p=2, n_max=15)
    H_sprung_minus2 = compute_H_zero_full(a_p=-2, p=2)
    candidates_minus2 = basis_change_sprung_to_sage(H_sprung_minus2, a_p=-2, p=2)

    print(f"\n--- EMPIRISCHE Sage-Werte zum Vergleich ---")
    print(f"u_Sage = (7/4, 3/4) bei a_p=-2 (V68/V69)")
    print(f"u_Sprung erwartet via ε^-T: (4/5, 9/10) gemäß SPRUNG_2012_U_IDENTIFIKATION.md")

    # ========== Fall a_p = 0 ==========
    print("\n" + "=" * 70)
    print("FALL: a_p = 0, p = 2")
    print("=" * 70)

    ref_zero = verify_n_independence(a_p=0, p=2, n_max=15)
    H_sprung_zero = compute_H_zero_full(a_p=0, p=2)
    candidates_zero = basis_change_sprung_to_sage(H_sprung_zero, a_p=0, p=2)

    print(f"\n--- EMPIRISCHE Sage-Werte zum Vergleich ---")
    print(f"v_Sage = (1, 2) bei a_p=0 (V67-V69)")
    print(f"v_Sprung erwartet via ε^-T: (-2/3, 4/3) gemäß SPRUNG_2012_U_IDENTIFIKATION.md")

    print("\n" + "=" * 70)
    print("ABSCHLUSS")
    print("=" * 70)
    print("""
Strukturelle Beobachtungen:

1. Bei X = 0 ist die H_n · A^(-(n+2))-Folge n-UNABHÄNGIG, sofort gleich A^(-2).
   Der Limes ist algebraisch trivial.

2. H(0) = A^(-2) · B liefert komplexe Werte in Q(alpha):
   - a_p = -2: in Q(i)
   - a_p = 0:  in Q(i·sqrt(2))

3. Der Bezug zur Sage-(omega, phi·omega)-Basis ist über (I - phi)^2 vermittelt.
   Welche der 4 Kandidaten-Formeln die richtige ist, klärt der Vergleich mit
   den empirischen Sage-Werten (7/4, 3/4) bzw. (1, 2).

4. NOCH ZU KLAEREN: Sprungs Halblogarithmen sind 1D-Vektoren (log^♯, log^♭),
   nicht 2x2-Matrizen. Die Identifikation, welcher SPALTEN- oder ZEILEN-Vektor
   von H(0) der ♯/♭-Halblogarithmus ist, kommt aus Sprung 2012 Theorem 1.1.
""")


if __name__ == "__main__":
    main()
