#!/usr/bin/env python3
"""
H_n-Limes-Konstruktor v2 — Systematische Suche der Brücken-Konvention.

Ausgangspunkt v1:
  H(0) = A^(-2) · B exakt bestimmt (n-unabhängig).
  Bei a_p=-2: H(0) ist 2x2 über Q(i) mit Eintraegen in {-i/2, i/2, 1/2±i/2}.
  Bei a_p=0:  H(0) ist 2x2 über Q(i·sqrt(2)).

  Empirisch (Sage): u_Sage = (7/4, 3/4) bei a_p=-2, v_Sage = (1, 2) bei a_p=0.

Aufgabe v2:
  Identifiziere die korrekte Brücken-Operation, die H(0) auf die Sage-Werte abbildet.

Strategie:
  Bei a_p=-2 ist H(0) komplex. Eine reale 2D-Reduktion von Q(i)^2 nach Q^2 macht
  jeder Eintrag c = re + im*i zu einem 2D-Vektor (re, im) oder zur Matrix
  [[re, -im], [im, re]] (regular representation of Q(i) on Q^2).

  Konkret: jede 2x2-Matrix M ueber Q(i) entspricht einer 4x4-Matrix ueber Q,
  wenn man Q(i) als 2D-Q-Vektorraum mit Basis (1, i) auffasst.

  Sage liefert L^♯, L^♭ in einem Q_p-2-dim Raum mit Basis (omega, phi*omega).
  Die Halblogarithmen log_alpha^♯, log_alpha^♭ leben in Q_p(alpha), so ist
  jeder zwei Q_p-Eintraege wert. u = (7/4, 3/4) ist genau dieser Q_p-Vektor.

Test 1:
  Reduziere H(0)[i,j] = re + im*alpha zur Q-Repraesentation (re, im).
  Vergleiche mit (7/4, 3/4) und permutationen.

Test 2:
  Welche Linearkombination der Eintraege von H(0) gibt (7/4, 3/4)?
"""

import sympy as sp
from sympy import I, Rational, Matrix, sqrt, simplify, expand, eye, Symbol, im, re


def make_A(a_p, p):
    return Matrix([[a_p, 1], [-p, 0]])


def make_B_concrete_ap_minus2():
    return Matrix([[-1, -1], [-1 - I, -1 + I]])


def make_B_concrete_ap_zero():
    s = I * sqrt(2)
    return Matrix([[-1, -1], [-s, s]])


def make_phi_sage(a_p, p):
    return Matrix([[0, Rational(-1, p)], [1, Rational(a_p, p)]])


def make_eps(a_p, p):
    return (eye(2) - make_phi_sage(a_p, p)) ** 2


def alpha_re_im(c, a_p):
    """Zerlege c = re + im*alpha in Q-Komponenten.

    Bei a_p=-2: alpha = -1+i. Schreibe c = u + v*i, dann c = (u+v) + v*alpha
    (weil i = alpha + 1).
    Bei a_p=0: alpha = i*sqrt(2). Schreibe c = u + v*i. Dann v = alpha*coeff/sqrt(2),
    also c = u + (v/sqrt(2))*alpha.
    """
    expanded = sp.expand(c)
    u = re(expanded)
    v = im(expanded)
    if a_p == -2:
        # i = alpha + 1, also c = u + v*i = u + v*(alpha+1) = (u+v) + v*alpha
        return (u + v, v)
    elif a_p == 0:
        # alpha = i*sqrt(2), also i = alpha/sqrt(2)
        # c = u + v*i = u + v*alpha/sqrt(2) = u + (v/sqrt(2))*alpha
        return (u, simplify(v / sqrt(2)))
    else:
        raise ValueError


def print_alpha_zerlegung(M, a_p, label):
    """Zerlege jeden Eintrag M[i,j] = re + im*alpha und drucke (re, im)."""
    print(f"\n--- alpha-Zerlegung {label} bei a_p={a_p} ---")
    print("Jeder Eintrag M[i,j] = c0 + c1*alpha mit c0, c1 in Q:")
    for i in range(2):
        for j in range(2):
            c0, c1 = alpha_re_im(M[i, j], a_p)
            print(f"  M[{i},{j}] = {str(M[i,j]):>20} =  {c0} + ({c1})*alpha")


def systematic_test(a_p, p, target_u):
    """Teste systematisch Linearkombinationen / Spalten / Zeilen von H(0) und Varianten."""
    A = make_A(a_p, p)
    A_inv = A.inv()
    H = A_inv * A_inv * (make_B_concrete_ap_minus2() if a_p == -2 else make_B_concrete_ap_zero())
    H = sp.simplify(H)
    eps = make_eps(a_p, p)
    eps_inv = eps.inv()

    print(f"\n{'='*70}\nSystematic test bei a_p={a_p}, target u = {target_u}\n{'='*70}")

    print(f"\nH(0) (in Sprung-Basis, Q(alpha)):")
    for row in H.tolist():
        print(f"  {row}")

    print_alpha_zerlegung(H, a_p, "H(0)")

    # Kandidat: H(0) interpretiert ueber Q via alpha-Zerlegung jeder Spalte.
    # Spalte 0 von H(0): (M[0,0], M[1,0]). In alpha-Basis: 2 Eintraege, jeweils (c0, c1).
    # Das ist also ein 4D-Q-Vektor (c0_00, c1_00, c0_10, c1_10).
    print(f"\nSpalten von H(0) als 4D-Q-Vektor (über alpha-Basis):")
    for col in range(2):
        v00 = alpha_re_im(H[0, col], a_p)
        v10 = alpha_re_im(H[1, col], a_p)
        full = (v00[0], v00[1], v10[0], v10[1])
        print(f"  Spalte {col}: {full}")

    # Test ob (4/5, 9/10) als eps oder eps^-1 angewandt auf eine Spalte herauskommt
    print(f"\nDirekttest: eps und eps^-1 auf jede Spalte von H(0):")
    for col in range(2):
        v = Matrix([H[0, col], H[1, col]])
        for name, M in [("eps", eps), ("eps^-1", eps_inv), ("eps^T", eps.T), ("eps^-T", eps_inv.T)]:
            res = sp.simplify(M * v)
            print(f"  Spalte {col} unter {name}: {res.T.tolist()}")

    # Direkter Test: ist (7/4, 3/4) in der u-Basis irgendwo direkt sichtbar?
    print(f"\nVergleich mit u_Sage = {target_u}:")
    print(f"  Faktor 7/4 = {Rational(7,4)} = {float(Rational(7,4))}")
    print(f"  Faktor 3/4 = {Rational(3,4)} = {float(Rational(3,4))}")
    print(f"  Faktor 4/5 = {Rational(4,5)} = {float(Rational(4,5))}")
    print(f"  Faktor 9/10 = {Rational(9,10)} = {float(Rational(9,10))}")

    # Reverse-Engineering: was muss eps mal etwas sein, damit target_u herauskommt?
    target_vec = Matrix(list(target_u))
    print(f"\nReverse: x = eps^-1 * {target_u} (Sage->Sprung):")
    x = sp.simplify(eps_inv * target_vec)
    print(f"  {x.T.tolist()}")

    print(f"\nReverse: x = eps * {target_u}:")
    x2 = sp.simplify(eps * target_vec)
    print(f"  {x2.T.tolist()}")

    print(f"\nReverse: x = eps^-T * {target_u}:")
    x3 = sp.simplify(eps_inv.T * target_vec)
    print(f"  {x3.T.tolist()}")

    print(f"\nReverse: x = eps^T * {target_u}:")
    x4 = sp.simplify(eps.T * target_vec)
    print(f"  {x4.T.tolist()}")


def main():
    print("=" * 70)
    print("H_n-LIMES-KONSTRUKTOR v2 — BRUECKEN-SUCHE")
    print("=" * 70)

    systematic_test(a_p=-2, p=2, target_u=(Rational(7, 4), Rational(3, 4)))
    systematic_test(a_p=0, p=2, target_u=(Rational(1), Rational(2)))


if __name__ == "__main__":
    main()
