#!/usr/bin/env python3
"""
SAGE-VERIFIKATION v2 an 17963c1 — KORRIGIERT.

Korrekturen gegenüber v1:
  (a) Rational(7,4) -> QQ(7)/4  (Sage's Rational akzeptiert nicht (numer, denom)
      wie sympy)
  (b) Statt naivem "Konstanten-Term von sp[i][0]" fragen wir nach
      DEN HALBLOGARITHMEN log_alpha^flat|_n und log_alpha^sharp|_n,
      extrahiert via Sprung-Hauptidentitaet aus den DPVAL-Reihen.

Methodik:
  Aus V58a-Linie (im Archiv): bei 17963c1 mit rank=2, V_0=3 gilt
    log^flat|_0  = -2 + (-7/4)*alpha   (G_3/L^flat_3, H_3/L^flat_3)
  wobei (G_k, H_k) = sage_vec_k * eps_inv_T.

  Das heisst: log^flat|_0 ist ein 2D-Q-Vektor (-2, -7/4) in der
  (1, alpha)-Basis des lokalen Koerpers Q_2(alpha) mit alpha^2 + 2alpha + 2 = 0.

  Frage: Welche Brücken-Operation auf welchen Sage-Vektor liefert (-2, -7/4)?

  Beobachtung: Das urspruengliche u_Sage = (7/4, 3/4) ist NICHT log^flat|_0,
  sondern u-als-Linearitaetsfaktor in log^flat|_n = c_n * u fuer n >= 1.

  Aus den Sage-DPVAL-Werten kann ich log^flat|_n iterativ bestimmen.
  Hier teste ich nur die nullte Stufe: ist (-2, -7/4) der erwartete Wert?

  Und dann: was ist (4/5, 9/10)?
"""

import sys
from sage.rings.padics import precision_error as _pe
sys.modules['precision_error'] = _pe

from sage.all import EllipticCurve, QQ, matrix, vector


P = 2
LABEL = "17963c1"
N_PARAM = 13
PREC = 30


def to_Q(x):
    if x == 0:
        return QQ(0)
    try:
        return QQ(x.lift())
    except Exception:
        pass
    try:
        return QQ(x)
    except Exception:
        return QQ(0)


def main():
    print("=" * 70)
    print(f"SAGE-VERIFIKATION v2: {LABEL}, p={P}")
    print("=" * 70)

    E = EllipticCurve(LABEL)
    a_p = int(E.ap(P))
    print(f"  Trigger: {LABEL}, a_2={a_p}, rank={E.rank()}, cond={E.conductor()}")
    print()

    p_q = QQ(P)
    ap_q = QQ(a_p)
    I2 = matrix.identity(QQ, 2)
    phi = matrix(QQ, [[0, -1 / p_q], [1, ap_q / p_q]])
    eps = (I2 - phi) ** 2
    eps_inv = eps.inverse()
    eps_inv_T = eps_inv.transpose()

    print("=== Sage-Basis-Matrizen ===")
    print(f"phi:     {phi.list()}")
    print(f"eps:     {eps.list()},  det = {eps.determinant()}")
    print(f"eps^-1:  {eps_inv.list()}")
    print(f"eps^-T:  {eps_inv_T.list()}")
    print()

    # DPVAL-Reihen
    Lp = E.padic_lseries(P)
    sp = Lp.Dp_valued_series(n=N_PARAM, prec=PREC)

    # Extrahiere (G_k, H_k) via eps_inv_T-Konvention (V58-Linie)
    print("=== Sage-DPVAL-Werte und Sprung-Konvertierung ===")
    print(f"{'k':>3} | {'sp[0][k]':>12} | {'sp[1][k]':>12} | "
          f"{'G_k = (sv*eps^-T)[0]':>22} | {'H_k = (sv*eps^-T)[1]':>22}")
    print("-" * 95)

    gh_list = []
    for k in range(N_PARAM):
        Ls = to_Q(sp[0][k])
        Lf = to_Q(sp[1][k])
        sage_vec = vector(QQ, [Ls, Lf])
        # V58-Konvention: lpv = sage_vec * eps_inv_T (Zeilen-Vektor * Matrix)
        lpv = sage_vec * eps_inv_T
        # In V58:  H_k = -lpv[1] / p_q,  G_k = lpv[0] - ap_q * H_k
        # Aber das ist nochmal eine weitere Transformation. Probieren wir
        # zuerst lpv direkt zu interpretieren.
        gh_list.append((lpv[0], lpv[1]))
        print(f"{k:>3} | {str(Ls):>12} | {str(Lf):>12} | "
              f"{str(lpv[0]):>22} | {str(lpv[1]):>22}")
    print()

    # V58-Originalformel
    print("=== V58-Original: (G,H) via H = -lpv[1]/p, G = lpv[0] - a_p*H ===")
    print(f"{'k':>3} | {'lpv[0]':>15} | {'lpv[1]':>15} | "
          f"{'G_k':>15} | {'H_k':>15}")
    print("-" * 80)
    GH_v58 = []
    for k in range(N_PARAM):
        Ls = to_Q(sp[0][k])
        Lf = to_Q(sp[1][k])
        sage_vec = vector(QQ, [Ls, Lf])
        lpv = sage_vec * eps_inv_T
        H_k = -lpv[1] / p_q
        G_k = lpv[0] - ap_q * H_k
        GH_v58.append((G_k, H_k))
        print(f"{k:>3} | {str(lpv[0]):>15} | {str(lpv[1]):>15} | "
              f"{str(G_k):>15} | {str(H_k):>15}")
    print()

    # log^flat|_0 = (G_3, H_3) / L_flat_3 in (1, alpha)-Basis
    print("=== log_alpha^flat|_0 via Sprung-Hauptidentitaet ===")
    print("V_0 = 3 (17963c1). log^flat|_0 = (G_3, H_3) / L^flat_3 in (1, alpha)-Basis")
    L_flat_3 = to_Q(sp[1][3])
    print(f"L^flat_3 = sp[1][3] = {L_flat_3}")
    G_3, H_3 = GH_v58[3]
    print(f"G_3 = {G_3},  H_3 = {H_3}")
    if L_flat_3 != 0:
        log_flat_0 = (G_3 / L_flat_3, H_3 / L_flat_3)
        print(f"log^flat|_0 = ({log_flat_0[0]}, {log_flat_0[1]})  in (1, alpha)-Basis")
        print(f"Archiv-Erwartung V58a: (-2, -7/4)")
        match = (log_flat_0[0] == QQ(-2)) and (log_flat_0[1] == QQ(-7) / 4)
        print(f"Match: {match}")
    else:
        print("L^flat_3 = 0, Division undefiniert.")
    print()

    # Universelle u-Einheit: log^flat|_n / c_n fuer n >= 1
    print("=== Test: u-Einheit (7/4, 3/4) vs. log^flat-Niederstufen ===")
    print("Hypothese V68/V69: log_alpha^flat|_n = c_n * (7/4, 3/4) fuer n >= 1")
    print("c_n trigger-spezifisch rational, (7/4, 3/4) universell.")
    print()
    # Wir extrahieren log^flat|_1, log^flat|_2, log^flat|_3 iterativ unter
    # der V58-Annahme log^sharp|_n = 0 fuer kleine n.
    # X^4: G_4 + H_4*alpha = L^flat_3 * log^flat|_1 + L^flat_4 * log^flat|_0
    # => log^flat|_1 = ((G_4, H_4) - L^flat_4 * log^flat|_0) / L^flat_3
    G_4, H_4 = GH_v58[4]
    L_flat_4 = to_Q(sp[1][4])
    log_flat_0_g = G_3 / L_flat_3
    log_flat_0_h = H_3 / L_flat_3
    log_flat_1_g = (G_4 - L_flat_4 * log_flat_0_g) / L_flat_3
    log_flat_1_h = (H_4 - L_flat_4 * log_flat_0_h) / L_flat_3
    print(f"log^flat|_1 = ({log_flat_1_g}, {log_flat_1_h})")
    # Wenn log^flat|_1 = c_1 * (7/4, 3/4), dann log^flat|_1[0]/log^flat|_1[1] = 7/3
    if log_flat_1_h != 0:
        ratio = log_flat_1_g / log_flat_1_h
        print(f"Verhaeltnis log^flat|_1[0]/log^flat|_1[1] = {ratio}")
        print(f"Erwartet (wenn u=(7/4,3/4)): 7/3 = {QQ(7)/3}")
    else:
        print("log^flat|_1[1] = 0, ratio undefiniert.")
    print()

    # Sprung-Brücken-Test mit korrekter Eingabe
    print("=" * 70)
    print("BRUECKEN-TEST: u_Sage -> u_Sprung via eps^-1")
    print("=" * 70)
    u_Sage = vector(QQ, [QQ(7) / 4, QQ(3) / 4])
    print(f"u_Sage = {tuple(u_Sage)} = ({float(u_Sage[0])}, {float(u_Sage[1])})")
    print()
    print(f"eps^-1 * u_Sage   = {tuple(eps_inv * u_Sage)}")
    print(f"eps^-T * u_Sage   = {tuple(eps_inv_T * u_Sage)}")
    print(f"u_Sage * eps^-1   = {tuple(u_Sage * eps_inv)}")
    print(f"u_Sage * eps^-T   = {tuple(u_Sage * eps_inv_T)}")
    print(f"eps * u_Sage      = {tuple(eps * u_Sage)}")
    print(f"u_Sage * eps      = {tuple(u_Sage * eps)}")
    print()
    print(f"Erwartet u_Sprung = (4/5, 9/10)")
    print()


if __name__ == "__main__":
    main()
