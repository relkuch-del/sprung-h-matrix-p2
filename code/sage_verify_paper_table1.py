#!/usr/bin/env python3
"""
sage_verify_paper_table1.py — Reproduzieres Tabelle 1 aus
sprung_halflog_bridge_p2.tex bit-genau.

Zielwerte (a_p = -2, am ersten nicht-trivialen Mazur-Tate-Niveau V_0):
  17963c1: V_0 = 3, (G, H) = (-2, -7/4)
  70449a1: V_0 = 4, (G, H) = (-2, -7/4)
  53461b1: V_0 = 5, (G, H) = (-1/4, -1)
  55935g1: V_0 = 5, (G, H) = (-2, -7/4)

Methode (Theorem 3.5 des Papers):
  1. sp = Lp.Dp_valued_series(n, prec)
  2. V_0 := erstes k mit sp[1][k] != 0  (♭-Vanishing-Niveau, Sprung-Konvention)
  3. phi = [[0, -1/p], [1, a_p/p]]
  4. lpv = sage_vec * ((1-phi)^2).transpose()
  5. H = -lpv[1] / p, G = lpv[0] - a_p * H
  6. Ausgabe: (V_0, G, H)

Aufruf:
  conda activate sage10 && cd /mnt/c/Temp/Forschung/paper_log_alpha && \
    python3 -u sage_verify_paper_table1.py 2>&1 | \
    tee sage_verify_paper_table1_log.txt
"""

import sys
from sage.rings.padics import precision_error as _pe
sys.modules['precision_error'] = _pe

from sage.all import EllipticCurve, QQ, matrix, vector

P = 2
N_PARAM = 13
PREC = 30

TRIGGERS = [
    # a_p = -2 (4 Trigger)
    "17963c1",  # rank=2, V_0=3 erwartet (G, H) = (-2, -7/4)
    "70449a1",  # V_0=4 erwartet (G, H) = (-2, -7/4)
    "53461b1",  # V_0=5 erwartet (G, H) = (-1/4, -1)
    "55935g1",  # V_0=5 erwartet (G, H) = (-2, -7/4)
    # a_p = 0 (3 Trigger — erwartete Werte beim ersten Lauf bestimmen)
    "30487a1",
    "60803a1",
    "58939a1",
]

EXPECTED = {
    # a_p = -2
    "17963c1": (3, QQ(-2), QQ(-7)/4),
    "70449a1": (4, QQ(-2), QQ(-7)/4),
    "53461b1": (5, QQ(-1)/4, QQ(-1)),
    "55935g1": (5, QQ(-2), QQ(-7)/4),
    # a_p = 0 (verifiziert 27.5.2026)
    "30487a1": (7, QQ(1), QQ(-1)/4),
    "60803a1": (7, QQ(3)/2, QQ(-3)/8),
    "58939a1": (3, QQ(1), QQ(-1)/4),
}


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


def extract_gh(label):
    """Extrahiere (V_0, G_{V_0}, H_{V_0}) für einen Trigger."""
    E = EllipticCurve(label)
    a_p = int(E.ap(P))
    p_q = QQ(P)
    ap_q = QQ(a_p)
    I2 = matrix.identity(QQ, 2)
    phi = matrix(QQ, [[0, -1 / p_q], [1, ap_q / p_q]])
    # (1 - phi)^2  als TRANSPONIERTE — wird auf den Sage-Zeilenvektor von rechts angewandt
    # (= "eps_inv_T" in V57-Schreibweise; eigentlich ((1-phi)^2)^T)
    bridge_T = ((I2 - phi) ** 2).transpose()

    Lp = E.padic_lseries(P)
    sp = Lp.Dp_valued_series(n=N_PARAM, prec=PREC)

    # V_0 := erstes k mit sp[1][k] != 0 (=♭-Vanishing-Niveau)
    # Sprung-Konvention: L_p^♭ verschwindet bis Stufe V_0-1, ist erstmals nicht-null bei V_0.
    V_0 = None
    for k in range(N_PARAM):
        Lf = to_Q(sp[1][k])
        if Lf != 0:
            V_0 = k
            break

    if V_0 is None:
        return (a_p, None, None, None)

    Ls = to_Q(sp[0][V_0])
    Lf = to_Q(sp[1][V_0])
    sage_vec = vector(QQ, [Ls, Lf])
    lpv = sage_vec * bridge_T
    H = -lpv[1] / p_q
    G = lpv[0] - ap_q * H
    return (a_p, V_0, G, H)


def main():
    print("=" * 70)
    print("Reproduzieren von Tabelle 1 aus sprung_halflog_bridge_p2.tex")
    print("=" * 70)
    print()
    print(f"{'Trigger':<12} | {'a_p':>3} | {'V_0':>3} | {'G':>10} | {'H':>10} | "
          f"{'Erwartet':>20} | {'Match':>5}")
    print("-" * 95)

    all_match = True
    n_match = 0
    n_discovery = 0
    for label in TRIGGERS:
        try:
            a_p, V_0, G, H = extract_gh(label)
            if V_0 is None:
                print(f"{label:<12} | {a_p:>3} | N/A | KEIN V_0 gefunden")
                all_match = False
                continue
            if label in EXPECTED:
                exp_V_0, exp_G, exp_H = EXPECTED[label]
                match = (V_0 == exp_V_0 and G == exp_G and H == exp_H)
                marker = "✓" if match else "✗"
                print(f"{label:<12} | {a_p:>3} | {V_0:>3} | {str(G):>10} | {str(H):>10} | "
                   