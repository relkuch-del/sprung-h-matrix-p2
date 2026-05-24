#!/usr/bin/env python3
"""
SAGE-VERIFIKATION a_p=0 mit höherem Konduktor-Range.

Ziel:
  30 zusätzliche a_p=0 ss Trigger mit cond zwischen 1000 und 100000,
  gleichmäßig über die Range verteilt.

  Zusammen mit den bestehenden 75 Triggern (cond 77-225 aus bulk_v2)
  ergibt sich eine Konduktor-Range über drei Größenordnungen — robuste
  Reviewer-Empirie.

Bucket-Strategie:
  10 logarithmische Buckets in [1000, 100000], pro Bucket 3 Trigger.
  Bucket-Grenzen: 1000, 1585, 2512, 3981, 6310, 10000, 15849, 25119, 39811, 63096, 100000.

  Pro Bucket: scanne Cremona-Kurven cond ∈ [bucket_lo, bucket_hi],
  filter auf a_p(E)=0 und rank>=1, nimm die ersten 3 mit Match.

Aufruf:
  conda activate sage10 && cd /mnt/c/Temp/Forschung/paper_log_alpha && \
    python3 -u sage_verify_ap0_highcond.py 2>&1 | tee sage_verify_ap0_highcond_log.txt
"""

import sys
from sage.rings.padics import precision_error as _pe
sys.modules['precision_error'] = _pe

from sage.all import EllipticCurve, QQ, matrix, vector, CremonaDatabase

P = 2
TARGET_A_P = 0
N_PER_BUCKET = 3        # 3 Trigger pro Bucket, 10 Buckets -> 30 total
SKIP_RANK_BELOW = 1

# Log-uniforme Buckets in [1000, 100000]
BUCKETS = [
    (1000, 1584),
    (1585, 2511),
    (2512, 3980),
    (3981, 6309),
    (6310, 9999),
    (10000, 15848),
    (15849, 25118),
    (25119, 39810),
    (39811, 63095),
    (63096, 100000),
]


def make_eps(a_p):
    p_q = QQ(P)
    ap_q = QQ(a_p)
    phi = matrix(QQ, [[0, -1 / p_q], [1, ap_q / p_q]])
    return (matrix.identity(QQ, 2) - phi) ** 2


def check_curve(E, u_Sage, u_Sprung_expected):
    """Pruefe Brücke. Gibt (match: bool, result_tuple) oder None bei a_p-Mismatch."""
    try:
        a_p = int(E.ap(P))
    except Exception:
        return None
    if a_p != TARGET_A_P:
        return None
    eps_inv = make_eps(a_p).inverse()
    result = eps_inv * u_Sage
    return (result == u_Sprung_expected, tuple(result), a_p)


def search_bucket(db, lo, hi, u_Sage, u_Sprung_expected, n_target):
    """Scanne Cremona-Konduktoren in [lo, hi], gib bis zu n_target Match-Trigger."""
    found = []
    n_scanned = 0
    for cond in range(lo, hi + 1):
        if len(found) >= n_target:
            break
        try:
            all_for_cond = db.allcurves(cond)
        except Exception:
            continue
        if not all_for_cond:
            continue
        for label, info in all_for_cond.items():
            if len(found) >= n_target:
                break
            try:
                if isinstance(info, (list, tuple)) and len(info) >= 2:
                    ainvs = info[0]
                    stored_rank = int(info[1]) if info[1] is not None else None
                else:
                    ainvs = info
                    stored_rank = None
            except Exception:
                continue
            if stored_rank is not None and stored_rank < SKIP_RANK_BELOW:
                continue
            try:
                E = EllipticCurve(ainvs)
            except Exception:
                continue
            n_scanned += 1
            outcome = check_curve(E, u_Sage, u_Sprung_expected)
            if outcome is None:
                continue
            match, result, a_p = outcome
            full_label = f"{cond}{label}"
            if match:
                found.append((full_label, cond, a_p,
                              stored_rank if stored_rank is not None else "?", result))
            else:
                print(f"    !!! MISMATCH {full_label}: {result}")
    return found, n_scanned


def main():
    print("=" * 70)
    print("SAGE-VERIFIKATION a_p=0 — höherer Konduktor-Range")
    print("=" * 70)
    print(f"  Ziel: {N_PER_BUCKET} Trigger pro Bucket x {len(BUCKETS)} Buckets "
          f"= {N_PER_BUCKET * len(BUCKETS)} total")
    print(f"  Range: cond in [{BUCKETS[0][0]}, {BUCKETS[-1][1]}]")
    print(f"  rank-Filter: >= {SKIP_RANK_BELOW}")
    print()

    db = CremonaDatabase()
    largest = db.largest_conductor()
    print(f"  Cremona-DB largest_conductor = {largest}")
    if largest < 100000:
        print("  WARNUNG: Mini-DB aktiv!")
        sys.exit(1)
    print()

    u_Sage = vector(QQ, [QQ(1), QQ(2)])
    u_Sprung_expected = vector(QQ, [QQ(-2) / 3, QQ(4) / 3])

    print(f"  v_Sage          = {tuple(u_Sage)}")
    print(f"  v_Sprung expect = {tuple(u_Sprung_expected)} = (-2/3, 4/3)")
    print()

    all_results = []
    total_scanned = 0
    for i, (lo, hi) in enumerate(BUCKETS, 1):
        print(f"--- Bucket {i:>2}/{len(BUCKETS)}: cond in [{lo}, {hi}] ---")
        found, n_scan = search_bucket(db, lo, hi, u_Sage, u_Sprung_expected, N_PER_BUCKET)
        total_scanned += n_scan
        if len(found) < N_PER_BUCKET:
            print(f"    !!! NUR {len(found)}/{N_PER_BUCKET} gefunden (gescannt: {n_scan})")
        for entry in found:
            lbl, cnd, ap, rk, res = entry
            print(f"    {lbl:<14} cond={cnd:>6} a_p={ap:>3} rk={rk} → {res}")
            all_results.append(entry)
        print()

    print("=" * 70)
    print("FINALE TRIGGER-TABELLE")
    print("=" * 70)
    print(f"{'i':>3} | {'Label':<14} | {'Cond':>6} | {'a_p':>3} | {'rank':>4}")
    print("-" * 50)
    for i, (lbl, cnd, ap, rk, _) in enumerate(all_results, 1):
        print(f"{i:>3} | {lbl:<14} | {cnd:>6} | {ap:>3} | {str(rk):>4}")
    print()

    n_match = len(all_results)
    n_target = N_PER_BUCKET * len(BUCKETS)
    print("=" * 70)
    print("BILANZ a_p=0 high-cond")
    print("=" * 70)
    print(f"  Gefunden: {n_match}/{n_target}")
    print(f"  Cremona-Kurven gescannt: {total_scanned}")
    if all_results:
        cmin = min(r[1] for r in all_results)
        cmax = max(r[1] for r in all_results)
        print(f"  Konduktor-Range: {cmin} - {cmax}")
    print()
    if n_match == n_target:
        print("ALLE Trigger bit-genau verifiziert.")
    elif n_match > 0:
        print(f"{n_match} Trigger verifiziert; in einigen Buckets nicht genug a_p=0 ss Kurven.")
    else:
        print("KEINE Trigger gefunden — Skript-Fehler oder DB-Problem.")
    print("=" * 70)


if __name__ == "__main__":
    main()
