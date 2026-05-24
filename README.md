# sprung-h-matrix-p2

Verification code and data accompanying the paper:

**"Explicit evaluation of Sprung's half-logarithm $H$-matrix at $X=0$ for $p=2$ supersingular reduction"**
by Raphael Elkuch (2026).

## Overview

This repository contains the SageMath verification code used to confirm
Theorem 1.2 of the paper, which provides an explicit bridge between
SageMath's Dieudonné basis and Sprung's canonical basis for the
half-logarithms of elliptic curves with supersingular reduction at $p = 2$.

The numerical verification covers **180 elliptic curves** drawn from the
Cremona database, spanning conductors from 37 to 63,096 — four orders of
magnitude. All 180 curves confirm the bridge formula bit-exactly.

## Repository structure

```
sprung-h-matrix-p2/
├── code/                          # SageMath verification scripts
│   ├── Hn_limes_konstruktor_v1.py        # H_n-limit constructor (Theorem 3.1)
│   ├── Hn_limes_konstruktor_v2.py        # epsilon-bridge identification (sympy)
│   ├── sage_verify_17963c1_v2.py         # main anchor verification (17963c1)
│   ├── sage_verify_30487a1.py            # a_p=0 anchor verification (30487a1)
│   ├── sage_verify_bulk.py               # manual pool (14 triggers)
│   ├── sage_verify_bulk_v2.py            # automatic Cremona search (150 triggers)
│   ├── sage_verify_ap0_highcond.py       # a_p=0 high-conductor (30 triggers)
│   ├── sage_verify_paper_table1.py       # reproducer for Table 1 of the paper
│   └── pollack_cross_check_v1.py         # Pollack-Sprung Remark 6.16 cross-check
├── data/                          # Verification log files
│   ├── sage_verify_17963c1_v2_log.txt
│   ├── sage_verify_bulk_v2_log.txt        # full 150-trigger log
│   ├── sage_verify_ap0_highcond_log.txt   # 30 high-conductor log
│   └── sage_verify_paper_table1_log.txt   # Table 1 reproduction log
├── paper/
│   ├── sprung_halflog_bridge_p2.tex       # LaTeX source (current version)
│   └── sprung_halflog_bridge_p2.pdf       # compiled paper
├── LICENSE                        # MIT (for code)
├── LICENSE-DATA                   # CC-BY 4.0 (for data and paper)
├── CITATION.cff                   # citation metadata
└── README.md                      # this file
```

## Reproducing the verification

### Prerequisites

* **SageMath 10.x** with the full Cremona elliptic-curve database
  (largest conductor 499,998). See
  [Cremona's elliptic curve database](https://github.com/JohnCremona/ecdata)
  for installation instructions.
* Python 3.10+ with `sympy` (for the algebraic constructor).

### Reproducing Table 1 of the paper

The four-curve sample of Table 1 in `paper/sprung_halflog_bridge_p2.pdf` is
reproduced bit-exactly by the dedicated script:

```bash
sage code/sage_verify_paper_table1.py
```

Expected output: four lines, each ending with ✓, and the final summary
`ALLE 4 TRIGGER bit-genau verifiziert`. Computation time under one minute.

### Running the algebraic part (sympy only)

The algebraic content of Theorem 3.1 (triviality of the limit at $X=0$)
and the bridge of Theorem 3.5 is verified without SageMath:

```bash
python3 code/Hn_limes_konstruktor_v1.py
python3 code/Hn_limes_konstruktor_v2.py
```

The first script verifies that $H_n(0) \cdot A^{-(n+2)} = A^{-2}$ for
$n = 1, \ldots, 15$. The second isolates the bridge between SageMath's
Dieudonné basis and Sprung's canonical $(1, \alpha)$-basis explicitly.

### Running the bulk SageMath verification

The 180-curve consistency check of SageMath's `Dp_valued_series` requires
SageMath:

```bash
# Main anchor (17963c1, a_p = -2)
sage code/sage_verify_17963c1_v2.py

# a_p = 0 anchor (30487a1)
sage code/sage_verify_30487a1.py

# Auto bulk over Cremona DB: 75 a_p=-2 + 75 a_p=0
sage code/sage_verify_bulk_v2.py

# Logarithmically-spaced high-conductor a_p=0: 30 curves
sage code/sage_verify_ap0_highcond.py
```

Total compute time is approximately five minutes on a standard workstation.

### Pollack consistency check (sympy only)

```bash
python3 code/pollack_cross_check_v1.py
```

This documents the structural consistency between Sprung's normalisation
$A^{-(n+2)}$ and the divergent Pollack constants at $T = 0$, in line with
[Sprung 2012, Remark 6.16].

## Bilanz of numerical verification

| Class                              | n   | Conductor range   | Bit-exact |
|------------------------------------|----:|-------------------|----------:|
| $a_p = -2$ supersingular           | 75  | $[37, 997]$        | 75/75     |
| $a_p = 0$  supersingular (low)     | 75  | $[77, 225]$        | 75/75     |
| $a_p = 0$  supersingular (high)    | 30  | $[1001, 63096]$    | 30/30     |
| **Total**                          | **180** | $[37, 63096]$  | **180/180** |

Zero mismatches, zero errors. Compute time approximately 5 minutes.

## Citation

If you use this code or data, please cite the paper. See `CITATION.cff` for
machine-readable citation metadata, or use:

```bibtex
@article{Elkuch2026sprung,
  author = {Elkuch, Raphael},
  title  = {Explicit half-logarithm constants and the {SageMath--Sprung} bridge
            for {$p=2$} supersingular reduction},
  year   = {2026},
  note   = {Preprint, arXiv:[to be added]},
}
```

## License

* **Code** (in `code/`): MIT License — see `LICENSE`.
* **Data, logs, paper PDF, and LaTeX source** (in `data/` and `paper/`):
  Creative Commons Attribution 4.0 (CC-BY 4.0) — see `LICENSE-DATA`.

## Contact

Raphael Elkuch — `relkuch@gmail.com`
