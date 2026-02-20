# verify-minimal-inverse-totients

Independent computational verification of the results in *"Minimal Inverse Totients via Assigned Valuation Covers, Shadow Prices, and Bounded-Width Structure"*.

## What this verifies

The paper develops an optimization framework for computing $N_{\min}(m) = \min\\{n \in \mathbb{N} : \varphi(n) = m\\}$ and presents Table 1 with 30 computed instances. This repository independently recomputes and checks **every value in all 9 columns of all 30 rows**:

| Column | Method |
|--------|--------|
| $m$ (factorization) | Verify prime factorization reconstructs to $m$ |
| $N_{\min}(m)$ | Sieve $\varphi(n)$ for $n \leq 10^7$; exhaustive minimality checks |
| $N_{\min}/m$ | Direct division, rounded to 3 decimal places |
| odd? | Parity check on $N_{\min}(m)$ |
| $\|\mathcal{A}\|$ | Build atom set per Definition 3.1 |
| $\|S\|$ | Count odd primes $s$ with $v_s(m) > 0$ |
| tw$\leq$ | Greedy min-degree elimination heuristic on $H(m)$ |
| LP | Solve both parity-correct LP relaxations (Theorem 4.1) via HiGHS |
| gap | Compute integrality gap per Section B.5 formula |

**Result: all 30 rows match across all 9 columns. No errors found.**

## Scripts

### `scripts/verify_table1.py`

The main verification script. Sieves $\varphi(n)$ for $n \leq 10^7$, builds atom sets $\mathcal{A}(m)$ and incidence graphs $H(m)$, estimates treewidth, solves both LP relaxations, and checks all columns against the paper's reported values.

```bash
python scripts/verify_table1.py
```

Runtime: ~3 minutes (dominated by the totient sieve).

### `scripts/verify_lp.py`

Focused verification of the LP bounds and gap values. Solves both parity-correct LP relaxations (even $\sum w \cdot x \leq B$ and odd $\sum w \cdot x = B$), computes the ratio lower bound and integrality gap, and prints a detailed atom-by-atom analysis for selected instances ($m = 2, 8, 12, 24, 32768, 5888, 86528$).

```bash
python scripts/verify_lp.py
```

### `scripts/verify_final.py`

Cross-checks and minimality proofs:
- Verifies $\varphi(N_{\min}(m)) = m$ directly for all 30 rows
- Confirms all $N_{\min}$ values are odd
- Validates all factorizations against the paper
- Checks $N_{\min}/m$ ratios to 3 decimal places
- Runs exhaustive minimality checks (no $n < N_{\min}$ with $\varphi(n) = m$)

```bash
python scripts/verify_final.py
```

## Requirements

```
numpy
scipy
networkx
```

Install with:

```bash
pip install -r requirements.txt
```

Tested with Python 3.10+ on Ubuntu 24.04. The LP solver used is SciPy's HiGHS interface (`scipy.optimize.linprog` with `method='highs'`).

## Key observations from verification

1. **Every minimizer is odd.** In all 30 instances, $W(x) = B$: the odd atoms exactly exhaust the 2-adic budget, so the odd preimage is valid and strictly smaller than the even-lift $2 \cdot n_{\text{odd}}$.

2. **LP gaps are small for structured instances.** When $|S| = 0$ (pure powers of 2) or treewidth is low, the LP bound is tight (gap = 0.0%). The largest gaps (~7.3%) occur for $m = 2^k \cdot p$ with large prime $p$.

3. **The $m = 86528$ rendering.** The compiled PDF shows "132" in the factorization column. This is a rendering artifact for $13^2$ — confirmed since $2^9 \times 13^2 = 512 \times 169 = 86528$.

## Implementation notes

The atom set $\mathcal{A}(m)$ is constructed per Definition 3.1: for each odd prime $p$ with $p - 1 \mid m$, atoms $(p, e)$ with $e \in E_p$ where $E_p = \\{1\\}$ if $p \notin S(m)$ and $E_p = \\{1, \ldots, b_p + 1\\}$ if $p \in S(m)$.

The incidence graph $H(m)$ is bipartite (Definition 6.1): left vertices are atoms, right vertices are constraints (valuation equalities $\sum A_s(a) x_a = b_s$ and per-prime inequalities $\sum_{e \in E_p} x_{p,e} \leq 1$). Treewidth is estimated by greedy min-degree elimination, which gives an upper bound.

Both LP relaxations are solved as stated in Theorem 4.1:
- **Even LP:** $\min \sum c(a) x_a$ subject to valuation equalities, per-prime constraints, $\sum w(a) x_a \leq B$, and $0 \leq x \leq 1$
- **Odd LP:** Same but with $\sum w(a) x_a = B$

The final bound is $\log N_{\min}(m) \geq \min\\{(B+1)\log 2 + \text{even\_opt},\; B \log 2 + \text{odd\_opt}\\}$.

## License

MIT
