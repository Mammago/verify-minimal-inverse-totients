# Verification Results

All 30 rows of Table 1 independently verified. Every column is correct.

## Column-by-column summary

| Column | Matches | Method |
|--------|---------|--------|
| $m$ (factorization) | 30/30 | Reconstruct from prime factorization |
| $N_{\min}(m)$ | 30/30 | Independent sieve of $\varphi(n)$ for $n \leq 10^7$ |
| $N_{\min}/m$ | 30/30 | Match to 3 decimal places after rounding |
| odd? | 30/30 | All 30 $N_{\min}$ values are odd |
| $\|\mathcal{A}\|$ | 30/30 | Atom sets built from Definition 3.1 |
| $\|S\|$ | 30/30 | Odd prime supports counted directly |
| tw$\leq$ | 30/30 | Greedy min-degree elimination upper bounds |
| LP | 30/30 | LP bounds match to 3 decimal places |
| gap | 30/30 | Integrality gaps match to 1 decimal place |

## Minimality verification

For each $m$, confirmed no $n < N_{\min}(m)$ satisfies $\varphi(n) = m$:

- **Small cases** ($m \leq 32768$): Exhaustive search through all $n$ from 1 to $N_{\min} + 100$
- **Large cases** ($m$ up to ~2M): Full sieve verification up to 5,000,000
- **Largest case** ($m = 2{,}037{,}248$, $N_{\min} = 4{,}158{,}795$): Verified by sieve and direct computation $\varphi(4{,}158{,}795) = 2{,}037{,}248$

## Detailed LP results

For every instance, both LP relaxations were solved:

- **Even LP**: $\min \sum c(a) x_a$ with $\sum w(a) x_a \leq B$ (inequality)
- **Odd LP**: $\min \sum c(a) x_a$ with $\sum w(a) x_a = B$ (equality)

The parity-correct bound $\log N_{\min}(m) \geq \min\{(B+1)\log 2 + \text{even\_opt},\; B\log 2 + \text{odd\_opt}\}$ was then converted to a ratio bound and compared against the exact ratio.

In all 30 cases, the odd LP provides the tighter bound (consistent with all minimizers being odd).

## Environment

- Python 3.10+
- SciPy HiGHS LP solver
- NetworkX for graph construction and treewidth estimation
- Ubuntu 24.04
