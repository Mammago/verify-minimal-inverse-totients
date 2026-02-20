#!/usr/bin/env python3
"""
Verify the LP bounds and gap values from Table 1.

For each m, we solve both parity-correct LP relaxations (Theorem 4.1):
  Even LP: min sum c(a)*x(a)  s.t. valuation eqs, per-prime ineqs, sum w*x <= B, 0<=x<=1
  Odd LP:  min sum c(a)*x(a)  s.t. valuation eqs, per-prime ineqs, sum w*x = B, 0<=x<=1

Then:
  log Nmin(m) >= min{ (B+1)log2 + even_opt, B*log2 + odd_opt }
  
LP bound on Nmin/m = exp(log_LP_bound - log(m))
gap = (log(Nmin/m) - log(LP_bound_on_ratio)) / log(Nmin)
"""

import math
from scipy.optimize import linprog
import numpy as np


def factorize(n):
    factors = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors


def is_prime(n):
    if n < 2: return False
    if n < 4: return True
    if n % 2 == 0 or n % 3 == 0: return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0: return False
        i += 6
    return True


def divisors(n):
    divs = []
    for i in range(1, int(math.isqrt(n)) + 1):
        if n % i == 0:
            divs.append(i)
            if i != n // i:
                divs.append(n // i)
    return sorted(divs)


def v_p(n, p):
    if n == 0: return float('inf')
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def euler_phi(n):
    result = n
    p = 2
    temp = n
    while p * p <= temp:
        if temp % p == 0:
            while temp % p == 0:
                temp //= p
            result -= result // p
        p += 1
    if temp > 1:
        result -= result // temp
    return result


def compute_nmin_sieve(limit):
    nmin = {}
    for n in range(1, limit + 1):
        m = euler_phi(n)
        if m not in nmin or n < nmin[m]:
            nmin[m] = n
    return nmin


def build_atoms(m):
    facts = factorize(m)
    B = facts.get(2, 0)
    S = {s for s in facts if s != 2}
    bs = {s: facts[s] for s in S}
    
    divs = divisors(m)
    P_set = set()
    for d in divs:
        p = d + 1
        if is_prime(p):
            P_set.add(p)
    
    atoms = []
    atom_data = []
    
    for p in sorted(P_set):
        if p == 2:
            continue
        if p in S:
            Ep = list(range(1, bs[p] + 2))
        else:
            Ep = [1]
        
        for e in Ep:
            As = {}
            for s in S:
                val = v_p(p - 1, s) + (e - 1 if s == p else 0)
                As[s] = val
            
            w = v_p(p - 1, 2)
            cost = e * math.log(p)
            c = cost - math.log(2) * w
            
            atoms.append((p, e))
            atom_data.append({
                'p': p, 'e': e,
                'As': As,
                'w': w,
                'cost': cost,
                'c': c
            })
    
    return atoms, S, bs, B, P_set, atom_data


def solve_lp(atoms, atom_data, S, bs, B, mode='even'):
    """
    Solve LP relaxation.
    mode='even': sum w*x <= B  (inequality)
    mode='odd':  sum w*x = B   (equality)
    """
    n_atoms = len(atoms)
    if n_atoms == 0:
        if len(S) == 0 and mode == 'even':
            return 0.0  # trivial: no odd constraints, no atoms needed
        return None
    
    c_vec = np.array([ad['c'] for ad in atom_data], dtype=float)
    S_list = sorted(S)
    
    # Equality constraints: valuation equalities
    A_eq_rows = []
    b_eq_vals = []
    for s in S_list:
        row = [float(ad['As'].get(s, 0)) for ad in atom_data]
        A_eq_rows.append(row)
        b_eq_vals.append(float(bs[s]))
    
    # Per-prime inequality constraints
    primes_in_atoms = sorted(set(a[0] for a in atoms))
    A_ub_rows = []
    b_ub_vals = []
    for p in primes_in_atoms:
        row = [1.0 if atoms[i][0] == p else 0.0 for i in range(n_atoms)]
        A_ub_rows.append(row)
        b_ub_vals.append(1.0)
    
    bounds = [(0.0, 1.0) for _ in range(n_atoms)]
    
    w_row = [float(ad['w']) for ad in atom_data]
    
    if mode == 'even':
        # Add sum w*x <= B as inequality
        A_ub_all = A_ub_rows + [w_row]
        b_ub_all = b_ub_vals + [float(B)]
        A_eq = np.array(A_eq_rows) if A_eq_rows else None
        b_eq = np.array(b_eq_vals) if b_eq_vals else None
        A_ub = np.array(A_ub_all)
        b_ub = np.array(b_ub_all)
    else:  # odd
        # Add sum w*x = B as equality
        A_eq_all = A_eq_rows + [w_row]
        b_eq_all = b_eq_vals + [float(B)]
        A_eq = np.array(A_eq_all)
        b_eq = np.array(b_eq_all)
        A_ub = np.array(A_ub_rows) if A_ub_rows else None
        b_ub = np.array(b_ub_vals) if b_ub_vals else None
    
    try:
        if A_ub is not None and len(A_ub) > 0:
            res = linprog(c_vec, A_ub=A_ub, b_ub=b_ub,
                          A_eq=A_eq, b_eq=b_eq,
                          bounds=bounds, method='highs')
        else:
            res = linprog(c_vec, A_eq=A_eq, b_eq=b_eq,
                          bounds=bounds, method='highs')
        if res.success:
            return res.fun
    except Exception as e:
        print(f"  LP solver error: {e}")
    return None


# Expected from Table 1
expected = {
    2:       {'Nmin': 3,       'LP': 1.500, 'gap': 0.0},
    4:       {'Nmin': 5,       'LP': 1.250, 'gap': 0.0},
    8:       {'Nmin': 15,      'LP': 1.875, 'gap': 0.0},
    12:      {'Nmin': 13,      'LP': 1.083, 'gap': 0.0},
    24:      {'Nmin': 35,      'LP': 1.211, 'gap': 5.2},
    48:      {'Nmin': 65,      'LP': 1.117, 'gap': 4.6},
    120:     {'Nmin': 143,     'LP': 1.047, 'gap': 2.6},
    480:     {'Nmin': 527,     'LP': 1.012, 'gap': 1.3},
    720:     {'Nmin': 779,     'LP': 1.009, 'gap': 1.1},
    5040:    {'Nmin': 5183,    'LP': 1.002, 'gap': 0.3},
    5888:    {'Nmin': 11985,   'LP': 1.025, 'gap': 7.3},
    10496:   {'Nmin': 21165,   'LP': 1.016, 'gap': 6.9},
    17408:   {'Nmin': 34935,   'LP': 1.011, 'gap': 6.6},
    32768:   {'Nmin': 65535,   'LP': 2.000, 'gap': 0.0},
    40080:   {'Nmin': 75165,   'LP': 1.013, 'gap': 5.5},
    40448:   {'Nmin': 80835,   'LP': 1.007, 'gap': 6.1},
    44288:   {'Nmin': 88485,   'LP': 1.006, 'gap': 6.0},
    55328:   {'Nmin': 103755,  'LP': 1.015, 'gap': 5.3},
    64256:   {'Nmin': 128265,  'LP': 1.005, 'gap': 5.8},
    65024:   {'Nmin': 129795,  'LP': 1.005, 'gap': 5.8},
    71168:   {'Nmin': 142035,  'LP': 1.005, 'gap': 5.8},
    72080:   {'Nmin': 135165,  'LP': 1.004, 'gap': 5.3},
    86528:   {'Nmin': 172635,  'LP': 1.001, 'gap': 5.7},
    88944:   {'Nmin': 166785,  'LP': 1.004, 'gap': 5.2},
    101888:  {'Nmin': 203235,  'LP': 1.005, 'gap': 5.6},
    107744:  {'Nmin': 202035,  'LP': 1.011, 'gap': 5.1},
    124928:  {'Nmin': 249135,  'LP': 1.004, 'gap': 5.5},
    140288:  {'Nmin': 279735,  'LP': 1.004, 'gap': 5.5},
    276736:  {'Nmin': 563295,  'LP': 1.025, 'gap': 5.2},
    2037248: {'Nmin': 4158795, 'LP': 1.006, 'gap': 4.6},
}

table_m_values = sorted(expected.keys())

print("=" * 130)
print("LP BOUND AND GAP VERIFICATION")
print("=" * 130)

# First, let's verify what the paper means by the factorization column
print("\nVerifying factorizations...")
factorization_check = {
    2: "2",
    4: "2^2",
    8: "2^3",
    12: "2^2 * 3",
    24: "2^3 * 3",
    48: "2^4 * 3",
    120: "2^3 * 3 * 5",
    480: "2^5 * 3 * 5",
    720: "2^4 * 3^2 * 5",
    5040: "2^4 * 3^2 * 5 * 7",
    5888: "2^8 * 23",
    10496: "2^8 * 41",
    17408: "2^10 * 17",
    32768: "2^15",
    40080: "2^4 * 3 * 5 * 167",
    40448: "2^9 * 79",
    44288: "2^8 * 173",
    55328: "2^5 * 7 * 13 * 19",
    64256: "2^8 * 251",
    65024: "2^9 * 127",
    71168: "2^9 * 139",
    72080: "2^4 * 5 * 17 * 53",
    86528: "2^9 * 169",  # Table says 13^2, let's check
    88944: "2^4 * 3 * 17 * 109",
    101888: "2^9 * 199",
    107744: "2^5 * 7 * 13 * 37",
    124928: "2^11 * 61",
    140288: "2^10 * 137",
    276736: "2^8 * 23 * 47",
    2037248: "2^9 * 23 * 173",
}

for m in table_m_values:
    facts = factorize(m)
    # Reconstruct
    reconstructed = 1
    for p, e in facts.items():
        reconstructed *= p ** e
    if reconstructed != m:
        print(f"  ERROR: factorization of {m} doesn't reconstruct!")

# Check m=86528 specifically since table says 2^9 * 13^2
m_check = 86528
facts_check = factorize(m_check)
print(f"\n  86528 factorization: {facts_check}")
print(f"  Table says: 2^9 * 13^2 = {2**9 * 13**2} (= {2**9 * 169})")
# 2^9 * 169 = 512 * 169 = 86528 ✓
# But 169 = 13^2, and the table writes "132" which in the PDF rendering means 13^2

print("\n" + "-" * 130)
print(f"{'m':>10} | {'B':>3} | {'even_opt':>12} | {'odd_opt':>12} | "
      f"{'log_lb_even':>12} | {'log_lb_odd':>12} | {'LP_ratio':>10} | "
      f"{'LP(exp)':>8} | {'LP_match':>8} | {'gap_comp':>8} | {'gap(exp)':>8} | {'gap_match':>9}")
print("-" * 130)

discrepancies = []

for m in table_m_values:
    exp = expected[m]
    Nmin = exp['Nmin']
    
    atoms, S, bs, B, P_set, atom_data = build_atoms(m)
    
    # Solve both LPs
    even_opt = solve_lp(atoms, atom_data, S, bs, B, mode='even')
    odd_opt = solve_lp(atoms, atom_data, S, bs, B, mode='odd')
    
    log2 = math.log(2)
    
    # Compute log lower bounds
    log_lb_even = (B + 1) * log2 + even_opt if even_opt is not None else None
    log_lb_odd = B * log2 + odd_opt if odd_opt is not None else None
    
    # Take minimum (best lower bound on log Nmin)
    candidates = []
    if log_lb_even is not None:
        candidates.append(log_lb_even)
    if log_lb_odd is not None:
        candidates.append(log_lb_odd)
    
    if candidates:
        log_lb = min(candidates)
        lp_nmin_lb = math.exp(log_lb)
        lp_ratio = lp_nmin_lb / m
    else:
        log_lb = None
        lp_ratio = None
    
    # Compute gap as defined in the paper:
    # gap = (log(Nmin/m) - log(LP)) / log(Nmin)
    if lp_ratio is not None and lp_ratio > 0:
        log_actual_ratio = math.log(Nmin / m)
        log_lp_ratio = math.log(lp_ratio)
        gap_computed = (log_actual_ratio - log_lp_ratio) / math.log(Nmin) * 100  # percentage
    else:
        gap_computed = None
    
    # Check match
    lp_match = ""
    if lp_ratio is not None:
        if abs(lp_ratio - exp['LP']) < 0.002:
            lp_match = "✓"
        else:
            lp_match = f"✗ ({lp_ratio:.4f})"
    else:
        lp_match = "N/A"
    
    gap_match = ""
    if gap_computed is not None:
        if abs(gap_computed - exp['gap']) < 0.3:
            gap_match = "✓"
        else:
            gap_match = f"✗"
    else:
        gap_match = "N/A"
    
    even_str = f"{even_opt:.6f}" if even_opt is not None else "N/A"
    odd_str = f"{odd_opt:.6f}" if odd_opt is not None else "N/A"
    lb_even_str = f"{log_lb_even:.6f}" if log_lb_even is not None else "N/A"
    lb_odd_str = f"{log_lb_odd:.6f}" if log_lb_odd is not None else "N/A"
    lp_ratio_str = f"{lp_ratio:.4f}" if lp_ratio is not None else "N/A"
    gap_str = f"{gap_computed:.1f}%" if gap_computed is not None else "N/A"
    
    print(f"{m:>10} | {B:>3} | {even_str:>12} | {odd_str:>12} | "
          f"{lb_even_str:>12} | {lb_odd_str:>12} | {lp_ratio_str:>10} | "
          f"{exp['LP']:>8.3f} | {lp_match:>8} | {gap_str:>8} | {exp['gap']:>7.1f}% | {gap_match:>9}")
    
    if lp_ratio is not None and abs(lp_ratio - exp['LP']) >= 0.002:
        discrepancies.append((m, 'LP', exp['LP'], lp_ratio))
    if gap_computed is not None and abs(gap_computed - exp['gap']) >= 0.3:
        discrepancies.append((m, 'gap', exp['gap'], gap_computed))


print("\n" + "=" * 130)
print("LP/GAP DISCREPANCIES")
print("=" * 130)
if discrepancies:
    for m, what, exp_val, comp_val in discrepancies:
        print(f"  m={m}: {what} expected {exp_val}, computed {comp_val:.4f}")
else:
    print("  None found!")


# Now let's also verify some specific cases in detail
print("\n" + "=" * 130)
print("DETAILED ATOM ANALYSIS FOR SELECTED CASES")
print("=" * 130)

for m in [2, 8, 12, 24, 32768, 5888, 86528]:
    atoms, S, bs, B, P_set, atom_data = build_atoms(m)
    print(f"\nm = {m}, B = {B}, S = {S}, bs = {bs}")
    print(f"  P(m) \\ {{2}} = {sorted(p for p in P_set if p != 2)}")
    print(f"  Atoms ({len(atoms)}):")
    for i, ((p, e), ad) in enumerate(zip(atoms, atom_data)):
        As_str = {s: ad['As'][s] for s in sorted(S) if ad['As'].get(s, 0) > 0}
        print(f"    ({p},{e}): w={ad['w']}, c={ad['c']:.4f}, As={As_str}")
    
    # Also verify: for Nmin, what atoms are selected?
    # Use sieve to find Nmin and factorize it
    nmin_val = expected[m]['Nmin']
    nmin_facts = factorize(nmin_val)
    print(f"  Nmin({m}) = {nmin_val} = {nmin_facts}")
    print(f"  φ({nmin_val}) = {euler_phi(nmin_val)} (should be {m})")
    
    # Check the atom selection for the minimizer
    odd_part_facts = {p: e for p, e in nmin_facts.items() if p != 2}
    print(f"  Odd prime powers in Nmin: {odd_part_facts}")
    
    # Verify W(x) for this selection
    W = sum(v_p(p - 1, 2) for p in odd_part_facts)
    is_odd = nmin_val % 2 == 1
    print(f"  W(x) = {W}, B = {B}, Nmin is {'odd' if is_odd else 'even'}")
    if is_odd:
        print(f"  W(x) == B? {W == B} (required for odd preimage)")
    else:
        print(f"  neven = 2^(B-W+1) * odd_part = 2^{B-W+1} * {nmin_val // (2**(B-W+1))}")
