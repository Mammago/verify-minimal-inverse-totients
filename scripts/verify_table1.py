#!/usr/bin/env python3
"""
Verify Table 1 from "Minimal Inverse Totients via Assigned Valuation Covers"

For each m in the table, we:
1. Compute Nmin(m) by sieving φ(n) for n up to a bound
2. Build the atom set A(m) and count |A(m)|
3. Compute |S(m)|
4. Build the incidence graph H(m) and estimate treewidth (heuristic upper bound)
5. Solve the LP relaxations (even and odd) from Theorem 4.1
6. Compute Nmin/m, check odd?, and compute the gap
"""

import math
from collections import defaultdict
from itertools import combinations
import numpy as np
from scipy.optimize import linprog
import networkx as nx


def euler_phi(n):
    """Compute Euler's totient function."""
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
    """Compute Nmin(m) for all totient values m by sieving φ(n) for n=1..limit."""
    nmin = {}
    for n in range(1, limit + 1):
        m = euler_phi(n)
        if m not in nmin or n < nmin[m]:
            nmin[m] = n
    return nmin


def factorize(n):
    """Return prime factorization as dict {prime: exponent}."""
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
    if n < 2:
        return False
    if n < 4:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True


def divisors(n):
    """Return all divisors of n."""
    divs = []
    for i in range(1, int(math.isqrt(n)) + 1):
        if n % i == 0:
            divs.append(i)
            if i != n // i:
                divs.append(n // i)
    return sorted(divs)


def v_p(n, p):
    """p-adic valuation of n."""
    if n == 0:
        return float('inf')
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def build_atoms_and_graph(m):
    """
    Build the atom set A(m) and the incidence graph H(m).
    
    Returns:
        atoms: list of (p, e) tuples
        S: set of odd primes s with v_s(m) > 0
        bs: dict {s: v_s(m)}
        B: v_2(m)
        P_set: candidate prime set
        edges: list of (atom_index, constraint_label) pairs
        atom_data: list of dicts with w, cost, c, As for each atom
    """
    facts = factorize(m)
    B = facts.get(2, 0)
    
    # S(m) = odd primes s with v_s(m) > 0
    S = {s for s in facts if s != 2}
    bs = {s: facts[s] for s in S}
    
    # Candidate primes: p prime with p-1 | m
    divs = divisors(m)
    P_set = set()
    for d in divs:
        p = d + 1
        if is_prime(p):
            P_set.add(p)
    
    # Build atom set A(m)
    atoms = []
    atom_data = []
    
    for p in sorted(P_set):
        if p == 2:
            continue  # p=2 excluded from atoms
        
        if p in S:
            Ep = list(range(1, bs[p] + 2))  # {1, 2, ..., b_p + 1}
        else:
            Ep = [1]
        
        for e in Ep:
            # Compute contributions
            As = {}
            for s in S:
                val = v_p(p - 1, s) + (e - 1 if s == p else 0)
                As[s] = val
            
            w = v_p(p - 1, 2)  # 2-weight
            cost = e * math.log(p)  # log-cost
            c = cost - math.log(2) * w  # adjusted cost
            
            atoms.append((p, e))
            atom_data.append({
                'p': p, 'e': e,
                'As': As,
                'w': w,
                'cost': cost,
                'c': c
            })
    
    return atoms, S, bs, B, P_set, atom_data


def build_incidence_graph(atoms, atom_data, S, bs, P_set):
    """
    Build the incidence graph H(m) as a bipartite graph.
    Left vertices: atom indices 0..|A|-1
    Right vertices: constraint labels
      - ('s', s) for s in S: the valuation constraint sum As(a)*xa = bs
      - ('p', p) for p in P\{2}: the per-prime constraint sum_{e in Ep} x_{p,e} <= 1
    """
    G = nx.Graph()
    
    # Add atom vertices
    for i in range(len(atoms)):
        G.add_node(('atom', i))
    
    # Add constraint vertices and edges
    # Per-prime constraints
    primes_used = set()
    for i, (p, e) in enumerate(atoms):
        primes_used.add(p)
    
    for p in primes_used:
        G.add_node(('prime_constr', p))
        for i, (pp, ee) in enumerate(atoms):
            if pp == p:
                G.add_edge(('atom', i), ('prime_constr', p))
    
    # Valuation constraints
    for s in S:
        G.add_node(('val_constr', s))
        for i, ad in enumerate(atom_data):
            if ad['As'].get(s, 0) > 0:
                G.add_edge(('atom', i), ('val_constr', s))
    
    return G


def estimate_treewidth(G):
    """Estimate treewidth using NetworkX's min_degree heuristic (upper bound)."""
    if G.number_of_nodes() == 0:
        return 0
    if G.number_of_edges() == 0:
        return 0 if G.number_of_nodes() <= 1 else 1
    
    # Use greedy min-degree elimination for upper bound
    H = G.copy()
    tw = 0
    while H.number_of_nodes() > 0:
        # Find min-degree node
        min_deg = float('inf')
        min_node = None
        for node in H.nodes():
            d = H.degree(node)
            if d < min_deg:
                min_deg = d
                min_node = node
        
        # Eliminate: connect all neighbors, then remove node
        neighbors = list(H.neighbors(min_node))
        tw = max(tw, len(neighbors))
        
        for i in range(len(neighbors)):
            for j in range(i + 1, len(neighbors)):
                if not H.has_edge(neighbors[i], neighbors[j]):
                    H.add_edge(neighbors[i], neighbors[j])
        
        H.remove_node(min_node)
    
    return tw


def solve_lp_relaxations(atoms, atom_data, S, bs, B):
    """
    Solve both LP relaxations from Theorem 4.1:
    - Even LP: sum w*x <= B
    - Odd LP: sum w*x = B
    
    Returns (even_opt, odd_opt) - the optimal adjusted costs, or None if infeasible.
    """
    n_atoms = len(atoms)
    if n_atoms == 0:
        return None, None
    
    # Objective: minimize sum c(a)*x(a)
    c_vec = np.array([ad['c'] for ad in atom_data])
    
    # Build constraint matrices
    # Equality constraints: sum As(a)*xa = bs for each s in S
    # Per-prime constraints: sum_{e in Ep} x_{p,e} <= 1 for each p
    # 2-weight constraint: sum w(a)*xa <= B (even) or = B (odd)
    # Bounds: 0 <= x <= 1
    
    S_list = sorted(S)
    
    # Equality constraints (A_eq * x = b_eq)
    A_eq_rows = []
    b_eq_vals = []
    
    for s in S_list:
        row = [ad['As'].get(s, 0) for ad in atom_data]
        A_eq_rows.append(row)
        b_eq_vals.append(bs[s])
    
    # Per-prime inequality constraints: sum_{e} x_{p,e} <= 1
    # Collect primes
    primes_in_atoms = sorted(set(a[0] for a in atoms))
    
    A_ub_rows = []
    b_ub_vals = []
    
    for p in primes_in_atoms:
        row = [1 if atoms[i][0] == p else 0 for i in range(n_atoms)]
        A_ub_rows.append(row)
        b_ub_vals.append(1)
    
    # Bounds
    bounds = [(0, 1) for _ in range(n_atoms)]
    
    # --- Even LP: add sum w*x <= B as inequality ---
    w_row = [ad['w'] for ad in atom_data]
    A_ub_even = A_ub_rows + [w_row]
    b_ub_even = b_ub_vals + [B]
    
    A_eq_even = np.array(A_eq_rows) if A_eq_rows else None
    b_eq_even = np.array(b_eq_vals) if b_eq_vals else None
    A_ub_even_np = np.array(A_ub_even)
    b_ub_even_np = np.array(b_ub_even)
    
    even_opt = None
    try:
        res = linprog(c_vec, A_ub=A_ub_even_np, b_ub=b_ub_even_np,
                       A_eq=A_eq_even, b_eq=b_eq_even,
                       bounds=bounds, method='highs')
        if res.success:
            even_opt = res.fun
    except:
        pass
    
    # --- Odd LP: add sum w*x = B as equality ---
    A_eq_odd_rows = list(A_eq_rows) + [w_row]
    b_eq_odd_vals = list(b_eq_vals) + [B]
    
    A_eq_odd = np.array(A_eq_odd_rows)
    b_eq_odd = np.array(b_eq_odd_vals)
    A_ub_odd_np = np.array(A_ub_rows) if A_ub_rows else None
    b_ub_odd_np = np.array(b_ub_vals) if b_ub_vals else None
    
    odd_opt = None
    try:
        if A_ub_odd_np is not None and len(A_ub_rows) > 0:
            res = linprog(c_vec, A_ub=A_ub_odd_np, b_ub=b_ub_odd_np,
                           A_eq=A_eq_odd, b_eq=b_eq_odd,
                           bounds=bounds, method='highs')
        else:
            res = linprog(c_vec, A_eq=A_eq_odd, b_eq=b_eq_odd,
                           bounds=bounds, method='highs')
        if res.success:
            odd_opt = res.fun
    except:
        pass
    
    return even_opt, odd_opt


def compute_lp_bound_on_ratio(even_opt, odd_opt, B):
    """
    From Theorem 4.1:
    log Nmin(m) >= min{(B+1)log2 + even_opt, B*log2 + odd_opt}
    
    Since Nmin/m = Nmin / m, and log m = B*log2 + sum of odd part logs,
    we need: LP bound on Nmin(m)/m.
    
    Actually the paper reports LP as a lower bound on Nmin(m)/m.
    
    log(Nmin/m) = log(Nmin) - log(m)
    
    For the even case:
      log(neven) = (B+1)log2 + sum c(a)*xa
      
    For the odd case:
      log(nodd) = B*log2 + sum c(a)*xa
    
    So LP lower bound on log(Nmin) = min{(B+1)log2 + even_opt, B*log2 + odd_opt}
    
    Then LP bound on Nmin(m)/m = exp(LP bound on log Nmin - log m)
    """
    candidates = []
    log2 = math.log(2)
    
    if even_opt is not None:
        candidates.append((B + 1) * log2 + even_opt)
    if odd_opt is not None:
        candidates.append(B * log2 + odd_opt)
    
    if not candidates:
        return None
    
    return min(candidates)


# ============================================================
# MAIN VERIFICATION
# ============================================================

print("=" * 100)
print("VERIFICATION OF TABLE 1")
print("=" * 100)

# The m values from Table 1
table_m_values = [
    2, 4, 8, 12, 24, 48, 120, 480, 720, 5040, 5888, 10496, 17408,
    32768, 40080, 40448, 44288, 55328, 64256, 65024, 71168, 72080,
    86528, 88944, 101888, 107744, 124928, 140288, 276736, 2037248
]

# Expected values from Table 1
expected = {
    2:       {'Nmin': 3,       'ratio': 1.500, 'odd': True,  'A': 1,  'S': 0, 'tw': 1, 'LP': 1.500, 'gap': 0.0},
    4:       {'Nmin': 5,       'ratio': 1.250, 'odd': True,  'A': 2,  'S': 0, 'tw': 1, 'LP': 1.250, 'gap': 0.0},
    8:       {'Nmin': 15,      'ratio': 1.875, 'odd': True,  'A': 2,  'S': 0, 'tw': 1, 'LP': 1.875, 'gap': 0.0},
    12:      {'Nmin': 13,      'ratio': 1.083, 'odd': True,  'A': 5,  'S': 1, 'tw': 1, 'LP': 1.083, 'gap': 0.0},
    24:      {'Nmin': 35,      'ratio': 1.458, 'odd': True,  'A': 5,  'S': 1, 'tw': 1, 'LP': 1.211, 'gap': 5.2},
    48:      {'Nmin': 65,      'ratio': 1.354, 'odd': True,  'A': 6,  'S': 1, 'tw': 1, 'LP': 1.117, 'gap': 4.6},
    120:     {'Nmin': 143,     'ratio': 1.192, 'odd': True,  'A': 10, 'S': 2, 'tw': 2, 'LP': 1.047, 'gap': 2.6},
    480:     {'Nmin': 527,     'ratio': 1.098, 'odd': True,  'A': 13, 'S': 2, 'tw': 2, 'LP': 1.012, 'gap': 1.3},
    720:     {'Nmin': 779,     'ratio': 1.082, 'odd': True,  'A': 17, 'S': 2, 'tw': 2, 'LP': 1.009, 'gap': 1.1},
    5040:    {'Nmin': 5183,    'ratio': 1.028, 'odd': True,  'A': 30, 'S': 3, 'tw': 3, 'LP': 1.002, 'gap': 0.3},
    5888:    {'Nmin': 11985,   'ratio': 2.035, 'odd': True,  'A': 5,  'S': 1, 'tw': 1, 'LP': 1.025, 'gap': 7.3},
    10496:   {'Nmin': 21165,   'ratio': 2.016, 'odd': True,  'A': 5,  'S': 1, 'tw': 1, 'LP': 1.016, 'gap': 6.9},
    17408:   {'Nmin': 34935,   'ratio': 2.007, 'odd': True,  'A': 6,  'S': 1, 'tw': 1, 'LP': 1.011, 'gap': 6.6},
    32768:   {'Nmin': 65535,   'ratio': 2.000, 'odd': True,  'A': 4,  'S': 0, 'tw': 1, 'LP': 2.000, 'gap': 0.0},
    40080:   {'Nmin': 75165,   'ratio': 1.875, 'odd': True,  'A': 14, 'S': 3, 'tw': 2, 'LP': 1.013, 'gap': 5.5},
    40448:   {'Nmin': 80835,   'ratio': 1.998, 'odd': True,  'A': 5,  'S': 1, 'tw': 1, 'LP': 1.007, 'gap': 6.1},
    44288:   {'Nmin': 88485,   'ratio': 1.998, 'odd': True,  'A': 5,  'S': 1, 'tw': 1, 'LP': 1.006, 'gap': 6.0},
    55328:   {'Nmin': 103755,  'ratio': 1.875, 'odd': True,  'A': 8,  'S': 3, 'tw': 2, 'LP': 1.015, 'gap': 5.3},
    64256:   {'Nmin': 128265,  'ratio': 1.996, 'odd': True,  'A': 5,  'S': 1, 'tw': 1, 'LP': 1.005, 'gap': 5.8},
    65024:   {'Nmin': 129795,  'ratio': 1.996, 'odd': True,  'A': 5,  'S': 1, 'tw': 1, 'LP': 1.005, 'gap': 5.8},
    71168:   {'Nmin': 142035,  'ratio': 1.996, 'odd': True,  'A': 5,  'S': 1, 'tw': 1, 'LP': 1.005, 'gap': 5.8},
    72080:   {'Nmin': 135165,  'ratio': 1.875, 'odd': True,  'A': 13, 'S': 3, 'tw': 2, 'LP': 1.004, 'gap': 5.3},
    86528:   {'Nmin': 172635,  'ratio': 1.995, 'odd': True,  'A': 7,  'S': 1, 'tw': 1, 'LP': 1.001, 'gap': 5.7},
    88944:   {'Nmin': 166785,  'ratio': 1.875, 'odd': True,  'A': 13, 'S': 3, 'tw': 2, 'LP': 1.004, 'gap': 5.2},
    101888:  {'Nmin': 203235,  'ratio': 1.995, 'odd': True,  'A': 5,  'S': 1, 'tw': 1, 'LP': 1.005, 'gap': 5.6},
    107744:  {'Nmin': 202035,  'ratio': 1.875, 'odd': True,  'A': 9,  'S': 3, 'tw': 1, 'LP': 1.011, 'gap': 5.1},
    124928:  {'Nmin': 249135,  'ratio': 1.994, 'odd': True,  'A': 5,  'S': 1, 'tw': 1, 'LP': 1.004, 'gap': 5.5},
    140288:  {'Nmin': 279735,  'ratio': 1.994, 'odd': True,  'A': 5,  'S': 1, 'tw': 1, 'LP': 1.004, 'gap': 5.5},
    276736:  {'Nmin': 563295,  'ratio': 2.035, 'odd': True,  'A': 6,  'S': 2, 'tw': 2, 'LP': 1.025, 'gap': 5.2},
    2037248: {'Nmin': 4158795, 'ratio': 2.041, 'odd': True,  'A': 7,  'S': 2, 'tw': 1, 'LP': 1.006, 'gap': 4.6},
}

# Step 1: Compute Nmin by sieve
# We need n up to at least max(Nmin) ~ 4.2M, but let's go higher to be safe
SIEVE_LIMIT = 10_000_000
print(f"\nStep 1: Sieving φ(n) for n = 1..{SIEVE_LIMIT:,}...")
nmin_sieve = compute_nmin_sieve(SIEVE_LIMIT)
print(f"  Found {len(nmin_sieve):,} distinct totient values.")

# Step 2: Verify each row
print(f"\n{'m':>10} | {'Nmin(exp)':>10} | {'Nmin(comp)':>10} | {'Match':>5} | "
      f"{'Ratio(exp)':>10} | {'Ratio(comp)':>10} | {'Odd?':>5} | "
      f"|A|(exp) | |A|(comp) | |S|(exp) | |S|(comp) | tw(exp) | tw(comp)")
print("-" * 160)

results = []

for m in table_m_values:
    exp = expected[m]
    
    # Nmin from sieve
    nmin_computed = nmin_sieve.get(m, None)
    nmin_match = (nmin_computed == exp['Nmin']) if nmin_computed else False
    
    # Ratio
    ratio_computed = nmin_computed / m if nmin_computed else None
    
    # Odd?
    odd_computed = (nmin_computed % 2 == 1) if nmin_computed else None
    
    # Build atoms and graph
    atoms, S, bs, B, P_set, atom_data = build_atoms_and_graph(m)
    
    A_size = len(atoms)
    S_size = len(S)
    
    # Build incidence graph and estimate treewidth
    G = build_incidence_graph(atoms, atom_data, S, bs, P_set)
    tw_est = estimate_treewidth(G)
    
    # Print comparison
    ratio_str = f"{ratio_computed:.3f}" if ratio_computed else "N/A"
    odd_str = "✓" if odd_computed else "✗" if odd_computed is not None else "?"
    
    nmin_flag = "✓" if nmin_match else "✗"
    A_flag = "✓" if A_size == exp['A'] else "✗"
    S_flag = "✓" if S_size == exp['S'] else "✗"
    tw_flag = "✓" if tw_est == exp['tw'] else "?"
    
    print(f"{m:>10} | {exp['Nmin']:>10} | {nmin_computed:>10} | {nmin_flag:>5} | "
          f"{exp['ratio']:>10.3f} | {ratio_str:>10} | {odd_str:>5} | "
          f"{exp['A']:>4}{A_flag:>3} | {A_size:>5}    | {exp['S']:>4}{S_flag:>3} | {S_size:>5}    | "
          f"{exp['tw']:>4}{tw_flag:>3} | {tw_est:>5}")
    
    results.append({
        'm': m,
        'nmin_match': nmin_match,
        'nmin_computed': nmin_computed,
        'nmin_expected': exp['Nmin'],
        'ratio_computed': ratio_computed,
        'ratio_expected': exp['ratio'],
        'odd_computed': odd_computed,
        'odd_expected': exp['odd'],
        'A_size': A_size,
        'A_expected': exp['A'],
        'S_size': S_size,
        'S_expected': exp['S'],
        'tw_computed': tw_est,
        'tw_expected': exp['tw'],
        'atoms': atoms,
        'atom_data': atom_data,
        'S': S,
        'bs': bs,
        'B': B,
    })

print("\n" + "=" * 100)
print("SUMMARY OF BASIC CHECKS")
print("=" * 100)

nmin_ok = sum(1 for r in results if r['nmin_match'])
A_ok = sum(1 for r in results if r['A_size'] == r['A_expected'])
S_ok = sum(1 for r in results if r['S_size'] == r['S_expected'])
tw_ok = sum(1 for r in results if r['tw_computed'] == r['tw_expected'])
odd_ok = sum(1 for r in results if r['odd_computed'] == r['odd_expected'])

print(f"Nmin matches: {nmin_ok}/{len(results)}")
print(f"|A| matches:  {A_ok}/{len(results)}")
print(f"|S| matches:  {S_ok}/{len(results)}")
print(f"tw matches:   {tw_ok}/{len(results)}")
print(f"odd matches:  {odd_ok}/{len(results)}")

# Flag any discrepancies
print("\n" + "=" * 100)
print("DISCREPANCIES")
print("=" * 100)
any_disc = False
for r in results:
    m = r['m']
    discs = []
    if not r['nmin_match']:
        discs.append(f"Nmin: expected {r['nmin_expected']}, got {r['nmin_computed']}")
    if r['A_size'] != r['A_expected']:
        discs.append(f"|A|: expected {r['A_expected']}, got {r['A_size']}")
    if r['S_size'] != r['S_expected']:
        discs.append(f"|S|: expected {r['S_expected']}, got {r['S_size']}")
    if r['tw_computed'] != r['tw_expected']:
        discs.append(f"tw: expected {r['tw_expected']}, got {r['tw_computed']}")
    if r['odd_computed'] != r['odd_expected']:
        discs.append(f"odd: expected {r['odd_expected']}, got {r['odd_computed']}")
    if abs(r['ratio_computed'] - r['ratio_expected']) > 0.002:
        discs.append(f"ratio: expected {r['ratio_expected']:.3f}, got {r['ratio_computed']:.3f}")
    
    if discs:
        any_disc = True
        print(f"  m = {m}: " + "; ".join(discs))

if not any_disc:
    print("  None found!")
