#!/usr/bin/env python3
"""
Final cross-checks:
1. Verify φ(Nmin) = m for every row
2. Verify that Nmin is truly minimal (no smaller n with φ(n)=m)
3. Check paper's factorization column matches actual factorization
4. Check the ratio column to 3 decimal places
5. Verify that all Nmin are odd as claimed
6. Check the "132" entry for m=86528 (rendered as 13^2 in the PDF)
7. Spot-check: verify no n < Nmin has φ(n) = m for a few cases
"""

import math

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

def factor_str(n):
    facts = factorize(n)
    parts = []
    for p in sorted(facts):
        if facts[p] == 1:
            parts.append(str(p))
        else:
            parts.append(f"{p}^{facts[p]}")
    return " · ".join(parts) if parts else "1"

# Table 1 data
rows = [
    (2,       3),
    (4,       5),
    (8,       15),
    (12,      13),
    (24,      35),
    (48,      65),
    (120,     143),
    (480,     527),
    (720,     779),
    (5040,    5183),
    (5888,    11985),
    (10496,   21165),
    (17408,   34935),
    (32768,   65535),
    (40080,   75165),
    (40448,   80835),
    (44288,   88485),
    (55328,   103755),
    (64256,   128265),
    (65024,   129795),
    (71168,   142035),
    (72080,   135165),
    (86528,   172635),
    (88944,   166785),
    (101888,  203235),
    (107744,  202035),
    (124928,  249135),
    (140288,  279735),
    (276736,  563295),
    (2037248, 4158795),
]

# Paper's factorizations (from the table)
paper_factorizations = {
    2:       {2: 1},
    4:       {2: 2},
    8:       {2: 3},
    12:      {2: 2, 3: 1},
    24:      {2: 3, 3: 1},
    48:      {2: 4, 3: 1},
    120:     {2: 3, 3: 1, 5: 1},
    480:     {2: 5, 3: 1, 5: 1},
    720:     {2: 4, 3: 2, 5: 1},
    5040:    {2: 4, 3: 2, 5: 1, 7: 1},
    5888:    {2: 8, 23: 1},
    10496:   {2: 8, 41: 1},
    17408:   {2: 10, 17: 1},
    32768:   {2: 15},
    40080:   {2: 4, 3: 1, 5: 1, 167: 1},
    40448:   {2: 9, 79: 1},
    44288:   {2: 8, 173: 1},
    55328:   {2: 5, 7: 1, 13: 1, 19: 1},
    64256:   {2: 8, 251: 1},
    65024:   {2: 9, 127: 1},
    71168:   {2: 9, 139: 1},
    72080:   {2: 4, 5: 1, 17: 1, 53: 1},
    86528:   {2: 9, 13: 2},  # "132" in the PDF means 13^2
    88944:   {2: 4, 3: 1, 17: 1, 109: 1},
    101888:  {2: 9, 199: 1},
    107744:  {2: 5, 7: 1, 13: 1, 37: 1},
    124928:  {2: 11, 61: 1},
    140288:  {2: 10, 137: 1},
    276736:  {2: 8, 23: 1, 47: 1},
    2037248: {2: 9, 23: 1, 173: 1},
}

print("=" * 120)
print("COMPREHENSIVE CROSS-CHECKS")
print("=" * 120)

all_ok = True

print(f"\n{'m':>10} | {'Nmin':>10} | φ(Nmin)=m? | odd? | {'ratio':>10} | {'factor_check':>12} | {'m factorization':>30} | {'Nmin factorization':>30}")
print("-" * 150)

for m, nmin in rows:
    # Check φ(Nmin) = m
    phi_val = euler_phi(nmin)
    phi_ok = (phi_val == m)
    
    # Check odd
    is_odd = (nmin % 2 == 1)
    
    # Ratio
    ratio = nmin / m
    
    # Check factorization
    actual_facts = factorize(m)
    paper_facts = paper_factorizations[m]
    fact_ok = (actual_facts == paper_facts)
    
    # Reconstruct from paper factorization
    reconstructed = 1
    for p, e in paper_facts.items():
        reconstructed *= p ** e
    
    phi_str = "✓" if phi_ok else f"✗ (got {phi_val})"
    odd_str = "✓" if is_odd else "✗"
    fact_str = "✓" if fact_ok else f"✗"
    
    if not phi_ok or not is_odd or not fact_ok:
        all_ok = False
    
    print(f"{m:>10} | {nmin:>10} | {phi_str:>10} | {odd_str:>4} | {ratio:>10.4f} | {fact_str:>12} | {factor_str(m):>30} | {factor_str(nmin):>30}")

print("\n" + "=" * 120)
print("MINIMALITY SPOT-CHECKS")
print("=" * 120)
print("Checking that NO n < Nmin has φ(n) = m for selected cases...")

# For smaller m values, exhaustively check
spot_check_cases = [
    (2, 3), (4, 5), (8, 15), (12, 13), (24, 35), (48, 65),
    (120, 143), (480, 527), (720, 779), (5040, 5183),
    (5888, 11985), (32768, 65535)
]

for m, claimed_nmin in spot_check_cases:
    actual_nmin = None
    for n in range(1, claimed_nmin + 100):
        if euler_phi(n) == m:
            actual_nmin = n
            break
    
    if actual_nmin == claimed_nmin:
        print(f"  m={m:>6}: Nmin = {claimed_nmin:>6} ✓ (verified exhaustively up to {claimed_nmin + 100})")
    else:
        print(f"  m={m:>6}: MISMATCH! Claimed {claimed_nmin}, found {actual_nmin}")
        all_ok = False

# For larger m values, use the sieve approach
print("\nUsing sieve for larger values...")
from collections import defaultdict

SIEVE_LIMIT = 5_000_000
print(f"Sieving φ(n) for n=1..{SIEVE_LIMIT:,}...")
nmin_sieve = {}
for n in range(1, SIEVE_LIMIT + 1):
    m = euler_phi(n)
    if m not in nmin_sieve or n < nmin_sieve[m]:
        nmin_sieve[m] = n

large_cases = [
    (86528, 172635), (101888, 203235), (276736, 563295),
    (2037248, 4158795)
]

for m, claimed_nmin in large_cases:
    if claimed_nmin > SIEVE_LIMIT:
        # Need larger sieve for this one
        continue
    sieved = nmin_sieve.get(m)
    if sieved == claimed_nmin:
        print(f"  m={m:>10}: Nmin = {claimed_nmin:>10} ✓ (verified by sieve up to {SIEVE_LIMIT:,})")
    elif sieved is not None:
        print(f"  m={m:>10}: MISMATCH! Claimed {claimed_nmin}, sieve found {sieved}")
        all_ok = False
    else:
        print(f"  m={m:>10}: m not found as totient value in range (need larger sieve)")

# Special sieve for m=2037248
print(f"\nExtended sieve for m=2037248 (checking n up to 5,000,000)...")
target_m = 2037248
claimed = 4158795
found_min = None
for n in range(1, 5_000_001):
    if euler_phi(n) == target_m:
        found_min = n
        break

if found_min is not None:
    if found_min == claimed:
        print(f"  m=2037248: Nmin = {claimed} ✓")
    else:
        print(f"  m=2037248: MISMATCH! Claimed {claimed}, found {found_min}")
        all_ok = False
else:
    print(f"  m=2037248: No preimage found up to 5M, need to go higher")
    # Try directly checking φ(4158795)
    phi_check = euler_phi(4158795)
    print(f"  φ(4158795) = {phi_check} (should be {target_m})")
    if phi_check == target_m:
        print(f"  Confirmed φ(4158795) = {target_m}")
        # Check a few smaller candidates
        # 4158795 = 3 * 5 * 17 * 47 * 347
        facts_n = factorize(4158795)
        print(f"  4158795 = {facts_n}")

print("\n" + "=" * 120)
print("RATIO PRECISION CHECK (to 3 decimal places)")
print("=" * 120)

paper_ratios = {
    2: 1.500, 4: 1.250, 8: 1.875, 12: 1.083, 24: 1.458, 48: 1.354,
    120: 1.192, 480: 1.098, 720: 1.082, 5040: 1.028,
    5888: 2.035, 10496: 2.016, 17408: 2.007, 32768: 2.000,
    40080: 1.875, 40448: 1.998, 44288: 1.998, 55328: 1.875,
    64256: 1.996, 65024: 1.996, 71168: 1.996, 72080: 1.875,
    86528: 1.995, 88944: 1.875, 101888: 1.995, 107744: 1.875,
    124928: 1.994, 140288: 1.994, 276736: 2.035, 2037248: 2.041
}

for m, nmin in rows:
    actual_ratio = nmin / m
    paper_ratio = paper_ratios[m]
    # Check if rounding to 3 decimal places matches
    rounded = round(actual_ratio, 3)
    match = (rounded == paper_ratio)
    if not match:
        print(f"  m={m}: actual ratio = {actual_ratio:.6f}, rounded = {rounded:.3f}, "
              f"paper says {paper_ratio:.3f} — {'✓' if match else '✗ MISMATCH'}")
        all_ok = False

# Check if any ratios mismatch
ratio_mismatches = []
for m, nmin in rows:
    actual = round(nmin / m, 3)
    if actual != paper_ratios[m]:
        ratio_mismatches.append((m, actual, paper_ratios[m]))

if not ratio_mismatches:
    print("  All ratios match to 3 decimal places ✓")
else:
    for m, actual, paper in ratio_mismatches:
        print(f"  m={m}: computed {actual:.3f}, paper says {paper:.3f}")

print("\n" + "=" * 120)
if all_ok:
    print("ALL CHECKS PASSED ✓")
else:
    print("SOME CHECKS FAILED ✗")
print("=" * 120)
