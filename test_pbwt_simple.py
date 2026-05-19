#!/usr/bin/env python3
import numpy as np
from recombtracer.core.pbwt import PBWT
from rich import inspect

def test_pbwt_simple():
    # 1. Define a small parental panel (P=4, N=10)
    # Rows are haplotypes, columns are SNP sites
    panel = np.array([
        [0, 0, 1, 1, 0, 0, 1, 1, 0, 0], # P0
        [0, 1, 0, 1, 0, 1, 0, 1, 0, 1], # P1
        [1, 1, 0, 0, 1, 1, 0, 0, 1, 1], # P2
        [1, 0, 1, 0, 1, 0, 1, 0, 1, 0]  # P3
    ], dtype=np.int8)
    
    print("Parental Panel (4x10):")
    print(panel)
    print("-" * 30)

    # 2. Initialize PBWT
    pbwt = PBWT(panel)

    inspect(pbwt)
    
    # 3. Define a query haplotype (e.g., a mosaic of P0 and P2)
    # Part 1 (0-4): 0, 0, 1, 1 (from P0)
    # Part 2 (4-10): 1, 1, 0, 0, 1, 1 (from P2)
    query = np.array([0, 0, 1, 1, 1, 1, 0, 0, 1, 1], dtype=np.int8)
    
    print("Query Haplotype:")
    print(query)
    print("-" * 30)

    # 4. Find matches with min_len=3
    min_len = 3
    matches = pbwt.match_query(query, min_len=min_len)

    inspect(matches)

    print(f"Matches found (min_len={min_len}):")
    for m in matches:
        match_seq = panel[m.j, m.start:m.end]
        print(f"Parent {m.j} matches on [{m.start}, {m.end}): {match_seq}")

if __name__ == "__main__":
    test_pbwt_simple()
