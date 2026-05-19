#!/usr/bin/env python3
"""
PBWT - Positional Burrows-Wheeler Transform implementation for haplotype matching.
Based on Durbin (2014), "Efficient haplotype matching and storage using 
the positional Burrows-Wheeler transform".
"""

import numpy as np
from dataclasses import dataclass
from typing import List


@dataclass
class MatchSegment:
    """
    Define a maximal match segment: query haplotype matches parent haplotype j on [start, end)
    Data Structure : MatchSegment(j=0, start=1000, end=4000)
    """
    j: int       # parent haplotype index
    start: int   # start position in the query
    end: int     # end position in the query


class PBWT:
    """
    Standard PBWT implementation for efficient haplotype matching.
    """
    
    def __init__(self, panel: np.ndarray):
        """
        Initialize PBWT with a parental haplotype panel.
        
        panel : np.ndarray of shape (P, N)
            P: number of haplotypes, N: number of sites
        """
        self.P, self.N = panel.shape
        self.panel = panel.astype(np.int8)
        
        # Precompute prefix arrays (a) and divergence arrays (d)
        # a[k] is the permutation of indices [0..P-1] that sorts haplotypes by suffixes ending at k-1
        # d[k][i] is the start position of the longest match between a[k][i] and a[k][i-1]
        self.a = np.zeros((self.N + 1, self.P), dtype=np.int32)
        self.d = np.zeros((self.N + 1, self.P), dtype=np.int32)
        
        self._build()
        
    def _build(self):
        """
        Build prefix and divergence arrays in O(PN) time using vectorized numpy operations.
        """
        P, N = self.P, self.N
        self.a[0] = np.arange(P)
        self.d[0] = 0
        
        for k in range(N):
            # Values at site k in the current sorted order
            vals = self.panel[self.a[k], k]
            
            # Indices of haplotypes with 0 and 1 at site k
            mask0 = (vals == 0)
            mask1 = ~mask0
            
            # Update prefix array for site k+1
            self.a[k+1] = np.concatenate([self.a[k, mask0], self.a[k, mask1]])
            
            # Update divergence array for site k+1
            # Standard PBWT divergence update:
            # We track the 'running' max divergence for 0s and 1s
            # This part is harder to fully vectorize without some tricks, 
            # but we can use a small loop or optimized logic.
            # For MAGIC panels (P is small), a simple loop over P is acceptable,
            # but we'll use a more efficient approach.
            
            p, q = k + 1, k + 1
            curr_d = self.d[k]
            
            # We still need a loop to maintain p and q correctly according to Algorithm 2
            # However, we only do this once during initialization.
            u_idx = 0
            v_idx = np.sum(mask0)
            
            for i in range(P):
                dist = curr_d[i]
                if dist > p: p = dist
                if dist > q: q = dist
                
                if vals[i] == 0:
                    self.d[k+1, u_idx] = p
                    u_idx += 1
                    p = 0
                else:
                    self.d[k+1, v_idx] = q
                    v_idx += 1
                    q = 0

    def match_query(self, query: np.ndarray, min_len: int = 3) -> List[MatchSegment]:
        """
        Find all maximal matches between query and panel haplotypes.
        
        While a "true" PBWT search (Algorithm 5) can find the longest match in O(N), 
        identifying ALL maximal matches for ALL parents is O(PN + #matches).
        
        This implementation uses a site-by-site scan which is O(PN), ensuring 
        exact compatibility with the painting weights required by MagicRecombiner.
        For very large panels, this can be further optimized using the precomputed 
        prefix and divergence arrays to jump between matching ranges.
        """
        matches = []
        P, N = self.P, self.N
        
        # f, g define the range [f, g) in the prefix array a[k] that matches the query
        f, g = 0, P
        e = 0 # match length
        
        # To find ALL maximal matches, we need to track when matches start and end.
        # However, for MagicRecombiner's "painting" purpose, we want maximal matches 
        # for each parent j.
        
        # We'll use a slightly different approach to find all maximal matches:
        # For each parent j, we keep track of its current match with the query.
        last_match_start = np.full(P, -1, dtype=np.int32)
        
        for k in range(N):
            y = query[k]
            
            # This inner loop is still O(P) per site, but it's more efficient than O(P*N) naive
            # because we only check the column k.
            # Actually, standard PBWT matching is O(N) to find SOME matches, 
            # but finding ALL maximal matches for ALL parents is naturally O(PN) 
            # if many matches exist.
            
            # Optimization: only check parents that were matching or could match
            # For simplicity and correctness in this first pass, we iterate P.
            # We can use the PBWT properties to optimize this if needed.
            for j in range(P):
                if self.panel[j, k] == y:
                    if last_match_start[j] == -1:
                        last_match_start[j] = k
                else:
                    if last_match_start[j] != -1:
                        length = k - last_match_start[j]
                        if length >= min_len:
                            matches.append(MatchSegment(j=j, start=last_match_start[j], end=k))
                        last_match_start[j] = -1
        
        # Close remaining matches at the end of the chromosome
        for j in range(P):
            if last_match_start[j] != -1:
                length = N - last_match_start[j]
                if length >= min_len:
                    matches.append(MatchSegment(j=j, start=last_match_start[j], end=N))
                    
        return matches
