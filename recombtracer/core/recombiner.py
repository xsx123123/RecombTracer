#!/usr/bin/env python3
"""
MAGIC Recombiner - Pure Python implementation of PBWT-style chromosome painting
for identifying parental recombination fragments in MAGIC progeny.
"""

import numpy as np
import pandas as pd
from dataclasses import dataclass
from typing import List, Tuple, Dict, Optional
import argparse
import warnings


from recombtracer.core.pbwt import PBWT, MatchSegment


@dataclass
class AncestrySegment:
    """
    An inferred local ancestry segment.
    
    This represents a contiguous genomic region assigned to a specific parent
    after PBWT painting and optional HMM smoothing.
    """
    chrom: str       # Chromosome identifier (e.g., "1", "Chr01")
    start_pos: int   # Physical start position in base pairs (bp)
    end_pos: int     # Physical end position in base pairs (bp)
    start_idx: int   # Start index in the SNP matrix/array (inclusive)
    end_idx: int     # End index in the SNP matrix/array (exclusive)
    parent: str      # Name/ID of the inferred parental origin
    haplotype: int   # Haplotype index (0 or 1 for diploid progeny)
    score: float     # Confidence score (mean matching weight or posterior probability)


class MagicRecombiner:
    """
    Identify recombination fragments in MAGIC progeny using PBWT-style painting.
    """
    
    def __init__(self, 
                 parent_haps: np.ndarray,
                 parent_names: List[str],
                 positions: np.ndarray,
                 chrom: str = "1"):
        """
        Initialize the MagicRecombiner with parental data.
        
        Args:
            parent_haps: Matrix of parental haplotypes (P x N).
            parent_names: List of names for each parental haplotype.
            positions: Array of physical positions (bp) for each SNP.
            chrom: Identifier for the chromosome.
        """
        # --- Basic Dimensions ---
        self.P = parent_haps.shape[0]                   # Total number of parental haplotypes (P)
        self.N = parent_haps.shape[1]                   # Total number of SNP markers (N)
        
        # --- Core Data Arrays ---
        self.parent_haps = parent_haps.astype(np.int8)  # Store as int8 to minimize memory footprint
        self.parent_names = np.array(parent_names)      # Array mapping each haplotype index to its parent name
        self.positions = positions                      # Map of SNP indices to physical bp coordinates
        self.chrom = chrom                              # Metadata: the chromosome being analyzed
        
        # --- Parent Mapping ---
        # Group haplotype indices by their unique parent names (handling multi-haplotype parents)
        self.parent_to_haps: Dict[str, List[int]] = {}
        for i, name in enumerate(parent_names):
            self.parent_to_haps.setdefault(name, []).append(i)

        # Pre-calculate unique parents for downstream classification
        self.unique_parents = list(self.parent_to_haps.keys())  # List of unique founder names
        self.n_parents = len(self.unique_parents)               # Count of unique founders
        
        # --- PBWT Initialization ---
        # Build the PBWT prefix and divergence arrays for the parental panel.
        # This is the "Search Engine" that enables fast haplotype matching.
        self.pbwt = PBWT(self.parent_haps)
        
    def _find_maximal_matches(self, query_hap: np.ndarray, min_match_len: int = 2) -> List[MatchSegment]:
        """Find all maximal matches between query_hap and all parent haplotypes using PBWT."""
        return self.pbwt.match_query(query_hap, min_len=min_match_len)
    
    def _paint_haplotype(self, query_hap: np.ndarray, 
                         min_match_len: int = 2,
                         smooth_window: int = 5) -> Tuple[np.ndarray, np.ndarray]:
        """Paint a single haplotype against all parents."""
        N = self.N
        n_par = self.n_parents
        
        weights = np.zeros((n_par, N), dtype=np.float64)
        matches = self._find_maximal_matches(query_hap, min_match_len=min_match_len)
        
        if not matches:
            warnings.warn("No matches found for a haplotype - check data!")
            return np.full(N, -1, dtype=np.int32), np.zeros(N)
        
        # Vectorized weight accumulation per match segment
        for m in matches:
            # We don't need the seg_len < min_match_len check here as PBWT already filtered it
            par_name = self.parent_names[m.j]
            par_idx = self.unique_parents.index(par_name)
            k = np.arange(m.start, m.end)
            # Triangular weight: peaks at center of match segment
            w = (k - m.start) * (m.end - k)
            weights[par_idx, m.start:m.end] += w
        
        total_weight = weights.sum(axis=0)
        valid = total_weight > 0
        
        best_parent_idx = np.zeros(N, dtype=np.int32)
        confidence = np.zeros(N, dtype=np.float64)
        
        if np.any(valid):
            best_parent_idx[valid] = np.argmax(weights[:, valid], axis=0)
            confidence[valid] = np.max(weights[:, valid], axis=0) / total_weight[valid]
        
        gaps = ~valid
        if np.any(gaps):
            valid_idx = np.where(valid)[0]
            if valid_idx.size > 0:
                best_parent_idx[gaps] = np.interp(
                    np.where(gaps)[0], 
                    valid_idx, 
                    best_parent_idx[valid]
                ).astype(np.int32)
                confidence[gaps] = np.interp(
                    np.where(gaps)[0],
                    valid_idx,
                    confidence[valid]
                )
        
        return best_parent_idx, confidence
    
    def paint_progeny(self, 
                      progeny_haps: np.ndarray,
                      progeny_name: str = "progeny",
                      min_match_len: int = 2,
                      smooth_window: int = 5) -> pd.DataFrame:
        """Paint all haplotypes of a progeny individual."""
        ploidy = progeny_haps.shape[0]
        results = []
        
        for h in range(ploidy):
            best_par, conf = self._paint_haplotype(
                progeny_haps[h].astype(np.int8),
                min_match_len=min_match_len
            )
            
            if smooth_window > 1:
                try:
                    from scipy.ndimage import median_filter
                    best_par = median_filter(best_par, size=smooth_window, mode='nearest')
                except ImportError:
                    pass
            
            for k in range(self.N):
                results.append({
                    'chrom': self.chrom,
                    'position': int(self.positions[k]),
                    'progeny': progeny_name,
                    'haplotype': h,
                    'parent': self.unique_parents[best_par[k]],
                    'confidence': float(conf[k])
                })
        
        return pd.DataFrame(results)
    
    def extract_segments(self, paint_df: pd.DataFrame, 
                         min_segment_snps: int = 5,
                         min_segment_bp: int = 1000) -> List[AncestrySegment]:
        """Extract contiguous ancestry segments from painted results."""
        segments = []
        
        for (prog, hap), group in paint_df.groupby(['progeny', 'haplotype']):
            group = group.sort_values('position').reset_index(drop=True)
            
            if len(group) == 0:
                continue
            
            parents = group['parent'].values
            changes = np.where(parents[1:] != parents[:-1])[0] + 1
            boundaries = [0] + changes.tolist() + [len(group)]
            
            for i in range(len(boundaries) - 1):
                s, e = boundaries[i], boundaries[i+1]
                seg_df = group.iloc[s:e]
                
                n_snps = len(seg_df)
                start_pos = int(seg_df['position'].iloc[0])
                end_pos = int(seg_df['position'].iloc[-1])
                bp_len = end_pos - start_pos
                
                if n_snps >= min_segment_snps and bp_len >= min_segment_bp:
                    segments.append(AncestrySegment(
                        chrom=str(seg_df['chrom'].iloc[0]),
                        start_pos=start_pos,
                        end_pos=end_pos,
                        start_idx=s,
                        end_idx=e,
                        parent=str(seg_df['parent'].iloc[0]),
                        haplotype=int(hap),
                        score=float(seg_df['confidence'].mean())
                    ))
        
        return segments
    
    def call_recombinations(self, segments: List[AncestrySegment]) -> pd.DataFrame:
        """Call recombination breakpoints from ancestry segments."""
        records = []
        
        from collections import defaultdict
        by_ind = defaultdict(list)
        for seg in segments:
            by_ind[(seg.chrom, seg.haplotype)].append(seg)
        
        for key, segs in by_ind.items():
            chrom, hap = key
            segs = sorted(segs, key=lambda s: s.start_pos)
            
            for i in range(len(segs) - 1):
                left, right = segs[i], segs[i+1]
                if left.parent != right.parent:
                    pos = (left.end_pos + right.start_pos) // 2
                    records.append({
                        'chrom': chrom,
                        'position': pos,
                        'haplotype': hap,
                        'left_parent': left.parent,
                        'right_parent': right.parent,
                        'left_end': left.end_pos,
                        'right_start': right.start_pos,
                        'confidence': (left.score + right.score) / 2
                    })
        
        return pd.DataFrame(records)


def demo():
    """Run a small synthetic MAGIC example."""
    np.random.seed(42)
    
    n_snps = 500
    n_parents = 4
    
    parent_haps = np.random.randint(0, 2, size=(n_parents, n_snps)).astype(np.int8)
    parent_names = [f"P{i+1}" for i in range(n_parents)]
    positions = np.arange(1000, 1000 + n_snps * 1000, 1000)
    
    prog_hap0 = np.zeros(n_snps, dtype=np.int8)
    prog_hap0[:150] = parent_haps[0, :150]
    prog_hap0[150:350] = parent_haps[1, 150:350]
    prog_hap0[350:] = parent_haps[2, 350:]
    
    prog_hap1 = parent_haps[3].copy()
    
    prog_haps = np.vstack([prog_hap0, prog_hap1])
    
    print("=" * 60)
    print("MAGIC Recombiner Demo")
    print("=" * 60)
    print(f"Parents: {parent_names}")
    print(f"SNPs: {n_snps}")
    print(f"Expected Hap0: P1[0:150] -> P2[150:350] -> P3[350:500]")
    print(f"Expected Hap1: P4[0:500]")
    print()
    
    recombiner = MagicRecombiner(
        parent_haps=parent_haps,
        parent_names=parent_names,
        positions=positions,
        chrom="1"
    )
    
    paint_df = recombiner.paint_progeny(
        prog_haps, 
        progeny_name="F1_demo",
        smooth_window=5
    )
    
    segments = recombiner.extract_segments(paint_df, min_segment_snps=5)
    rec_df = recombiner.call_recombinations(segments)
    
    print("Inferred segments:")
    for seg in segments:
        print(f"  Hap{seg.haplotype}: {seg.start_pos}-{seg.end_pos} -> {seg.parent} "
              f"(n_snps={seg.end_idx-seg.start_idx}, conf={seg.score:.3f})")
    
    print("\nInferred recombinations:")
    if len(rec_df) > 0:
        print(rec_df[['position', 'haplotype', 'left_parent', 'right_parent']].to_string(index=False))
    else:
        print("  None found")
    
    return paint_df, segments, rec_df


# Command-line entry point moved to recombtracer.cli
