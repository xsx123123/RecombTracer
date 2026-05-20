#!/usr/bin/env python3
"""
MAGIC HMM - Link PBWT chromosome painting with Hidden Markov Model
for refined local ancestry and recombination breakpoint inference.

Key idea:
---------
PBWT paint gives us P(parent=j | data at SNP k) as emission probabilities.
HMM adds transition probabilities based on genetic distance between SNPs,
smoothing out noise and giving principled recombination calls.

Model specification:
--------------------
Hidden states:    z_k ∈ {0, 1, ..., N_parent-1}  (which parent haplotype)
Emissions:        x_k ∈ {0, 1}                   (progeny allele)

Transition: P(z_k = j | z_{k-1} = i) = 
    if i == j:  exp(-r * d_k)           # no recombination
    if i != j:  (1 - exp(-r * d_k)) / (N-1)   # recombination to another parent

    where r = recombination rate per bp (~1e-8 per gen per bp)
          d_k = physical distance from SNP k-1 to k

Emission (from PBWT paint + error model):
    P(x_k | z_k = j) ∝ w_jk * (1-e) + (1-w_jk) * e
    
    where w_jk = normalized PBWT paint weight for parent j at SNP k
          e    = genotyping error rate (~0.001-0.01)

We use the forward-backward algorithm for posterior decoding,
and Viterbi for most likely path.
"""

import numpy as np
import pandas as pd
from dataclasses import dataclass
from typing import List, Tuple, Optional, Dict
from scipy.special import logsumexp
import warnings
from ..utils.log_utils import logger


@dataclass
class HMMParams:
    """HMM hyperparameters for MAGIC recombination inference."""
    rec_rate_per_bp: float = 1e-8      # Recombination rate per bp per generation
    genotyping_error: float = 0.005    # Genotyping error rate (0.5% typical for SNP chip)
    pbwt_weight_floor: float = 0.01    # Minimum weight to avoid log(0)
    effective_pop_size: float = 100.0  # For MAGIC, effective number of founders


class MagicHMM:
    """
    HMM smoother for MAGIC local ancestry, taking PBWT paint weights as input.
    
    Usage:
        1. Run MagicRecombiner.paint_progeny() to get per-SNP weights
        2. Feed weights into MagicHMM.decode() 
        3. Get Viterbi path + posterior probabilities
    """
    
    def __init__(self, 
                 parent_names: List[str],
                 positions: np.ndarray,
                 params: HMMParams = None):
        """
        Parameters
        ----------
        parent_names : list of str
            Unique parent names (hidden states)
        positions : np.ndarray
            SNP positions in bp
        params : HMMParams
            Model hyperparameters
        """
        self.parent_names = parent_names
        self.n_states = len(parent_names)
        self.positions = np.asarray(positions)
        self.N = len(positions)
        self.params = params or HMMParams()
        
        # Precompute genetic distances (cM) between adjacent SNPs
        self.distances = np.diff(self.positions).astype(np.float64)
        self.distances = np.concatenate([[self.distances[0]], self.distances])
        # Cap huge distances (e.g., centromere gaps) to avoid numerical issues
        self.distances = np.minimum(self.distances, 10_000_000)
    
    def _build_transition_matrix(self, dist: float) -> np.ndarray:
        """
        Build N x N transition matrix for a given physical distance.
        
        Uses the standard Haldane mapping function approximation:
        P(no recombination) = exp(-r * d)
        P(recombination)    = 1 - exp(-r * d)
        
        For MAGIC with N founders, if recombination happens,
        we can jump to any of the N founders (including self, but that's
        handled by the diagonal).
        """
        n = self.n_states
        r = self.params.rec_rate_per_bp
        
        # Probability of at least one recombination event in this interval
        # For MAGIC advanced generations, the effective recombination rate
        # is higher. We can approximate by scaling:
        # After g generations, effective rate ≈ g * r (for small intervals)
        # For MAGIC, g is typically 6-8. We'll allow user to tune via rec_rate.
        p_rec = 1.0 - np.exp(-r * dist)
        
        # Transition matrix
        if n == 1:
            T = np.ones((1, 1), dtype=np.float64)
        else:
            T = np.full((n, n), p_rec / max(n - 1, 1), dtype=np.float64)
            np.fill_diagonal(T, 1.0 - p_rec)
        
        # Ensure valid probability matrix
        T = np.clip(T, 1e-12, 1.0)
        T = T / T.sum(axis=1, keepdims=True)
        
        return T
    
    def _build_emission_probs(self, 
                              pbwt_weights: np.ndarray,
                              progeny_alleles: Optional[np.ndarray] = None,
                              parent_haps: Optional[np.ndarray] = None) -> np.ndarray:
        """
        Build emission probabilities from PBWT paint weights.
        
        Parameters
        ----------
        pbwt_weights : np.ndarray, shape (n_states, N)
            Normalized PBWT paint weights for each parent at each SNP
        progeny_alleles : np.ndarray, shape (N,)
            Progeny allele calls (optional, for error model)
        parent_haps : np.ndarray, shape (n_states, N)
            Parent haplotypes (optional, for explicit error model)
            
        Returns
        -------
        emission : np.ndarray, shape (N, n_states)
            P(observation | state=j) for each SNP and state
        """
        N, n = self.N, self.n_states
        e = self.params.genotyping_error
        floor = self.params.pbwt_weight_floor
        
        # weights: (n_states, N) -> transpose to (N, n_states)
        w = pbwt_weights.T.copy()
        w = np.clip(w, floor, 1.0)
        w = w / w.sum(axis=1, keepdims=True)
        
        if progeny_alleles is not None and parent_haps is not None:
            # Explicit error model: P(x_k | z_k=j) based on whether parent_j[k] == x_k
            # If match: (1-e) * w_jk + e * (1-w_jk)  ... but actually simpler:
            # We use PBWT weight as prior, and genotype match as likelihood
            emission = np.zeros((N, n), dtype=np.float64)
            for j in range(n):
                match = (parent_haps[j] == progeny_alleles).astype(float)
                # If parent allele matches progeny: high emission
                # If not: low emission (could be error or recombination)
                emission[:, j] = w[:, j] * (match * (1 - e) + (1 - match) * e)
        else:
            # Simplified model: PBWT weight directly modulates emission
            # Add small noise term so even low-weight parents have non-zero prob
            emission = w * (1 - e) + (1 - w) * e / (n - 1)
        
        # Normalize
        emission = np.clip(emission, 1e-12, 1.0)
        emission = emission / emission.sum(axis=1, keepdims=True)
        
        return emission
    
    def viterbi(self, 
                pbwt_weights: np.ndarray,
                progeny_alleles: Optional[np.ndarray] = None,
                parent_haps: Optional[np.ndarray] = None) -> Tuple[np.ndarray, float]:
        """
        Viterbi decoding: find most likely parent-of-origin path.
        
        Parameters
        ----------
        pbwt_weights : np.ndarray, shape (n_states, N)
        progeny_alleles : np.ndarray, shape (N,)
        parent_haps : np.ndarray, shape (n_states, N)
        
        Returns
        -------
        path : np.ndarray, shape (N,)
            Most likely state sequence
        log_prob : float
            Log probability of the path
        """
        N, n = self.N, self.n_states
        logger.debug(f"    - Starting Viterbi: N={N}, states={n}")
        
        emission = self._build_emission_probs(pbwt_weights, progeny_alleles, parent_haps)
        
        # Log-space Viterbi for numerical stability
        log_emission = np.log(emission)
        log_V = np.zeros((N, n), dtype=np.float64)
        backptr = np.zeros((N, n), dtype=np.int32)
        
        # Initialize
        log_V[0, :] = log_emission[0, :] + np.log(1.0 / n)
        
        # Recursion (vectorized over states)
        for k in range(1, N):
            T = self._build_transition_matrix(self.distances[k])
            log_T = np.log(T)
            # scores[i, j] = log_V[k-1, i] + log_T[i, j]
            scores = log_V[k-1, :][:, None] + log_T
            backptr[k, :] = np.argmax(scores, axis=0)
            log_V[k, :] = np.max(scores, axis=0) + log_emission[k, :]
        
        # Backtrack
        path = np.zeros(N, dtype=np.int32)
        path[-1] = np.argmax(log_V[-1, :])
        log_prob = log_V[-1, path[-1]]
        
        for k in range(N-2, -1, -1):
            path[k] = backptr[k+1, path[k+1]]
        
        logger.debug(f"    - Viterbi completed. Log-prob: {log_prob:.2f}")
        return path, log_prob
    
    def forward_backward(self,
                         pbwt_weights: np.ndarray,
                         progeny_alleles: Optional[np.ndarray] = None,
                         parent_haps: Optional[np.ndarray] = None) -> Tuple[np.ndarray, np.ndarray]:
        """
        Forward-backward algorithm: compute posterior state probabilities.
        
        Returns
        -------
        posterior : np.ndarray, shape (N, n_states)
            P(z_k = j | x_1:N) for each SNP k and state j
        log_likelihood : float
            Log likelihood of the data
        """
        N, n = self.N, self.n_states
        logger.debug(f"    - Starting Forward-Backward: N={N}, states={n}")
        
        emission = self._build_emission_probs(pbwt_weights, progeny_alleles, parent_haps)
        log_emission = np.log(emission)
        
        # Forward pass
        log_alpha = np.zeros((N, n), dtype=np.float64)
        log_alpha[0, :] = log_emission[0, :] + np.log(1.0 / n)
        
        for k in range(1, N):
            T = self._build_transition_matrix(self.distances[k])
            log_T = np.log(T)
            # log_alpha[k, j] = logsumexp_i [log_alpha[k-1, i] + log_T[i,j]] + log_emission[k,j]
            log_alpha[k, :] = logsumexp(log_alpha[k-1, :][:, None] + log_T, axis=0) + log_emission[k, :]
        
        log_likelihood = logsumexp(log_alpha[-1, :])
        
        # Backward pass
        log_beta = np.zeros((N, n), dtype=np.float64)
        log_beta[-1, :] = 0.0  # log(1)
        
        for k in range(N-2, -1, -1):
            T = self._build_transition_matrix(self.distances[k+1])
            log_T = np.log(T)
            # log_beta[k, i] = logsumexp_j [log_T[i,j] + log_emission[k+1,j] + log_beta[k+1,j]]
            log_beta[k, :] = logsumexp(log_T + log_emission[k+1, :][None, :] + log_beta[k+1, :][None, :], axis=1)
        
        # Posterior
        log_posterior = log_alpha + log_beta - log_likelihood
        posterior = np.exp(log_posterior)
        posterior = posterior / posterior.sum(axis=1, keepdims=True)
        
        logger.debug(f"    - Forward-Backward completed. Log-likelihood: {log_likelihood:.2f}")
        return posterior, log_likelihood
    
    def decode(self,
               pbwt_weights: np.ndarray,
               progeny_alleles: Optional[np.ndarray] = None,
               parent_haps: Optional[np.ndarray] = None,
               method: str = "viterbi") -> pd.DataFrame:
        """
        Main decoding interface.
        
        Parameters
        ----------
        pbwt_weights : np.ndarray, shape (n_states, N)
        method : str
            "viterbi" for hard path, "posterior" for soft probabilities
            
        Returns
        -------
        df : pd.DataFrame
            Columns: position, viterbi_parent (or best_posterior_parent), 
                     posterior_prob, all_parent_probs...
        """
        if method == "viterbi":
            path, log_prob = self.viterbi(pbwt_weights, progeny_alleles, parent_haps)
            posterior, _ = self.forward_backward(pbwt_weights, progeny_alleles, parent_haps)
            
            records = []
            for k in range(self.N):
                rec = {
                    'position': int(self.positions[k]),
                    'viterbi_parent': self.parent_names[path[k]],
                    'viterbi_posterior': float(posterior[k, path[k]]),
                }
                for j, name in enumerate(self.parent_names):
                    rec[f'p_{name}'] = float(posterior[k, j])
                records.append(rec)
            
            return pd.DataFrame(records)
        
        elif method == "posterior":
            posterior, log_likelihood = self.forward_backward(
                pbwt_weights, progeny_alleles, parent_haps
            )
            
            records = []
            for k in range(self.N):
                best_j = np.argmax(posterior[k])
                rec = {
                    'position': int(self.positions[k]),
                    'best_parent': self.parent_names[best_j],
                    'best_posterior': float(posterior[k, best_j]),
                }
                for j, name in enumerate(self.parent_names):
                    rec[f'p_{name}'] = float(posterior[k, j])
                records.append(rec)
            
            return pd.DataFrame(records)
        
        else:
            raise ValueError(f"Unknown method: {method}")
    
    def call_recombinations_hmm(self, 
                                viterbi_path: np.ndarray,
                                posterior: np.ndarray,
                                min_posterior: float = 0.8,
                                min_segment_snps: int = 3) -> pd.DataFrame:
        """
        Call recombination breakpoints from Viterbi path with confidence filtering.
        
        Parameters
        ----------
        viterbi_path : np.ndarray, shape (N,)
        posterior : np.ndarray, shape (N, n_states)
        min_posterior : float
            Minimum posterior probability at breakpoint to keep
        min_segment_snps : int
            Minimum SNPs per segment (filter spurious jumps)
            
        Returns
        -------
        df : pd.DataFrame
            Recombination events with confidence scores
        """
        N = len(viterbi_path)
        
        # Find breakpoints
        changes = np.where(viterbi_path[1:] != viterbi_path[:-1])[0]
        
        records = []
        prev_idx = 0
        
        for cp in changes:
            left_end = cp
            right_start = cp + 1
            
            # Segment lengths
            left_len = left_end - prev_idx + 1
            right_len = N - right_start  # approximate, will be refined next iter
            
            if left_len < min_segment_snps:
                continue
            
            # Confidence at breakpoint: average posterior of left and right states
            left_post = posterior[left_end, viterbi_path[left_end]]
            right_post = posterior[right_start, viterbi_path[right_start]]
            conf = (left_post + right_post) / 2
            
            if conf >= min_posterior:
                pos = (self.positions[left_end] + self.positions[right_start]) // 2
                records.append({
                    'position': int(pos),
                    'left_parent': self.parent_names[viterbi_path[left_end]],
                    'right_parent': self.parent_names[viterbi_path[right_start]],
                    'left_posterior': float(left_post),
                    'right_posterior': float(right_post),
                    'confidence': float(conf),
                    'left_snps': int(left_len),
                })
            
            prev_idx = right_start
        
        return pd.DataFrame(records)


# =============================================================================
# Integration with MagicRecombiner
# =============================================================================

def run_hmm_refinement(paint_df: pd.DataFrame,
                       parent_haps: np.ndarray,
                       parent_names: List[str],
                       progeny_hap: np.ndarray,
                       progeny_name: str = "",
                       haplotype: int = 0,
                       params: HMMParams = None,
                       chrom: str = "1") -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Run complete PBWT -> HMM pipeline for one haplotype.
    
    Parameters
    ----------
    paint_df : pd.DataFrame
        Output from MagicRecombiner.paint_progeny() for one haplotype
    parent_haps : np.ndarray, shape (n_parents, N)
    parent_names : list of str
    progeny_hap : np.ndarray, shape (N,)
    
    Returns
    -------
    viterbi_df : pd.DataFrame
        Per-SNP Viterbi decoding
    segments_df : pd.DataFrame
        Inferred ancestry segments
    rec_df : pd.DataFrame
        Recombination breakpoints
    """
    positions = paint_df['position'].values
    N = len(paint_df)
    
    logger.debug(f"  - Refinement for {progeny_name} hap {haplotype}: N={N} SNPs")
    
    # Extract PBWT weights from paint_df (need p_parent columns or reconstruct)
    # For now, reconstruct from parent assignments (simplified)
    # In practice, modify MagicRecombiner to output raw weights
    unique_parents = sorted(paint_df['parent'].unique())
    n_par = len(unique_parents)
    logger.debug(f"  - Building HMM weight matrix for {n_par} unique parents found in paint")
    
    # Build weight matrix from paint_df confidence
    # This is a simplification; ideally we pass raw weights from recombiner
    weights = np.zeros((n_par, N), dtype=np.float64)
    for k in range(N):
        par = paint_df.iloc[k]['parent']
        conf = paint_df.iloc[k]['confidence']
        j = unique_parents.index(par)
        weights[j, k] = conf
        # Spread small weight to others
        for other in range(n_par):
            if other != j:
                weights[other, k] = (1 - conf) / (n_par - 1)
    
    # Make sure parent_haps aligns with unique_parents
    parent_name_to_idx = {name: i for i, name in enumerate(parent_names)}
    aligned_haps = np.zeros((n_par, N), dtype=np.int8)
    for j, name in enumerate(unique_parents):
        if name in parent_name_to_idx:
            aligned_haps[j] = parent_haps[parent_name_to_idx[name]]
    
    # Run HMM
    hmm = MagicHMM(
        parent_names=unique_parents,
        positions=positions,
        params=params
    )
    
    # Viterbi + forward-backward (compute posterior only once)
    path, _ = hmm.viterbi(weights, progeny_hap, aligned_haps)
    posterior, _ = hmm.forward_backward(weights, progeny_hap, aligned_haps)
    
    # Build Viterbi dataframe
    records = []
    for k in range(N):
        rec = {
            'position': int(positions[k]),
            'viterbi_parent': unique_parents[path[k]],
            'viterbi_posterior': float(posterior[k, path[k]]),
        }
        for j, name in enumerate(unique_parents):
            rec[f'p_{name}'] = float(posterior[k, j])
        records.append(rec)
    viterbi_df = pd.DataFrame(records)
    
    rec_df = hmm.call_recombinations_hmm(path, posterior)
    logger.debug(f"  - HMM called {len(rec_df)} recombinations")
    
    # Build segments from Viterbi path
    segments = []
    if len(path) > 0:
        changes = np.where(path[1:] != path[:-1])[0]
        boundaries = [0] + (changes + 1).tolist() + [N]
        
        for i in range(len(boundaries) - 1):
            s, e = boundaries[i], boundaries[i+1]
            seg_post = posterior[s:e, path[s]].mean()
            segments.append({
                'chrom': chrom,
                'start': int(positions[s]),
                'end': int(positions[e-1]),
                'parent': unique_parents[path[s]],
                'n_snps': e - s,
                'mean_posterior': float(seg_post)
            })
    
    segments_df = pd.DataFrame(segments)
    logger.debug(f"  - HMM identified {len(segments_df)} ancestry segments")
    
    return viterbi_df, segments_df, rec_df


def demo_hmm():
    """Demonstrate PBWT + HMM integration with synthetic data."""
    np.random.seed(42)
    
    n_snps = 300
    n_parents = 4
    positions = np.cumsum(np.random.randint(10000, 50000, size=n_snps))
    
    # Create 4 independent founder haplotypes
    parent_haps = np.random.randint(0, 2, size=(n_parents, n_snps)).astype(np.int8)
    parent_names = [f"P{i+1}" for i in range(n_parents)]
    
    # Synthetic progeny with recombination at SNP 100 and 200
    prog_hap = np.zeros(n_snps, dtype=np.int8)
    prog_hap[:100] = parent_haps[0, :100]
    prog_hap[100:200] = parent_haps[1, 100:200]
    prog_hap[200:] = parent_haps[2, 200:]
    
    # Add genotyping errors (~2%)
    error_mask = np.random.random(n_snps) < 0.02
    prog_hap[error_mask] = 1 - prog_hap[error_mask]
    
    print("=" * 70)
    print("PBWT + HMM Integration Demo")
    print("=" * 70)
    print(f"Setup: {n_parents} parents, {n_snps} SNPs")
    print(f"True path: P1[0:100] -> P2[100:200] -> P3[200:300]")
    print(f"Genotyping error rate: 2% ({error_mask.sum()} errors)")
    print()
    
    from .recombiner import MagicRecombiner
    
    recombiner = MagicRecombiner(parent_haps, parent_names, positions, "1")
    paint_df = recombiner.paint_progeny(
        prog_hap.reshape(1, -1),
        progeny_name="demo",
        smooth_window=1
    )
    
    viterbi_df, seg_df, rec_df = run_hmm_refinement(
        paint_df[paint_df['haplotype'] == 0],
        parent_haps,
        parent_names,
        prog_hap,
        params=HMMParams(rec_rate_per_bp=2e-8, genotyping_error=0.02)
    )
    
    print("Viterbi path (transition regions):")
    show_idx = list(range(95, 105)) + list(range(195, 205))
    print(viterbi_df.iloc[show_idx][['position', 'viterbi_parent', 'viterbi_posterior']].to_string(index=False))
    print()
    
    print("Inferred segments:")
    print(seg_df.to_string(index=False))
    print()
    
    print("Inferred recombinations:")
    if len(rec_df) > 0:
        print(rec_df[['position', 'left_parent', 'right_parent', 'confidence']].to_string(index=False))
    else:
        print("  None")
    print()
    
    raw_parents = paint_df[paint_df['haplotype'] == 0]['parent'].values
    viterbi_parents = viterbi_df['viterbi_parent'].values
    true_parents = np.array(['P1']*100 + ['P2']*100 + ['P3']*100)
    
    raw_acc = (raw_parents == true_parents).mean()
    hmm_acc = (viterbi_parents == true_parents).mean()
    
    print(f"Accuracy comparison (per-SNP):")
    print(f"  Raw PBWT paint:  {raw_acc:.3f}")
    print(f"  HMM Viterbi:     {hmm_acc:.3f}")
    print(f"  Improvement:     {hmm_acc - raw_acc:+.3f}")
    
    if hmm_acc > raw_acc:
        print(f"\nHMM corrected {int((hmm_acc - raw_acc) * n_snps)} SNP calls!")
    print()
    print("Key insight: PBWT finds local matches; HMM enforces continuity.")


# Command-line entry point moved to recombtracer.cli
