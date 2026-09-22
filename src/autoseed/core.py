#!/usr/bin/env python3

# Generative python version of autoseed script that can extract motifs from sequences and generate new sequences containing the motifs (with option -distill)
# depends on: localmax-motif, seedextender, motifsimilarity and genint-PWM
# requires free tmp folder (file handling will overwrite files), and background PWM if used in generative mode (position-specific background for the generated sequences, width of PWM defines seq length)
# ./autoseed_v2.2.py [background.seq] [signal.seq] [min_kmer_length] [max_initial_kmer_length] [number_of_seeds] [max_kmer_extension_above_initial] [count_cutoff] [multinomial] [motif_similarity_cutoff]
# options: -distill                : generates new sequences with the motifs discovered, runs 10 iterations
#          -newonly motif_folder   : outputs motifs to specified folder, only adds motifs that are different from already mined ones
#          -debug                  : debug output
# example: ./autoseed_v2.2.py 80bpshuffled.seq STARRseq_promoter_80bp.seq 5 8 100 4 100 0.35 1 0.20 -distill

import sys
import os
import subprocess
import re
from dataclasses import dataclass
from typing import List, Tuple, Optional
import glob
import math
import numpy as np
import io
import contextlib

# Global settings and paths
debug = False
USE_IUPAC = True
_PKG_DIR = os.path.dirname(os.path.abspath(__file__))  # Installed package folder (holds bin/ and data/)
LOCALMAX_PATH = os.path.join(_PKG_DIR, "bin", "localmax-motif")           # Path to localmax-motif (local max kmer counter and motif generator)
EXTENDER_PATH = os.path.join(_PKG_DIR, "bin", "seedextender")             # Path to seedextender (program to extend seeds beyond what can be counted in memory)
MOTIFSIMILARITY_PATH = os.path.join(_PKG_DIR, "bin", "motifsimilarity")   # Path to motifsimilarity (program to measure similarity of motifs)
GENINT_PATH = os.path.join(_PKG_DIR, "bin", "genint-PWM")                 # Path to genint-PWM (generative part)
BACKGROUND_PWM = os.path.join(_PKG_DIR, "data", "80bak.pfm")               # Path to background PWM (position-specific background for the generated sequences, width of PWM defines seq length)

seed_history = []  # Global list for convergence checking
newonly_mode = False  # Global flag for -newonly mode
motif_folder = None   # Global folder path for existing motifs
genint_options = ""   # Global extra options passed through to genint-PWM (-genint)
keep_mode = False     # Global flag for -keep mode (genint-PWM -keepfeatures replaces the generate step)
keep_pvalue = "3.5"   # Global -log10 p-value for -keep mode (-keeppvalue)
keep_nofilter = False # Global flag: -keep mode without rescan filter (-nofilter, passed to genint-PWM)
keep_flanks = 0       # Global: bases protected on each side of kept features in -keep mode (-protect_flanks=N)
quiet = False         # Global flag for -quiet mode (no per-seed refinement output)

# TF characteristic k-mers and names (kmer, name)
tf_data = [
    ("AAAAAA", "NONE"), ("AATCAA", "PBX"), ("ACATGT", "P53"), ("ACCGCA", "RUNX"),
    ("ACCGGA", "ETSI"), ("AGATAA", "GATA"), ("AGGTCA", "NucRes1"), ("ATCGAT", "CUT"),
    ("CAACAG", "SCRT"), ("CACCTG", "TCF4"), ("CACGCA", "EGR"), ("CACGTC", "CREB3"),
    ("CACGTG", "EboxHLH"), ("CACTTA", "NKX"), ("CAGCTG", "NHLH"), ("CATAAA", "HOX9to12"),
    ("CATATG", "bHLHcat"), ("CATTCC", "TEAD"), ("CATATA", "SRF"), ("CCATTA", "PITX"),
    ("CCCGCC", "HighE2F"), ("CCCGGA", "ETSII"), ("CCGGAT", "SPDEF"), ("CCGTTA", "MYB"),
    ("CGAAAC", "IRF"), ("CGTAAA", "HOX13"), ("CTAGTG", "CTCF"), ("CTGTCA", "MEIS"),
    ("GAAACC", "IRF2"), ("GAACAA", "SOX"), ("GACCAC", "GLI"), ("GACGTC", "cre_bZIP"),
    ("GAGGAA", "ETSIII"), ("GCCACG", "SP_KLF"), ("GGCAAC", "RFX"), ("GGCGCC", "LowE2F"),
    ("GGGGAA", "NFKB1"), ("GGTACA", "AR"), ("GGTGTG", "Tbox"), ("GTCACG", "PAX"),
    ("TAAACA", "FOX"), ("TAATTA", "Homeo"), ("TACGTA", "TEF_GMEB"), ("TATGCA", "POU"),
    ("TGACAG", "MEIS"), ("TGCATA", "CUT2"), ("TGCCAA", "NFI"), ("TGCGCA", "CEBP"),
    ("TGCGGG", "GCM"), ("TGCTGA", "MAF"), ("TGGAAA", "NFAT"), ("TTCTAG", "HSF"),
    ("ATGCCC", "HIC"), ("GTCCGC", "HINFP"), ("GTGAAA", "PRD_IRF"), ("CGCCAT", "YY"),
    ("TATCGC", "ZBED"), ("CGACCA", "ZBTB7"), ("CGCTGT", "ZIC"), ("ACCCAC", "ZNF143"),
    ("ACCGGT", "GRHL"), ("CCATGG", "EBF"), ("ATCAAA", "TCF7_LEF"), ("AACGAT", "HOMEZ"),
    ("TTCGAA", "HSFY"), ("AAATAG", "MEF"), ("TCTAGA", "SMAD3"), ("TGCCCT", "TFAP"),
    ("CACGCC", "SREBF"), ("GATGCA", "CRE_CEBP"), ("TGACTC", "prebZIP1"), ("TGAGTC", "prebZIP2"),
    ("TGACAC", "oddFox_a"), ("TCCCCA", "FOXO"), ("TAAACG", "BARHLa"), ("TAATTG", "BARHLb"),
    ("CAATAA", "HOX9to12b"), ("ACATGA", "IRX"), ("CAAGGT", "ESRR"), ("GTCCAA", "HNF4"),
    ("AAGTCA", "NRE"), ("AGTTCA", "VDR"), ("GCATGC", "NRF1"), ("TCTCGC", "BANP")
]

# Extract separate lists for compatibility with existing code
tf_kmers = [kmer for kmer, name in tf_data]
tf_names = [name for kmer, name in tf_data]

# DNA encoding
dnaforward = "ACGTN"


@dataclass
class LocalMaxResult:
    kmer: str
    gap: int
    position: int
    background: int
    signal: int
    is_localmax: bool
    shift_ratio: float
    ic: float
    avg_var: float
    exp_var: float

@dataclass
class ExtendedSeedResult:
    seed: str
    bg_count: int
    sig_count: int
    parent_count: int
    parent_kmer: str
    extension: str
    shift_ratio: float = 0.0
    ic: float = 0.0
    avg_var: float = 0.0
    exp_var: float = 0.0
    max_type: str = ""
    confidence: str = ""
    is_localmax: bool = True

class OrientedMatch:
    def __init__(self):
        self.position = 0
        self.strand = 0
        self.score = 0.0
        self.id = 0
        
def read_pfm(filename: str) -> np.ndarray:
    """Read a PFM file with tab-separated counts for A, C, G, T in rows."""
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    counts = []
    for line in lines:
        if line.strip():
            row = [float(x) for x in line.strip().split('\t')]
            counts.append(row)
    
    return np.array(counts)

def normalize_pfm_to_pwm(counts: np.ndarray) -> np.ndarray:
    """Normalize PFM counts to frequencies."""
    column_sums = np.sum(counts, axis=0)
    column_sums = np.where(column_sums == 0, 1, column_sums)
    frequencies = counts / column_sums
    return frequencies

def svg_logo(pfm_file: str, output_file: str, title: str = "Sequence Logo"):
   """Generate SVG logo from PFM file using font rendering to match C code exactly."""
   
   # Read and normalize the PFM
   counts = read_pfm(pfm_file)
   n_fraction = normalize_pfm_to_pwm(counts)  # This is n[current_pwm].fraction in C
   
   # C code exact variables
   nucleotide_char = ['A', 'C', 'G', 'T']  # dnaforward in C
   colors = ['green', 'blue', 'orange', 'red', 'black', 'lightgray', 'white',
             'cornflowerblue', 'none', 'orchid', 'midnightblue', 'aliceblue',
             'saddlebrown', 'moccasin']  # exact colors array from C
   
   width = n_fraction.shape[1]  # n[current_pwm].width
   current_pwm = 0
   top_position = 20  # exact from C
   offset = 0
   paths = 0  # we're using paths=0 mode (font rendering)
   font = "Courier"
   
   # SVG dimensions from C code
   svg_width = width * 20  # n[0].width * 20
   svg_height = 100 + 130 * 0  # contacts = 0, so 100 + 130 * 0
   
   # Start building SVG - exact from C code
   svg_lines = []
   svg_lines.append('<?xml version="1.0" standalone="no"?>')
   svg_lines.append('<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN"')
   svg_lines.append('"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">')
   svg_lines.append(f'<svg width="{svg_width}" height="{svg_height}" x="0" y="0" version="1.1" xmlns="http://www.w3.org/2000/svg">')
   svg_lines.append(f'<title>{title}</title>')
   
   # No path definitions needed for font-based rendering
   
   # Group - exact from C code
   svg_lines.append(f'<g id="group{current_pwm}" transform="translate({0}, {current_pwm * 95 + offset})" >')
   
   # Main loop - exact translation from C
   for pwm_position in range(width):
       
       # DETERMINES ORDER BY BUBBLE SORT - exact from C
       positive_sum = 0
       negative_sum = 0
       
       # order array - exact from C
       order = [[0, 0, 0, 0], [0.0, 0.0, 0.0, 0.0]]
       
       for counter in range(4):
           order[0][counter] = counter
           order[1][counter] = n_fraction[counter, pwm_position]
           if order[1][counter] > 0:
               positive_sum = positive_sum + order[1][counter]
           else:
               negative_sum = negative_sum + order[1][counter]
       
       # Bubble sort - exact from C
       for counter in range(3):
           for nucleotide_value in range(counter, 4):
               # Handle NaN/Inf like C code
               if math.isnan(order[1][counter]) or math.isinf(order[1][counter]):
                   order[1][counter] = 0
               if order[1][counter] < order[1][nucleotide_value]:
                   # Swap
                   swap = order[0][counter]
                   order[0][counter] = order[0][nucleotide_value]
                   order[0][nucleotide_value] = swap
                   swap = order[1][counter]
                   order[1][counter] = order[1][nucleotide_value]
                   order[1][nucleotide_value] = swap
       
       # Print nucleotides using font rendering (C code paths=0 section exactly):
       font_position = 0
       for nucleotide_value in range(4):
           if order[1][nucleotide_value] > 0:
               nuc_idx = int(order[0][nucleotide_value])
               nuc_char = nucleotide_char[nuc_idx]
               value = order[1][nucleotide_value]
               color = colors[nuc_idx]
               
               # EXACT translation from C paths=0 section:
               # fprintf(outfile, "<text x=\"%i\" y=\"%f\" fill=\"%s\" font-size=\"30\" font-family=\"%s\"
               #         transform=\"scale(1, %f)\">%c</text>",
               #         pwm_position * 20, font_position / (order[1][nucleotide_value] * 4.5) + top_position - order[1][0] * 0.9,
               #         colors[...], font, order[1][nucleotide_value] * 4.5, nucleotide_char[...]);
               # font_position += (order[1][nucleotide_value] * 90);
               
               x_pos = pwm_position * 20
               y_pos = font_position / (value * 4.5) + top_position - order[1][0] * 0.9
               y_scale = value * 4.5
               
               svg_lines.append(f'<text x="{x_pos}" y="{y_pos}" fill="{color}" font-size="30" font-family="{font}" transform="scale(1, {y_scale})">{nuc_char}</text>')
               
               font_position += (value * 90)
   
   # Close group
   svg_lines.append('</g>')
   
   # Close SVG
   svg_lines.append('</svg>')
   
   # Write to file
   with open(output_file, 'w') as f:
       f.write('\n'.join(svg_lines))
   
   return '\n'.join(svg_lines)
   
def generate_logo_from_pfm(pfm_file: str, output_file: str = None, title: str = "Sequence Logo"):
    """Generate SVG logo from PFM file."""
    if output_file is None:
        output_file = pfm_file.replace('.pfm', '.svg')
        if output_file == pfm_file:
            output_file = pfm_file + '.svg'
    
    svg_logo(pfm_file, output_file, title)
    return f"SVG logo generated: {output_file}"
    
def parse_localmax_output(output: str) -> List[LocalMaxResult]:
    results = []
    for line in output.splitlines():  # Don't skip first line
        if not line or line.startswith("Total") or line.startswith("kmer"):  # Skip actual headers
            continue
        fields = line.split()
        if len(fields) < 10:
            continue
        try:
            results.append(LocalMaxResult(
                kmer=fields[0],
                gap=int(fields[1]),
                position=int(fields[2]),
                background=int(fields[3]),
                signal=int(fields[4]),
                is_localmax="Localmax" in fields[5],
                shift_ratio=float(fields[6]) if fields[6] != "NA" else 0.0,
                ic=float(fields[7]),
                avg_var=float(fields[8].strip('%'))/100,
                exp_var=float(fields[9].strip('%'))/100
            ))
        except (ValueError, IndexError):
            continue
    return results

def parse_extender_output(output: str) -> List[ExtendedSeedResult]:
    results = []
    grouped_results = {}
    for line in output.splitlines():
        if "Full_kmer" in line or not line:
            continue
        fields = line.split()
        if len(fields) < 9:
            continue
        try:
            # Handle 'NA' values for numeric fields
            shift_ratio = 0.0 if fields[6] == 'NA' else float(fields[6])
            ic = 0.0 if fields[7] == 'NA' else float(fields[7].strip('%'))/100
            avg_var = 0.0 if fields[8] == 'NA' else float(fields[8].strip('%'))/100
            # Fixed the logic for exp_var
            exp_var = 0.0 if (len(fields) <= 9 or fields[9] == 'NA') else float(fields[9].strip('%'))/100

            res = ExtendedSeedResult(
                seed=fields[0],
                bg_count=int(fields[1]),
                sig_count=int(fields[2]),
                parent_count=int(fields[3]),
                parent_kmer=fields[4],
                extension=fields[5],
                shift_ratio=shift_ratio,
                ic=ic,
                avg_var=avg_var,
                exp_var=exp_var,
                max_type=fields[10] if len(fields) > 10 else "",
                confidence=fields[11] if len(fields) > 11 else "",
                is_localmax=True
            )
            if res.parent_kmer not in grouped_results or res.sig_count > grouped_results[res.parent_kmer].sig_count:
                grouped_results[res.parent_kmer] = res
        except (ValueError, IndexError) as e:
            if debug:
                print(f"Failed to parse line: {line}")
                print(f"Error: {e}")
            continue
    return sorted(grouped_results.values(), key=lambda x: x.sig_count, reverse=True)

def parse_motif_output(output: str) -> str:
    if debug:
        print("\nParsing PFM from motif output:")
        print(output)
    IUPAC_MAP = {1: 'A', 2: 'C', 3: 'M', 4: 'G', 5: 'R', 6: 'S', 7: 'V', 8: 'T', 9: 'W', 10: 'Y', 11: 'H', 12: 'K', 13: 'D', 14: 'B', 15: 'N'}
    matrix = []
    skipped = 0
    started = False
    for line in output.splitlines():
        if "Motif from all matches" in line:
            started = True
            if debug: print("\nFound motif start")
            continue
        if "Match statistics:" in line:
            break
        if started and line and not line.startswith("Total"):
            nums = line.split()
            if len(nums) > 1:
                if skipped == 0 and nums[0].isdigit() and int(nums[0]) == 1:
                    skipped = 1
                    continue
                values = [float(x) for x in nums[1:]]
                matrix.append(values)
    if len(matrix) != 4:
        if debug: print(f"Error: Got {len(matrix)} rows instead of 4 ACGT rows")
        return ""
    if debug:
        print("\nFinal matrix:")
        for row in matrix:
            print(row)
    
    # Start with reasonable cutoffs
    primary_cutoff = 0.4    # For single base
    secondary_cutoff = 0.5  # For two base redundant
    tertiary_cutoff = 0.6   # For three base redundant
    
    iter = 0
    while primary_cutoff <= 0.95 and primary_cutoff > 0 and iter < 50:
        if debug: print(f"\nTrying cutoffs: primary={primary_cutoff}, secondary={secondary_cutoff}, tertiary={tertiary_cutoff}")
        consensus = []
        for col in range(len(matrix[0])):
            col_vals = [row[col] for row in matrix]
            total = sum(col_vals)
            if total <= 0:
                consensus.append('N')
                if debug: print("All zeros -> N")
                continue
            
            # Normalize to get frequencies
            freqs = [val/total for val in col_vals]
            max_freq = max(freqs)
            max_idx = freqs.index(max_freq)
            
            # Sort frequencies in descending order
            sorted_freqs = sorted(freqs, reverse=True)
            
            # Calculate IUPAC based on new rules
            if max_freq >= primary_cutoff:
                # Use single base (most abundant > 50%)
                nuc_value = 1 << max_idx
            elif sorted_freqs[1] >= secondary_cutoff * max_freq:
                # Second most abundant base is >50% of max base, use two-base code
                indices = [i for i, freq in enumerate(freqs) if freq >= secondary_cutoff * max_freq]
                nuc_value = sum(1 << i for i in indices)
            elif sorted_freqs[2] >= tertiary_cutoff * max_freq:
                # Third most abundant base is >50% of max base, use three-base code
                indices = [i for i, freq in enumerate(freqs) if freq >= tertiary_cutoff * max_freq]
                nuc_value = sum(1 << i for i in indices)
            else:
                # No clear pattern, use N
                nuc_value = 15
            
            base = IUPAC_MAP.get(nuc_value, 'N')
            consensus.append(base)
        
        consensus = ''.join(consensus)
        # Define the low information characters
        low_info = set('BDHVN')

        # Trim left end until max 2 low info characters remain
        while len(consensus) >= 3:
            if (consensus[0] in low_info and
                consensus[1] in low_info and
                consensus[2] in low_info):
                consensus = consensus[1:]  # Remove just one character at a time
            else:
                break

        # Trim right end until max 2 low info characters remain
        while len(consensus) >= 3:
            if (consensus[-3] in low_info and
                consensus[-2] in low_info and
                consensus[-1] in low_info):
                consensus = consensus[:-1]  # Remove just one character at a time
            else:
                break
        
        ic = sum(2 if c in 'ACGT' else 1 if c in 'RYMKWS' else 0.415 if c in 'BDHV' else 0 for c in consensus)
        if debug: print(f"\nConsensus: {consensus} (IC: {ic:.2f}, Length excluding N: {len(consensus)-consensus.upper().count('N')})")
        
        if 8 <= ic <= 20:
            if debug: print("IC in range 8-20")
            if ic/(len(consensus)-consensus.upper().count('N')+0.001) >= 1:
                if debug: print("IC per base also acceptable, returning consensus")
                return consensus
        
        iter += 1
        if ic > 20:
            # Too high information content, need MORE redundancy
            primary_cutoff += 0.025  # Make it harder to use single bases
            secondary_cutoff -= 0.025  # Make it easier to use two-base codes
            tertiary_cutoff -= 0.025  # Make it easier to use three-base codes
        else:
            # Too low information content, need LESS redundancy
            primary_cutoff -= 0.025  # Make it easier to use single bases
            secondary_cutoff += 0.025  # Make it harder to use two-base codes
            tertiary_cutoff += 0.025  # Make it harder to use three-base codes
    
    if debug: print("Failed to find consensus with IC 8-20")
    return ""
   
def run_localmax_kmer(bg_file: str, signal_file: str, min_len: int, max_len: int,
                      count_threshold: int, length_cutoff: float, top_n: int) -> List[LocalMaxResult]:
    cmd = [LOCALMAX_PATH, bg_file, signal_file, str(min_len), str(max_len), str(count_threshold), str(length_cutoff)]
    full_cmd = f"{' '.join(cmd)} | grep Localmax | sort -n -r -k5 | head -n {top_n}"
    if debug: print(f"\nRunning command: {full_cmd}")
    result = subprocess.run(full_cmd, shell=True, capture_output=True, text=True)
    with open("short_seeds_filtered.tmp", 'w') as f:
        f.write(result.stdout)
    return parse_localmax_output(result.stdout)

def run_seedextender(bg_file: str, signal_file: str, seeds_file: str, max_extension: int,
                    min_count: int, length_cutoff: float, robust: bool = True) -> List[ExtendedSeedResult]:
   base_cmd = [EXTENDER_PATH, bg_file, signal_file, seeds_file, str(max_extension), str(min_count), str(length_cutoff)]
   if robust:
       base_cmd.append("robust")
       
   cmd_str = f"{' '.join(base_cmd)} | grep Localmax"
   
   if debug:
       print(f"\nRunning seedextender command: {cmd_str}")
       print(f"Input seeds file: {seeds_file}")
       os.system(f"cat {seeds_file}")

   result = subprocess.run(cmd_str, shell=True, capture_output=True, text=True, check=True)
   
   with open("extended_seeds.tmp", 'w') as f:
       f.write(result.stdout)

   return parse_extender_output(result.stdout)

def run_localmax_motif(bg_file: str, signal_file: str, seed: str, multinomial: int,
                       lambda_val: Optional[float] = None) -> Tuple[str, str]:
    cmd = [LOCALMAX_PATH, bg_file, signal_file, seed, str(multinomial)]
    if lambda_val is not None:
        cmd.append(str(lambda_val))
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    with open("motif.tmp", 'w') as f:
        f.write(result.stdout)
    return result.stdout, parse_motif_output(result.stdout)

def seeds_converged(old_seed: str, new_seed: str) -> bool:
    global seed_history
    if new_seed == old_seed:
        return True
    if new_seed in seed_history:
        return True
    seed_history.append(old_seed)
    return False

def generate_sequence_value(searchstring):
    """Convert DNA sequence to numerical representation"""
    query_sequence_length = len(searchstring)
    query_sequence_value = 0
    position_value = 4 ** (query_sequence_length - 1)
    
    for position in range(query_sequence_length):
        nucleotide_value = 0
        for nuc_idx in range(4):
            if searchstring[position] == dnaforward[nuc_idx]:
                nucleotide_value = nuc_idx
                break
        else:
            raise ValueError(f"Invalid nucleotide in sequence: {searchstring[position]}")
        
        query_sequence_value += position_value * nucleotide_value
        position_value //= 4
    
    return query_sequence_value

def kmer_pos_and_score(pwm, pwm_width, kmer_sequence_value, kmer_length, o, o2):
    """Find best and second best score and position for kmer in PWM"""
    best_score = -100.0
    second_best_score = -100.0
    
    # Generate reverse complement
    reverse_complement_kmer_sequence_value = 0
    current_kmer_sequence_value = kmer_sequence_value
    
    for counter in range(kmer_length):
        reverse_complement_kmer_sequence_value <<= 2
        reverse_complement_kmer_sequence_value |= (3 ^ (current_kmer_sequence_value & 3))
        current_kmer_sequence_value >>= 2
    
    # Scan all positions
    for pwm_position in range(pwm_width - kmer_length + 1):
        current_position = pwm_position + kmer_length - 1  # Convert to rightmost position
        
        score = 0.0
        revscore = 0.0
        current_kmer_sequence_value = kmer_sequence_value
        current_reverse_complement_kmer_sequence_value = reverse_complement_kmer_sequence_value
        
        for kmer_position in range(kmer_length):
            nucleotide = current_kmer_sequence_value & 3
            rev_nucleotide = current_reverse_complement_kmer_sequence_value & 3
            
            score += pwm[nucleotide][current_position - kmer_position]
            revscore += pwm[rev_nucleotide][current_position - kmer_position]
            
            current_kmer_sequence_value >>= 2
            current_reverse_complement_kmer_sequence_value >>= 2
        
        current_strand = 1
        if revscore == score:
            current_strand = 0
        if revscore > score:
            score = revscore
            current_strand = -1
        
        if score >= best_score:
            second_best_score = best_score
            best_score = score
            o2.position = o.position
            o2.strand = o.strand
            o.position = pwm_position  # Store leftmost position (0-indexed)
            o.strand = current_strand
        elif score >= second_best_score:
            second_best_score = score
            o2.position = pwm_position  # Store leftmost position (0-indexed)
            o2.strand = current_strand
    
    o.score = best_score
    o2.score = second_best_score
    return 0

def best_characteristic_kmer_match_to_pwm(pwm_matrix):
    """Find best TF match for a PWM matrix"""
    pwm = pwm_matrix
    pwm_width = len(pwm[0])
    kmer_length = len(tf_kmers[0])
    
    # Initialize TF kmer sequence values
    number_of_tf_kmers = len(tf_kmers)
    tf_kmer_values = []
    
    for counter in range(number_of_tf_kmers):
        tf_kmer_values.append(generate_sequence_value(tf_kmers[counter]))
    
    # Create match results
    matches = []
    for current_kmer_number in range(number_of_tf_kmers):
        match1 = OrientedMatch()
        match2 = OrientedMatch()
        match1.id = current_kmer_number
        match2.id = current_kmer_number
        
        kmer_pos_and_score(pwm, pwm_width, tf_kmer_values[current_kmer_number],
                          kmer_length, match1, match2)
        
        matches.append(match1)
        matches.append(match2)
    
    # Sort by score (descending)
    matches.sort(key=lambda x: x.score, reverse=True)
    
    # Return best match
    best_match = matches[0]
    return {
        'tf_name': tf_names[best_match.id],
        'kmer': tf_kmers[best_match.id],
        'strand': 'top' if best_match.strand == 1 else ('bottom' if best_match.strand == -1 else 'both'),
        'position': best_match.position + 1,  # Convert to 1-indexed
        'score': best_match.score
    }

def predict_tf_from_pfm(pfm_file):
    """Load PFM file and predict TF"""
    try:
        counts = read_pfm(pfm_file)
        pwm = normalize_pfm_to_pwm(counts)
        return best_characteristic_kmer_match_to_pwm(pwm)
    except Exception as e:
        print(f"Error predicting TF from {pfm_file}: {e}")
        return None

def create_merged_svg(final_seeds, seed_counts, iteration):
    """Create merged SVG with TF predictions shown as arrows with labels above logos"""
    
    def get_text_width(text: str, font_size: int, font_family: str = "Courier") -> int:
        """Estimate text width in pixels. For monospace fonts only."""
        if font_family == "Courier":
            return len(text) * (font_size * 0.6)  # Courier character width is ~0.6 of font size
        raise ValueError("Width calculation only supported for Courier font")

    def svg_table(data: List[List[str]], column_widths: List[int], row_height: int,
                 start_x: int, start_y: int, font_size: int = 20,
                 font_family: str = "Courier", font_color: str = "black",
                 show_grid: bool = False, grid_color: str = "#CCCCCC",
                 h_padding: int = 15, alignments: List[str] = None) -> str:
        if not data or not data[0]:
            return ""
        
        n_cols = len(data[0])
        n_rows = len(data)
        if len(column_widths) != n_cols:
            raise ValueError("Number of column widths must match number of columns in data")
        
        if alignments is None:
            alignments = ["left"] + ["right"] * (n_cols - 1)
        
        x_positions = []
        curr_x = start_x
        for width in column_widths:
            curr_x += width
            x_positions.append(curr_x)
        
        svg = []
        
        if show_grid:
            curr_x = start_x
            for width in column_widths:
                curr_x += width
                svg.append(
                    f'<line x1="{curr_x}" y1="{start_y-row_height/2}" '
                    f'x2="{curr_x}" y2="{start_y + (n_rows-0.5)*row_height}" '
                    f'stroke="{grid_color}" stroke-width="1"/>'
                )
            
            for row_idx in range(n_rows+1):
                y = start_y + (row_idx * row_height) - row_height/2
                svg.append(
                    f'<line x1="{start_x}" y1="{y}" '
                    f'x2="{x_positions[-1]}" y2="{y}" '
                    f'stroke="{grid_color}" stroke-width="1"/>'
                )
        
        for row_idx, row in enumerate(data):
            y = start_y + (row_idx * row_height)
            for col_idx, cell in enumerate(row):
                if alignments[col_idx] == "left":
                    x = x_positions[col_idx] - column_widths[col_idx] + h_padding
                    anchor = "start"
                elif alignments[col_idx] == "center":
                    x = x_positions[col_idx] - column_widths[col_idx]/2
                    anchor = "middle"
                else:  # right
                    x = x_positions[col_idx] - h_padding
                    anchor = "end"
                
                if row_idx == 0:
                    x = x_positions[col_idx] - column_widths[col_idx]/2
                    anchor = "middle"
                
                svg.append(
                    f'<text x="{x}" y="{y}" '
                    f'font-size="{font_size}" font-family="{font_family}" '
                    f'text-anchor="{anchor}" fill="{font_color}">{cell}</text>'
                )
        
        return "\n".join(svg)

# Font settings
    title_font_size = 24
    table_font_size = 20
    filename_font_size = 16
    tf_label_font_size = 16
    font_family = "Courier"
    title_color = "black"
    table_color = "black"
    filename_color = "gray"
    tf_color = "black"

    # Vertical spacing and dimensions
    initial_offset = 50     # More space at top for TF labels
    logo_height = 120       # Height for each motif logo
    extra_spacing = 30      # Additional space between motif groups
    group_height = logo_height + extra_spacing
    merged_width = 3000     # Total SVG width

    # Horizontal positioning
    left_margin = 10        # Left edge margin
    table_position = 550    # Where table text begins
    
    # Vertical positioning within each group
    title_offset = 20       # Distance from group top to title (moved down to make room for TF)
    header_offset = 40      # Distance from group top to table header (moved down)
    filename_offset = logo_height + 5  # Distance from group top to filename

    # Initialize SVG
    merged_svg_filename = f"merged_logos_{iteration}.svg"
    current_y = initial_offset

    merged_svg = (
        '<?xml version="1.0" standalone="no"?>\n'
        '<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" '
        '"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">\n'
        f'<svg version="1.1" xmlns="http://www.w3.org/2000/svg" '
        f'width="{merged_width}" height="{len(final_seeds) * group_height + initial_offset}">\n'
    )

    # Process each motif and create its group
    for i, (orig_seed, refined_seed, pfm_filename) in enumerate(final_seeds, start=1):
        svg_filename = pfm_filename.replace(".pfm", ".svg")
        
        try:
            with open(svg_filename, "r") as f:
                content = f.read()
        except Exception as e:
            print(f"Error reading {svg_filename}: {e}")
            continue

        # Clean SVG tags
        content = re.sub(r'<\?xml.*?\?>', '', content, flags=re.DOTALL)
        content = re.sub(r'<!DOCTYPE.*?>', '', content, flags=re.DOTALL)
        content = re.sub(r'<svg[^>]*>', '', content, flags=re.DOTALL)
        content = re.sub(r'</svg>', '', content, flags=re.DOTALL)
        content = content.strip()

        # Get TF prediction
        tf_prediction = predict_tf_from_pfm(pfm_filename)
        if tf_prediction:
            tf_name = tf_prediction['tf_name']
            tf_strand = tf_prediction['strand']
            tf_position = tf_prediction['position']
            tf_score = tf_prediction['score']
        else:
            tf_name = "Unknown"
            tf_strand = "both"
            tf_position = 1
            tf_score = 0.0

        # Get statistics
        orig_sig = seed_counts.get(orig_seed, {}).get('signal', 0)
        orig_bg = seed_counts.get(orig_seed, {}).get('background', 0)
        orig_sig_pct = seed_counts.get(orig_seed, {}).get('signal_percent', 0)
        orig_bg_pct = seed_counts.get(orig_seed, {}).get('background_percent', 0)
        
        refined_sig = seed_counts.get(refined_seed, {}).get('signal', 0)
        refined_bg = seed_counts.get(refined_seed, {}).get('background', 0)
        refined_sig_pct = seed_counts.get(refined_seed, {}).get('signal_percent', 0)
        refined_bg_pct = seed_counts.get(refined_seed, {}).get('background_percent', 0)

        # Calculate column widths based on content
        seed_width = max(
            get_text_width("Seed", table_font_size, font_family),
            get_text_width(f"Original: {orig_seed}", table_font_size, font_family),
            get_text_width(f"Refined:  {refined_seed}", table_font_size, font_family)
        )

        bg_width = max(
            get_text_width("Background", table_font_size, font_family),
            get_text_width(f"{orig_bg} ({orig_bg_pct:.2f}%)", table_font_size, font_family),
            get_text_width(f"{refined_bg} ({refined_bg_pct:.2f}%)", table_font_size, font_family)
        )

        sig_width = max(
            get_text_width("Signal", table_font_size, font_family),
            get_text_width(f"{orig_sig} ({orig_sig_pct:.2f}%)", table_font_size, font_family),
            get_text_width(f"{refined_sig} ({refined_sig_pct:.2f}%)", table_font_size, font_family)
        )

        extra_width = 20
        column_widths = [seed_width + extra_width, bg_width + extra_width, sig_width + extra_width]

        # Format data for table (no title text above logos anymore)
        data = [
            ["Seed", "Background", "Signal"],
            [f"Original: {orig_seed}",
             f"{orig_bg} ({orig_bg_pct:.2f}%)",
             f"{orig_sig} ({orig_sig_pct:.2f}%)"],
            [f"Refined : {refined_seed}",
             f"{refined_bg} ({refined_bg_pct:.2f}%)",
             f"{refined_sig} ({refined_sig_pct:.2f}%)"]
        ]

        # Add group to SVG
        merged_svg += f'  <g transform="translate({left_margin}, {current_y})">\n'
        
        # Add TF prediction arrow and label REPLACING the title
        # Calculate arrow position: k-mer is always 6 bases wide, spans positions tf_position to tf_position+5
        arrow_left = (tf_position - 1) * 20  # Start of k-mer (0-indexed in pixels)
        arrow_right = arrow_left + 6 * 20  # End of k-mer (6 bases * 20 pixels each)
        arrow_y = title_offset-2  # Where the title used to be
        
        # Draw arrow line spanning exactly 6 bases
        merged_svg += f'    <polyline points="{arrow_left},{arrow_y} {arrow_right},{arrow_y}" fill="none" stroke="{tf_color}" stroke-width="2"/>\n'
        
        # Draw arrowhead based on strand
        if tf_strand == "top":
            # Forward strand arrowhead (pointing right)
            merged_svg += f'    <polyline points="{arrow_right-4},{arrow_y-3} {arrow_right},{arrow_y} {arrow_right-4},{arrow_y+3}" fill="none" stroke="{tf_color}" stroke-width="2"/>\n'
        elif tf_strand == "bottom":
            # Reverse strand arrowhead (pointing left)
            merged_svg += f'    <polyline points="{arrow_left+4},{arrow_y-3} {arrow_left},{arrow_y} {arrow_left+4},{arrow_y+3}" fill="none" stroke="{tf_color}" stroke-width="2"/>\n'
        # For "both" strand, no arrowhead
        
        # Add TF name label centered above the arrow
        label_x = (arrow_left + arrow_right) // 2 - len(tf_name) * 5  # Center the text over the arrow
        label_y = arrow_y - 8  # Above the arrow
        merged_svg += f'    <text x="{label_x}" y="{label_y}" fill="{tf_color}" font-size="{tf_label_font_size}" font-family="{font_family}" font-style="normal">{tf_name}</text>\n'
        
        # Add the logo
        merged_svg += f'    <g transform="translate(0, 20)">\n{content}\n    </g>\n'

        # Add the table
        merged_svg += svg_table(
            data=data,
            column_widths=column_widths,
            row_height=25,
            start_x=table_position,
            start_y=header_offset,
            font_size=table_font_size,
            font_family=font_family,
            font_color=table_color,
            show_grid=False,
            h_padding=10,
            alignments=["left", "right", "right"]
        ) + "\n"

        # Add filename only (no TF details)
        merged_svg += (
            f'    <text x="{left_margin}" y="{filename_offset}" '
            f'font-size="{filename_font_size}" font-family="{font_family}" '
            f'fill="{filename_color}">{pfm_filename}</text>\n'
        )

        merged_svg += "  </g>\n"
        current_y += group_height

    merged_svg += "</svg>\n"

    # Save the SVG
    with open(merged_svg_filename, "w") as out:
        out.write(merged_svg)

    if debug:
        print("Merged SVG generated:", merged_svg_filename)
        print(f"Total height: {current_y}")
        
def check_motif_similarity(pfm1: str, pfm2: str, spacing: int = 6) -> float:
   """Check similarity between two PFM files using motifsimilarity.
   Returns similarity score or 0.0 if error."""
   try:
       cmd = f"{MOTIFSIMILARITY_PATH} {pfm1} {pfm2} gapped {spacing}"
       if debug:
           print(f"Running similarity command: {cmd}")
       result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
       if debug:
           print("Motifsimilarity output:", result.stdout)
           if result.stderr:
               print("Motifsimilarity stderr:", result.stderr)
       for line in result.stdout.splitlines():
           if "gapped" in line:
               # Extract similarity score
               score = line.split('\t')[1]  # Use tab as separator
               if debug:
                   print(f"Found similarity line: {line}")
                   print(f"Extracted score: {score}")
               return 0.0 if 'e' in score.lower() else float(score)
   except Exception as e:
       if debug:
           print(f"Error checking motif similarity: {e}")
   return 0.0
    
    

def load_existing_motifs(folder_path: str) -> List[str]:
    """Load all existing .pfm files from the specified folder"""
    if not os.path.exists(folder_path):
        if debug:
            print(f"Motif folder does not exist: {folder_path}")
        return []
    
    pfm_files = glob.glob(os.path.join(folder_path, "*.pfm"))
    if debug:
        print(f"Found {len(pfm_files)} existing motif files in {folder_path}:")
        for pfm in pfm_files:
            print(f"  - {pfm}")
    
    return pfm_files

def write_refined_models_list(final_seeds: List[Tuple[str, str, str]], output_file: str = "refined_models.txt"):
    """Write the list of refined models to a flat file in genint-PWM format"""
    global newonly_mode, motif_folder
    
    try:
        # Determine output file path and mode
        if newonly_mode and motif_folder:
            # In newonly mode, append to the file in the motif folder
            output_path = os.path.join(motif_folder, "refined_models.txt")
            mode = 'a'  # Append mode
            print(f"Appending {len(final_seeds)} new models to: {output_path}")
        else:
            # Normal mode, create/overwrite in current directory
            output_path = output_file
            mode = 'w'  # Write mode (overwrite)
            print(f"Writing {len(final_seeds)} models to: {output_path}")
        
        with open(output_path, mode) as f:
            for original_seed, refined_seed, pfm_file in final_seeds:
                # Extract just the filename from the full path
                pfm_filename = os.path.basename(pfm_file)
                f.write(f"{pfm_filename}\t-1\n")
        
        if debug:
            print(f"  - {len(final_seeds)} models listed in genint-PWM format")
    except Exception as e:
        print(f"Error writing refined models list: {e}")

def parse_arguments(args: List[str]) -> Tuple[List[str], bool]:
    """Parse command line arguments and return filtered args and distill mode flag"""
    global debug, newonly_mode, motif_folder, genint_options, keep_mode, keep_pvalue, keep_nofilter, keep_flanks, quiet
    
    # Reset option globals so repeated calls (e.g. Jupyter) do not inherit previous settings
    debug = False
    newonly_mode = False
    motif_folder = None
    genint_options = ""
    keep_mode = False
    keep_pvalue = "3.5"
    keep_nofilter = False
    keep_flanks = 0
    quiet = False
    
    filtered_args = []
    distill_mode = False
    i = 0
    while i < len(args):
        if args[i] == "-debug":
            debug = True
        elif args[i] == "-newonly":
            newonly_mode = True
            if i + 1 < len(args) and not args[i + 1].startswith("-"):
                motif_folder = args[i + 1]
                i += 1  # Skip the folder path argument
            else:
                print("Error: -newonly requires a folder path")
                return [], False
        elif args[i] == "-distill":
            distill_mode = True
        elif args[i] == "-keep":
            keep_mode = True
        elif args[i] == "-quiet":
            quiet = True
        elif args[i] == "-nofilter":
            keep_nofilter = True
        elif args[i].startswith("-protect_flanks="):
            value = args[i][len("-protect_flanks="):]
            if not value.isdigit():
                print("Error: -protect_flanks=N requires a non-negative integer")
                return [], False
            keep_flanks = int(value)
        elif args[i] == "-keeppvalue":
            if i + 1 < len(args):
                try:
                    float(args[i + 1])
                except ValueError:
                    print("Error: -keeppvalue requires a number")
                    return [], False
                keep_pvalue = args[i + 1]
                i += 1  # Skip the value argument
            else:
                print("Error: -keeppvalue requires a number")
                return [], False
        elif args[i] == "-genint":
            if i + 1 < len(args):
                genint_options = args[i + 1]
                i += 1  # Skip the options string argument
            else:
                print("Error: -genint requires an options string")
                return [], False
        else:
            filtered_args.append(args[i])
        i += 1
    
    return filtered_args, distill_mode

def validate_arguments(filtered_args: List[str]) -> bool:
    """Validate command line arguments"""
    if len(filtered_args) < 6:
        print("Usage: ./autoseed.py <background_file> <signal_file> <shortest_kmer_length> "
              "<longest_kmer_length> <top_N_seeds> <max_extension> [count_threshold (default: 100)] "
              "[length_diff_cutoff (default: 0.35)] [multinomial (default: 1)] "
              "[similarity_cutoff (default: 0.20)] [-debug] [-newonly <motif_folder>] [-distill]")
        print("\nOptions:")
        print("  -debug                Enable debug output")
        print("  -newonly <folder>     Only accept motifs that are not similar to existing .pfm files in <folder>")
        print("  -distill              Run distillation mode")
        print("  -quiet                Do not print per-seed refinement output, only main steps and results")
        print("  -genint \"<options>\"   Extra options passed to both genint-PWM calls")
        print("  -keep                 With -distill: next iteration keeps matched features of current sequences (genint-PWM -keepfeatures)")
        print("  -keeppvalue <value>   -log10 p-value used in -keep mode (default: 3.5)")
        print("  -nofilter             -keep mode: re-sample once, no rescan for new matches (genint-PWM -nofilter)")
        print("  -protect_flanks=N     -keep mode: also keep N bases on each side of each feature (default: 0)")
        return False
    return True

def extract_parameters(filtered_args: List[str]) -> dict:
    """Extract and return parameters from filtered arguments"""
    return {
        'bg_file': filtered_args[0],
        'signal_file': filtered_args[1],
        'min_len': int(filtered_args[2]),
        'max_len': int(filtered_args[3]),
        'top_n': int(filtered_args[4]),
        'max_extension': int(filtered_args[5]),
        'count_threshold': int(filtered_args[6]) if len(filtered_args) > 6 else 100,
        'length_cutoff': float(filtered_args[7]) if len(filtered_args) > 7 else 0.35,
        'multinomial': int(filtered_args[8]) if len(filtered_args) > 8 else 1,
        'similarity_cutoff': float(filtered_args[9]) if len(filtered_args) > 9 else 0.20
    }

def setup_existing_motifs() -> List[str]:
    """Setup existing motifs if in newonly mode"""
    existing_motifs = []
    if newonly_mode and motif_folder:
        existing_motifs = load_existing_motifs(motif_folder)
        print(f"Running in -newonly mode with {len(existing_motifs)} existing motifs from {motif_folder}")
    return existing_motifs

def process_seed_refinement(result, idx, total_seeds, params, seed_counts, accepted_pwms) -> Optional[Tuple[str, str, str]]:
    """Process a single seed through refinement iterations"""
    global seed_history
    seed_history = [result.seed]  # Reset convergence history per initial seed, starting with the initial seed
    print(f"\nProcessing seed {idx}/{total_seeds}: {result.seed}")
    
    current_seed = result.seed
    current_seed_history = []
    iter_count = 0
    final_motif_output = ""
    max_iterations = 20
    
    # Capture original seed match counts
    if iter_count == 0:
        orig_motif_output, _ = run_localmax_motif(params['bg_file'], params['signal_file'],
                                                 result.seed, params['multinomial'])
        for line in orig_motif_output.splitlines():
            if line.startswith("Background:"):
                parts = line.split()
                bg_count = int(parts[1])
                bg_percent = float(parts[6].strip('()%'))
                seed_counts[result.seed] = {'signal': 0, 'background': bg_count, 'background_percent': bg_percent}
            elif line.startswith("Signal:"):
                parts = line.split()
                sig_count = int(parts[1])
                sig_percent = float(parts[6].strip('()%'))
                seed_counts[result.seed].update({'signal': sig_count, 'signal_percent': sig_percent})

    # Refinement loop
    while iter_count <= max_iterations:
        motif_output, new_seed = run_localmax_motif(params['bg_file'], params['signal_file'],
                                                   current_seed, params['multinomial'])

        # Extract counts and percentages
        for line in motif_output.splitlines():
            if line.startswith("Background:"):
                parts = line.split()
                bg_count = int(parts[1])
                bg_percent = float(parts[6].strip('()%'))
                seed_counts[new_seed] = {'signal': 0, 'background': bg_count, 'background_percent': bg_percent}
            elif line.startswith("Signal:"):
                parts = line.split()
                sig_count = int(parts[1])
                sig_percent = float(parts[6].strip('()%'))
                seed_counts[new_seed].update({'signal': sig_count, 'signal_percent': sig_percent})
                
        print(f"Iteration {iter_count+1}: Seed change {current_seed} -> {new_seed}", flush=True)
        
        if iter_count == max_iterations:
            print(f"Seed error: failed to converge in {iter_count+1} iterations")
            return None
        if not new_seed:
            print("Seed error: Empty seed")
            return None
        if iter_count > 0 and seeds_converged(current_seed, new_seed):
            print(f"Converged after {iter_count+1} iterations")
            final_motif_output = motif_output
            break
            
        current_seed_history.append(current_seed)
        current_seed = new_seed
        iter_count += 1
        final_motif_output = motif_output
    
    if not final_motif_output or not new_seed:
        print(f"Seed error for seed {result.seed}: no valid consensus found.")
        return None
    
    # Extract PWM matrix
    pwm_matrix = extract_pwm_matrix(final_motif_output)
    if len(pwm_matrix) != 4:
        print(f"PWM matrix extraction failed for seed {result.seed}")
        return None
    
    # Save and validate motif
    original_seed = current_seed_history[0] if current_seed_history else result.seed
    
    filename = f"refined_model_{idx}_for_seed_{original_seed}_wide.pfm"
    output_filename = filename
    if newonly_mode and motif_folder:
        output_filename = os.path.join(motif_folder, os.path.basename(filename))
    # Write full PWM file
    with open(output_filename, "w") as f:
        for row in pwm_matrix:
            f.write("\t".join(row) + "\n")
    filename = f"refined_model_{idx}_for_seed_{original_seed}.pfm"
    output_filename = filename
    if newonly_mode and motif_folder:
        output_filename = os.path.join(motif_folder, os.path.basename(filename))
    # Write seed-width PWM file
    with open(output_filename, "w") as f:
        for row in pwm_matrix:
            # Remove columns that are not part of seed
            seed_length = len(new_seed)
            # Always cut 5 from start
            start_trim = 5
            # Calculate end trim based on seed length
            if seed_length > 20:
                end_trim = 5 - (seed_length - 20)
                end_trim = max(0, end_trim)  # Don't go negative
            else:
                end_trim = 5

            # Apply trimming
            if end_trim > 0:
                trimmed_row = row[start_trim:-end_trim]
            else:
                trimmed_row = row[start_trim:]  # No end trimming if end_trim is 0

            f.write("\t".join(trimmed_row) + "\n")
            
    # Check similarity
    if is_motif_similar(output_filename, accepted_pwms, params['similarity_cutoff']):
        print(f"Seed {original_seed} REJECTED due to similarity with previous model")
        os.remove(output_filename)
        return None
    
    # Accept motif
    print(f"Seed {original_seed} -> {new_seed} ACCEPTED")
    svg_filename = output_filename.replace(".pfm", ".svg")
    generate_logo_from_pfm(output_filename, svg_filename, "Sequence Logo")
    print(f"Generated SVG logo: {svg_filename}")
    
    return (original_seed, new_seed, output_filename)

def extract_pwm_matrix(motif_output: str) -> List[List[str]]:
    """Extract PWM matrix from motif output"""
    pwm_matrix = []
    skipped = 0
    started = False
    for line in motif_output.splitlines():
        if "Motif from all matches" in line:
            started = True
            continue
        if "Match statistics:" in line:
            break
        if started and line and not line.startswith("Total"):
            nums = line.split()
            if len(nums) > 1:
                # Check if first column is A, C, G, or T (nucleotide labels)
                if nums[0].upper() in ['A', 'C', 'G', 'T']:
                    pwm_matrix.append(nums[1:])  # Skip nucleotide label
                elif skipped == 0 and nums[0].isdigit() and int(nums[0]) == 1:
                    skipped = 1
                    continue  # Skip position header row
                else:
                    pwm_matrix.append(nums)  # Keep all columns including first
    return pwm_matrix
    
    
def is_motif_similar(output_filename: str, accepted_pwms: List[str], similarity_cutoff: float) -> bool:
    """Check if motif is similar to any existing motifs"""
    if not accepted_pwms:
        return False
    
    for prev_pwm in accepted_pwms:
        if prev_pwm == output_filename:
            continue
        similarity = check_motif_similarity(prev_pwm, output_filename)
        if debug:
            print(f"Comparing {output_filename} with {prev_pwm}: similarity = {similarity}")
        if similarity > similarity_cutoff:
            if debug:
                print(f"Rejecting motif - similarity {similarity:.3f} to {prev_pwm}")
            return True
    return False

def run_autoseed_iteration(params: dict) -> List[Tuple[str, str, str]]:
    """Run a single autoseed iteration and return final seeds"""
    # Setup
    existing_motifs = setup_existing_motifs()
    seed_counts = {}
    
    if debug:
        print("\n=== Initial Parameters ===")
        for key, value in params.items():
            print(f"{key}: {value}")
        print("\n=== Starting Analysis ===")
    
    # Step 1: Find initial seeds
    print("\n1. Running localmax-motif to find short seeds...")
    initial_results = run_localmax_kmer(params['bg_file'], params['signal_file'],
                                       params['min_len'], params['max_len'],
                                       params['count_threshold'], params['length_cutoff'],
                                       params['top_n']+1)
    
    print(f"Found {len(initial_results)} initial seeds")
    if not initial_results:
        print("No initial seeds found. Exiting.")
        return [], seed_counts
    
    # Step 2: Extend seeds
    filtered_results = sorted(initial_results, key=lambda x: x.signal, reverse=True)[:params['top_n']]
    with open("short_seeds_filtered.tmp", 'w') as f:
        for result in filtered_results:
            f.write(f"{result.kmer}\t{result.gap}\t{result.position}\t{result.background}\t{result.signal}\t{result.is_localmax}\n")
    
    print(f"\n2. Running seedextender...")
    extended_results = run_seedextender(params['bg_file'], params['signal_file'],
                                       "short_seeds_filtered.tmp", params['max_extension'],
                                       params['count_threshold'], params['length_cutoff'], False)
    
    print(f"Found {len(extended_results)} extended seeds")
    if not extended_results:
        print("No extended seeds found. Exiting.")
        return [], seed_counts
    
    # Step 3: Refine seeds
    print("\n3. Starting PWM refinement for each extended seed...")
    final_seeds = []
    accepted_pwms = existing_motifs.copy()
    
    for idx, result in enumerate(extended_results, start=1):
        if quiet:
            # Suppress per-seed refinement output (-quiet); the result is unchanged
            with contextlib.redirect_stdout(io.StringIO()):
                refined_seed = process_seed_refinement(result, idx, len(extended_results), params, seed_counts, accepted_pwms)
        else:
            print("\n-----")
            refined_seed = process_seed_refinement(result, idx, len(extended_results), params, seed_counts, accepted_pwms)
        if refined_seed:
            final_seeds.append(refined_seed)
            accepted_pwms.append(refined_seed[2])  # Add the filename
    
    return final_seeds, seed_counts

def generate_random_background(pfm_file: str, n_sequences: int, output_file: str):
    """Write n_sequences random sequences sampled position-specifically from a background PFM
    (the same background model genint-PWM uses to generate / re-sample bases)"""
    freqs = normalize_pfm_to_pwm(read_pfm(pfm_file))  # 4 x width (A, C, G, T rows), columns sum to 1
    width = freqs.shape[1]
    bases = np.array(list("ACGT"))
    rng = np.random.default_rng()
    seqs = np.empty((n_sequences, width), dtype='<U1')
    for pos in range(width):
        seqs[:, pos] = rng.choice(bases, size=n_sequences, p=freqs[:, pos])
    with open(output_file, 'w') as f:
        for seq in seqs:
            f.write(''.join(seq) + "\n")

def seed_lists_identical(current_seeds: List[Tuple[str, str, str]],
                        previous_seeds: List[Tuple[str, str, str]]) -> bool:
    """Check if the current seed list is identical to previous iteration"""
    # Extract just the refined seeds (second element of each tuple)
    current_refined = [refined for _, refined, _ in current_seeds]
    previous_refined = [refined for _, refined, _ in previous_seeds]
    
    # Sort and compare
    return sorted(current_refined) == sorted(previous_refined)


def run(args: List[str]):
    global seed_history
    
    # Parse arguments
    filtered_args, distill_mode = parse_arguments(args)
    if not validate_arguments(filtered_args):
        return
    
    params = extract_parameters(filtered_args)
    seed_history = []  # Reset seed history
    
    # Distillation loop (breaks early if not in distill mode)
    current_signal_file = params['signal_file']
    iteration = 1
    max_iterations = 10
    
    while iteration <= max_iterations:
        print(f"\n{'='*60}")
        if distill_mode:
            print(f"DISTILLATION ITERATION {iteration}")
        else:
            print("RUNNING AUTOSEED")
        print(f"{'='*60}")
        
        # Update params with current signal file
        current_params = params.copy()
        current_params['signal_file'] = current_signal_file
        
        # Remember seeds from previous generation
        if iteration > 1:
            previous_iteration_seeds = final_seeds;
        else:
            previous_iteration_seeds = []
            frequency_tracker = {}
            
        # Run autoseed iteration
        final_seeds, seed_counts = run_autoseed_iteration(current_params)
        
        # Write results
        write_refined_models_list(final_seeds)
        
        # Display results
        print("\n=== Final Results ===")
        for original, refined, fname in final_seeds:
            print(f"{original}\t{refined}\t{fname}")
        
        if debug and seed_counts:
            print("\n=== Debug: All Seed Counts ===")
            for seed, counts in seed_counts.items():
                print(f"{seed}: signal={counts['signal']} ({counts.get('signal_percent', 0):.4f}%), "
                      f"background={counts['background']} ({counts.get('background_percent', 0):.4f}%)")

        print("\n----------------")
        if final_seeds:
            print(f"A TOTAL OF {len(final_seeds)} MODELS FOUND")
            if newonly_mode:
                print(f"New models saved to: {motif_folder}")
            print(f"Model list saved to: refined_models.txt")
            print("\nCreating merged SVG file...")
            create_merged_svg(final_seeds, seed_counts, iteration)
        else:
            print("NO MODELS FOUND")
        
        # Break early if not in distillation mode
        if not distill_mode:
            break
        
        # If no motifs found, stop distillation
        if not final_seeds:
            print("No motifs found. Stopping distillation.")
            break
        
        # Generate new sequences for next iteration using genint-PWM
        print(f"\nGenerating sequences for iteration {iteration + 1}...")
        
        # Run genint-PWM to count matches
        cmd_count = f"{GENINT_PATH} {BACKGROUND_PWM} refined_models.txt 200000 3.5 -file {current_signal_file} {genint_options}"
        print(f"Running: {cmd_count}")
        try:
            result = subprocess.run(cmd_count, shell=True, capture_output=True, text=True, check=True)
            #with open("match_results.txt", 'w') as f:
            #    f.write(result.stdout)
        except subprocess.CalledProcessError as e:
            print(f"Error running genint count: {e}")
            break
        
        # Generate new sequences
        next_signal_file = f"genint_iteration{iteration + 1}.seq"
        cmd_generate = f"{GENINT_PATH} {BACKGROUND_PWM} match_results.txt 200000 3.2 {genint_options} | cut -f1 > {next_signal_file}"
        if keep_mode:
            # Keep matched features of current sequences, re-sample all other bases (sequence lines only)
            cmd_generate = f"{GENINT_PATH} {BACKGROUND_PWM} refined_models.txt 200000 {keep_pvalue} -file {current_signal_file} -keepfeatures{' -nofilter' if keep_nofilter else ''}{f' -protect_flanks={keep_flanks}' if keep_flanks else ''} {genint_options} | cut -f1 | grep -E '^[ACGTNacgtn]+$' > {next_signal_file}"
        print(f"Running: {cmd_generate}")
        try:
            subprocess.run(cmd_generate, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error generating sequences: {e}")
            break
        
        # Check if generated file has sequences
        try:
            with open(next_signal_file, 'r') as f:
                lines = f.readlines()
            if len(lines) < 10:  # Minimum threshold
                print(f"Too few sequences generated ({len(lines)}). Stopping distillation.")
                break
            print(f"Generated {len(lines)} sequences for next iteration")
        except FileNotFoundError:
            print(f"Failed to create {next_signal_file}. Stopping distillation.")
            break
        
        # Check if seed lists are identical to previous iteration
        if seed_lists_identical(final_seeds, previous_iteration_seeds):
            print(f"\nCONVERGED! Identical seed list as previous iteration.")
            print("Stopping distillation due to convergence.")
            break
        
        # After iteration 1 (detection against the original background file), switch the background to
        # random sequences from the same background PFM that genint-PWM uses, so signal and background match
        if iteration == 1:
            with open(params['bg_file'], 'r') as f:
                n_background = sum(1 for line in f if line.strip() and set(line.strip().upper()) <= set("ACGTN"))
            random_bg_file = "random_background.seq"
            generate_random_background(BACKGROUND_PWM, n_background, random_bg_file)
            print(f"Background for next iterations: {random_bg_file} ({n_background} sequences from {BACKGROUND_PWM})")
            params['bg_file'] = random_bg_file
        
        # Update for next iteration
        current_signal_file = next_signal_file
        iteration += 1
    
    # Print frequency summary at the end
    if distill_mode and frequency_tracker:
        print_frequency_summary(frequency_tracker, iteration)
        print(f"\nDistillation completed after {iteration} iterations.")
    elif distill_mode:
        print(f"\nDistillation completed after {iteration} iterations.")
        
if __name__ == "__main__":
    run(sys.argv[1:])
