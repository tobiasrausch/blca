#!/usr/bin/env python3

import argparse
import gzip
from collections import defaultdict

def parse_cna_file(cna_file):
    breakpoints = []    
    with open(cna_file, 'r') as f:
        header = f.readline()
        prev_line = None        
        for line in f:
            if prev_line is None:
                prev_line = line.strip().split('\t')
                continue
            curr_line = line.strip().split('\t')
            prev_chrom = prev_line[0]
            curr_chrom = curr_line[0]
            prev_end = int(prev_line[2])
            curr_start = int(curr_line[1])
            try:
                prev_total_cn = float(prev_line[3]) if prev_line[3] != 'NA' else None
                curr_total_cn = float(curr_line[3]) if curr_line[3] != 'NA' else None
            except (ValueError, IndexError):
                prev_line = curr_line
                continue            
            if prev_chrom == curr_chrom and prev_total_cn is not None and curr_total_cn is not None and prev_total_cn != curr_total_cn:
                breakpoint_pos = prev_end
                breakpoints.append((prev_chrom, breakpoint_pos))            
            prev_line = curr_line
    return breakpoints


def parse_sv_file(sv_file):
    sv_breakpoints = []
    open_func = gzip.open if sv_file.endswith('.gz') else open
    mode = 'rt' if sv_file.endswith('.gz') else 'r'
    with open_func(sv_file, mode) as f:
        header = f.readline()        
        for line in f:
            fields = line.strip().split('\t')
            chrom1 = fields[0]
            start1 = int(fields[1])
            end1 = int(fields[2])
            chrom2 = fields[3]
            start2 = int(fields[4])
            end2 = int(fields[5])
            sv_breakpoints.append((chrom1, (start1 + end1) // 2))
            sv_breakpoints.append((chrom2, (start2 + end2) // 2))    
    return sv_breakpoints


def parse_l1_file(l1_file):
    l1_positions = []    
    with open(l1_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            pos = (start + end) // 2
            l1_positions.append((chrom, pos))
    return l1_positions


def is_near(pos1, chrom1, pos2, chrom2, distance):
    if chrom1 != chrom2:
        return False
    return abs(pos1 - pos2) <= distance


def classify_cna_breakpoints(cna_breakpoints, sv_breakpoints, distance):
    sv_explained = []
    loose_ends = []    
    for cna_chrom, cna_pos in cna_breakpoints:
        has_sv = False
        for sv_chrom, sv_pos in sv_breakpoints:
            if is_near(cna_pos, cna_chrom, sv_pos, sv_chrom, distance):
                has_sv = True
                break        
        if has_sv:
            sv_explained.append((cna_chrom, cna_pos))
        else:
            loose_ends.append((cna_chrom, cna_pos))    
    return sv_explained, loose_ends


def count_l1_near_breakpoints(breakpoints, l1_positions, distance):
    count = 0
    
    for bp_chrom, bp_pos in breakpoints:
        has_l1 = False
        
        for l1_chrom, l1_pos in l1_positions:
            if is_near(bp_pos, bp_chrom, l1_pos, l1_chrom, distance):
                has_l1 = True
                break
        
        if has_l1:
            count += 1
    
    return count


def main():
    parser = argparse.ArgumentParser(description='Analyze somatic L1s near CNV breakpoints')
    parser.add_argument('--cna', required=True, help='CNV file')
    parser.add_argument('--sv', required=True, help='SV BEDPE file')
    parser.add_argument('--l1', required=True, help='L1 BED file')
    parser.add_argument('--sv-distance', type=int, default=5000, help='Distance for SV-CNA overlap (default: 5000)')
    parser.add_argument('--l1-distance', type=int, default=5000, help='Distance for L1-CNA overlap (default: 5000)')
    args = parser.parse_args()

    ## Parse input files
    cna_breakpoints = parse_cna_file(args.cna)
    print(f"Found {len(cna_breakpoints)} CNA breakpoints")
    sv_breakpoints = parse_sv_file(args.sv)
    print(f"Found {len(sv_breakpoints)} SV breakpoints")
    l1_positions = parse_l1_file(args.l1)
    print(f"Found {len(l1_positions)} L1 insertions")

    ## Get SV-explained CNAs and loose ends
    sv_explained, loose_ends = classify_cna_breakpoints(cna_breakpoints, sv_breakpoints, args.sv_distance)
    print(f"SV-explained breakpoints: {len(sv_explained)}")
    print(f"Loose ends: {len(loose_ends)}")

    ## Count L1s
    l1_near_sv_explained = count_l1_near_breakpoints(sv_explained, l1_positions, args.l1_distance)
    l1_near_loose_ends = count_l1_near_breakpoints(loose_ends, l1_positions, args.l1_distance)
    frac_sv_explained = None
    if len(sv_explained) > 0:
        frac_sv_explained = l1_near_sv_explained / len(sv_explained)
        print(f"SV-explained breakpoints with nearby L1: {l1_near_sv_explained}/{len(sv_explained)} ({frac_sv_explained:.2%})")
    else:
        print(f"SV-explained breakpoints with nearby L1: 0/0 (N/A)")
    frac_loose_ends = None
    if len(loose_ends) > 0:
        frac_loose_ends = l1_near_loose_ends / len(loose_ends)
        print(f"Loose ends with nearby L1: {l1_near_loose_ends}/{len(loose_ends)} ({frac_loose_ends:.2%})")
    else:
        print(f"Loose ends with nearby L1: 0/0 (N/A)")
    fold_enrichment = None
    if len(sv_explained) > 0 and len(loose_ends) > 0:
        fold_enrichment = (frac_loose_ends / frac_sv_explained) if frac_sv_explained > 0 else float('inf')
        print(f"Fold enrichment at loose ends: {fold_enrichment:.2f}x")
    print(f"Stats\t{args.l1}\t{args.sv_distance}\t{args.l1_distance}\t{len(l1_positions)}\t{frac_sv_explained}\t{frac_loose_ends}\t{fold_enrichment}")

if __name__ == '__main__':
    main()
