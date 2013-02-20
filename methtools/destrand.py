#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os, sys
import argparse

def merge_sites(c1, c2, keep_positions = False):
    """
        Merge two lines of a bed file together into one line if they only differ in one position and are located on different strands.
        In that case these two lines represent one symetric CpG and we can merge the coverage and methylation rate.
    """
    c1_chrom, c1_start, c1_end, c1_cov, c1_meth, c1_strand = c1.strip().split('\t')
    c2_chrom, c2_start, c2_end, c2_cov, c2_meth, c2_strand = c2.strip().split('\t')
    if c1_strand != c2_strand and int(c1_end) == (int(c2_end) - 1):
        # merge with weigthed mean
        c1_cov, c2_cov, c1_meth, c2_meth = map(float, [c1_cov, c2_cov, c1_meth, c2_meth])
        merged_cov = c1_cov + c2_cov
        merged_meth = '%.02f' % ( ((c1_meth * c1_cov) + (c2_meth * c2_cov)) / merged_cov )
        
        if keep_positions:
            start = min(c1_start, c2_start)
            end = max(c1_end, c2_end)

        if c1_strand == '+':
            if keep_positions:
                cm = c1_chrom, start, end, str(merged_cov), merged_meth, c1_strand
            else:
                cm = c1_chrom, c1_start, c1_end, str(merged_cov), merged_meth, c1_strand
        else:
            if keep_positions:
                cm = c2_chrom, start, end, str(merged_cov), merged_meth, c2_strand
            else:
                cm = c2_chrom, c2_start, c2_end, str(merged_cov), merged_meth, c2_strand

        return True, cm
    else:
        c1 = (c1_chrom, c1_start, c1_end, c1_cov, c1_meth, c1_strand)
        return False,c1


def merge(sample, outfile, keep_positions = False):
    previouse_side = None
    for sample_line in sample:
        if not previouse_side:
            previouse_side = sample_line
            continue
        else:
            merged, merged_sample = merge_sites (previouse_side, sample_line, keep_positions)

            if merged:
                previouse_side = None
                outfile.write( '\t'.join(merged_sample) + '\n' )
            else:
                previouse_side = sample_line
                outfile.write( '\t'.join(merged_sample) + '\n' )
    if previouse_side:
        outfile.write(previouse_side)
    outfile.close()


def main():
    parser = argparse.ArgumentParser(
        description='Merge CpGs together that are located next to each other. Methylation is symetric, so we can use that trick to enhance the coverage.')

    parser.add_argument("-i", "--infile", required=True, 
        type=argparse.FileType('r'),
        help="Path to the sample file.")

    parser.add_argument('-o', '--outfile', required=True, 
        type=argparse.FileType('w'),
        default=sys.stdout)

    parser.add_argument('-k', '--keep-positions', dest="keep_positions", action='store_true', default=False, help='Keep the position from both sites.')

    options = parser.parse_args()
    merge(options.infile, options.outfile, options.keep_positions)


if __name__ == '__main__':
    main()






