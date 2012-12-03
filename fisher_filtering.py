#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os, sys
import argparse
import numpy as np
from itertools import izip
from scipy import stats
try:
    import fisher as fisher_exact
except:
    print 'The much faster fisher library is not installed. Fallback to scipy.'


def fisher_filtering(control_file, affected_file, filtered_control_file, filtered_affected_file, cutoff):
    for control_line, affected_line in izip(open(control_file), open(affected_file)):
        c_chrom, c_start, c_end, c_cov, c_meth, c_strand = control_line.strip().split('\t')
        a_chrom, a_start, a_end, a_cov, a_meth, a_strand = affected_line.strip().split('\t')
        try:
            assert( c_chrom == a_chrom )
            assert( c_start == a_start )
            assert( c_end == a_end )
            assert( c_strand == a_strand )
        except AssertionError:
            raise
        control = (float(c_cov), float(c_meth))
        affected = (float(a_cov), float(a_meth))
        try:
            #Try to use the much faster fisher module from http://pypi.python.org/pypi/fisher/
            p = fisher_exact.pvalue(control[0], control[1], affected[0], affected[1])
            pvalue = p.two_tail
        except:
            #raise
            oddsratio, pvalue = stats.fisher_exact([control, affected], alternative='two-sided')

        if pvalue <= cutoff:
            filtered_control_file.write(control_line)
            filtered_affected_file.write(affected_line)
    filtered_affected_file.close()
    filtered_control_file.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extracting all reads that are differentialy methylated, according to the fisher exact test.')

    parser.add_argument("--control", required=True,
                    help="Path to the control file.")

    parser.add_argument("--affected", required=True,
                    help="Path to the affected file.")

    parser.add_argument('--ocontrol', required=True, type=argparse.FileType('w+'),
                     default=sys.stdout)
    parser.add_argument('--oaffected', required=True, type=argparse.FileType('w+'),
                     default=sys.stdout)

    parser.add_argument("--pvalue", default=0.05, type=float,
                    help="maximal pvalue to classify a site as differential methylated (default: 0.05 )")

    options = parser.parse_args()
    fisher_filtering(options.control, options.affected, options.ocontrol, options.oaffected, options.pvalue)






