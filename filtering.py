#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os, sys
import argparse
import numpy as np
from itertools import izip
from scipy.stats.mstats import mquantiles
from scipy import stats
try:
    import fisher as fisher_exact
except:
    print 'The much faster fisher library is not installed. Fallback to scipy.'


"""
    That file needs intersected inputfiles, so that each site is present in both files, affected and control.
"""

def filtering(control_file, affected_file, filtered_control_file, filtered_affected_file, max_pvalue = None, min_cov = None, max_cov = None, min_delta_methylation = None, filter_quantil = None):

    control_quantil = None
    affected_quantil = None
    if filter_quantil:
        control_quantil = mquantiles( np.loadtxt(control_file, delimiter='\t', usecols=(3,)), prob = [filter_quantil])[0]
        affected_quantil = mquantiles( np.loadtxt(affected_file, delimiter='\t', usecols=(3,)), prob = [filter_quantil])[0]
        print control_quantil, affected_quantil

    for control_line, affected_line in izip(open(control_file), open(affected_file)):
        c_chrom, c_start, c_end, c_cov, c_meth, c_strand = control_line.strip().split('\t')
        a_chrom, a_start, a_end, a_cov, a_meth, a_strand = affected_line.strip().split('\t')
        try:
            assert( c_chrom == a_chrom )
            assert( c_start == a_start )
            assert( c_end == a_end )
            assert( c_strand == a_strand )
        except AssertionError:
            sys.exit('That file needs intersected inputfiles, so that each site is present in both files, affected and control.')

        c_cov, c_meth, a_cov, a_meth = map(float, [c_cov, c_meth, a_cov, a_meth])
        if min_cov != None and (a_cov < min_cov or c_cov < min_cov):
            continue
        if max_cov != None and (a_cov > max_cov or c_cov > max_cov):
            continue
        if min_delta_methylation != None and abs(a_meth - c_meth) < min_delta_methylation:
            continue
        if filter_quantil and (c_cov > control_quantil or a_cov > affected_quantil):
            continue

        if max_pvalue != None:
            control_methylated = c_cov * c_meth / 100
            control_unmethylated = c_cov - control_methylated
            affected_methylated = a_cov * a_meth / 100
            affected_unmethylated = a_cov - affected_methylated
            try:
                #Try to use the much faster fisher module from http://pypi.python.org/pypi/fisher/
                p = fisher_exact.pvalue(control_methylated, control_unmethylated, affected_methylated, affected_unmethylated)
                pvalue = p.two_tail
            except:
                oddsratio, pvalue = stats.fisher_exact([(control_methylated, control_unmethylated), (affected_methylated, affected_unmethylated)], alternative='two-sided')

            if pvalue > max_pvalue:
                continue

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

    parser.add_argument("--pvalue", default=None, type=float,
                    help="maximal pvalue to classify a site as differential methylated")

    parser.add_argument("--min-coverage", dest="min_cov", default=None, type=int,
                    help="minimal allowed coverage")

    parser.add_argument("--max-coverage", dest="max_cov", default=None, type=int,
                    help="maximal allowed coverage")

    parser.add_argument("--min-delta-methylation", dest="min_delta_methylation", default=None, type=float,
                    help="Minimal delta between the two mehylation states")

    parser.add_argument("--quantil", dest="filter_quantil", default=None, type=float,
                    help="coverage quantil filter, example for the 99.9 quantil: 0.999")


    options = parser.parse_args()
    if [options.pvalue, options.min_cov, options.max_cov, options.filter_quantil].count(None) == 4:
        sys.exit('You need to specify at least one filter parameter: --pvalue, --min-coverage, --quantil or --max-coverage')
    filtering(options.control, options.affected, options.ocontrol, options.oaffected, options.pvalue, options.min_cov, options.max_cov, options.min_delta_methylation, options.filter_quantil)






