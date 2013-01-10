#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os, sys
import argparse
import numpy as np
import tempfile
from scipy.stats.mstats import mquantiles
from itertools import izip
from scipy import stats
from destrand import merge
try:
    import fisher as fisher_exact
except:
    print 'The much faster fisher library is not installed. Fallback to scipy.'


class CpG():
    def __init__(self, chrom, start, end, strand, min_cov = 0, max_cov = 50):
        # cpg properties
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.cov_control = 0.0
        self.meth_control = 0.0
        self.cov_affected = 0.0
        self.meth_affected = 0.0
        self.delta = 0.0
        self.weighted_methylation_control = 0.0
        self.weighted_methylation_affected = 0.0

    def __repr__(self):
        return '%s\t%s\t%s control(c:%s, m:%s) affected(c:%s, m:%s) delta:%s' % (self.chrom, 
            self.start, self.strand, 
            self.cov_control, self.meth_control, self.cov_affected, 
            self.meth_affected, self.delta)

    def add_control(self, cov, meth):
        # if the cpg site does not fullfil the requirements, we add a zero one and increase the skip counter
        # that is necessary to not decrease the mean methylation state but be able to count these site as stopping criteria
        # if skip is already set, than the other sample failed already and we should set all values to null
        self.cov_control = cov
        self.meth_control = meth
        self._calculate_delta()
        self.weighted_methylation_control = meth * cov

    def add_affected(self, cov, meth):
        self.cov_affected = cov
        self.meth_affected = meth
        self._calculate_delta()
        self.weighted_methylation_affected = meth * cov

    def _calculate_delta(self):
        self.delta = self.meth_affected - self.meth_control


class Window():
    """
        Defines a window of differential methylated CpG sites.
        Such a window will start with a min_size and grow dynamically until a 
        treshold is undercuted.
    """
    def __init__(self, min_window_length = 4, max_cpg_distance = None, min_delta_methylation = 25, last_n = 4, allow_failed = None):
        # constraints
        self.min_window_length = min_window_length
        self.max_cpg_distance = max_cpg_distance
        self.min_delta_methylation = min_delta_methylation
        self.last_n = last_n # check last n postitions ot the cpgs if they are over min_delta_methylation
        self.allow_failed = allow_failed
        # window properties
        self.chrom = None
        self.start = None
        self.end = None
        self.cpgs = list()
        self.window_length = 0
        self.delta_sum = 0.0

        # save intermediate results for methylation calculation
        self.delta = 0.0

    def __len__(self):
        return self.window_length

    def __repr__(self):
        return 'Window: %s - %s, %s methylation sites and %s%% methylation.' % (self.start, self.end, self.window_length, self.delta)

    def _check_cpg(self, cpg):
        """
            Ordering of the evalutation is crucial
        """
        # if that is 0, than wen add the first site and we should skip the distance check
        last_site_position_end = self.get_last_cpg()
        if last_site_position_end != []:
            last_site_position_end = last_site_position_end.end
            if self.max_cpg_distance and (cpg.end - last_site_position_end) > self.max_cpg_distance:
                #print 'kick: max_cpg_distance', (cpg.end - self.get_last_cpg().end), cpg.end, - self.get_last_cpg().end, self.max_cpg_distance
                return False
        if len(self) + 1 < self.min_window_length:
            return True
        #print '---',self.calculate_window_methylation( self.get_last_n_cpgs( self.last_n ) ), len(self.get_last_n_cpgs( self.last_n )), self.last_n
        # get the last_n cpg site from the current window, n-1 because we want to check with the new cpg site
        if self.allow_failed:
            last_n_cpgs = self.get_last_n_cpgs( self.allow_failed -1 )
            # add new cpg site to the end
            last_n_cpgs.append(cpg)

            failed_cpgs = 0
            for test_cpg in last_n_cpgs:
                if abs(test_cpg.delta) < self.min_delta_methylation:
                    failed_cpgs += 1
            if failed_cpgs == self.allow_failed:
                #print '#########failed_cpg_in a row', last_n_cpgs
                return False

        if self.last_n:
            last_n_cpgs = self.get_last_n_cpgs( self.last_n -1 )
            last_n_cpgs.append(cpg)
            # check the last n cpg sites if they are in the mean over the min_delta_methylation
            if abs(self.calculate_window_methylation( last_n_cpgs )) < self.min_delta_methylation:
                return False

        return True

    def get_last_cpg(self):
        """
            returns the last cpg site
        """
        try:
            return self.cpgs[ -1 ]
        except IndexError:
            return []

    def get_last_n_cpgs(self, n):
        index = min( n, len(self) )
        return self.cpgs[ len(self) - index : ]

    def add_cpg(self, cpg):
        if self._check_cpg( cpg ):
            # only set if self.start is None, that is during the first cpg insert
            if not self.start:
                self.start = cpg.start
                self.chrom = cpg.chrom
                assert(self.start+1 == cpg.end)
            self.end = cpg.end

            self.cpgs.append(cpg)
            self.window_length += 1
            self.delta_sum += cpg.delta
            self.delta = self.delta_sum / len(self)
            return True
        else:
            return False

    def calculate_window_methylation(self, cpgs):
        try:
            return sum([cpg.delta for cpg in cpgs]) / len(cpgs)
        except ZeroDivisionError:
            return 0.0

    def calculate_differential_methylation_fisher_exact(self, weighted = False):
        sum_meth_control = 0
        sum_meth_affected = 0
        sum_cov_control = 0
        sum_cov_affected = 0
        for cpg in self.cpgs:
            if weighted:
                sum_meth_control += cpg.weighted_methylation_control
                sum_meth_affected += cpg.weighted_methylation_affected
                sum_cov_control += cpg.cov_control
                sum_cov_affected += cpg.cov_affected
            else:
                sum_meth_control += cpg.meth_control
                sum_meth_affected += cpg.meth_affected
                sum_cov_control += cpg.cov_control
                sum_cov_affected += cpg.cov_affected

        control = sum_meth_control / sum_cov_control
        affected = sum_meth_affected / sum_cov_affected
        control_methylated = sum_cov_control * control / 100
        control_unmethylated = sum_cov_control - control_methylated
        affected_methylated = sum_cov_affected * affected / 100
        affected_unmethylated = sum_cov_affected - affected_methylated
        try:
            #Try to use the much faster fisher module from http://pypi.python.org/pypi/fisher/
            p = fisher_exact.pvalue(control_methylated, control_unmethylated, affected_methylated, affected_unmethylated)
            pvalue = p.two_tail
        except:
            oddsratio, pvalue = stats.fisher_exact([(control_methylated, control_unmethylated), (affected_methylated, affected_unmethylated)], alternative='two-sided')
        return pvalue


    def strip_cpgs(self):
        # left strip
        lstrip = []
        while self.cpgs:
            # remove from the left end of the list all bad methylation sites
            # if the first good site is encounterd, leave the loop
            if abs(self.cpgs[0].delta) < self.min_delta_methylation:
                #print 'lstrip', self.cpgs[0]
                lstrip.append( self.cpgs.pop(0) )
            else:
                break
        # right strip
        rstrip = []
        while self.cpgs:
            # remove from the right end of the list all bad methylation sites
            # if the first good site is encounterd, leave the loop
            if abs(self.cpgs[-1].delta) < self.min_delta_methylation:
                #print 'rstrip', self.cpgs[-1]
                rstrip.append( self.cpgs.pop() )
            else:
                break
        if rstrip or lstrip:
            self.delta = self.calculate_window_methylation(self.cpgs)
            # it meight be faster to substract the lstrip and rstrip objects from skip_sites
            self.window_length = len(self.cpgs)
            self.start = self.cpgs[0].start
            self.end = self.cpgs[-1].end
        return lstrip, rstrip

    def write_to_bed_string(self, fisher = False, hyper = False, hypo = False):
        self.strip_cpgs()
        if len(self) < self.min_window_length:
            return False

        if abs( self.delta ) < self.min_delta_methylation:
            return False

        if hypo and self.delta > 0:
            return False
        if hyper and self.delta < 0:
            return False

        ret = "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (self.chrom, self.start, self.end, self.window_length, self.delta, 0, self.end - self.start)
        if fisher:
            pvalue = self.calculate_differential_methylation_fisher_exact()
            ret += '\t%e' % pvalue
            pvalue = self.calculate_differential_methylation_fisher_exact( weighted = True )
            ret += '\t%e\n' % pvalue
        else:
            ret += '\n'
        return ret

    def write_to_bed_file(self, handle, fisher = False, hyper = False, hypo = False):
        text = self.write_to_bed_string(fisher, hyper, hypo)
        if text:
            handle.write( text )


def dmr(options):

    control_quantil = None
    affected_quantil = None
    if options.filter_quantil:
        control_quantil = mquantiles( np.loadtxt(options.control, delimiter='\t', usecols=(3,)), prob = [options.filter_quantil])[0]
        affected_quantil = mquantiles( np.loadtxt(options.affected, delimiter='\t', usecols=(3,)), prob = [options.filter_quantil])[0]

    win = Window(options.min_window_length, options.max_cpg_distance, options.min_delta_methylation, options.check_last_n, options.allow_failed)
    old_chrom = False
    if options.destrand:
        tmp_control = tempfile.NamedTemporaryFile(delete=False, prefix='/home/bag/projects/')
        merge(open(options.control), open(tmp_control.name, 'w+'))
        tmp_control.close()
        tmp_affected = tempfile.NamedTemporaryFile(delete=False, prefix='/home/bag/projects/')
        tmp_affected.close()
        merge(open(options.affected), open(tmp_affected.name,'w+'))
        options.control = tmp_control.name
        options.affected = tmp_affected.name

    for control, affected in izip(open(options.control), open(options.affected)):
        c_chrom, c_start, c_end, c_cov, c_meth, c_strand = control.strip().split('\t')
        a_chrom, a_start, a_end, a_cov, a_meth, a_strand = affected.strip().split('\t')
        try:
            assert( c_chrom == a_chrom )
            assert( c_start == a_start )
            assert( c_end == a_end )
            assert( c_strand == a_strand )
        except AssertionError:
            raise

        c_cov, c_meth, a_cov, a_meth = map(float, [c_cov, c_meth, a_cov, a_meth])

        if old_chrom and old_chrom != c_chrom:
            # write window to the file if they fullfil the requirements
            if len(win) >= options.min_window_length:
                win.write_to_bed_file( options.outfile, options.fisher, options.hyper, options.hypo)
            # init a new window
            win = Window(options.min_window_length, options.max_cpg_distance, options.min_delta_methylation, options.check_last_n, options.allow_failed)

        if c_cov < options.min_cov or a_cov < options.min_cov:
            continue
        if c_cov > options.max_cov or a_cov > options.max_cov:
            continue
        if options.filter_quantil and (c_cov > control_quantil or a_cov > affected_quantil):
            continue

        cpg = CpG( c_chrom, int(c_start), int(c_end), c_strand )
        cpg.add_control( float(c_cov), float(c_meth) )
        cpg.add_affected( float(a_cov), float(a_meth) )
        if not win.add_cpg( cpg ):
            if len(win) >= options.min_window_length:
                win.write_to_bed_file( options.outfile, options.fisher, options.hyper, options.hypo)

            # create a new window with, if we have remainings frome the previouse window, add these at start cpgs
            win = Window(options.min_window_length, options.max_cpg_distance, options.min_delta_methylation, options.check_last_n, options.allow_failed)
        old_chrom = c_chrom

    if options.destrand:
        # remove tempfiles
        os.remove(tmp_affected.name)
        os.remove(tmp_control.name)


def main():
    parser = argparse.ArgumentParser(description='Extract differential methylated regions.')

    parser.add_argument("--control", required=True,
                    help="Path to the control file.")

    parser.add_argument("--affected", required=True,
                    help="Path to the affected file.")

    parser.add_argument('-o', '--outfile', type=argparse.FileType('w'),
                     default=sys.stdout)

    parser.add_argument("--mincov", dest="min_cov", default=0, type=int,
                    help="min coverage (default:0)")

    parser.add_argument("--maxcov", dest="max_cov", default=50, type=int,
                    help="max coverage (default:50)")

    parser.add_argument("--filter-quantil", dest="filter_quantil", default=None, type=float,
                    help="coverage quantil filter (default:None), example for the 99.9 quantil: 0.999")

    parser.add_argument("--min-window-length", dest="min_window_length", default=4, type=int,
                    help="minimal window length (default:4)")

    parser.add_argument("--min-delta-methylation", dest="min_delta_methylation", default=25, type=int,
                    help="minimal delta between the two mehylation states (default:25)")

    parser.add_argument("--check-last-n", dest="check_last_n", metavar = 'N', default=4, type=int,
                    help="check last N CpG sites if they fullfill all constraints (default:4)")

    parser.add_argument("--allow-failed", dest="allow_failed", metavar = 'N', default=None, type=int,
                    help="how many sites are allowed to fail in the checked last-n sites (default: )")

    parser.add_argument("--max-cpg-distance", dest="max_cpg_distance", default=None, type=int,
                    help="maximal CpG distance (default:None)")

    parser.add_argument('--fisher', action='store_true', default=False, help='Calculate the pvalue of each window with a fisher-exact-test')

    parser.add_argument('--destrand', action='store_true', default=False, help='Combine CpGs.')

    parser.add_argument('--hyper', action='store_true', default=False, help='Output only hyper methylated DMRs.')
    parser.add_argument('--hypo', action='store_true', default=False, help='Output only hypo methylated DMRs.')


    options = parser.parse_args()
    dmr(options)

if __name__ == '__main__':
    main()
