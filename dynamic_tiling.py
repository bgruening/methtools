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
    def __init__(self, chrom, start, end, strand, min_cov = 0, max_cov = 50, control_filter_quantil = None, affected_filter_quantil = None):
        # constraints
        self.min_cov = min_cov
        self.max_cov = max_cov
        self.control_filter_quantil = control_filter_quantil
        self.affected_filter_quantil = affected_filter_quantil
        # cpg properties
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.cov_control = 0.0
        self.cov_control_original = 0.0
        self.meth_control = 0.0
        self.meth_control_original = 0.0
        self.cov_affected = 0.0
        self.cov_affected_original = 0.0
        self.meth_affected = 0.0
        self.meth_affected_original = 0.0
        self.weighted_methylation_control = 0.0
        self.weighted_methylation_affected = 0.0
        self.skip = 0

    def __repr__(self):
        return '%s\t%s\t%s control(c:%s, m:%s) affected(c:%s, m:%s) skip:%s\n' % (self.chrom, 
            self.start, self.strand, 
            self.cov_control, self.meth_control, self.cov_affected, self.meth_affected, self.skip)

    """
    @property
    def meth_control(self):
        return self.meth_control

    @meth_control.setter
    def meth_control(self, meth):
        self.meth_control = meth
        self.meth_control_original = meth
        print 'hhhhhhhhhhhhhhhhhhhhhhh'
    """
    def add_control(self, cov, meth):
        # if the cpg site does not fullfil the requirements, we add a zero one and increase the skip counter
        # that is necessary to not decrease the mean methylation state but be able to count these site as stopping criteria
        # if skip is already set, than the other sample failed already and we should set all values to null
        if (self.control_filter_quantil and cov >= self.control_filter_quantil) or \
            cov < self.min_cov or cov > self.max_cov or \
            self.skip == 1:
            self.cov_control = 0.0
            self.cov_control_original = cov
            self.meth_control = 0.0
            self.meth_control_original = meth
            self.cov_affected = 0.0
            self.meth_affected = 0.0
            self.skip = 1
        else:
            self.skip = 0
            self.cov_control = cov
            self.cov_control_original = cov
            self.meth_control = meth
            self.meth_control_original = meth

    def add_affected(self, cov, meth):
        if (self.affected_filter_quantil and cov >= self.affected_filter_quantil) or \
            cov < self.min_cov or cov > self.max_cov or self.skip == 1:
            self.cov_control = 0.0
            self.meth_control = 0.0
            self.cov_affected = 0.0
            self.cov_affected_original = cov
            self.meth_affected = 0.0
            self.meth_affected_original = meth
            self.skip = 1
        else:
            self.skip = 0
            self.cov_affected = cov
            self.cov_affected_original = cov
            self.meth_affected = meth
            self.meth_affected_original = meth

    def calculate_weighted_methylation(self):
        self.weighted_methylation_control = self.meth_control * self.cov_control
        self.weighted_methylation_affected = self.meth_affected * self.cov_affected



class Window():
    """
        Defines a window of differential methylated CpG sites.
        Such a window will start with a min_size and grow dynamically until a 
        treshold is undercuted.
    """
    def __init__(self, min_window_length = 4, max_cpg_distance = None, min_delta_methylation = 25, last = 4):
        # constraints
        self.min_window_length = min_window_length
        self.max_cpg_distance = max_cpg_distance
        self.min_delta_methylation = min_delta_methylation
        self.last = last # check last n postitions ot the cpgs if they are over min_delta_methylation
        # window properties
        self.chrom = None
        self.start = None
        self.end = None
        self.cpgs = list()
        self.window_length = 0
        self.skip_sites = 0

        # save intermediate results for methylation calculation
        self.sum_meth_affected = 0.0
        self.sum_meth_control = 0.0
        self.sum_cov_control = 0.0
        self.sum_cov_affected = 0.0
        self.delta = 0.0

    def __len__(self):
        return self.window_length - self.skip_sites

    def __repr__(self):
        return 'Window: %s - %s, %s methylation sites and %s%% methylation.' % (self.start, self.end, self.window_length, self.delta)

    def _check_cpg(self, cpg):
        """
            Ordering of the evalutation is crucial
        """
        # if that is 0, than wen add the first site and we should skip the distance check
        last_site_position_end = self.get_last_cpg()
        if last_site_position_end != 0:
            last_site_position_end = last_site_position_end.end
            if self.max_cpg_distance and (cpg.end - last_site_position_end) > self.max_cpg_distance:
                #print 'kick: max_cpg_distance', (cpg.end - self.get_last_cpg().end), cpg.end, - self.get_last_cpg().end, self.max_cpg_distance
                return False

        if self.window_length + 1 < self.min_window_length:
            return True
        
        #print 'adrgg', self.window_length, cpg.differential_methylation
        #print self.methylation
        #print cpg.differential_methylation
        try:
            con = (self.sum_meth_control + cpg.weighted_methylation_control) / (self.sum_cov_control + cpg.cov_control)
            aff = (self.sum_meth_affected + cpg.weighted_methylation_affected) / (self.sum_cov_affected + cpg.cov_affected)
            delta = abs( con - aff )
        except ZeroDivisionError:
            delta = 0.0

        if  delta < self.min_delta_methylation:
            #print 'kick: differential meth to low', self.delta, delta, self.min_delta_methylation, cpg.start
            return False

        if abs( self.calculate_window_methylation( self.get_last_n_cpgs( self.last ) ) ) < self.min_delta_methylation:
            #print self.get_last_n_cpgs( self.last )
            #print self.calculate_window_methylation( self.get_last_n_cpgs( self.last ) )
            #print 'kick: last_n', self.calculate_window_methylation( self.get_last_n_cpgs( self.last ) )
            return False
        return True

    def get_last_cpg(self):
        """
            returns the last cpg site that is not marked as skiped
        """
        rev_iterator = -1
        while True:
            try:
                if self.cpgs[ rev_iterator ].skip == 0:
                    return self.cpgs[ rev_iterator ]
                else:
                    rev_iterator -= 1
            except IndexError:
                return 0.0
    """
    def cpg_pop(self):

            #Delete the last CpG postition and decrease all variables

        temp_cpg = self.cpgs.pop()
        self.skip_sites -= temp_cpg.skip
        self.window_length -= temp_cpg.skip
        if self.cpgs:
            self.end = self.cpgs[-1].end
        else:
            self.end = 0
        if temp_cpg.skip == 0:
            self.sum_meth_affected -= temp_cpg.weighted_methylation_affected
            self.sum_meth_control -= temp_cpg.weighted_methylation_control
            self.sum_cov_control -= temp_cpg.cov_control
            self.sum_cov_affected -= temp_cpg.cov_affected
            try:
                con = self.sum_meth_control / self.sum_cov_control
                aff = self.sum_meth_affected / self.sum_cov_affected
                self.delta = abs( con - aff )
            except ZeroDivisionError:
                self.delta = 0
        
        return temp_cpg
    """
    def add_cpg(self, cpg, destrand = False):
        if self._check_cpg(cpg):
            if cpg.skip:
                return True
            # only set if self.start is None, that is during the first cpg insert
            if not self.start:
                self.start = cpg.start
                assert(self.start+1 == cpg.end)
            self.end = cpg.end

            self.chrom = cpg.chrom
            self.cpgs.append(cpg)
            self.window_length += 1
            self.skip_sites += cpg.skip
            if cpg.skip == 0:
                self.sum_meth_affected += cpg.weighted_methylation_affected
                self.sum_meth_control += cpg.weighted_methylation_control
                self.sum_cov_control += cpg.cov_control
                self.sum_cov_affected += cpg.cov_affected
                con = self.sum_meth_control / self.sum_cov_control
                aff = self.sum_meth_affected / self.sum_cov_affected
                self.delta = con - aff

            return True
        else:
            return False

    def calculate_window_methylation(self, cpgs):
        sum_meth_control = 0.0
        sum_meth_affected = 0.0
        sum_cov_control = 0.0
        sum_cov_affected = 0.0
        for cpg in cpgs:
            if cpg.skip == 0:
                sum_meth_affected += cpg.weighted_methylation_affected
                sum_meth_control += cpg.weighted_methylation_control
                sum_cov_control += cpg.cov_control
                sum_cov_affected += cpg.cov_affected
        try:
            control = sum_meth_control / sum_cov_control
            affected = sum_meth_affected / sum_cov_affected
        except ZeroDivisionError:
            return 0.0
        return control - affected


    def calculate_differential_methylation_fisher_exact(self):
        """
        sum_meth_control = 0
        sum_meth_affected = 0
        for cpg in self.cpgs:
            if cpg.skip == 0:
                sum_meth_control += cpg.meth_control
                sum_meth_affected += cpg.meth_affected
        control_methylated = self.sum_cov_control / len(self) * sum_meth_control / len(self) / 100
        control_unmethylated = self.sum_cov_control - control_methylated
        affected_methylated = self.sum_cov_affected / len(self) * sum_meth_affected / len(self) / 100
        affected_unmethylated = self.sum_cov_affected - affected_methylated
        """

        control = self.sum_meth_control / self.sum_cov_control
        affected = self.sum_meth_affected / self.sum_cov_affected
        control_methylated = self.sum_cov_control * control / 100
        control_unmethylated = self.sum_cov_control - control_methylated
        affected_methylated = self.sum_cov_affected * affected / 100
        affected_unmethylated = self.sum_cov_affected - affected_methylated
        try:
            #Try to use the much faster fisher module from http://pypi.python.org/pypi/fisher/
            p = fisher_exact.pvalue(control_methylated, control_unmethylated, affected_methylated, affected_unmethylated)
            pvalue = p.two_tail
        except:
            oddsratio, pvalue = stats.fisher_exact([(control_methylated, control_unmethylated), (affected_methylated, affected_unmethylated)], alternative='two-sided')
        """print pvalue, len(self)

        print control_methylated, control_unmethylated
        for c in self.cpgs:
            print c
        sys.exit()"""
        return pvalue


    def strip_cpgs(self):
        # left strip
        lstrip = []
        while self.cpgs:
            # remove from the left end of the list all bad methylation sites
            # if the first good site is encounterd, leave the loop
            if self.cpgs[0].skip:
                lstrip.append( self.cpgs.pop(0) )
            else:
                break
        # right strip
        rstrip = []
        while self.cpgs:
            # remove from the right end of the list all bad methylation sites
            # if the first good site is encounterd, leave the loop
            if self.cpgs[-1].skip:
                rstrip.append( self.cpgs.pop() )
            else:
                break
        if rstrip or lstrip:
            self.delta = self.calculate_window_methylation(self.cpgs)
            # it meight be faster to substract the lstrip and rstrip objects from skip_sites
            self.skip_sites = sum([s.skip for s in self.cpgs])
            self.window_length = len(self.cpgs)
            self.start = self.cpgs[0].start
            self.end = self.cpgs[-1].end
        return lstrip, rstrip

    def get_last_n_cpgs(self, n):
        index = max(n, self.window_length)
        return self.cpgs[ self.window_length - index:]

    def write_to_bed_string(self, fisher = False):
        ret = "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (self.chrom, self.start, self.end, self.window_length, self.delta, 0, self.end - self.start)
        if fisher:
            pvalue = self.calculate_differential_methylation_fisher_exact()
            ret += '\t%e\n' % pvalue
        else:
            ret += '\n'
        return ret

    def write_to_bed_file(self, handle, fisher = False):
        handle.write( self.write_to_bed_string(fisher) )
        #handle.flush()


def main(options):

    control_quantil = None
    affected_quantil = None
    if options.filter_quantil:
        control_quantil = mquantiles( np.loadtxt(options.control, delimiter='\t', usecols=(3,)), prob = [options.filter_quantil])[0]
        affected_quantil = mquantiles( np.loadtxt(options.affected, delimiter='\t', usecols=(3,)), prob = [options.filter_quantil])[0]

    win = Window(options.min_window_length, options.max_cpg_distance, options.min_delta_methylation, options.check_last_n)
    old_chrom = False
    if options.destrand:
        tmp_control = tempfile.NamedTemporaryFile(delete=False, prefix='/home/bag/projects/')
        merge(open(options.control), open(tmp_control.name, 'w+'))
        tmp_control.close()
        tmp_affected = tempfile.NamedTemporaryFile(delete=False, prefix='/home/bag/projects/')
        tmp_affected.close()
        merge(open(options.affected), open(tmp_affected.name,'w+'))
        print 'merging done'
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

        if old_chrom and old_chrom != c_chrom:
            # write window to the file if they fullfil the requirements
            if len(win) >= options.min_window_length and abs(win.delta) >= options.min_delta_methylation:
                win.strip_cpgs()
                if len(win) >= options.min_window_length:
                    win.write_to_bed_file( options.outfile, options.fisher )
            win = Window(options.min_window_length, options.max_cpg_distance, options.min_delta_methylation, options.check_last_n)

        cpg = CpG( c_chrom, int(c_start), int(c_end), c_strand, options.min_cov, options.max_cov, control_quantil, affected_quantil )

        cpg.add_control( float(c_cov), float(c_meth) )
        cpg.add_affected( float(a_cov), float(a_meth) )

        cpg.calculate_weighted_methylation()

        if not win.add_cpg( cpg, options.destrand):
            if len(win) >= options.min_window_length and abs(win.delta) >= options.min_delta_methylation:
                #if counter > 30:
                #    sys.exit()
                #print win.delta
                #for a in win.cpgs:
                #    print 'k',a
                # strip bad cpgs from both ends and remember how many are striped from the right side
                #print win, win.delta
                left, right_side_cpgs = win.strip_cpgs()

                """if win.start == 42293256 and win.chrom == '12':
                    for b in [(a.start,a.end,a.differential_methylation, a.skip,win.window_length) for a in win.cpgs]:
                        print b
                    print win.write_to_bed_string(  )
                    sys.exit()
                """
                a  = False
                if win.start > 9997837:
                    a = True
                # if after stripping the window is still larger than the min_window_length size, write it out
                if len(win) >= options.min_window_length:
                    win.write_to_bed_file( options.outfile, options.fisher )
                if len(win) == 4 and abs(win.delta) < options.min_delta_methylation:
                    print 'wat?', abs(win.delta)
                    for a in win.cpgs:
                        print a
                    sys.exit()
                # create a new window with, if we have remainings frome the previouse window, add these at start cpgs
                win = Window(options.min_window_length, options.max_cpg_distance, options.min_delta_methylation, options.check_last_n)

                if right_side_cpgs:
                    #print 'right_side_cpgs'
                    for cpg in right_side_cpgs:
                        win.add_cpg( cpg )

            else:
                win = Window(options.min_window_length, options.max_cpg_distance, options.min_delta_methylation, options.check_last_n)

        old_chrom = c_chrom
    if options.destrand:
        # remove tempfiles
        os.remove(tmp_affected.name)
        os.remove(tmp_control.name)

if __name__ == '__main__':
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
                    help="quantil filter (default:None), example for the 99.9 quantil: 0.999")

    parser.add_argument("--min-window-length", dest="min_window_length", default=4, type=int,
                    help="minimal window length (default:4)")

    parser.add_argument("--min-delta-methylation", dest="min_delta_methylation", default=25, type=int,
                    help="minimal delta between the two mehylation states (default:25)")

    parser.add_argument("--check-last-n", dest="check_last_n", metavar = 'N', default=4, type=int,
                    help="check last N CpG sites if they fullfill all constraints (default:4)")

    parser.add_argument("--max-cpg-distance", dest="max_cpg_distance", default=None, type=int,
                    help="maximal CpG distance (default:None)")

    parser.add_argument('--fisher', action='store_true', default=False, help='Calculate the pvalue of each window with a fisher-exact-test')

    parser.add_argument('--destrand', action='store_true', default=False, help='Combine CpGs.')

    options = parser.parse_args()
    main(options)
