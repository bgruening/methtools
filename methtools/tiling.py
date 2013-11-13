#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import argparse
import sys
import sqlalchemy
import StringIO
from contextlib import closing

__doc__ = """
    Example methylation call bed file (it needs to be sorted):

    chr7	3295868	3295868	methylation with a coverage of 14	85.71	-
    chr7	4469971	4469971	methylation with a coverage of 30	6.67	+

    Tip:
        You can sort your BED file with bedtools, for example:
        bedtools sort -i data/control.bed > data/control_sort.bed

    The genome file should tab delimited and structured as follows:
         <chromName><TAB><chromSize>

        For example, Human (hg19):
        chr1    249250621
        chr2    243199373
        ...
        chr18**gl000207**random 4262

    Tip: 
        One can use the UCSC Genome Browser's MySQL database to extract
        chromosome sizes. For example, H. sapiens:

        mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from mm10.chromInfo" > mm10.genome

    There are two calculation modes:
        mean methylation -> sum( all methylation sites ) / count( methylation sites )
        methylation density -> sum( all methylation sites ) / window_length, where window_length is count( nucleotides ) in a window

        Default mode is mean methylation and if --density is set methylation density mode is on.

    Example program call for mouse genome mm10 and a control.bed file with window length 1000 and step-size 500:
        ./tiling.py -g mm10.genome -i control_sorted.bed -w 1000 -s 500 -o control_w1000_s500.bed
"""

def create_window_broders(window_length, step_size, chrom_size):
    """
        Arguments:
            window_length -- size of the window
            step_size -- size of the shift between to adjacent windows
            chrom_size -- size of a chromosome

        Return:
            iteratable of window borders -> [(0,1000), (500, 1500) ... ]

        Creates a iteratable with tuples.
        Each tuple represents the borders of a window.
    """
    start = None
    for window in zip(xrange(0, chrom_size, step_size), xrange(window_length, chrom_size, step_size)):
        start = window[1]
        yield window
    # return the last window with the exact end of the chromosome, that window is not window_length long but shorter
    if start:
        yield (start, chrom_size)

    # if no accurate mapping to chrom_size is found we end up here and need to raise an error
    yield (None, None)


class Windows():
    """
        Main class to handle the sliding windows and corresponding functions.
    """
    def __init__(self,options, chrom, genome_size):
        self.start = None
        self.stop = None
        self.chrom = chrom
        self._windows = create_window_broders( options.window_length, options.step_size, genome_size.get(chrom, -1))
        self.density = {
                '+': {},
                '-': {},
                }

    def new_window(self):
        self.start, self.stop = self._windows.next()
        if self.start == None and self.stop == None:
            """ the chromosome tag is probably not correct, so the end of one 
            chromosome given in the genome file does not fit with the 
            coordinates from the input file
            """
            sys.stderr.write('Given chromosome tag does not fit with the given coordinates.\n')
            sys.exit()

        # delete scores that do not belong to our window
        for base_position in self.density['+'].keys():
            if base_position < self.start:
                del self.density['+'][ base_position ]
        for base_position in self.density['-'].keys():
            if base_position < self.start:
                del self.density['-'][ base_position ]

    def update_desity(self, score, strand, base_position):
        try:
            assert not self.density[ strand ].has_key( base_position )
        except AssertionError:
            sys.stderr.write( 'Found douplicate base position: %s %s\n' % (self.chrom, base_position) )
            self.density[ strand ][ base_position ] += score
            return
        self.density[ strand ][ base_position ] = score

    def density_sum(self, strand = False):
        if not strand:
            return self._sum_density( self.density['+'] ) + self._sum_density( self.density['-'] )
        else:
            return self._sum_density( self.density[strand] )

    def _sum_density(self, density_dict):
        sum = 0.0
        for value in density_dict.values():
            sum += value
        return sum

    def get_methylation_sites_count(self, strand = False):
        if not strand:
            return len( self.density['+'].keys() ) + len( self.density['-'].keys() )
        else:
            return len( self.density[strand].keys() )


def read_genome_file(genome_path):
    """
        Arguments:
            genome_path -- file like object

        Return:
            dictionary with chromosome name <-> size mapping

        Reads a genome file into a dictionary. 
        The genome file in needed to know the borders of each chromosom to 
        adjust the windows properly.
        closing() is needed, because if the user specified an organism_tag the
        genome_path will be a StringIO Stream, which has no close definition.
    """
    genome_dict = dict()
    with closing(genome_path) as handle:
        for line in handle:
            line = line.strip()
            if line:
                chrom, size = line.split('\t')
                try:
                    size = int(size)
                except:
                    continue
                genome_dict[chrom] = size
    return genome_dict


def write_to_bedfile( options, chrom, name, window ):
    """
        Write one line to the options.output result file.
        Depending on the options.all_windows the result file is either in BED or in bedgraph format.

        There are two calculation modes:
            mean methylation -> sum over all methylation sites / #methylation sites
            methylation density -> sum over all methylation sites / window_length, where window_length is #nucleotides

            Default mode is mean methylation and if --density is set methylation density mode is on.
    """
    window_length = (window.stop - window.start)
    if options.merge_strands:

        density = window.density_sum()
        if options.density:
            # methylation density
            denominator = window_length
        else:
            # mean methylation
            denominator = window.get_methylation_sites_count()

        if density == 0.0 and not options.all_windows:
            # if we have no methylation state in that window and the option to write such states into the result file is False we do not wirte it
            return True
        elif density == 0.0:
            options.outfile.write( '%s\t%s\t%s\t%s\n' % (chrom, window.start, window.stop, '0.0') )
        else:
            options.outfile.write( '%s\t%s\t%s\t%s\n' % (chrom, window.start, window.stop, density / denominator) )

    else:
        density_forward = window.density_sum('+')
        density_reverse = window.density_sum('-')

        if options.density:
            # methylation density
            denominator_forward = window_length
            denominator_reverse = window_length
        else:
            # mean methylation
            denominator_forward = window.get_methylation_sites_count('+')
            denominator_reverse = window.get_methylation_sites_count('-')

        if density_forward != 0.0 or options.all_windows:
            options.outfile.write( '%s\t%s\t%s\t%s\t%s\t%s\n' % (chrom, window.start, window.stop, name, density_forward / denominator_forward, '+') )
        if density_reverse != 0.0 or options.all_windows:
            options.outfile.write( '%s\t%s\t%s\t%s\t%s\t%s\n' % (chrom, window.start, window.stop, name, density_reverse / denominator_reverse, '-') )


def tiling( options ):
    if options.genome_file:
        genome_size = read_genome_file( open(options.genome_file) )
    elif options.organism_tag:
        options.organism_tag = options.organism_tag.strip()
        try:
            engine = sqlalchemy.create_engine('mysql://genome@genome-mysql.cse.ucsc.edu')
            output = engine.execute("select chrom, size from %s.chromInfo" % options.organism_tag)
            output = ''.join(['%s\t%s\n' % (r[0],r[1]) for r in output])
            genome_file = StringIO.StringIO( output )
            genome_size = read_genome_file( genome_file )
        except:
            sys.exit('Fetching of the genome file failed! You need a working internet connection and mysql client installed.\nPlease download manually and specify it ith --genome_file.')

    else:
        sys.exit('Please specify a genome file or an organism tag.')

    window_start = 0
    window_counter = 0
    old_chrom = None
    blacklist_chrom = False
    with options.infile as meth_call:
        # for every bedfile line
        for line_counter, line in enumerate( meth_call ):
            line = line.strip()
            if line:
                chrom, base_start, base_end, name, score, strand = line.split('\t')
                base_start = int(base_start)
                score = float(score)
                if chrom == blacklist_chrom:
                    continue
                if chrom != old_chrom:
                    # leave one window and create a new one, in fact leave the whole chromosome
                    if old_chrom != None:
                        # during the first iteration the old_chr is null, in that case do not write out the results
                        write_to_bedfile( options, old_chrom, name, windows )
                    windows = Windows( options, chrom, genome_size ) #create_window_broders( options.window_length, options.step_size, genome_size.get(chrom, -1))
                    old_chrom = chrom

                    blacklist_chrom = False

                    windows.new_window()
                    window_start = windows.start
                    window_end = windows.stop
                    """
                        We get (None, None) if no mapping with to the chromosome
                        size is found in genome_size, that happens when obscure 
                        chromosome ids are present in the input files and not 
                        in the genome-file
                    """
                    if window_start == None:
                        blacklist_chrom = chrom
                        old_chrom = None
                        continue

                if base_start >= window_start and base_start < window_end:
                    windows.update_desity(score, strand, base_start)
                else:
                    # leave one window and create a new one
                    write_to_bedfile( options, chrom, name, windows )

                    windows.new_window()
                    window_start = windows.start
                    window_end = windows.stop

                    #if window_start < 0:
                    #    # happens when no chrom_size is present in the mapping file specified with the -g option
                    #    continue
                    while True:
                        """
                            iterate over all windows until the next methylation
                            site pops up
                        """
                        if base_start >= window_start and base_start < window_end:
                            windows.update_desity(score, strand, base_start)
                            # leave the while loop
                            break
                        else:
                            """
                                in that window no mehtylation site exist
                                write to output if options.all_windows is set
                            """
                            write_to_bedfile( options, chrom, name, windows )

                        windows.new_window()
                        window_start = windows.start
                        window_end = windows.stop

        if chrom != blacklist_chrom:
            write_to_bedfile( options, chrom, name, windows )

    options.outfile.close()

def main():
    parser = argparse.ArgumentParser(
        description='Calcualtes methylation desity for a given sequence window.',
        epilog = __doc__,
        formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--infile', type=argparse.FileType('r'),
        default=sys.stdin)
    parser.add_argument('-o', '--outfile', type=argparse.FileType('w'),
                     default=sys.stdout)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-g', '--genome_file', help="Genome file should tab delimited and structured as follows: <chrom name><TAB><chrom size>")
    group.add_argument('--organism-tag', dest='organism_tag' ,help="If no genome file is specified we try to download the file with the unique organism identifier, e.g. mm10 or hg18")
    parser.add_argument('-m','--merge-strands', dest='merge_strands', action='store_true', default=False, help='Sum up all methylation sites independent from the strand. In that case the output will be a BED-graph file.')
    parser.add_argument('--all-windows', dest="all_windows", action='store_true', default=False, help='Write also windows with no methylation sites to the result file - default: False')
    parser.add_argument('-w', '--window-length', dest="window_length", type=int, default=1000, help='Length of the sliding window - default: 1000')
    parser.add_argument('-s', '--step-size', dest="step_size", type=int, default=500, help='Step size - default: 500')
    parser.add_argument('--density', action='store_true', default=False, 
        help='Calculate the methylation density: Sum over all methylation sites / nucleotides (window_length). Default calculation mode is the mean methylation: Sum over all methylation sites / methylated sites')

    options = parser.parse_args()

    tiling(options)

if __name__ == '__main__':
    main()
