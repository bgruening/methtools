#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import argparse
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import multiprocessing
import copy

__doc__ = """

    To enable patch processing from the commandline you can specify a
    configfile with --configfile. That file should have all options TAB
    separated in one line. Every line is a new experiment.
    The configfile _must_ have the following colums , TAB separated:

        input file
        chromosome
        position
        window length
        gene name
        image path
        text path
        scale   -- leave empty if not desired otherwise insert True or scale
        bed     -- leave empty if not desired otherwise insert True or bed
        yaxis max   -- leave empty if not desired otherwise insert a number


    Examples:
    ./plot_neighborhood.py -i ControlPCM_rmdup_methcall_fixed_w1000_s200.bedgraph -c chrX -p 103460373 -w 4000 --scale --gene_name XIST --image xist.svg
    ./plot_neighborhood.py -i ControlPCM_rmdup_methcall_fixed_sort.bed -c chrX -p 103460373 -w 4000 --bed --gene_name XIST --image xist_total.svg

    ./plot_neighborhood.py --configfile example.cfg

"""

def parse_configfile( options ):
    """
        To enable patch processing from the commandline you can specify a
        configfile with --configfile. That file should have all options TAB
        separated in one line. Every line is a new experiment.
    """
    with open(options.configfile) as conf:
        for line in conf:
            if line.strip() == '':
                continue
            if not line.startswith('#'):
                try:
                    ifile, chrom, pos, winlen, reverse, gname, impath, txtpath, scale, bed, yaxis_max = line.split('\t')
                    options.infile = ifile
                    options.chromosome = chrom
                    options.position = int(pos)
                    options.window_length = int(winlen)
                    options.gene_name = gname
                    options.image = impath
                    options.text = None
                    if txtpath.strip():
                        options.text = txtpath
                    options.scale = False
                    options.bed = False
                    options.reverse = False
                    options.yaxis_max = 105
                    if scale.strip():
                        options.scale = True
                    if bed.strip():
                        options.bed = True
                    if reverse.strip():
                        options.reverse = True
                    if yaxis_max.strip():
                        try:
                            options.yaxis_max = int(yaxis_max)
                        except:
                            sys.stdout.write('WARNING:\n\tThe yaxis value in the following configfile line is not an integer. \n\t%s\n' % line )
                except:
                    raise
                    sys.stdout.write('WARNING:\n\tThe following line in your configfile has an error. \n\t%s\n\tProbably a column is missing.\n' % line )
                    continue
                check_options( options, exit = False )
                yield copy.deepcopy( options )

def check_options( options, exit = True ):
    """
        Add a few checks to the provided user options.
    """
    if not options.image and not options.text:
        if not options.gene_name:
            if exit:
                sys.exit('Please specify a path to a resulting image and/or a coordinate file with --image or --text.')
            else:
                sys.stdout.write('Warning for %s %s:\n\tPlease specify a path to a resulting image and/or a coordinate file with --image or --text.\n' % (options.chromosome, options.position))
        else:
            options.image = '%s.svg' % options.gene_name
            options.text = '%s.text' % options.gene_name
            sys.stdout.write('No --image and --text is specified, we will create both files with the given --gene_name (%s).\n' % options.gene_name)

def plot( scores, positions, options ):
    fig = plt.figure(figsize = (20,6))
    ax = fig.add_subplot(111)
    ax.plot(positions, scores,'-o', ms=15, lw=2, alpha=0.7, mfc='orange')
    ax.grid()
    # default for options.yaxis_max is 105
    ax.vlines(0, -5, options.yaxis_max, color='r', linewidth=5.0, linestyle=':', alpha=0.7)

    if options.gene_name:
        tss = 'TSS - %s' % options.gene_name
    else:
        tss = 'TSS'

    el = Ellipse((2, -1), 0.5, 0.5)
    ax.annotate(tss, xy=(0, 0),  xycoords='data',
                xytext=(-100, -100), textcoords='offset points',
                size=20,
                arrowprops=dict(arrowstyle='fancy',#"wedge,tail_width=0.7",
                                fc="0.6", ec="none",
                                patchB=el,
                                connectionstyle="arc3,rad=-0.3"),
                )

    # default is 105
    ax.set_ylim([-5, options.yaxis_max])
    extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(options.image, bbox_inches=extent.expanded(1.1, 1.6))


def main( options ):
    """
        Main plotting function.
        The functions assumes that the bed- or bedgraph file is sorted.
        
        If we find a base that fits our fiter criteria we save that state,
        so we can leave the bedline loop as soon as the creteria does not 
        match anymore. In a orted file that means that we can safely leave 
        the loop without missing bases.
    """
    if options.text:
        options.text = open(options.text, 'w+')
    
    distance = options.window_length / 2
    scores = list()
    positions = list()
    hit_found = False
    for bed_line in open(options.infile):
        if bed_line.startswith(options.chromosome):
            try:
                if options.bed:
                    chrom, start, stop, text, score, strand = bed_line.strip().split()
                else:
                    chrom, start, stop, score = bed_line.strip().split()
            except:
                sys.exit('Wrong input format. Have a look at the --bed option.')

            if chrom.strip() != options.chromosome and hit_found:
                # if the chromosome did not match but we already encouterd a hit
                # in previous iterations we can safely leave the loop
                break
            if chrom.strip() != options.chromosome:
                continue

            start = int(start)
            score = float(score)
            if start >= options.position - distance and start <= options.position + distance:
                hit_found = True
                if options.text:
                    options.text.write( bed_line )
                if options.scale:
                    scores.append( score * 100 )
                else:
                    scores.append( score )
                
                if options.reverse:
                    positions.append( (start - options.position)*-1 )
                else:
                    positions.append( start - options.position )
        elif hit_found:
            break

    if options.image:
        plot( scores, positions, options )

    if options.text:
        options.text.close()

    return (scores, positions)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Plots the methylation surrounding of one DNA spot.',
        epilog = __doc__,
        formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--infile',
        help='Input file with all methylation sites aggregated into windows. In BED or bedgraph format.')
    parser.add_argument('--image', help='Path to the resulting image file. If not specified no image will be created.')
    parser.add_argument('--text', help='Path to the resulting coordinate file. If not specified no file will be created.')

    parser.add_argument('-w', '--window_length', type=int, default=1000, help='Length of the surrounding window - default: 1000')
    parser.add_argument('-c', '--chromosome', help='Chromosome name')
    parser.add_argument('-r', '--reverse', action='store_true', default=False, help='Invert the x-axis to reflect the reverse strand of the gene.')
    parser.add_argument('-p', '--position', type=int, help='Chromosome position')
    parser.add_argument('--scale', action='store_true', default=False, help='scale the scores to be between 0 and 100')
    parser.add_argument('--gene_name', help='If given use that name in the plot.')
    parser.add_argument('--bed', action='store_true', default=False, help='Input is in BED format, rather than in bedgraph format.')
    parser.add_argument('--configfile', help='Specify your input parameters in a config file. See more information below.')
    parser.add_argument('--yaxis_max', type=int, default=105, help='Set the maximum on the yaxis.')

    parser.add_argument('--processors', type=int, default=multiprocessing.cpu_count())

    options = parser.parse_args()

    if not options.configfile and not options.position and not options.chromosome:
        sys.exit('If you do not use --configfile you need to specify --chromosome and --positions.')

    if options.configfile:
        scores = list()
        positions = list()

        p = multiprocessing.Pool( options.processors )
        results = p.map_async(main, parse_configfile(options))
        results_iterator = results.get()

        for s, p in results_iterator:
            scores += s
            positions += p

        positions, scores = zip(*sorted(zip(positions, scores)))
        options.image = 'all_merged.png'
        options.gene_name = 'Superimpose'
        plot( scores, positions, options )

        with open('all_merged.txt', 'w+') as handle:
            for pos, score in zip(positions, scores):
                handle.write('%s\t%s\n' % (pos, score))
    else:
        check_options( options )
        main(options)


    test = """
chr1	17907070	17907070	cov:23	100.00	+
chr1	17907071	17907071	cov:25	96.00	-
chr1	17907114	17907114	cov:15	26.67	+
chr1	17907115	17907115	cov:14	35.71	-
chr1	17907178	17907178	cov:10	70.00	+
chr1	17907307	17907307	cov:10	40.00	+
chr1	17907326	17907326	cov:11	36.36	+
chr1	17907384	17907384	cov:17	100.00	+
chr1	17907393	17907393	cov:15	100.00	+
"""
    test = """
chr18	61414351	61414351	cov:19	100.00	+
chr18	61414352	61414352	cov:19	100.00	-
chr18	61414432	61414432	cov:22	63.64	+
chr18	61414433	61414433	cov:18	100.00	-
chr18	61414439	61414439	cov:20	95.00	+
chr18	61414440	61414440	cov:15	100.00	-
chr18	61414451	61414451	cov:20	100.00	+
chr18	61414460	61414460	cov:18	100.00	+
chr18	61414978	61414978	cov:10	70.00	+
chr18	61415032	61415032	cov:11	100.00	-
chr18	61415060	61415060	cov:12	100.00	+
chr18	61415475	61415475	cov:14	100.00	-
chr18	61415748	61415748	cov:16	75.00	+
chr18	61415749	61415749	cov:28	100.00	-
chr18	61415785	61415785	cov:17	94.12	-
chr18	61415816	61415816	cov:13	84.62	+
chr18	61415817	61415817	cov:14	57.14	-
chr18	61416205	61416205	cov:19	100.00	+
chr18	61416206	61416206	cov:15	93.33	-
chr18	61416277	61416277	cov:16	50.00	-
chr18	61416351	61416351	cov:16	0.00	-
chr18	61416401	61416401	cov:17	5.88	+
chr18	61416402	61416402	cov:43	2.33	-
chr18	61416702	61416702	cov:26	0.00	+
chr18	61416703	61416703	cov:18	5.56	-
chr18	61416794	61416794	cov:20	10.00	+
chr18	61416795	61416795	cov:15	0.00	-
chr18	61416809	61416809	cov:18	0.00	+
chr18	61416810	61416810	cov:17	0.00	-
chr18	61416826	61416826	cov:21	4.76	+
chr18	61416827	61416827	cov:20	0.00	-
chr18	61416883	61416883	cov:17	5.88	+
chr18	61416884	61416884	cov:20	15.00	-
chr18	61416959	61416959	cov:12	58.33	+
chr18	61416960	61416960	cov:12	16.67	-
chr18	61416993	61416993	cov:13	30.77	+
chr18	61417066	61417066	cov:10	30.00	-
chr18	61417084	61417084	cov:16	6.25	-
chr18	61417175	61417175	cov:13	100.00	+
chr18	61417176	61417176	cov:16	75.00	-
chr18	61417202	61417202	cov:14	85.71	+
chr18	61417230	61417230	cov:15	100.00	+
chr18	61417288	61417288	cov:11	90.91	-
chr18	61417313	61417313	cov:14	92.86	-
chr18	61417596	61417596	cov:14	100.00	+
chr18	61417597	61417597	cov:25	52.00	-
chr18	61417822	61417822	cov:10	60.00	-
chr18	61418073	61418073	cov:10	50.00	-
chr18	61418173	61418173	cov:17	94.12	+
"""


