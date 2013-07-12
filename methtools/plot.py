#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import argparse
import os, sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import multiprocessing
import copy
from collections import defaultdict
import math
import signal
import tempfile

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
    ./plot_neighborhood.py -i ControlPCM_rmdup_methcall_fixed_w1000_s200.bedgraph -c chrX -p 103460373 -w 4000 --scale --gene-name XIST --image xist.svg
    ./plot_neighborhood.py -i ControlPCM_rmdup_methcall_fixed_sort.bed -c chrX -p 103460373 -w 4000 --bed --gene-name XIST --image xist_total.svg

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
                    options.infiles = ifile.split()
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
                    options.yaxis_max = 110
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

def parse_coordinatefile( options ):
    """
    """
    if options.coordinatefile == '-':
        coordinatefile = sys.stdin
    else:
        coordinatefile = open(options.coordinatefile)

    for line in coordinatefile:
        if line.strip() == '':
            continue
        if not line.startswith('#'):
            try:
                chrom, start, end, foo = line.strip().split('\t', 3)
                options.chromosome = chrom
                dist = int(end) - int(start)   #  144311276-144311255 = 21
                rounded_win_len = dist + 1000#(int(math.ceil( dist / 100.0)) * 100) + 1000 # -> 200
                options.position = int(start) + (dist / 2 ) # 144311255 + 11
                options.window_length = rounded_win_len
                #options.gene_name = ""options.position
                options.image = os.path.join(os.getcwd(),'chr%s_%s_%s.png' % (chrom, start, end))
                options.text = os.path.join(os.getcwd(),'chr%s_%s_%s.txt' % (chrom, start, end))
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


def moving_average(x, n, type='simple'):
    """
    compute an n period moving average.

    type is 'simple' | 'exponential'
    """
    x = np.asarray(x)
    if type=='simple':
        weights = np.ones(n)
    else:
        weights = np.exp(np.linspace(-1., 0., n))

    weights /= weights.sum()

    a =  np.convolve(x, weights, mode='full')[:len(x)]
    a[:n] = a[n]
    return a


def plot( scores, positions, options, overlay=False, no_dots=False):
    """
        overlay is True if many dots are plotted, in that case we adjust the settings
    """
    colors = ['orange', 'green', 'blue', 'red']
    fig = plt.figure(figsize = (20,6))
    ax = fig.add_subplot(111)

    # the score and positions are a list ot lists, the user inputs more than one file and wants to have two different plots in on figure
    if not overlay:#type(scores[0]) == list:
        #print 'multiple plotting mode'
        index = 0
        # save the plots for the legend afterwarts
        plots = list()
        labels = list()
        for sco, pos in zip(scores, positions):
            color = colors[index]
            pl = ax.plot(pos, sco, '-o', ms=15, lw=2, alpha=0.7, color=color)
            plots.append(pl)
            labels.append( os.path.basename(options.infiles[index]) )
            index += 1
        #ax.legend( plots, labels )
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])
        ax.legend(labels, loc='upper center', bbox_to_anchor=(0.5, 1.20),
          ncol=1, fancybox=True, shadow=True, numpoints=1, markerscale=10)
        # workaround for broken markerscale
        #for lines in l.get_lines():
        #    lines._legmarker.set_ms(10)
    else:
        #if overlay:
        """
            creating data for the moving average
        """
        unique_positions = set(positions)
        result = []
        pos_counter = 0
        positions = np.asarray(positions)
        scores = np.asarray(scores)

        for pos in unique_positions:
            temp = scores[ positions==pos ]
            temp = np.average(temp)
            result.append((pos, temp))

        ma_pos, ma_scores = zip(*sorted(result))
        if not no_dots:
            ax.plot(positions, scores, '.', ms=3, lw=2, alpha=0.7, color='orange')
        #ax.plot(ma_pos, moving_average(ma_scores, 7), color='blue', lw=2, label='MA (7)')
        ax.plot(ma_pos, moving_average(ma_scores, 20), color='red', lw=2, label='MA (20)')
        ax.legend()

    ax.grid()
    # default for options.yaxis_max is 110

    if overlay:
        yaxis_min = -105
        yaxis_annotation = -104
    else:
        yaxis_min = -5
        yaxis_annotation = -4

    ax.vlines(0, yaxis_min, options.yaxis_max, color='r', linewidth=5.0, linestyle=':', alpha=0.7)

    if options.gene_name:
        tss = 'chr: %s; pos: %s; %s' % (options.chromosome, options.position, options.gene_name)
    elif overlay:
        tss = 'Overlay point'
    else:
        tss = 'chr: %s; pos: %s' % (options.chromosome, options.position)

    el = Ellipse((2, -1), 0.5, 0.5)
    ax.annotate(tss, xy=(0, yaxis_annotation),  xycoords='data',
                xytext=(-100, -90), textcoords='offset points',
                size=20,
                arrowprops=dict(arrowstyle='fancy',#"wedge,tail_width=0.7",
                                fc="0.6", ec="none",
                                patchB=el,
                                connectionstyle="arc3,rad=-0.3"),
                )

    # default for yaxis is 105
    ax.set_ylim([yaxis_min, options.yaxis_max])
    extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(options.image, bbox_inches=extent.expanded(1.1, 1.6))


def plot_regions( options ):
    """
        Main plotting function.
        The functions assumes that the bed- or bedgraph file is sorted.

        If we find a base that fits our fiter criteria we save that state,
        so we can leave the bedline loop as soon as the creteria does not 
        match anymore. In a orted file that means that we can safely leave 
        the loop without missing bases.
    """
    if options.text:
        outfiles = list()
        if len(options.infiles) > 1:
            for counter, foo in enumerate(options.infiles):
                temp = list(os.path.splitext( options.text ))
                temp.insert(1, '_%s' % (counter+1))
                outfiles.append( ''.join(temp) )
        else:
            outfiles.append( options.text )

    distance = options.window_length / 2
    scores = list()
    positions = list()
    coverage = list()

    for infile, outfile in zip( options.infiles, outfiles ):
        if options.text and not options.overlay_only:
            outfile = open( outfile, 'w+' )
        hit_found = False
        scores_temp = list()
        positions_temp = list()
        coverage_temp = list()

        """
        INSTALL
        
        http://pysam.googlecode.com/files/pysam-0.7.4.tar.gz
        python ./setup.py install --home $INSTALL_DIR
        
        http://wwwfgu.anat.ox.ac.uk/~andreas/documentation/samtools/api.html

        Index Bedfile beforehand:
        pysam.tabix_compress() tabix_compress(filename_in, filename_out, force=False)
        pysam.tabix_index()

        tabix_index(filename, force=False, seq_col=None, start_col=None, end_col=None, preset=None, meta_char=’#’, zerobased=False)
        If preset is provided, the column coordinates are taken from a preset. Valid values for preset are “gff”, “bed”, “sam”, “vcf”, psltbl”, “pileup”.

        pysam.Tabixfile(filename, mode=’r’)
        opens a tabix file for reading. A missing index (filename + ”.tbi”) will raise an exception.

        contigs
        chromosome names

        fetch
        Tabixfile.fetch(self, reference=None, start=None, end=None, region=None, parser=None)

        fetch one or more rows in a region using 0-based indexing. The region is specified by reference, start and end. Alternatively, a samtools region string can be supplied.

        Without reference or region all entries will be fetched.

        If only reference is set, all reads matching on reference will be fetched.

        If parser is None, the results are returned as an unparsed string. Otherwise, parser is assumed to be a functor that will return parsed data (see for example asTuple() and asGTF()).

        """
    """
    for infile, outfile in zip( options.infiles, outfiles ):
        hit_found = False
        scores_temp = list()
        positions_temp = list()
        coverage_temp = list()

        sort -k1,1 -k2,2n Galaxy4.bed | bgzip > sorted_galaxy4.bed.gz
        compressed_sorted_bedfile = tempfile.NamedTemporaryFile( suffix='.gz' ).name
        tmp_out = open( compressed_sorted_bedfile, 'wb' )

        cmd = 'sort -k1,1 -k2,2n %s' % ( infile )

        p1 = subprocess.Popen( args=shlex.split( cmd ), stdout=subprocess.PIPE )
        proc = subprocess.Popen( ['bgzip'], stdin=p1.stdout, stdout=tmp_out )
        returncode = proc.wait()
        if returncode != 0:
            raise Exception, 'sorting and creating tabix index failed'
        tmp_out.close()

        cmd = "tabix %s %s:%s-%s" % (compressed_sorted_bedfile, options.chromosome, (options.position - distance), (options.position + distance))
        proc = subprocess.Popen( args=shlex.split( cmd ), stdout=outfile )
        returncode = proc.wait()

        if returncode != 0:
            raise Exception, 'tabix reading failed'

        for bed_line in open( outfile ):
            try:
                if options.bed:
                    chrom, start, stop, cov, score, strand = bed_line.strip().split()
                    # im fall von 2 files, schmeist er nur einen raus, evt beides implementieren? TODO
                    if float(cov) < options.min_cov:
                        continue
                else:
                    chrom, start, stop, score = bed_line.strip().split()
            except:
                sys.exit('Wrong input format. Have a look at the --bed option.')
            start = int(start)
            score = float(score)

            if options.scale:
                scores_temp.append( score * 100 )
            else:
                scores_temp.append( score )

            if options.reverse:
                positions_temp.append( (start - options.position)*-1 )
            else:
                positions_temp.append( start - options.position )
            
            coverage_temp.append( cov )

        scores.append(scores_temp)
        positions.append(positions_temp)
        coverage.append(coverage_temp)

        if not options.text or options.overlay_only:
            os.remove( outfile )
    """



        ########################################################################## old method ###################################

        for bed_line in open(infile):
            if bed_line.startswith(options.chromosome):
                try:
                    if options.bed:
                        chrom, start, stop, cov, score, strand = bed_line.strip().split()
                        # im fall von 2 files, schmeist er nur einnen raus, evt beides implementieren? TODO
                        if float(cov) < options.min_cov:
                            continue
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
                if start >= (options.position - distance) and start <= (options.position + distance):
                    hit_found = True
                    if options.text and not options.overlay_only:
                        outfile.write( bed_line )
                    if options.scale:
                        scores_temp.append( score * 100 )
                    else:
                        scores_temp.append( score )

                    if options.reverse:
                        positions_temp.append( (start - options.position)*-1 )
                    else:
                        positions_temp.append( start - options.position )

                    coverage_temp.append( cov )
            elif hit_found:
                break
        scores.append(scores_temp)
        positions.append(positions_temp)
        coverage.append(coverage_temp)
        if options.text and not options.overlay_only:
            outfile.close()

    ######################### old zu ende ############

    if options.image and not options.overlay_only:
        plot( scores, positions, options )

    return (scores, positions, coverage)

positions_scores = list()
merged_handle = open('all_merged.txt', 'w+')
position_dict = defaultdict(dict)

def log_results( results ):
    """
        log all results generated by the multiprocessing pool workers
        for each gene region we get a list of scores (methylation levels), gene positions and the corresponding coverage
    """
    s, p, c = results
    # if you habe specified two file for one gene region, iterate over it
    position_dict = defaultdict(dict)
    for pos_list, score_list, cov_list in zip(p, s, c):
        for pos, score, cov in zip(pos_list, score_list, cov_list):
            if pos in position_dict:
                # second value for that position
                position_dict[pos]['delta'] -= score
                position_dict[pos]['meth2'] = score
                position_dict[pos]['cov2'] = cov
            else:
                # first value for that position
                position_dict[pos]['delta'] = score
                position_dict[pos]['meth1'] = score
                position_dict[pos]['cov1'] = cov

    for pos, value_dict in position_dict.items():
        # write position into the dictionary to use it in the output string
        value_dict.update({'pos': pos})
        if not 'meth2' in value_dict:
            # only one file is processed
            merged_handle.write('%(pos)s\t%(meth1)s\t%(cov1)s\t%(delta)s\n' % (value_dict))
        else:
            merged_handle.write('%(pos)s\t%(meth1)s\t%(cov1)s\t%(meth2)s\t%(cov2)s\t%(delta)s\n' % (value_dict))
        positions_scores.append( [pos, value_dict['delta']] )


def main():
    parser = argparse.ArgumentParser(
        description='Plots the methylation surrounding of one DNA spot.',
        epilog = __doc__,
        formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--infiles', nargs='+',
        help='Input file with all methylation sites aggregated into windows. In BED or bedgraph format.')
    parser.add_argument('--image', help='Path to the resulting image file. If not specified no image will be created.')
    parser.add_argument('--text', help='Path to the resulting coordinate file. If not specified no file will be created.')

    parser.add_argument('-w', '--window-length', dest='window_length', type=int, default=1000, help='Length of the surrounding window - default: 1000')
    parser.add_argument('-c', '--chromosome', help='Chromosome name', default=None)
    parser.add_argument('-r', '--reverse', action='store_true', default=False, help='Invert the x-axis to reflect the reverse strand of the gene.')
    parser.add_argument('-p', '--position', type=int, help='Chromosome position', default=None)
    parser.add_argument('--scale', action='store_true', default=False, help='scale the scores to be between 0 and 100')
    parser.add_argument('--gene-name', dest='gene_name', help='If given use that name in the plot.')
    parser.add_argument('--bed', action='store_true', default=False, help='Input is in BED format, rather than in bedgraph format.')
    parser.add_argument('--configfile', default=None, help='Specify your input parameters in a config file. See more information below.')
    parser.add_argument('--coordinatefile', default=None, help='Extract the coordinates from a bed file.')
    parser.add_argument('--y-max', dest='yaxis_max', type=int, default=110, help='Set the maximum on the yaxis.')
    parser.add_argument('--overlay-only', dest='overlay_only', action='store_true', default=False, help='In configfile mode, only create the overlay image')
    parser.add_argument('--no-dots', dest='no_dots', action='store_true', default=False, help='In overlay image, omit the dots only print the smooth line')
    parser.add_argument('--min-coverage', dest='min_cov', type=int, default=0, help='Minimal allowed coverage to plot a methylation site. Default: 0')

    parser.add_argument('--processors', type=int, default=multiprocessing.cpu_count())

    options = parser.parse_args()

    if [options.coordinatefile, options.configfile].count(None) == 2 and (options.position == None or options.chromosome == None):
        sys.exit('If you do not use --configfile you need to specify the chromosome (-c), the position (-p) and the inputfile (-i).')

    if options.configfile or options.coordinatefile:
        """
            The config file mode is used to execute plottings in batches. Each line in such a configfile is used to plot one gene region.
            If two input files are given in one line, its also possible to plot control vs. affected per gene.
            In the end all gene will be plotted in one plot using the TSS as overlay-point.
        """

        pool = multiprocessing.Pool( options.processors )
        if options.configfile:
            pool.map_async(plot_regions, parse_configfile(options), callback=log_results)
        else:
            # coordinatefile
            pool.map_async(plot_regions, parse_coordinatefile(options), callback=log_results)
            
            # for debugging!
            
            #for opt in parse_coordinatefile(options):
            #    print opt
            #    #pool.apply_async( plot_regions, args=opt, callback=log_results)
            #    log_results( plot_regions(opt) )

        pool.close()
        pool.join()

        positions, scores = zip(*sorted( positions_scores ))
        options.image = 'all_merged.png'
        options.gene_name = None
        plot( scores, positions, options, overlay=True, no_dots=options.no_dots)

        #with open('all_merged.txt', 'w+') as handle:
        #    for pos, score, cov in zip(positions, scores, coverage):
        #        handle.write('%s\t%s\t%s\n' % (pos, score, cov))

    else:
        check_options( options )
        plot_regions(options)


if __name__ == '__main__':
    main()

