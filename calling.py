#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os, sys
import argparse
import pysam
from collections import defaultdict
import numpy
import multiprocessing
import tempfile
import shutil

def get_reading_mode( options ):
    if options.is_bam:
        reading_mode = 'rb'
    else:
        ext = os.path.splitext( options.input_path )[-1]
        if ext == '.bam':
            reading_mode = 'rb'
        else:
            # here we assume a SAM file
            reading_mode = 'r'
    return reading_mode


def process_call_string( letter, key, CGmethHash, nonCGmethHash,  CHHmethHash,  CHGmethHash ):

    if letter.upper() == "Z":
        # if genomic base is CpG
        if not CGmethHash.has_key(key):
            CGmethHash[ key ] = [ 0, 0, 0 ]

        if letter == "Z":
            # update Cs
            CGmethHash[key][0] += 1
        elif letter == "z":
            # update Ts
            CGmethHash[key][1] += 1
        else:
            # update other bases
            CGmethHash[key][2] += 1
    else:
        #if genomic base is non-CpG
        if not nonCGmethHash.has_key(key):
            nonCGmethHash[ key ] = [ 0, 0, 0 ]

        if letter.upper() == "X" and not CHGmethHash.has_key(key):
            CHGmethHash[key] = [0,0,0]
        elif letter.upper() == "H" and not CHHmethHash.has_key(key):
            CHHmethHash[key] = [0,0,0]

        if letter == "X":
            nonCGmethHash[key][0] += 1
            CHGmethHash[key][0] += 1
        elif letter == "H":
            nonCGmethHash[key][0] += 1
            CHHmethHash[key][0] += 1
        elif letter == "x":
            nonCGmethHash[key][1] += 1
            CHGmethHash[key][1] += 1
        elif letter == "h":
            nonCGmethHash[key][1] += 1
            CHHmethHash[key][1] += 1
        else:
           # this condition will never be used
            nonCGmethHash[key][2] += 1
            if  letter.upper() == "X":
                CHGmethHash[key][2] += 1
            else:
                CHHmethHash[key][2] += 1
    return (CGmethHash, nonCGmethHash,  CHHmethHash,  CHGmethHash)


# process a given CG methlation hash
# writes the filter passing CGs to output file
def processCGmethHash( CGmethHash, out, min_cov ):
    for key in CGmethHash.keys():
        strand, chr, loc = key.split('|')
        noCs,noTs,noOs = CGmethHash[key]
        temp_sum = (noTs + noCs + noOs)
        Cperc = "%.2f" % ( 100.0 * noCs / temp_sum )
        Tperc = "%.2f" % ( 100.0 * noTs / temp_sum )
        if float( noTs + noCs ) / temp_sum > 0.9 and temp_sum >= min_cov:
            out.write( "\t".join( [chr+"."+loc, chr, loc, strand, str(temp_sum), Cperc, Tperc] ) + "\n")

def processCHmethHash( CGmethHash, out, min_cov ):
    for key in CGmethHash.keys():
        strand, chr, loc = key.split('|')
        noCs,noTs,noOs = CGmethHash[key]
        temp_sum = (noTs + noCs + noOs)
        Cperc = "%.2f" % ( 100.0 * noCs / temp_sum )
        Tperc = "%.2f" % ( 100.0 * noTs / temp_sum )
        if float( noTs + noCs ) / temp_sum > 0.9 and temp_sum >= min_cov:
            out.write( "\t".join( [chr+"."+loc, chr, loc, strand, str( temp_sum ), Cperc, Tperc] ) + "\n")


# process a given non CG methlation hash
# writes the filter passing Cs to a hash, that hash will be used to calculate conversion rate later on
def processnonCGmethHash( nonCGmethHash, CTconvArray, min_cov ):
    for key in nonCGmethHash.keys():
        strand, chr, loc = key.split('|')
        noCs,noTs,noOs = nonCGmethHash[key]
        temp_sum = noTs + noCs + noOs
        Cperc = "%.2f" % (100.0 * noCs / temp_sum )
        Tperc = "%.2f" % (100.0 * noTs / temp_sum )
        if float( noTs + noCs ) / temp_sum > 0.95 and temp_sum >= min_cov:
            CTconvArray[strand].append( (noTs * 100.0) / temp_sum )
    return CTconvArray


def process_sam(options, chromosome, temp_dir):

    min_cov = options.min_cov
    min_qual = options.min_qual

    # create temp files
    if options.CpG:
        if temp_dir != 'temp_dir':
            # multiprocessing mode
            out_temp = tempfile.NamedTemporaryFile(dir=temp_dir, prefix='CpG', delete=False)

    if options.CHH:
        if temp_dir != 'temp_dir':
            # multiprocessing mode
            CHH_out_temp = tempfile.NamedTemporaryFile(dir=temp_dir, prefix='CHH', delete=False)

    if options.CHG:
        if temp_dir != 'temp_dir':
            # multiprocessing mode
            CHG_out_temp =  tempfile.NamedTemporaryFile(dir=temp_dir, prefix='CHG', delete=False)


    reading_mode = get_reading_mode( options )
    samfile = pysam.Samfile( options.input_path, reading_mode )

    offset = 33
    if options.phred64:
        offset = 64

    nonCGmethHash = dict()
    pMeth_nonCG = defaultdict(list)
    CGmethHash = dict()
    CHHmethHash = dict()
    CHGmethHash = dict()

    start_pre = -1
    chr_pre = ''
    last_pos  =-1
    last_chrom = None

    # if the reading morde is 'r' and not 'rb', we have a sam file and multiprocessing is disabled
    # and we iterate over the whole genome
    if reading_mode == 'r' or options.processors == 1:
        samfile_iterator = samfile
    else:
        samfile_iterator = samfile.fetch(chromosome)

    # for every read in the sam file
    for iteration, read in enumerate(samfile_iterator):
        if iteration % 500000 == 0 and iteration != 0:
            # write out the collected results and clean the dictionary
            # otherwise the dictionary would be to large and slowing the process down dramatically
            if temp_dir != 'temp_dir':
                if options.CpG:
                    processCGmethHash( CGmethHash, out_temp, min_cov )
                if options.CHH:
                    processCHmethHash( CHHmethHash, CHH_out_temp, min_cov)
                if options.CHG:
                    processCHmethHash( CHGmethHash, CHG_out_temp, min_cov)
                CGmethHash = dict()
                CHHmethHash = dict()
                CHGmethHash = dict()

        start = read.pos + 1
        end = read.aend #read.rlen + start + 1 # or len(read.seq)+start+1
        chr = samfile.getrname( read.tid )
        methc = read.opt('XM')
        mcalls = methc
        quals = read.qual
        mrnm = read.mrnm
        mpos = read.mpos + 1
        isize = read.isize #read.tlen
        slen = read.rlen

        if read.opt('XR') == 'CT' and read.opt('XG') == 'CT':
            # original top strand
            strand = '+'
        elif read.opt('XR') == 'CT' and read.opt('XG') == 'GA':
            # original bottom strand
            strand = '-'
        elif read.opt('XR') == 'GA' and read.opt('XG') == 'CT':
            # complementary to original top strand, bismark says - strand to this
            strand = '+'
        elif read.opt('XR') == 'GA' and read.opt('XG') == 'GA':
            # complementary to original bottom strand, bismark says + strand to this
            strand = '-'

        # Check if the file is sorted
        if chr == chr_pre:
            if start_pre > start:
                sys.exit("The sam file is not sorted properly you can sort the file in unix-like machines using:\n grep -v '^[[:space:]]*\@' test.sam | sort -k3,3 -k4,4n  > test.sorted.sam \n")
        chr_pre = chr
        start_pre = start

        """
            if there is no_overlap, trim the mcalls and quals
            adjust the start
        """
        """
        if options.no_overlap and mrnm == 0 and paired:
            if (start + slen -1) > mpos:
                if (mpos - start)

        if($nolap && ( ($mrnm eq "=") && $paired ) ){

            if( ($start+$slen-1)>$mpos){
                if(($mpos-$start)>0)
                    { #{continue;}
                splice @mcalls,($mpos-$start);
                splice @quals,($mpos-$start);
            }
          }
        }
        """

        if (start - last_pos > 100 and last_pos != -1) or (chr != last_chrom and last_chrom != None):
            processnonCGmethHash( nonCGmethHash, pMeth_nonCG, min_cov )
            nonCGmethHash = dict()

        for index, letter in enumerate(quals):
            if ord(letter) - offset < min_qual or mcalls[index] == '.':
                continue
            if strand == '+':
                key = '|'.join( ["F",chr,str(start + index)] )
            else:
                key = '|'.join( ["R",chr,str(start + index)] )

            process_call_string(mcalls[index], key, CGmethHash, nonCGmethHash, CHHmethHash, CHGmethHash)

        last_pos = end
        last_chrom = chr

    samfile.close()
    processnonCGmethHash( nonCGmethHash, pMeth_nonCG, min_cov )

    if temp_dir != 'temp_dir':
        # in multiprocessing mode
        if options.CpG:
            processCGmethHash( CGmethHash, out_temp, min_cov )
        if options.CHH:
            processCHmethHash( CHHmethHash, CHH_out_temp, min_cov)
        if options.CHG:
            processCHmethHash( CHGmethHash, CHG_out_temp, min_cov)

        # create temp files
        if options.CpG:
            out_temp.close()

        if options.CHH:
            CHH_out_temp.close()

        if options.CHG:
            CHG_out_temp.close()

    return (CGmethHash, CHHmethHash, CHGmethHash, pMeth_nonCG)



def run_calc( args ):
    """
        Multiprocessing helper function.
        Getting Jobs out of the in_queue, calculate the percentage for that chromosome and returns the results.
    """
    temp_dir, options, chromosome = args
    return process_sam( options, chromosome, temp_dir )


def main( options ):
    """
        Read the files. First file is for now a special file. From that we only take the first gene and calculate
        that against the other files.
    """
    reading_mode = get_reading_mode( options )

    # check the file status produce flags 
    out = None
    CHH_out = None
    CHG_out = None
    pMeth_nonCG = defaultdict(list)
    min_cov = options.min_cov

    if options.CpG:
        out = open( options.CpG, 'w+')
        out.write( "chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT\n" )

    if options.CHH:
        CHH_out = open( options.CHH, 'w+' )
        CHH_out.write( "chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT\n" )

    if options.CHG:
        CHG_out =  open( options.CHG, 'w+' )
        CHG_out.write( "chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT\n" )

    # if the reading morde is 'r' and not 'rb', we have a sam file and multiprocessing is disabled.
    if reading_mode == 'rb' and options.processors > 1:
        multiproc = True
        print "Multiprocessing mode started with %s" % options.processors
        
        # creating temp dir and and store all temporary results there
        temp_dir = tempfile.mkdtemp()
        samfile = pysam.Samfile( options.input_path, reading_mode )
        # building a triple for each multiprocessing run -> (temp_dir, options, one chromosome)
        references = zip([temp_dir]*len(samfile.references), [options]*len(samfile.references), samfile.references)
        samfile.close()

        p = multiprocessing.Pool( options.processors )
        results = p.map_async(run_calc, references)
        results_iterator = results.get()
    else:
        multiproc = False
        print 'Single Mode activated. You can use multiple processors if you use a indexed BAM file.'
        results_iterator = [ run_calc( ('temp_dir', options, 'all_chromosomes_at_once') ) ]


    for result in results_iterator:
        (CGmethHash, CHHmethHash, CHGmethHash, pMeth_nonCG_temp) = result
        pMeth_nonCG["F"] = pMeth_nonCG.get( "F", []) + pMeth_nonCG_temp.get( "F", [])
        pMeth_nonCG["R"] = pMeth_nonCG.get( "R", []) + pMeth_nonCG_temp.get( "R", [])

        # if not in multiprocessing mode, than write down all collected results
        if not multiproc:
            if options.CpG:
                processCGmethHash( CGmethHash, out, min_cov )
            if options.CHH:
                processCHmethHash( CHHmethHash, CHH_out, min_cov)
            if options.CHG:
                processCHmethHash( CHGmethHash, CHG_out, min_cov)

    # if in multiprocessing mode, than collect all temporary files and write them down to the output file
    if options.CpG and multiproc:
        for filename in os.listdir(temp_dir):
            if filename.startswith('CpG'):
                shutil.copyfileobj( open(os.path.join(temp_dir, filename)), out )
    if options.CHH and multiproc:
        for filename in os.listdir(temp_dir):
            if filename.startswith('CHH'):
                shutil.copyfileobj( open(os.path.join(temp_dir, filename)), CHH_out )
    if options.CHG and multiproc:
        for filename in os.listdir(temp_dir):
            if filename.startswith('CHG'):
                shutil.copyfileobj( open(os.path.join(temp_dir, filename)), CHG_out )

    # cleaning temporary working directory
    if multiproc:
        if options.CHG:
            CHG_out.close()
        if options.CHH:
            CHH_out.close()
        if options.CpG:
            out.close()
        shutil.rmtree( temp_dir )

    # get the conversion rate and write it out!!
    numF = len( pMeth_nonCG.get( "F", []) )
    numR = len( pMeth_nonCG.get( "R", []) )

    if numF == 0 and numR == 0:
        sys.exit("\nnot enough alignments that pass coverage and phred score thresholds to calculate conversion rates\n EXITING....\n")

    both_strands = pMeth_nonCG.get( "F", []) + pMeth_nonCG.get( "R", [])

    medFconvRate = 0
    medRconvRate = 0
    AvFconvRate = float( sum( pMeth_nonCG["F"] ) ) / numF
    AvRconvRate = float( sum( pMeth_nonCG["R"] ) ) / numR
    AvconvRate = float( sum( both_strands ) ) / (numF + numR)
    print sum( both_strands ), numF, numR

    if numF > 0:
        medFconvRate = numpy.median( pMeth_nonCG["F"] )

    if numR > 0:
        medRconvRate = numpy.median( pMeth_nonCG[ "R" ] )

    medconvRate = numpy.median( both_strands )

    totCpG = len( both_strands )

    res = "total otherC considered (>95%% C+T): %s\n" % totCpG
    res += "average conversion rate = %s\n" % AvconvRate
    res += "median conversion rate = %s\n\n" % medconvRate 

    res += "total otherC considered (Forward) (>95%% C+T): %s\n" % numF
    res += "average conversion rate (Forward) = %s\n" % AvFconvRate
    res += "median conversion rate (Forward) = %s\n\n" % medFconvRate

    res += "total otherC considered (Reverse) (>95%% C+T): %s\n" % numR
    res += "average conversion rate (Reverse) = %s\n" % AvRconvRate
    res += "median conversion rate (Reverse) = %s\n" % medRconvRate
    print res



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Base methylation calling from Bismark SAM files.')

    parser.add_argument("-i", "--input", dest="input_path",
                    required=True,
                    help="Path to the input file.")

    parser.add_argument("--CpG", dest="CpG",
                    help="output filename for CpG methylation scores (if not specified no file is written out)")

    parser.add_argument("--CHH", dest="CHH",
                    help="output filename for CHH methylation scores (if not specified no file is written out)")

    parser.add_argument("--CHG", dest="CHG",
                    help="output filename for CHG methylation scores (if not specified no file is written out)")

    parser.add_argument("--phred64", dest="phred64", action="store_true", default=False,
                    help="quality scores phred64 scale used otherwise phred33 is the default")

    parser.add_argument("--mincov", dest="min_cov", default=10, type=int,
                    help="min coverage (default:10)")

    parser.add_argument("--minqual", dest="min_qual", default=20, type=int,
                    help="minquality   (default:20)")

    parser.add_argument("--bam", dest="is_bam", action="store_true", default=False,
                    help="If specified the input file is in BAM format, otherwise the type is guessed from the filename extension.")

    parser.add_argument('-p', '--processors', type=int, 
        default=multiprocessing.cpu_count())

    options = parser.parse_args()

    main(options)
