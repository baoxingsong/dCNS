#!python
import re
import numpy as np
from argparse import ArgumentParser
import sys
import random
# bs674@cornell.edu
# Baoxing Song
# 11 Aug, 2019, updated

from utils import readFastaFile
from utils import getSubSequence
from utils import allThemRnaInterval
from utils import Interval

def gcContent(seq):
    gclength = 0
    for i in range(0, len(seq)):
        if seq[i] == 'g' or seq[i] == 'c' or seq[i] == 'G' or seq[i] == 'C':
            gclength = gclength + 1
    return gclength/len(seq)

if __name__ == '__main__':
    parser = ArgumentParser(description='extract CNS sequence for ML purpose')

    parser.add_argument("-g", "--ref_gff",
                        dest="ref_gff",
                        type=str,
                        default="",
                        help="gff file of reference genome")

    parser.add_argument("-r", "--ref_geno",
                        dest="ref_geno",
                        type=str,
                        default="",
                        help="reference genome sequence file in fasta format")

    parser.add_argument("-l", "--length",
                        dest="seq_length",
                        type=int,
                        default=200,
                        help="output sequence fragment length")

    parser.add_argument("-s", "--sam",
                        dest="sam",
                        type=str,
                        default="",
                        help="input sam files")

    parser.add_argument("-o", "--output",
                        dest="output",
                        type=str,
                        default="",
                        help="output file")

    args = parser.parse_args()

    if args.ref_gff == "":
        print("Error: please specify --ref_gff", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    if args.ref_geno == "":
        print("Error: please specify --ref_geno", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    if args.sam == "":
        print("Error: please specify --sam", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    if args.output == "":
        print("Error: please specify --output", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    ref_genome = readFastaFile(args.ref_geno)  # maize genome
    print("genome reading done")
    chromosome_gene_list_ref = allThemRnaInterval(args.ref_gff) #maize gff file
    print("reference gff reading done")

    ref_black_list_interval = dict()

    outputsequences = dict()

    with open( args.sam) as f:
        for line in f:
            # print(line)
            m = re.search('^@', line)
            if( m == None ):
                m = re.search('^(.*?)\s+(\d+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\w+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\w+)\s+', line)
                if( m != None and m.group(3) != "*" ):
                    chromosome_name_ref = m.group(3)

                    start = int(m.group(4))
                    end = int(m.group(4))
                    ms = re.findall('(\d+)[MDN=X]', m.group(6))
                    for mm in ms:
                        end += int(mm)
                    if( end > start ):
                        length = end - start
                        end -= 1
                        if length < args.seq_length :
                            lacking_length = args.seq_length - length
                            start_lacking = int(lacking_length/2)
                            if start_lacking >= start:
                                start_lacking = start - 1
                            end_lacking = lacking_length - start_lacking
                            start = start - start_lacking
                            end = end + end_lacking

                        everoverlap = False
                        if chromosome_name_ref in chromosome_gene_list_ref:
                            for g in chromosome_gene_list_ref[chromosome_name_ref]:
                                if g.start <  end and g.start > start:
                                    everoverlap = True
                                if start <  g.end and start > g.start:
                                    everoverlap = True

                        if (not chromosome_name_ref in ref_black_list_interval):
                            ref_black_list_interval[chromosome_name_ref] = np.empty([0, 1], Interval)
                        gRef = Interval (start, end)
                        ref_black_list_interval[chromosome_name_ref] = np.append(ref_black_list_interval[chromosome_name_ref], [gRef])
                        if (not everoverlap):
                            seq = getSubSequence(ref_genome, chromosome_name_ref, start, end)
                            if len(seq) > args.seq_length:
                                distance = len(seq) - args.seq_length
                                s_distance = int(distance/2)
                                seq = seq[s_distance:(s_distance+args.seq_length)]
                            outputString = seq
                            if outputString not in outputsequences:
                                if (not 'n' in outputString) and (not 'N' in outputString):
                                    outputsequences[outputString] = 1

    print ("sam file reading done")
    ref_chrs = []
    for chr in chromosome_gene_list_ref:
        ref_chrs.append(chr)
        for g in chromosome_gene_list_ref[chr]:
            if (not chr in ref_black_list_interval):
                ref_black_list_interval[chr] = np.empty([0, 1], Interval)
            gRef = Interval (g.start, g.end)
            ref_black_list_interval[chr] = np.append(ref_black_list_interval[chr], [gRef])
    print ("black list region setting done")
    random.seed( 30 )

    output = open(args.output, 'w')
    for o in outputsequences:
        ogc = gcContent(o)
        output.write(o + "\t")

        refNotOutPuted = True
        ref_seq = ""
        while refNotOutPuted:
            ref_chr = ref_chrs[random.randint(0, len(ref_chrs)-1)]
            if len(ref_genome[ref_chr]) > len(o):
                ref_start = random.randint(1, len(ref_genome[ref_chr])-args.seq_length+1 )
                ref_end = ref_start + args.seq_length - 1
                everoverlap = False
                if ref_chr in ref_black_list_interval:
                    for g in ref_black_list_interval[ref_chr]:
                        if g.start <  ref_end and g.start > ref_start:
                            everoverlap = True
                        if ref_start <  g.end and ref_start > g.start:
                            everoverlap = True
                ref_seq = getSubSequence(ref_genome, ref_chr, ref_start, ref_end)
                if (not everoverlap) and ('n' not in ref_seq) and ('N' not in ref_seq):
                    if ( abs(gcContent(ref_seq) - ogc)<0.15 ):
                        refNotOutPuted = False

        output.write(ref_seq + "\n")

    output.close()
    sys.exit(0)
