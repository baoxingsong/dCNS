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
if __name__ == '__main__':
    parser = ArgumentParser(description='extract CNS sequence for ML purpose')

    parser.add_argument("-g", "--ref_gff",
                        dest="ref_gff",
                        type=str,
                        default="",
                        help="gff file of reference genome")


    parser.add_argument("-f", "--query_gff",
                        dest="query_gff",
                        type=str,
                        default="",
                        help="gff file of query genome")

    parser.add_argument("-r", "--ref_geno",
                        dest="ref_geno",
                        type=str,
                        default="",
                        help="reference genome sequence file in fasta format")

    parser.add_argument("-q", "--que_geno",
                        dest="que_geno",
                        type=str,
                        default="",
                        help="query genome sequence file in fasta format")

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

    if args.query_gff == "":
        print("Error: please specify --query_gff", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    if args.ref_geno == "":
        print("Error: please specify --ref_geno", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    if args.que_geno == "":
        print("Error: please specify --que_geno", file=sys.stderr)
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

    query_genome = readFastaFile(args.que_geno)  # sorghum genome
    ref_genome = readFastaFile(args.ref_geno)  # maize genome
    chromosome_gene_list_ref = allThemRnaInterval(args.ref_gff) #maize gff file
    chromosome_gene_list_query = allThemRnaInterval(args.query_gff) #query gff file



    ref_black_list_interval = dict()
    query_black_list_interval = dict()  # map<chr, vector<gene>>

    outputsequences = dict()

    with open( args.sam) as f:
        for line in f:
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
                        end -= 1

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


                        ms2 = re.search('^(\d+)[SH]', m.group(6))
                        if ms2 != None:
                            query_start = int(ms2.group(1))
                            query_end = query_start + len(m.group(10)) - 1
                            query_chr = m.group(1)
                            everoverlap_query = False
                            if query_chr in chromosome_gene_list_query:
                                for g in chromosome_gene_list_query[query_chr]:
                                    if g.start <  query_end and g.start > query_start:
                                        everoverlap_query = True
                                    if query_start <  g.end and query_start > g.start:
                                        everoverlap_query = True

                            if (not everoverlap) and (not everoverlap_query):
                                seq = getSubSequence(ref_genome, chromosome_name_ref, start, end)
                                outputString = seq + "\t" + m.group(10)
                                if outputString not in outputsequences:
                                    if (not 'n' in outputString) and (not 'N' in outputString):
                                        outputsequences[outputString] = 1

                            if (not m.group(1) in query_black_list_interval):
                                query_black_list_interval[m.group(1)] = np.empty([0, 1], Interval)
                            gQue = Interval (query_start, query_end)
                            query_black_list_interval[ m.group(1)] = np.append(query_black_list_interval[ query_chr], [gQue])

    ref_chrs = []
    for chr in chromosome_gene_list_ref:
        ref_chrs.append(chr)
        for g in chromosome_gene_list_ref[chr]:
            if (not chr in ref_black_list_interval):
                ref_black_list_interval[chr] = np.empty([0, 1], Interval)
            gRef = Interval (g.start, g.end)
            ref_black_list_interval[chr] = np.append(ref_black_list_interval[chr], [gRef])

    query_chrs = []
    for chr in chromosome_gene_list_query:
        query_chrs.append(chr)
        for g in chromosome_gene_list_query[chr]:
            if (not chr in query_black_list_interval):
                query_black_list_interval[chr] = np.empty([0, 1], Interval)
            gQue = Interval (g.start, g.end)
            query_black_list_interval[ chr] = np.append(query_black_list_interval[chr], [gQue])

    random.seed( 30 )

    output1 = open(args.output+".positive", 'w')
    output2 = open(args.output+".negative", 'w')
    for o in outputsequences:
        output1.write(o + "\n")


        output = o.split("\t")
        referLength = len(output[0])
        queryLength = len(output[1])

        refNotOutPuted = True
        ref_seq = ""
        while refNotOutPuted:
            ref_chr = ref_chrs[random.randint(0, len(ref_chrs)-1)]
            if len(ref_genome[ref_chr]) > referLength:
                ref_start = random.randint(1, len(ref_genome[ref_chr])-referLength+1 )
                ref_end = ref_start + referLength - 1
                everoverlap = False
                if ref_chr in ref_black_list_interval:
                    for g in ref_black_list_interval[ref_chr]:
                        if g.start <  ref_end and g.start > ref_start:
                            everoverlap = True
                        if ref_start <  g.end and ref_start > g.start:
                            everoverlap = True
                ref_seq = getSubSequence(ref_genome, ref_chr, ref_start, ref_end)
                if (not everoverlap) and ('n' not in ref_seq) and ('N' not in ref_seq):
                    refNotOutPuted = False

        queryNotOutPuted = True
        query_seq = ""
        while queryNotOutPuted:
            query_chr = query_chrs[random.randint(0, len(query_chrs)-1)]
            if len(query_genome[query_chr]) > queryLength:
                query_start = random.randint(1, len(query_genome[query_chr])-queryLength+1 )
                query_end = query_start + queryLength - 1
                everoverlap = False
                if query_chr in query_black_list_interval:
                    for g in query_black_list_interval[query_chr]:
                        if g.start <  query_end and g.start > query_start:
                            everoverlap = True
                        if query_start <  g.end and query_start > g.start:
                            everoverlap = True
                query_seq = getSubSequence(query_genome, query_chr, query_start, query_end)
                if (not everoverlap) and ('n' not in query_seq) and ('N' not in query_seq):
                    queryNotOutPuted = False
        output2.write(ref_seq + "\t" + query_seq + "\n")

    output1.close()
    output2.close()
    sys.exit(0)
