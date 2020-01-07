#!python
import re
import os
import numpy as np
from argparse import ArgumentParser
import sys


#import self defined functions and classes

from utils import readFastaFile
from utils import readGff
from utils import readSam
from utils import getSubSequence
from utils import SimilarFragment
from operator import attrgetter

# bs674@cornell.edu
# Baoxing Song
# 18 June, 2019, updated

# set a minimum output length

if __name__ == '__main__':
    parser = ArgumentParser(description='prepare sequence for and_CNS software')
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

    parser.add_argument("-c", "--refcds",
                        dest="refcds",
                        type=str,
                        default="",
                        help="the reference CDS sequence file in fasta format")

    parser.add_argument("-m", "--similarity",
                        dest="similarity",
                        type=float,
                        default=0.6,
                        help="the similarity of alignment to define as a gene, (default:0.6 )")

    parser.add_argument("-s", "--sam",
                        dest="sam",
                        type=str,
                        default="",
                        help="input sam files")

    parser.add_argument("-d", "--cisdistance",
                        dest="cisdistance",
                        type=int,
                        default=100000,
                        help="maximum cis distance for the contig edge region (default: 100000)")

    parser.add_argument("-x", "--maxOutPutLength",
                        dest="maxOutPutLength",
                        type=int,
                        default=150000,
                        help="maximum output sequence length (default: 150000)")

    parser.add_argument("-l", "--minimumOutPutLength",
                        dest="minimumOutPutLength",
                        type=int,
                        default=20,
                        help="minimum output sequence length (default: 20)")

    parser.add_argument("-q", "--query",
                        dest="query",
                        type=str,
                        default="",
                        help="query genome sequence file")

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

    if args.refcds == "":
        print("Error: please specify --refcds", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    if args.sam == "":
        print("Error: please specify --sam", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    if args.query == "":
        print("Error: please specify --query", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    if args.output == "":
        print("Error: please specify --output", file=sys.stderr)
        parser.print_help()
        sys.exit(1)


    ref_genome = readFastaFile(args.ref_geno)
    fastas = readFastaFile(args.refcds)
    chromosome_gene_list1 = readGff(args.ref_gff)

    chromosome_gene_list2, gene_name_dict = readSam(args.sam, fastas, args.similarity )
    query_genome = readFastaFile(args.query)

    try:
        os.rmdir(args.output)
    except OSError:
        1 == 1
    os.mkdir(args.output)


    #get all the similar fragments
    similarFragments_dict = dict()
    for chr in chromosome_gene_list1:
        similarFragments = np.empty([0, 1], SimilarFragment)
        for g in chromosome_gene_list1[chr]:
            if g.name in gene_name_dict:
                for g2 in gene_name_dict[g.name]:
                    ref_chr = chr
                    query_chr = g2.chr

                    ref_strand = g.strand
                    query_strand = g2.strand
                    if ref_strand == query_strand:
                        query_strand = "+"
                    else:
                        query_strand = "-"

                    ref_strand = "+"

                    if query_strand == "+":
                        ref_start = g.start-args.cisdistance
                        ref_end = g.start - 1
                        if ref_start < 1:
                            ref_start=1
                        if ref_end > len(ref_genome[g.chr]):
                            ref_end = len(ref_genome[g.chr])

                        query_start = g2.start - args.cisdistance
                        query_end = g2.start - 1
                        if query_start < 1:
                            query_start = 1
                        if query_end > len(query_genome[g2.chr]):
                            query_end = len(query_genome[g2.chr])

                        similarFragment = SimilarFragment (ref_chr, ref_strand, ref_start, ref_end, query_chr, query_strand, query_start, query_end, g.name+"up")
                        similarFragments = np.append(similarFragments, [similarFragment])

                        ref_start = g.start
                        ref_end = g.end
                        if ref_start < 1:
                            ref_start=1
                        if ref_end > len(ref_genome[g.chr]):
                            ref_end = len(ref_genome[g.chr])

                        query_start = g2.start
                        query_end = g2.end
                        if query_start < 1:
                            query_start = 1
                        if query_end > len(query_genome[g2.chr]):
                            query_end = len(query_genome[g2.chr])

                        similarFragment = SimilarFragment (ref_chr, ref_strand, ref_start, ref_end, query_chr, query_strand, query_start, query_end, g.name+"in")
                        similarFragments = np.append(similarFragments, [similarFragment])


                        ref_start = g.end + 1
                        ref_end = g.end + args.cisdistance
                        if ref_start < 1:
                            ref_start=1
                        if ref_end > len(ref_genome[g.chr]):
                            ref_end = len(ref_genome[g.chr])

                        query_start = g2.end + 1
                        query_end = g2.end + args.cisdistance
                        if query_start < 1:
                            query_start = 1
                        if query_end > len(query_genome[g2.chr]):
                            query_end = len(query_genome[g2.chr])

                        similarFragment = SimilarFragment (ref_chr, ref_strand, ref_start, ref_end, query_chr, query_strand, query_start, query_end, g.name+"do")
                        similarFragments = np.append(similarFragments, [similarFragment])
                    else:
                        ref_start = g.start-args.cisdistance
                        ref_end = g.start - 1
                        if ref_start < 1:
                            ref_start=1
                        if ref_end > len(ref_genome[g.chr]):
                            ref_end = len(ref_genome[g.chr])

                        query_start = g2.end + 1
                        query_end = g2.end + args.cisdistance
                        if query_start < 1:
                            query_start = 1
                        if query_end > len(query_genome[g2.chr]):
                            query_end = len(query_genome[g2.chr])

                        similarFragment = SimilarFragment (ref_chr, ref_strand, ref_start, ref_end, query_chr, query_strand, query_start, query_end, g.name+"up")
                        similarFragments = np.append(similarFragments, [similarFragment])

                        ref_start = g.start
                        ref_end = g.end
                        if ref_start < 1:
                            ref_start=1
                        if ref_end > len(ref_genome[g.chr]):
                            ref_end = len(ref_genome[g.chr])

                        query_start = g2.start
                        query_end = g2.end
                        if query_start < 1:
                            query_start = 1
                        if query_end > len(query_genome[g2.chr]):
                            query_end = len(query_genome[g2.chr])

                        similarFragment = SimilarFragment (ref_chr, ref_strand, ref_start, ref_end, query_chr, query_strand, query_start, query_end, g.name+"in")
                        similarFragments = np.append(similarFragments, [similarFragment])

                        ref_start = g.end + 1
                        ref_end = g.end + args.cisdistance
                        if ref_start < 1:
                            ref_start=1
                        if ref_end > len(ref_genome[g.chr]):
                            ref_end = len(ref_genome[g.chr])

                        query_start = g2.start - args.cisdistance
                        query_end = g2.start - 1
                        if query_start < 1:
                            query_start = 1
                        if query_end > len(query_genome[g2.chr]):
                            query_end = len(query_genome[g2.chr])

                        similarFragment = SimilarFragment (ref_chr, ref_strand, ref_start, ref_end, query_chr, query_strand, query_start, query_end, g.name+"do")
                        similarFragments = np.append(similarFragments, [similarFragment])

        similarFragments_dict[chr] = similarFragments


    for chr in chromosome_gene_list1:
        similarFragments = similarFragments_dict[chr]
        # merge similar fragments
        similarFragments = sorted(similarFragments, key=attrgetter('ref_start', 'ref_end'))
        haveChange = True
        while haveChange:
            haveChange = False
            usedList = dict()
            newSimilarFragments = np.empty([0, 1], SimilarFragment)
            for i in range(0, len(similarFragments)):
                if i not in usedList:
                    for j in range(i+1, len(similarFragments)):
                        if (similarFragments[j].ref_end - similarFragments[i].ref_start) < args.maxOutPutLength:
                            if (j not in usedList) and similarFragments[j].query_chr == similarFragments[i].query_chr:
                                if similarFragments[j].ref_start < (similarFragments[i].ref_end): #if they overlap with each other by 1 bp, do not combine
                                    if(similarFragments[i].query_strand==similarFragments[j].query_strand and "+" == similarFragments[i].query_strand):
                                        if similarFragments[j].query_start >= similarFragments[i].query_start and similarFragments[j].query_start < (similarFragments[i].query_end):
                                            if similarFragments[j].ref_end > similarFragments[i].ref_end:
                                                similarFragments[i].ref_end = similarFragments[j].ref_end
                                            if similarFragments[j].query_end > similarFragments[i].query_end:
                                                similarFragments[i].query_end = similarFragments[j].query_end
                                            usedList[j] = 1
                                            haveChange = True
                                    elif (similarFragments[i].query_strand==similarFragments[j].query_strand and "-" == similarFragments[i].query_strand): #here is difficult to understant
                                        if ((similarFragments[j].query_end) > similarFragments[i].query_start and similarFragments[j].query_end <= similarFragments[i].query_end): #if they overlap with each other by 1 bp, do not combine
                                            if similarFragments[j].ref_end > similarFragments[i].ref_end:
                                                similarFragments[i].ref_end = similarFragments[j].ref_end
                                            if similarFragments[j].query_start < similarFragments[i].query_start:
                                                similarFragments[i].query_satrt = similarFragments[j].query_start
                                            usedList[j] = 1
                                            haveChange = True
                                else:
                                    break
                        else:
                            break
                    newSimilarFragments = np.append(newSimilarFragments, [similarFragments[i]])

            similarFragments = newSimilarFragments
            similarFragments = sorted(similarFragments, key=attrgetter('ref_start', 'ref_end'))
        similarFragments_dict[chr] = similarFragments

    minimumOutPutLength = args.minimumOutPutLength - 1
    #output
    index = 1
    for chr in similarFragments_dict:
        for similarFragment in similarFragments_dict[chr]:
            if (similarFragment.ref_end-similarFragment.ref_start) >= minimumOutPutLength and (similarFragment.query_end-similarFragment.query_start) >= minimumOutPutLength:
                strand = 1
                if similarFragment.ref_strand == "-":
                    strand = 0
                output = open(args.output+"/" + str(index), 'w')
                sequence = getSubSequence(ref_genome, similarFragment.ref_chr, similarFragment.ref_start, similarFragment.ref_end, similarFragment.ref_strand)
                output.write(">reference" + "\t" + similarFragment.name + "\tspecies:reference"+"\tchr:"+similarFragment.ref_chr+"\tstrand:"+ str(strand) +"\tstart:"+str(similarFragment.ref_start)+"\tend:"+str(similarFragment.ref_end)+"\n"+sequence+"\n")

                if similarFragment.query_strand == "-":
                    strand = 0
                else:
                    strand = 1
                sequence = getSubSequence(query_genome, similarFragment.query_chr, similarFragment.query_start, similarFragment.query_end, similarFragment.query_strand)
                output.write(">query" + "\tspecies:query"+"\tchr:"+similarFragment.query_chr+"\tstrand:"+ str(strand) +"\tstart:"+str(similarFragment.query_start) +"\tend:"+str(similarFragment.query_end)+"\n"+sequence+"\n")

                output.close()
                index = index + 1

    sys.exit(0)
