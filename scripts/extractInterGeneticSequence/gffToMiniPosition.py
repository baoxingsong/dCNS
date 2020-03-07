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
# 26 Jan, 2020, updated

# set a minimum output length
# the update is just remove the merging function

if __name__ == '__main__':
    parser = ArgumentParser(description='minimum flanking distance for CNS analysis')
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

    args = parser.parse_args()

    if args.ref_gff == "":
        print("Error: please specify --ref_gff", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    if args.ref_geno == "":
        print("Error: please specify --ref_geno", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    ref_genome = readFastaFile(args.ref_geno)
    chromosome_gene_list1 = readGff(args.ref_gff)


    for chr in chromosome_gene_list1:
        for g in chromosome_gene_list1[chr]:
            miniDistance = min(g.start, len(ref_genome[chr]) - g.end)
            for g2 in chromosome_gene_list1[chr]:
                if g.name != g2.name:
                    miniDistance = min(miniDistance, abs(g2.start-g.start))
                    miniDistance = min(miniDistance, abs(g2.start-g.end))
                    miniDistance = min(miniDistance, abs(g2.end-g.start))
                    miniDistance = min(miniDistance, abs(g2.end-g.end))
            print(g.name + "\t" + str(miniDistance))
    sys.exit(0)
