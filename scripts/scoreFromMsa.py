#!python
import re
import os
import numpy as np
from argparse import ArgumentParser
import sys

# bs674@cornell.edu
# Baoxing Song
# 26 April, 2019, updated


from utils import readFastaFile

if __name__ == '__main__':
    parser = ArgumentParser(description='score each basepair from MSA file')
    parser.add_argument("-s", "--seq",
                        dest="seq",
                        type=str,
                        default="",
                        help=" multiple sequence alignment file")
    # parser.add_argument("-m", "--miniscore",
    #                     dest="miniscore",
    #                     type=float,
    #                     default=0.7,
    #                     help="minimum score to take as a CNS site (default: 0.7)")
    # parser.add_argument("-d", "--maxdis",
    #                     dest="maxdis",
    #                     type=int,
    #                     default=3,
    #                     help="maximum distance to merger two adjacent consevered fragments (default: 3)")
    # parser.add_argument("-l", "--minilen",
    #                     dest="minilen",
    #                     type=int,
    #                     default=6,
    #                     help="minimum CNS fragment length to output (default: 6)")
    parser.add_argument("-o", "--output",
                        dest="output",
                        type=str,
                        default="",
                        help="output file")

    args = parser.parse_args()
    if args.seq == "":
        print("Error: please specify --seq", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    if args.output == "":
        print("Error: please specify --output", file=sys.stderr)
        parser.print_help()
        sys.exit(1)
    seqs = readFastaFile(args.seq)
    seqLength = 0
    numSeq = 0
    for seqEntry in seqs:
        seqLength = len(seqs[seqEntry])
        numSeq = numSeq + 1
    for i in range(0, seqLength):
        frequence = dict()
        for seqEntry in seqs:
            if seqs[seqEntry][i] != '-':
                if seqs[seqEntry][i] in frequence:
                    frequence[seqs[seqEntry][i]] = frequence[seqs[seqEntry][i]] + 1
                else:
                    frequence[seqs[seqEntry][i]] = 1
        maxFrq = 0
        maxChar = '-'
        for k in frequence:
            if frequence[k] > maxFrq:
                maxFrq = frequence[k]
                maxChar = k
        score = maxFrq/numSeq
        print(str(score), end='')
    print()
    sys.exit(0)

