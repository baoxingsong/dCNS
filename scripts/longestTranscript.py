#!python
import re

import NucleotideCodeSubstitution
import FastaFile
import GffFile
import sys
import numpy as np
from argparse import ArgumentParser
from  MyUtil import str2bool2
from  MyUtil import str2bool

# bs674@cornell.edu
def read_data(gffFile, fastaFile):
    gene_dict = dict()
    chromosome_gene_dict, chromosome_gene_list = GffFile.readGff( gffFile )
    chromosome_names, fastas = FastaFile.readFastaFile( fastaFile )
    GffFile.update_sequence_information(fastas, chromosome_gene_dict)
    #keep only one transcript for each gene
    #delete gene without transcript
    for chromosome_name in chromosome_names:
        if chromosome_name in chromosome_gene_dict:
            gene_dict[chromosome_name] = dict()
            for gene_name in chromosome_gene_dict[chromosome_name]:
                ##delete transcript so that the there is only one transcript left for each gene
                if len(chromosome_gene_dict[chromosome_name][gene_name].transcripts)>0:
                    longest = 0
                    index = 0
                    for i in range(0, len(chromosome_gene_dict[chromosome_name][gene_name].transcripts)):
                        if len(chromosome_gene_dict[chromosome_name][gene_name].transcripts[i].cds_sequence) > longest:
                            longest = len(chromosome_gene_dict[chromosome_name][gene_name].transcripts[i].cds_sequence)
                            index = i
                    g = GffFile.Gene(chromosome_gene_dict[chromosome_name][gene_name].name, chromosome_gene_dict[chromosome_name][gene_name].strand)
                    g.add_transcript(chromosome_gene_dict[chromosome_name][gene_name].transcripts[index])
                    gene_dict[chromosome_name][gene_name] = g
    return gene_dict
# print ("begin to run")


if __name__ == '__main__':
    parser = ArgumentParser(description='get the primary transcript sequence for each gene')
    parser.add_argument("-g", "--genome",
                        dest="genome_file",
                        type=str,
                        default="",
                        help="genome sequence file in fasta format")
    parser.add_argument("-f", "--gff",
                        dest="gff",
                        type=str,
                        default="",
                        help="gff file")
    parser.add_argument("-t", "--transcript_name",
                        dest="transcript_name",
                        type=str2bool,
                        default=True,
                        help="output transcript name (default True)")
    parser.add_argument("-n", "--output_genome_seq",
                        dest="output_genome_seq",
                        type=str2bool2,
                        default=False,
                        help="output genome Sequence (default False, which output CDS sequence)")

    parser.add_argument("-o", "--output",
                        dest="output",
                        type=str,
                        default="",
                        help="output file")
    args = parser.parse_args()
    if args.genome_file == "":
        print("Error: please specify --genome", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    if args.gff == "":
        print("Error: please specify --gff", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    if args.output == "":
        print("Error: please specify --output", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    output = open(args.output, 'w')
    chromosome_gene_dict = read_data(args.gff, args.genome_file)
    for chromosome_name in chromosome_gene_dict:
        for gene_name in chromosome_gene_dict[chromosome_name]:
            for transcript in chromosome_gene_dict[chromosome_name][gene_name].transcripts:
                if args.transcript_name:
                    output.write(">" + transcript.name.replace("transcript:", "")+"\n")
                else:
                    output.write(">" + gene_name.replace("gene:", "") + "\n")
                if args.output_genome_seq:
                    output.write(transcript.genome_sequence + "\n")
                else:
                    output.write(transcript.cds_sequence+"\n")
    output.close()
