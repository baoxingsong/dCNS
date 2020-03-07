#!python
import numpy as np

import NucleotideCodeSubstitution
import FastaFile
import GffFile

# bs674@cornell.edu


def overlap_with_certain_gene(position, chromosome_name, strand, chromosome_gene_dict):
    for gene_name in chromosome_gene_dict[chromosome_name]:
        if chromosome_gene_dict[chromosome_name][gene_name].start > position:
            return None
        if chromosome_gene_dict[chromosome_name][gene_name].end < position:
            continue
        if (chromosome_gene_dict[chromosome_name][gene_name].start < position) & (chromosome_gene_dict[chromosome_name][gene_name].end > position) & (chromosome_gene_dict[chromosome_name][gene_name].strand == strand):
            return gene_name
    return None

def get_genetic_region_states(chromosome_name, strand, chromosome_gene_dict, fastas):
    states = []
    for gene_name in chromosome_gene_dict[chromosome_name]:
        # positive strand
        if (strand == "+") & (strand == chromosome_gene_dict[chromosome_name][gene_name].strand):
            start =  chromosome_gene_dict[chromosome_name][gene_name].transcripts[0].start
            while len(states) < start-1:
                states.append("I") # intergenetic
            states.append("S")  # start codon
            states.append("S")  # start codon
            states.append("S")  # start codon
            while len(states) < (chromosome_gene_dict[chromosome_name][gene_name].transcripts[0].Cds[0][1]):
                states.append("C") # cds this is the first cds
            cds_size = len(chromosome_gene_dict[chromosome_name][gene_name].transcripts[0].Cds)
            if cds_size > 1:
                for i in range(1, cds_size):
                    # add donor splice sites information between the last cds and this cds
                    if (cds_size >1) & (i != cds_size):
                        states.append("D")  # dornor splice sites
                        states.append("D")  # dornor splice sites
                    cds = chromosome_gene_dict[chromosome_name][gene_name].transcripts[0].Cds[i]
                    while len(states) < (cds[0]-1):
                        states.append("N") # intron
                    states[len(states) - 1] = "A"  # acceptor splice sites
                    states[len(states) - 2] = "A"  # acceptor splice sites
                    while len(states) < (cds[1]):
                        states.append("C")  # cds
            states[len(states) - 1] = "T"  # stop codon
            states[len(states) - 2] = "T"  # stop codon
            states[len(states) - 3] = "T"  # stop codon
        # negative strand
        if (strand == "-") & (strand == chromosome_gene_dict[chromosome_name][gene_name].strand):
            start = chromosome_gene_dict[chromosome_name][gene_name].transcripts[0].start
            while len(states) < start-1:
                states.append("I") # intergenetic
            states.append("T")  # stop codon
            states.append("T")  # stop codon
            states.append("T")  # stop codon
            cds_size = len(chromosome_gene_dict[chromosome_name][gene_name].transcripts[0].Cds)
            while len(states) < (chromosome_gene_dict[chromosome_name][gene_name].transcripts[0].Cds[cds_size - 1][1]):
                states.append("C") # cds this is the first cds
            if cds_size > 1:
                for i in reversed(range(0, cds_size - 1)):
                    states.append("A")  # acceptor splice sites
                    states.append("A")  # acceptor splice sites
                    cds = chromosome_gene_dict[chromosome_name][gene_name].transcripts[0].Cds[i]
                    while len(states) < (cds[0] - 1):
                        states.append("N")  # intron
                    states[len(states) - 1] = "D"  # dornor splice sites
                    states[len(states) - 2] = "D"  # dornor splice sites
                    while len(states) < (cds[1]):
                        states.append("C")  # cds
            states[len(states) - 1] = "S"  # start codon
            states[len(states) - 2] = "S"  # start codon
            states[len(states) - 3] = "S"  # start codon

    while len(states) < (len(fastas[chromosome_name].seq)):
        states.append("I")  # intergenetic
    return states


def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        return True

def str2bool2(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        return False

