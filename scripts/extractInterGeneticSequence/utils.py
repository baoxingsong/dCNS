#!python
import re
import numpy as np
from argparse import ArgumentParser
import sys

# bs674@cornell.edu
# Baoxing Song
# 26 April, 2019, updated

#

# this is to define what is a gene
# not matter what the strand is, the start is always smaller than end
class Gene:
    def __init__(self, name, strand, start, end, chr):
        self.name = name
        self.strand = strand
        self.start = start
        self.end = end
        self.chr = chr

    # this is for sorting purpose
    def __lt__(self, other):
        if (self.start < other.start):
            return True
        if (self.start == other.start) & (self.end < other.end):
            return True
        return False
    # this is for sorting purpose
    def __gt__(self, other):
        if (self.start < other.start):
            return False
        if (self.start == other.start) & (self.end < other.end):
            return False
        if (self.start == other.start) & (self.end == other.end):
            return False
        return True
    # this is for sorting purpose
    def __eq__(self, other):
        if (self.start == other.start) & (self.end == other.end):
            return True
        return False
#end of gene class

class Interval:
    def __init__(self, start, end):
        self.start = start
        self.end = end
    # this is for sorting purpose
    def __lt__(self, other):
        if (self.start < other.start):
            return True
        if (self.start == other.start) & (self.end < other.end):
            return True
        return False
    # this is for sorting purpose
    def __gt__(self, other):
        if (self.start < other.start):
            return False
        if (self.start == other.start) & (self.end < other.end):
            return False
        if (self.start == other.start) & (self.end == other.end):
            return False
        return True
    # this is for sorting purpose
    def __eq__(self, other):
        if (self.start == other.start) & (self.end == other.end):
            return True
        return False
#end of gene class

#
# class GenePair:
#     def __init__(self, geneA, geneB, strandA, strandB, basePairDistance, numberOfGeneDistance):
#         self.geneA = geneA
#         self.geneB = geneB
#         self.strandA = strandA
#         self.strandB = strandB
#         self.basePairDistance = basePairDistance
#         self.numberOfGeneDistance = numberOfGeneDistance
#
#     # the __hash__ and __eq__ functions are for dictionary data structure purpose
#     def __hash__(self):
#         return hash(self.geneA.name)+hash(self.geneB.name)  # the order has no effect on the return value
#
#     def __eq__(self, other):
#         if (self.geneA.name, self.geneB.name, self.strandA, self.strandB) == (other.geneA.name, other.geneB.name, other.strandA, other.strandB): #or (self.geneA, self.geneB) == (other.geneB, other.geneA)
#             # a+,b+,A+,B+   a-,b-,A-,B-   a-,b+,A-,B+  a+,b-,A+,B-
#             return True
#         elif ( self.geneA.name==other.geneB.name and self.geneB.name==other.geneA.name and self.strandA!=other.strandA and self.strandB!=other.strandB):
#             # a+,b+,B-,A-   a-,b-,B+,A+    a+,b-,B+,A-   a-,b+,B-,A+
#             return True
#         return False
#
#     def __ne__(self, other):
#         # Not strictly necessary, but to avoid having both x==y and x!=y
#         # True at the same time
#         return not (self == other)
# # end of GenePair class


#read a fasta file and return a dictionary, the key is entry id and the value is the sequence in upcase
def readFastaFile(fastaFile):
    fastas = {}
    name = ""
    seq = []
    with open(fastaFile) as f:
        for line in f:
            m = re.search('^>(\S+)', line)
            if (m != None):
                if (len(name) > 0) & (len(seq) > 0):
                    s = ''.join(seq)
                    s = re.sub("\\s", "", s)
                    s = s.upper()
                    fastas[name] = s
                name = m.group(1)
                seq = []
            else:
                seq.append(line)
        if (len(name) > 0) & (len(seq) > 0):
            s = ''.join(seq)
            s = re.sub("\\s", "", s)
            s = s.upper()
            fastas[name] = s
    return fastas


def getReverseComplementary(sequence):
    return sequence.translate(str.maketrans('ATUCGRYKMBVDH', 'TAAGCYRMKVBHD'))[::-1]

def getSubSequence(fastas, name, start, end, strand="+"):
    # get a sequence fragment from fasta records
    start = start - 1
    if start < 0:
        start = 0
    if start > len(fastas[name]):
        return ""
    if end  >  len(fastas[name]):
        end = len(fastas[name])

    seq = (fastas[name])[start:end]
    if( "+" == strand ):
        return seq
    else:
        return getReverseComplementary(seq)

#read gff file
def readGff(gffFilePath):
    # the input if the gff file
    # return a dictionary [key(chromosome_name)]=genes  genes is a np.empty list
    chromosome_gene_list = dict()  # map<chr, vector<gene>>
    transcript_gene = dict() # key is transcript id and value is gene id
    gene_start = dict() # the gene start position, which is defined as the start codon position of the primary transcript
    gene_end = dict() # the gene end position, which is defined as the end codon position of the primary transcript
    gene_strand = dict()
    gene_name = ""
    chromosome_name = ""
    with open(gffFilePath) as f:
        for line in f:
            m = re.search('^#', line)
            if( m == None ):
                m = re.search('^(\S+)\s+(\S+)\s+gene\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+.*ID=(.*?);', line)
                if (m != None): # The begaining of the new gene entry suggests the end of information of last gene entry. So here we could organize the data structure of the last gene
                    if len(gene_name)>0 and gene_name in gene_start:
                        g = Gene(gene_name, gene_strand[gene_name], gene_start[gene_name], gene_end[gene_name], chromosome_name)
                        chromosome_gene_list[chromosome_name] = np.append(chromosome_gene_list[chromosome_name], [g])

                m1 = re.search('^(\S+)\s+(\S+)\s+mRNA\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+.*ID=(.*?);.*?Parent=(.*?);', line)
                m2 = re.search('^(\S+)\s+(\S+)\s+mRNA\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+.*ID=(.*?);.*?Parent=(.*?)$', line)
                if m1 != None:
                    m = m1
                else:
                    m = m2
                if (m != None) and "curator_summary=pseudogene" not in line:
                    chromosome_name = m.group(1)
                    if (not chromosome_name in chromosome_gene_list):
                        chromosome_gene_list[chromosome_name] = np.empty([0, 1], Gene)
                    strand = m.group(6)
                    name = m.group(8)
                    name = name.replace("transcript:", "")
                    pname = m.group(9)
                    pname = pname.replace("gene:", "")
                    gene_name=pname
                    if pname not in gene_strand: # if there is no gene stand information, suggest this is the first transcript of this gene, and we define the first transcript as the primary transcript
                        gene_strand[pname] = strand
                        transcript_gene[name] = pname
                else:
                    m = re.search('^(\S+)\s+(\S+)\s+CDS\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+.*ID=(.*?);Parent=(.*?)[;\n]', line) #update the start and end information of the gene
                    # and no matter the strand information the start is always smaller than end
                    if (m != None):
                        pname = m.group(9)
                        pname = pname.replace("transcript:", "")
                        if pname in transcript_gene:
                            start = int(m.group(3))
                            end = int(m.group(4))
                            if start > end:
                                temp = start
                                start = end
                                end = temp
                            if transcript_gene[pname] in gene_start:
                                if gene_start[transcript_gene[pname]] > start:
                                    gene_start[transcript_gene[pname]] = start
                                if gene_end[transcript_gene[pname]] < end:
                                    gene_end[transcript_gene[pname]] = end
                            else:
                                gene_start[transcript_gene[pname]] = start
                                gene_end[transcript_gene[pname]] = end
        if len(gene_name)>0:
            g = Gene(gene_name, gene_strand[gene_name], gene_start[gene_name], gene_end[gene_name], chromosome_name)
            chromosome_gene_list[chromosome_name] = np.append(chromosome_gene_list[chromosome_name], [g])
    for chromosome_name in chromosome_gene_list:
        chromosome_gene_list[chromosome_name] = sorted(chromosome_gene_list[chromosome_name]) # sort the gene using the start position
    return chromosome_gene_list



#read gff file
def allTheCodingInterval(gffFilePath):
    # the input if the gff file
    # return a dictionary [key(chromosome_name)]=genes  genes is a np.empty list
    chromosome_gene_list = dict()  # map<chr, vector<gene>>
    transcript_gene = dict() # key is transcript id and value is gene id
    gene_start = dict() # the gene start position, which is defined as the start codon position of the primary transcript
    gene_end = dict() # the gene end position, which is defined as the end codon position of the primary transcript
    gene_strand = dict()
    gene_name = ""
    chromosome_name = ""
    with open(gffFilePath) as f:
        for line in f:
            m = re.search('^#', line)
            if( m == None ):
                #m = re.search('^(\S+)\s+(\S+)\s+gene\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+.*ID=(.*?);.*?biotype=protein_coding;', line)
                m = re.search('^(\S+)\s+(\S+)\s+gene\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+.*ID=(.*?);', line)
                if (m != None): # The begaining of the new gene entry suggests the end of information of last gene entry. So here we could organize the data structure of the last gene
                    if len(gene_name)>0:
                        g = Gene(gene_name, gene_strand[gene_name], gene_start[gene_name], gene_end[gene_name], chromosome_name)
                        chromosome_gene_list[chromosome_name] = np.append(chromosome_gene_list[chromosome_name], [g])

                #m = re.search('^(\S+)\s+(\S+)\s+mRNA\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+.*ID=(.*?);Parent=(.*?);.*?biotype=protein_coding;', line)
                m = re.search('^(\S+)\s+(\S+)\s+mRNA\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+.*ID=(.*?);Parent=(.*?);', line)
                if (m != None):
                    chromosome_name = m.group(1)
                    if (not chromosome_name in chromosome_gene_list):
                        chromosome_gene_list[chromosome_name] = np.empty([0, 1], Gene)
                    strand = m.group(6)
                    name = m.group(8)
                    name = name.replace("transcript:", "")
                    pname = m.group(9)
                    pname = pname.replace("gene:", "")
                    gene_name=pname
                    if pname not in gene_strand: # is there is gene stand information, suggest this is the first transcript of this gene, and we define the first transcript as the primary transcript
                        gene_strand[pname] = strand
                    transcript_gene[name] = pname
                else:
                    m = re.search('^(\S+)\s+(\S+)\s+CDS\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+.*ID=(.*?);Parent=(.*?);', line) #update the start and end information of the gene
                    # and no matter the strand information the start is always smaller than end
                    if (m != None):
                        pname = m.group(9)
                        pname = pname.replace("transcript:", "")
                        if pname in transcript_gene:
                            start = int(m.group(3))
                            end = int(m.group(4))
                            if start > end:
                                temp = start
                                start = end
                                end = temp
                            if transcript_gene[pname] in gene_start:
                                if gene_start[transcript_gene[pname]] > start:
                                    gene_start[transcript_gene[pname]] = start
                                if gene_end[transcript_gene[pname]] < end:
                                    gene_end[transcript_gene[pname]] = end
                            else:
                                gene_start[transcript_gene[pname]] = start
                                gene_end[transcript_gene[pname]] = end
        if len(gene_name)>0:
            g = Gene(gene_name, gene_strand[gene_name], gene_start[gene_name], gene_end[gene_name], chromosome_name)
            chromosome_gene_list[chromosome_name] = np.append(chromosome_gene_list[chromosome_name], [g])

    for chromosome_name in chromosome_gene_list:
        chromosome_gene_list[chromosome_name] = sorted(chromosome_gene_list[chromosome_name]) # sort the gene using the start position
    return chromosome_gene_list


#read gff file
def allThemRnaInterval(gffFilePath):
    # the input if the gff file
    # return a dictionary [key(chromosome_name)]=genes  genes is a np.empty list
    chromosome_gene_list = dict()  # map<chr, vector<gene>>

    with open(gffFilePath) as f:
        for line in f:
            m = re.search('^#', line)
            if( m == None ):
                m = re.search('^(\S+)\s+(\S+)\s+mRNA\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+.*ID=(.*?);Parent=(.*?);.*?biotype=protein_coding;', line)
                if (m != None):
                    chromosome_name = m.group(1)
                    if (not chromosome_name in chromosome_gene_list):
                        chromosome_gene_list[chromosome_name] = np.empty([0, 1], Interval)
                    strand = m.group(6)
                    name = m.group(8)
                    name = name.replace("transcript:", "")
                    start = int(m.group(3))
                    end = int(m.group(4))
                    if start > end:
                        temp = start
                        start = end
                        end = temp
                    g = Interval(start, end)
                    chromosome_gene_list[chromosome_name] = np.append(chromosome_gene_list[chromosome_name], [g])


    for chromosome_name in chromosome_gene_list:
        chromosome_gene_list[chromosome_name] = sorted(chromosome_gene_list[chromosome_name]) # sort the gene using the start position
    return chromosome_gene_list


# read the sam file to define gene positon on the new assemblyPolishing
# if the ratio of identical to CDS sequence length is larger than similarPercentage, we think this is gene
def readSam(minimap2SamFile, fastas, similarPercentage):
    gene_name_dict = dict() # key is the gene name, and value is a list of gene
    chromosome_gene_list = dict()
    with open(minimap2SamFile) as f:
        for line in f:
            m = re.search('^#', line)
            if( m == None ):
                m = re.search('^(.*?)\s+(\d+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\w+)\t.*NM:i:(\d+).*cs:(\S+)', line)
                if( m != None and m.group(3) != "*" ):
                    chromosome_name = m.group(3)
                    if (not chromosome_name in chromosome_gene_list):
                        chromosome_gene_list[chromosome_name] = np.empty([0, 1], Gene)

                    match=0
                    ms = re.findall(':(\d+)', m.group(8))
                    for mm in ms:
                        match += int(mm)

                    start = int(m.group(4))
                    end = int(m.group(4))
                    ms = re.findall('(\d+)[MDN=X]', m.group(6))
                    for mm in ms:
                        end += int(mm)
                    tag = int(m.group(2))
                    if (m.group(1) in fastas) and ((float(match)/len(fastas[m.group(1)])) >= similarPercentage) : ## only use the alignment pass thredshole
                        strand = ""
                        if( 0 == tag % 32 ):
                            strand = "+"
                        else:
                            strand = "-"
                        if end > start:   # minute by one
                            end -= 1
                        name = m.group(1)
                        g = Gene(name, strand, start, end, m.group(3))
                        if (not name in gene_name_dict):
                            gene_name_dict[name] = np.empty([0, 1], Gene)
                        gene_name_dict[name] = np.append(gene_name_dict[name], [g]) # a gene could have multiple copies, so here use a list
                        chromosome_gene_list[chromosome_name] = np.append(chromosome_gene_list[chromosome_name], [g])

    for chromosome_name in chromosome_gene_list:
        chromosome_gene_list[chromosome_name] = sorted(chromosome_gene_list[chromosome_name])# sort the gene using the start position

    return chromosome_gene_list, gene_name_dict


def chromosome_gene_list_2_genepairst(chromosome_gene_list, windowSize, onlySameStrand, maxDistance, couldGenePairOverlap):
    genepairs = dict()
    if onlySameStrand:
        for chromosome_name in chromosome_gene_list:
            gene_list = chromosome_gene_list[chromosome_name]
            for j in range(1, len(gene_list)):
                if gene_list[j].strand == "+":
                    i = j
                    count = 0
                    while i > 0:
                        i = i - 1
                        if gene_list[i].strand == "+":
                            distance = gene_list[j].start - gene_list[i].end - 1
                            if (couldGenePairOverlap or distance > 0):
                                count = count + 1
                                if count < windowSize:
                                    if (maxDistance<0 or distance<=maxDistance):
                                        gp = GenePair(gene_list[i], gene_list[j], gene_list[i].strand, gene_list[j].strand, distance, count) # here we have i<j and always put the gene_list[i].name before gene_list[j].name
                                        if gp in genepairs:
                                            if gp.basePairDistance < genepairs[gp]:
                                                genepairs.pop('gp', None)
                                                genepairs[gp] = distance
                                        else:
                                            genepairs[gp] = distance

                                else:
                                    break
                else:
                    i = j
                    count = 0
                    while i > 0:
                        i = i - 1
                        if gene_list[i].strand == "-":
                            distance = gene_list[j].start - gene_list[i].end - 1
                            if (couldGenePairOverlap or distance > 0):
                                count = count + 1
                                if count < windowSize:
                                    if (maxDistance<0 or distance<=maxDistance):
                                        gp = GenePair(gene_list[i], gene_list[j], gene_list[i].strand, gene_list[j].strand, distance, count) ## the order of gene name is important here
                                        if gp in genepairs:
                                            if gp.basePairDistance < genepairs[gp]:
                                                genepairs.pop('gp', None)
                                                genepairs[gp] = distance
                                        else:
                                            genepairs[gp] = distance
                                else:
                                    break
    else:
        for chromosome_name in chromosome_gene_list:
            gene_list = chromosome_gene_list[chromosome_name]
            for j in range(1, len(gene_list)):
                i = j
                count = 0
                while i > 0:
                    i = i - 1
                    distance = gene_list[j].start - gene_list[i].end - 1
                    if (couldGenePairOverlap or distance > 0):
                        count = count + 1
                        if count < windowSize:
                            if (maxDistance<0 or distance<=maxDistance):
                                gp = GenePair(gene_list[i], gene_list[j], gene_list[i].strand, gene_list[j].strand, distance, count)
                                if gp in genepairs:
                                    if gp.basePairDistance < genepairs[gp]:
                                        genepairs.pop('gp', None)
                                        genepairs[gp] = distance
                                else:
                                    genepairs[gp] = distance

                        else:
                            break
    return genepairs

def two_genepairs_overlaps(genepairs1, genepairs2):
    genepairs = dict()
    for genepair in genepairs1:
        if genepair in genepairs2:
            genepairs[genepair] = 1
    return genepairs

# for parameter parsing begin
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
# for parameter parsing end


class SimilarFragment:
    def __init__(self, ref_chr, ref_strand, ref_start, ref_end, query_chr, query_strand, query_start, query_end, name):
        self.ref_chr=ref_chr
        self.query_chr=query_chr
        if ref_strand == query_strand:
            self.ref_strand = "+"
            self.query_strand = "+"
        else:
            self.ref_strand = "+"
            self.query_strand = "-"
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.query_start = query_start
        self.query_end = query_end
        self.name = name

    def __lt__(self, other):
        if (self.ref_start < other.ref_start):
            return True
        if (self.ref_start == other.ref_start) & (self.ref_end < other.ref_end):
            return True
        return False
    # this is for sorting purpose
    def __gt__(self, other):
        if (self.ref_start < other.ref_start):
            return False
        if (self.ref_start == other.ref_start) & (self.ref_end < other.ref_end):
            return False
        if (self.ref_start == other.ref_start) & (self.ref_end == other.ref_end):
            return False
        return True
    # this is for sorting purpose
    def __eq__(self, other):
        if (self.ref_start == other.ref_start) & (self.ref_end == other.ref_end):
            return True
        return False
