#!python
import re
import numpy as np
import FastaFile

# bs674@cornell.edu

class Transcript:
    chromosome_name = ""
    if_orf_conserved = True
    meta_informaiton = ""
    name = ""
    start = 1000000000
    end = 0
    strand = ""
    genome_sequence = ""
    cds_sequence = ""
    exons = np.empty([0,2], int) #start and end position
    Cds = np.empty([0,2], int)
    def __init__(self, name, strand, chromosome_name):
        self.name = name
        self.strand = strand
        self.chromosome_name = chromosome_name

    def add_cds(self, start, end):
        self.Cds = np.append(self.Cds, np.array([[start, end]]), axis=0)
        # self.Cds.append([start, end])

    def add_exon(self, start, end):
        self.exons = np.append(self.exons, np.array([[start, end]]), axis=0)

    def updateCordinate(self):
        if ("+" == self.strand):
            self.Cds=np.sort(self.Cds, axis=0 )
        else:
            self.Cds=-np.sort(-self.Cds, axis=0)
        if( len(self.Cds) >0 ):
            if("+" == self.strand):
                self.start = self.Cds[0][0]
                self.end = self.Cds[len(self.Cds)-1][1]
            else:
                self.start = self.Cds[len(self.Cds) - 1][0]
                self.end = self.Cds[0][1]

    def __lt__(self, other):
        if (self.start < other.start):
            return True
        if (self.start == other.start) & (self.end > other.end):
            return True
        if (self.start == other.start) & (self.end == other.end):
            return False
        return False

    def __gt__(self, other):
        if (self.start < other.start):
            return False
        if (self.start == other.start) & (self.end > other.end):
            return False
        if (self.start == other.start) & (self.end == other.end):
            return False
        return True

    def __eq__(self, other):
        if (self.start == other.start) & (self.end == other.end):
            return True
        return False


class Gene:
    name = ""
    start = 100000000
    end = 0
    transcripts = np.empty([0, 1], Transcript)  # start and end position
    length = 0
    strand = ""
    def __init__(self, name, strand):
        self.name = name
        self.strand = strand
    def add_transcript(self, transcript):
        self.transcripts=np.append(self.transcripts,[transcript])

    def updateCordinate(self):
        self.transcripts = np.sort(self.transcripts)
        for transcript in self.transcripts:
            self.start = self.transcripts[0].start
            if(transcript.end > self.end):
                self.end = transcript.end
        self.length = self.end - self.start

    def __lt__(self, other):
        if (self.start < other.start):
            return True
        if (self.start == other.start) & (self.end > other.end):
            return True
        if (self.start == other.start) & (self.end == other.end):
            return False
        return False

    def __gt__(self, other):
        if (self.start < other.start):
            return False
        if (self.start == other.start) & (self.end > other.end):
            return False
        if (self.start == other.start) & (self.end == other.end):
            return False
        return True

    def __eq__(self, other):
        if (self.start == other.start) & (self.end == other.end):
            return True
        return False
#end of gene class


def readGff(gffFilePath):
    # the input if the gff file
    # return a dictionary [key(chromosome_name)]=genes  genes is a np.empty list
    chromosome_gene_dict=dict()
    chromosome_transcript_dict = dict()
    chromosome_gene_list = dict()  # start and end position
    transcript_gene_map=dict()
    with open(gffFilePath) as f:
        for line in f:
            m = re.search('^#', line)
            if( m == None ):
                m = re.search('^(\S+)\s+(\S+)\s+CDS\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+.*Parent=(\S+?)(;|$)', line)
                if( m != None ):
                    chromosome_name = m.group(1)
                    if (not chromosome_name in chromosome_transcript_dict):
                        chromosome_transcript_dict[chromosome_name]=dict()

                    start = int(m.group(3))
                    end = int(m.group(4))
                    if( start > end ):
                        temp = start
                        start = end
                        end = temp

                    transcript_name = m.group(8)
                    if (not transcript_name in chromosome_transcript_dict[chromosome_name]):
                        strand = m.group(6)
                        transcript = Transcript(transcript_name, strand, chromosome_name)
                        (chromosome_transcript_dict[chromosome_name])[transcript_name] = transcript

                    chromosome_transcript_dict[chromosome_name][transcript_name].add_cds(start, end)

                m = re.search('^(\S+)\s+(\S+)\s+mRNA\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+.*ID=(\S+?);.*Parent=(\S+?)(;|$)', line)
                if( m != None ):
                    transcript_name = m.group(8)
                    gene_name = m.group(9)
                    # print("transcript_name: "+transcript_name+" gene_name: "+gene_name)
                    transcript_gene_map[transcript_name] = gene_name

    #update transcript, get the start end information
    for chromosome_name in chromosome_transcript_dict:
        if (not chromosome_name in chromosome_gene_dict):
            chromosome_gene_dict[chromosome_name] = dict()
        for transcript_name in chromosome_transcript_dict[chromosome_name]:
            transcript = chromosome_transcript_dict[chromosome_name][transcript_name]
            transcript.updateCordinate()
            if transcript_name not in transcript_gene_map:
                transcript_gene_map[transcript_name] = transcript_name

            gene_name = transcript_gene_map[transcript_name]
            if not gene_name in chromosome_gene_dict[chromosome_name]:
                chromosome_gene_dict[chromosome_name][gene_name] = Gene(gene_name, transcript.strand)

            chromosome_gene_dict[chromosome_name][gene_name].add_transcript(transcript)
    #print("transcript update done")
    for chromosome_name in chromosome_gene_dict:
        gene_list = np.empty([0, 1], Gene)
        for gene_name in chromosome_gene_dict[chromosome_name]:
            chromosome_gene_dict[chromosome_name][gene_name].updateCordinate()
            gene_list = np.append(gene_list, [chromosome_gene_dict[chromosome_name][gene_name]])

        gene_list = np.sort(gene_list)
        chromosome_gene_list[chromosome_name] = gene_list

    return chromosome_gene_dict, chromosome_gene_list

def update_sequence_information_onechromosome(fastas, chromosome_gene_dict, chromosome_name):
    if chromosome_name in fastas:
        if chromosome_name in chromosome_gene_dict:
            for gene_name in chromosome_gene_dict[chromosome_name]:
                for transcript in chromosome_gene_dict[chromosome_name][gene_name].transcripts:
                    transcript.genome_sequence = FastaFile.getSubSequence(fastas, chromosome_name, transcript.start, transcript.end, transcript.strand)
                    cds_seq_list = []
                    for cds in transcript.Cds:
                        cds_seq_list.append(
                            FastaFile.getSubSequence(fastas, chromosome_name, cds[0], cds[1], transcript.strand))
                        transcript.cds_sequence = "".join(cds_seq_list)


def update_sequence_information(fastas, chromosome_gene_dict):
    for chromosome_name in fastas:
        update_sequence_information_onechromosome(fastas, chromosome_gene_dict, chromosome_name)



#test readGff function

# chromosome_gene_dict = readGff("/Users/song/tair11/Araport11_GFF3_genes_transposons.201606.gff")
# for chromosome_name in chromosome_gene_dict:
#     for gene_name in  chromosome_gene_dict[chromosome_name]:
#         print(chromosome_name + " " + gene_name)

#
# chromosome_gene_dict, chromosome_gene_list = readGff("/Users/song/tair11/TAIR10_GFF3_genes_no_UM.gff")
# for chromosome_name in chromosome_gene_list:
#     for gene in chromosome_gene_list[chromosome_name]:
#         print(gene.name)
#         for transcript in gene.transcripts:
#             print(transcript.name)
#
# # #
# chromosome_gene_dict, chromosome_gene_list = readGff("/Users/song/tair11/TAIR10_GFF3_genes_no_UM.gff")
# chromosome_names, fastas = FastaFile.readFastaFile("/Users/song/Col.fa")
# update_sequence_information_onechromosome(fastas, chromosome_gene_dict, "Chr1")
# print (chromosome_gene_dict["Chr1"]["AT1G01010"].transcripts[0].cds_sequence)
# print (chromosome_gene_dict["Chr1"]["AT1G01010"].transcripts[0].genome_sequence)
# assert(chromosome_gene_dict["Chr1"]["AT1G01010"].transcripts[0].cds_sequence == "ATGGAGGATCAAGTTGGGTTTGGGTTCCGTCCGAACGACGAGGAGCTCGTTGGTCACTATCTCCGTAACAAAATCGAAGGAAACACTAGCCGCGACGTTGAAGTAGCCATCAGCGAGGTCAACATCTGTAGCTACGATCCTTGGAACTTGCGCTTCCAGTCAAAGTACAAATCGAGAGATGCTATGTGGTACTTCTTCTCTCGTAGAGAAAACAACAAAGGGAATCGACAGAGCAGGACAACGGTTTCTGGTAAATGGAAGCTTACCGGAGAATCTGTTGAGGTCAAGGACCAGTGGGGATTTTGTAGTGAGGGCTTTCGTGGTAAGATTGGTCATAAAAGGGTTTTGGTGTTCCTCGATGGAAGATACCCTGACAAAACCAAATCTGATTGGGTTATCCACGAGTTCCACTACGACCTCTTACCAGAACATCAGAGGACATATGTCATCTGCAGACTTGAGTACAAGGGTGATGATGCGGACATTCTATCTGCTTATGCAATAGATCCCACTCCCGCTTTTGTCCCCAATATGACTAGTAGTGCAGGTTCTGTGGTCAACCAATCACGTCAACGAAATTCAGGATCTTACAACACTTACTCTGAGTATGATTCAGCAAATCATGGCCAGCAGTTTAATGAAAACTCTAACATTATGCAGCAGCAACCACTTCAAGGATCATTCAACCCTCTCCTTGAGTATGATTTTGCAAATCACGGCGGTCAGTGGCTGAGTGACTATATCGACCTGCAACAGCAAGTTCCTTACTTGGCACCTTATGAAAATGAGTCGGAGATGATTTGGAAGCATGTGATTGAAGAAAATTTTGAGTTTTTGGTAGATGAAAGGACATCTATGCAACAGCATTACAGTGATCACCGGCCCAAAAAACCTGTGTCTGGGGTTTTGCCTGATGATAGCAGTGATACTGAAACTGGATCAATGATTTTCGAAGACACTTCGAGCTCCACTGATAGTGTTGGTAGTTCAGATGAACCGGGCCATACTCGTATAGATGATATTCCATCATTGAACATTATTGAGCCTTTGCACAATTATAAGGCACAAGAGCAACCAAAGCAGCAGAGCAAAGAAAAGGTGATAAGTTCGCAGAAAAGCGAATGCGAGTGGAAAATGGCTGAAGACTCGATCAAGATACCTCCATCCACCAACACGGTGAAGCAGAGCTGGATTGTTTTGGAGAATGCACAGTGGAACTATCTCAAGAACATGATCATTGGTGTCTTGTTGTTCATCTCCGTCATTAGTTGGATCATTCTTGTTGGTTAA")
# assert(chromosome_gene_dict["Chr1"]["AT1G01010"].transcripts[0].genome_sequence == "ATGGAGGATCAAGTTGGGTTTGGGTTCCGTCCGAACGACGAGGAGCTCGTTGGTCACTATCTCCGTAACAAAATCGAAGGAAACACTAGCCGCGACGTTGAAGTAGCCATCAGCGAGGTCAACATCTGTAGCTACGATCCTTGGAACTTGCGCTGTAAGTTCCGAATTTTCTGAATTTCATTTGCAAGTAATCGATTTAGGTTTTTGATTTTAGGGTTTTTTTTTGTTTTGAACAGTCCAGTCAAAGTACAAATCGAGAGATGCTATGTGGTACTTCTTCTCTCGTAGAGAAAACAACAAAGGGAATCGACAGAGCAGGACAACGGTTTCTGGTAAATGGAAGCTTACCGGAGAATCTGTTGAGGTCAAGGACCAGTGGGGATTTTGTAGTGAGGGCTTTCGTGGTAAGATTGGTCATAAAAGGGTTTTGGTGTTCCTCGATGGAAGATACCCTGACAAAACCAAATCTGATTGGGTTATCCACGAGTTCCACTACGACCTCTTACCAGAACATCAGGTTTTCTTCTATTCATATATATATATATATATATATGTGGATATATATATATGTGGTTTCTGCTGATTCATAGTTAGAATTTGAGTTATGCAAATTAGAAACTATGTAATGTAACTCTATTTAGGTTCAGCAGCTATTTTAGGCTTAGCTTACTCTCACCAATGTTTTATACTGATGAACTTATGTGCTTACCTCCGGAAATTTTACAGAGGACATATGTCATCTGCAGACTTGAGTACAAGGGTGATGATGCGGACATTCTATCTGCTTATGCAATAGATCCCACTCCCGCTTTTGTCCCCAATATGACTAGTAGTGCAGGTTCTGTGGTGAGTCTTTCTCCATATACACTTAGCTTTGAGTAGGCAGATCAAAAAAGAGCTTGTGTCTACTGATTTGATGTTTTCCTAAACTGTTGATTCGTTTCAGGTCAACCAATCACGTCAACGAAATTCAGGATCTTACAACACTTACTCTGAGTATGATTCAGCAAATCATGGCCAGCAGTTTAATGAAAACTCTAACATTATGCAGCAGCAACCACTTCAAGGATCATTCAACCCTCTCCTTGAGTATGATTTTGCAAATCACGGCGGTCAGTGGCTGAGTGACTATATCGACCTGCAACAGCAAGTTCCTTACTTGGCACCTTATGAAAATGAGTCGGAGATGATTTGGAAGCATGTGATTGAAGAAAATTTTGAGTTTTTGGTAGATGAAAGGACATCTATGCAACAGCATTACAGTGATCACCGGCCCAAAAAACCTGTGTCTGGGGTTTTGCCTGATGATAGCAGTGATACTGAAACTGGATCAATGGTAAGCTTTTTTTACTCATATATAATCACAACCTATATCGCTTCTATATCTCACACGCTGAATTTTGGCTTTTAACAGATTTTCGAAGACACTTCGAGCTCCACTGATAGTGTTGGTAGTTCAGATGAACCGGGCCATACTCGTATAGATGATATTCCATCATTGAACATTATTGAGCCTTTGCACAATTATAAGGCACAAGAGCAACCAAAGCAGCAGAGCAAAGAAAAGGTTTAACACTCTCACTGAGAAACATGACTTTGATACGAAATCTGAATCAACATTTCATCAAAAAGATTTAGTCAAATGACCTCTAAATTATGAGCTATGGGTCTGCTTTCAGGTGATAAGTTCGCAGAAAAGCGAATGCGAGTGGAAAATGGCTGAAGACTCGATCAAGATACCTCCATCCACCAACACGGTGAAGCAGAGCTGGATTGTTTTGGAGAATGCACAGTGGAACTATCTCAAGAACATGATCATTGGTGTCTTGTTGTTCATCTCCGTCATTAGTTGGATCATTCTTGTTGGTTAA")
#
# print (chromosome_gene_dict["Chr1"]["AT1G01020"].transcripts[0].name)
# print (chromosome_gene_dict["Chr1"]["AT1G01020"].transcripts[0].cds_sequence)
# print (chromosome_gene_dict["Chr1"]["AT1G01020"].transcripts[0].genome_sequence)
#
# assert(chromosome_gene_dict["Chr1"]["AT1G01020"].transcripts[0].cds_sequence == "ATGGCGGCGAGTGAACACAGATGCGTGGGATGTGGTTTTAGGGTAAAGTCATTGTTCATTCAATACTCTCCGGGTAACATTCGTCTCATGAAATGCGGAAATTGCAAGGAAGTAGCAGATGAGTACATCGAGTGTGAACGCATGATTATTTTCATCGATTTAATCCTTCACAGACCAAAGGTATATAGACACGTCCTCTACAATGCAATTAATCCAGCAACTGTCAATATTCAGCATCTGTTGTGGAAGTTGGTCTTCGCCTATCTTCTTCTAGACTGTTATAGAAGCTTGCTACTGAGAAAAAGTGATGAAGAATCGAGCTTTTCTGATAGCCCTGTTCTTCTATCTATAAAGGTTCTGATTGGTGTCTTATCTGCAAACGCTGCATTTATCATCTCTTTTGCCATTGCGACTAAGGGTTTGCTAAATGAAGTTTCCAGAAGAAGAGAGATTATGTTGGGGATATTCATCTCTAGTTACTTCAAGATATTTCTGCTTGCGATGTTGGTATGGGAATTCCCAATGTCAGTGATTTTTTTTGTCGATATACTTCTCTTAACATCAAACTCCATGGCTCTTAAAGTGATGACTGAATCAACAATGACCAGATGCATAGCCGTATGCTTAATCGCGCACTTGATTAGATTCTTGGTGGGTCAGATTTTTGAGCCGACAATATTTTTGATACAAATTGGATCTCTGTTGCAATATATGTCTTATTTTTTCAGAATCGTATGA")
# assert(chromosome_gene_dict["Chr1"]["AT1G01020"].transcripts[0].genome_sequence == "ATGGCGGCGAGTGAACACAGATGCGTGGGATGTGGTTTTAGGGTAAAGTCATTGTTCATTCAATACTCTCCGGGTAACATTCGTCTCATGAAATGCGTAAGATTTTTGTATTAATGGATCTTCTTTGACCCAATTAGGGTTTATAAGATTCTCTATTAAATTTCTTCAAAGTTTCCATCTTTGATCATATTGTTATATTTAGGGAAATTGCAAGGAAGTAGCAGATGAGTACATCGAGTGTGAACGCATGGTCTGTTTTAATCACTTTCTTTCCCTTTTTGGGTGTTTCCCTATTTCTTTCCAATTTTGGTTTTTGGTTCTTTAATTGTTTTCTTTGGCAGATTATTTTCATCGATTTAATCCTTCACAGACCAAAGGTATATAGACACGTCCTCTACAATGCAATTAATCCAGCAACTGTCAATATTCAGGTAAAAAACTATTTCAATTCAAGTTCTTGATTTCAAGCATTGCTAGTTCCTGTCTGATTCTAGCTCTCTTATCTTGAATTGATATTACAAACTTGTGGAAATATCTGGCTTTGAAAATCTTTTTCCATATATAGAAGTATCTTAGACCTCAAGCAACAAACTAATTTCCAAGTTTGGGGATATTGAGTAGTCTCTGATCCAATAAGAATTGTTGGGTTCTTTTTTTTTTTTTGTGTGTGTGTGTGCAGCATCTGTTGTGGAAGTTGGTCTTCGCCTATCTTCTTCTAGACTGTTGTATCCTTTTGTTTTATAATCTTTTACTATCAAAGTGTAGCTAGCAGATGTGTAGTTGGTGAGATTCATCGGCAAACTTATTTGAAAAGCCCTGAACTCATCATAGATAGAAGCTTGCTACTGAGAAAAAGTGATGAAGAATCGAGCTTTTCTGATAGCCCTGTTCTTCTATCTATAAAGGTTAGATCGTTTCTTTTTAATGGACTTAATTAGATTTTTTTTCCTTTAAGGTTTAAATGATTTAAGTCATTGAGATTGTTTAGAAAGGGTAATAATTTCGCATTTGAAACAGGTTCTGATTGGTGTCTTATCTGCAAACGCTGCATTTATCATCTCTTTTGCCATTGCGACTAAGGGTTTGCTAAATGAAGTTTCCAGGTGATTACTGATTCTTCTTCTTCTTAGTTCTTATTATGAGTTTTTCATCAAAGATCTCAGACTTCTGTCACCTGAAGATTGTTTATGTATAATTTTTTTTAAAAACTTTGCAGAAGAAGAGAGATTATGTTGGGGATATTCATCTCTAGTTACTTCAAGATATTTCTGCTTGCGATGTTGGTATGTTGTAGCTTTACCTCTCACTTAATTCCTAATATTGAAGTTCCAAACTTCTTAAGCATTCCATAAATCATGCACTAAAGTTTCTTGTATTGATTAAACATGGTGTTATGTCTCTTTGCTCAAAAATAGAGAATCTTTGTTATTACAGGTATGGGAATTCCCAATGTCAGTGATTTTTTTTGTCGATATACTTCTCTTAACATCAAACTCCATGGCTCTTAAAGGTCAGACCTTCAAAATGTTCTCAATGCAAATTGTTTTTTGTTGTTGTTATTTTGGAATCTCACAATGCAAGTTTGTCTTTAAACCAGTGATGACTGAATCAACAATGACCAGATGCATAGCCGTATGCTTAATCGCGCACTTGATTAGATTCTTGGTGGGTCAGATTTTTGAGCCGACAATATTTTTGATACAAATTGGATCTCTGTTGCAATATATGTCTTATTTTTTCAGAATCGTATGA")

# ATGATTATTTTCATCGATTTAATCCTTCACAGACCAAAGGTATATAGACACGTCCTCTACAATGCAATTAATCCAGCAACTGTCAATATTCAGCATCTGTTGTGGAAGTTGGTCTTCGCCTATCTTCTTCTAGACTGTTATAGAAGCTTGCTACTGAGAAAAAGTGATGAAGAATCGAGCTTTTCTGATAGCCCTGTTCTTCTATCTATAAAGGTTCTGATTGGTGTCTTATCTGCAAACGCTGCATTTATCATCTCTTTTGCCATTGCGACTAAGGGTTTGCTAAATGAAGTTTCCAGAAGAAGAGAGATTATGTTGGGGATATTCATCTCTAGTTACTTCAAGATATTTCTGCTTGCGATGTTGGTATGGGAATTCCCAATGTCAGTGATTTTTTTTGTCGATATACTTCTCTTAACATCAAACTCCATGGCTCTTAAAGTGATGACTGAATCAACAATGACCAGATGCATAGCCGTATGCTTAATCGCGCACTTGATTAGATTCTTGGTGGGTCAGATTTTTGAGCCGACAATATTTTTGATACAAATTGGATCTCTGTTGCAATATATGTCTTATTTTTTCAGAATCGTATGA
# ATGGCGGCGAGTGAACACAGATGCGTGGGATGTGGTTTTAGGGTAAAGTCATTGTTCATTCAATACTCTCCGGGTAACATTCGTCTCATGAAATGCGGAAATTGCAAGGAAGTAGCAGATGAGTACATCGAGTGTGAACGCATGATTATTTTCATCGATTTAATCCTTCACAGACCAAAGGTATATAGACACGTCCTCTACAATGCAATTAATCCAGCAACTGTCAATATTCAGCATCTGTTGTGGAAGTTGGTCTTCGCCTATCTTCTTCTAGACTGTTATAGAAGCTTGCTACTGAGAAAAAGTGATGAAGAATCGAGCTTTTCTGATAGCCCTGTTCTTCTATCTATAAAGGTTCTGATTGGTGTCTTATCTGCAAACGCTGCATTTATCATCTCTTTTGCCATTGCGACTAAGGGTTTGCTAAATGAAGTTTCCAGAAGAAGAGAGATTATGTTGGGGATATTCATCTCTAGTTACTTCAAGATATTTCTGCTTGCGATGTTGGTATGGGAATTCCCAATGTCAGTGATTTTTTTTGTCGATATACTTCTCTTAACATCAAACTCCATGGCTCTTAAAGTGATGACTGAATCAACAATGACCAGATGCATAGCCGTATGCTTAATCGCGCACTTGATTAGATTCTTGGTGGGTCAGATTTTTGAGCCGACAATATTTTTGATACAAATTGGATCTCTGTTGCAATATATGTCTTATTTTTTCAGAATCGTATGA

# ATGGCGGCGAGTGAACACAGATGCGTGGGATGTGGTTTTAGGGTAAAGTCATTGTTCATTCAATACTCTCCGGGTAACATTCGTCTCATGAAATGCGGAAATTGCAAGGAAGTAGCAGATGAGTACATCGAGTGTGAACGCATGATTATTTTCATCGATTTAATCCTTCACAGACCAAAGGTATATAGACACGTCCTCTACAATGCAATTAATCCAGCAACTGTCAATATTCAGCATCTGTTGTGGAAGTTGGTCTTCGCCTATCTTCTTCTAGACTGTTATAGAAGCTTGCTACTGAGAAAAAGTGATGAAGAATCGAGCTTTTCTGATAGCCCTGTTCTTCTATCTATAAAGGTTCTGATTGGTGTCTTATCTGCAAACGCTGCATTTATCATCTCTTTTGCCATTGCGACTAAGGGTTTGCTAAATGAAGTTTCCAGAAGAAGAGAGATTATGTTGGGGATATTCATCTCTAGTTACTTCAAGATATTTCTGCTTGCGATGTTGGTATGTTGTAGCTTTACCTCTCACTTAATTCCTAATATTGAAGTTCCAAACTTCTTAAGCATTCCATAA
# ATGGCGGCGAGTGAACACAGATGCGTGGGATGTGGTTTTAGGGTAAAGTCATTGTTCATTCAATACTCTCCGGGTAACATTCGTCTCATGAAATGCGGAAATTGCAAGGAAGTAGCAGATGAGTACATCGAGTGTGAACGCATGATTATTTTCATCGATTTAATCCTTCACAGACCAAAGGTATATAGACACGTCCTCTACAATGCAATTAATCCAGCAACTGTCAATATTCAGCATCTGTTGTGGAAGTTGGTCTTCGCCTATCTTCTTCTAGACTGTTATAGAAGCTTGCTACTGAGAAAAAGTGATGAAGAATCGAGCTTTTCTGATAGCCCTGTTCTTCTATCTATAAAGGTTCTGATTGGTGTCTTATCTGCAAACGCTGCATTTATCATCTCTTTTGCCATTGCGACTAAGGGTTTGCTAAATGAAGTTTCCAGAAGAAGAGAGATTATGTTGGGGATATTCATCTCTAGTTACTTCAAGATATTTCTGCTTGCGATGTTGGTATGTTGTAGCTTTACCTCTCACTTAATTCCTAATATTGAAGTTCCAAACTTCTTAAGCATTCCATAA

