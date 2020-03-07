#!python
import re

# bs674@cornell.edu

class Fasta:
    name = ""
    seq = ""
    def __init__(self, name, seq):
        self.name = name
        self.seq = seq

def readFastaFile(fastaFile):
    fastas = {}
    chromosome_names = []
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
                    fasta = Fasta(name, s)
                    fastas[name]=(fasta)
                    chromosome_names.append(name)
                name = m.group(1)
                seq = []
            else:
                seq.append(line)
        if (len(name) > 0) & (len(seq) > 0):
            s = ''.join(seq)
            s = re.sub("\\s", "", s)
            s = s.upper()
            fasta = Fasta(name, s)
            fastas[name] = (fasta)
            chromosome_names.append(name)

    return chromosome_names, fastas

def getReverseComplementary(sequence):
    return sequence.translate(str.maketrans('ATUCGRYKMBVDH', 'TAAGCYRMKVBHD'))[::-1]
    # reversecomplementary=[]
    # for c in  sequence[::-1]:
    #     if ('A' == c):
    #         c = 'T'
    #     elif ('T' == c):
    #         c = 'A'
    #     elif ('U' == c):
    #         c = 'A'
    #     elif ('C' == c):
    #         c = 'G'
    #     elif ('G' == c):
    #         c = 'C'
    #     elif ('R' == c):
    #         c = 'Y'
    #     elif ('Y' == c):
    #         c = 'R'
    #     elif ('K' == c):
    #         c = 'M'
    #     elif ('M' == c):
    #         c = 'K'
    #     elif ('B' == c):
    #         c = 'V'
    #     elif ('V' == c):
    #         c = 'B'
    #     elif ('D' == c):
    #         c = 'H'
    #     elif ('H' == c):
    #         c = 'D'
    #     reversecomplementary.append(c)
    #
    # return ''.join(reversecomplementary)


def getSubSequence(fastas, name, start, end, strand="+"):
    # get a sequence fragment from fasta records
    start = start -1
    if start > len(fastas[name].seq):
        return ""
    if end  >  len(fastas[name].seq):
        end = len(fastas[name].seq)

    seq = (fastas[name].seq)[start:end]
    if( "+" == strand ):
        return seq
    else:
        return getReverseComplementary(seq)
#
# # ##test
# print(time.time())
# chromosome_names, fastas = readFastaFile("/Users/song/Col.fa")
# seq = getSubSequence(fastas, "Chr1", 3, 11, "+")
# print(seq)
# print (getReverseComplementary(seq))
# print(getSubSequence(fastas, "Chr1", 3, 11, "-"))
# print(time.time())
# #
# # print("substr"[2:4])
