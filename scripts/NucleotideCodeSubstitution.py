#!python
import os.path
import re
import sys

import numpy as np

import FastaFile

# bs674@cornell.edu


location = os.path.dirname(sys.argv[0])
location = os.path.join(location,"SpliceSites")

dornors   = []
acceptors = []
with open(location) as f:
    for line in f:
        line = re.sub("\#.*$", "", line)
        m = re.search("^(\w+)\s+(\w+)", line)
        if m != None:
            dornor = m.group(1)
            acceptor = m.group(2)
            dornor = dornor.upper()
            acceptor = acceptor.upper()
            dornors.append(dornor)
            acceptors.append(acceptor)

#iupac code begin
dnaIupacCode = {}
dnaIupacCode["A"]={}
dnaIupacCode["A"]["A"]=1

dnaIupacCode["C"]={}
dnaIupacCode["C"]["C"]=1

dnaIupacCode["G"]={}
dnaIupacCode["G"]["G"]=1

dnaIupacCode["T"]={}
dnaIupacCode["T"]["T"]=1
dnaIupacCode["T"]["U"]=1

dnaIupacCode["U"]={}
dnaIupacCode["U"]["T"]=1
dnaIupacCode["U"]["U"]=1

dnaIupacCode["R"]={}
dnaIupacCode["R"]["A"]=1
dnaIupacCode["R"]["G"]=1


dnaIupacCode["Y"] = {}
dnaIupacCode["Y"]["C"] = 1
dnaIupacCode["Y"]["T"] = 1
dnaIupacCode["Y"]["U"] = 1

dnaIupacCode["S"] = {}
dnaIupacCode["S"]["G"] = 1
dnaIupacCode["S"]["C"] = 1

dnaIupacCode["W"] = {}
dnaIupacCode["W"]["A"] = 1
dnaIupacCode["W"]["T"] = 1
dnaIupacCode["W"]["U"] = 1

dnaIupacCode["K"] = {}
dnaIupacCode["K"]["G"] = 1
dnaIupacCode["K"]["T"] = 1
dnaIupacCode["K"]["U"] = 1

dnaIupacCode["M"] = {}
dnaIupacCode["M"]["A"] = 1
dnaIupacCode["M"]["C"] = 1

dnaIupacCode["B"] = {}
dnaIupacCode["B"]["C"] = 1
dnaIupacCode["B"]["T"] = 1
dnaIupacCode["B"]["U"] = 1
dnaIupacCode["B"]["G"] = 1

dnaIupacCode["D"] = {}
dnaIupacCode["D"]["A"] = 1
dnaIupacCode["D"]["T"] = 1
dnaIupacCode["D"]["U"] = 1
dnaIupacCode["D"]["G"] = 1

dnaIupacCode["H"] = {}
dnaIupacCode["H"]["A"] = 1
dnaIupacCode["H"]["T"] = 1
dnaIupacCode["H"]["U"] = 1
dnaIupacCode["H"]["C"] = 1

dnaIupacCode["V"] = {}
dnaIupacCode["V"]["A"] = 1
dnaIupacCode["V"]["G"] = 1
dnaIupacCode["V"]["C"] = 1

dnaIupacCode["N"] = {}
dnaIupacCode["N"]["A"] = 1
dnaIupacCode["N"]["G"] = 1
dnaIupacCode["N"]["C"] = 1
dnaIupacCode["N"]["T"] = 1
dnaIupacCode["N"]["U"] = 1


#start codon and stop codon
mustStartCodons = dict()
mustStartCodons["TTG"] = 1
mustStartCodons["CTG"] = 1
mustStartCodons["ATG"] = 1
mustStartCodons["YTG"] = 1
mustStartCodons["WTG"] = 1
mustStartCodons["MTG"] = 1
mustStartCodons["HTG"] = 1

mustStopCodons = dict()
mustStopCodons["TAA"]=1
mustStopCodons["TAG"]=1
mustStopCodons["TGA"]=1
mustStopCodons["TAR"]=1
mustStopCodons["TRA"]=1

basicStartCodons = dict()
basicStartCodons["ATG"] = 1
basicStartCodons["TTG"] = 1
basicStartCodons["CTG"] = 1

basicStopCodons = dict()
basicStopCodons["TAA"] = 1
basicStopCodons["TAG"] = 1
basicStopCodons["TGA"] = 1

possibleStartCodons = {}
for basicStartCodon in  basicStartCodons:
        a = basicStartCodon[0]
        b = basicStartCodon[1]
        c = basicStartCodon[2]

        for key1 in dnaIupacCode:
            if a in dnaIupacCode[key1]:
                for key2 in dnaIupacCode:
                    if b in dnaIupacCode[key2]:
                        for key3 in dnaIupacCode:
                            if c in dnaIupacCode[key3]:
                                threeNa = key1 + key2 + key3
                                possibleStartCodons[threeNa] = 1

possibleStopCodons = {}
for basicStopCodon in  basicStopCodons:
        a = basicStopCodon[0]
        b = basicStopCodon[1]
        c = basicStopCodon[2]

        for key1 in dnaIupacCode:
            if a in dnaIupacCode[key1]:
                for key2 in dnaIupacCode:
                    if b in dnaIupacCode[key2]:
                        for key3 in dnaIupacCode:
                            if c in dnaIupacCode[key3]:
                                threeNa = key1 + key2 + key3
                                possibleStopCodons[threeNa] = 1

# genetic codon
middleStandardGeneticCode={}

middleStandardGeneticCode["TAA"] = '*'
middleStandardGeneticCode["TAG"] = '*'
middleStandardGeneticCode["TGA"] = '*'

middleStandardGeneticCode["TRA"] = '*'
middleStandardGeneticCode["TAR"] = '*'

middleStandardGeneticCode["TTT"] = 'F'
middleStandardGeneticCode["TTC"] = 'F'

middleStandardGeneticCode["TTY"] = 'F'


middleStandardGeneticCode["TTA"] = 'L'
middleStandardGeneticCode["TTG"] = 'L'

middleStandardGeneticCode["CTC"] = 'L'
middleStandardGeneticCode["CTT"] = 'L'
middleStandardGeneticCode["CTA"] = 'L'
middleStandardGeneticCode["CTG"] = 'L'

middleStandardGeneticCode["TTR"] = 'L'

middleStandardGeneticCode["YTA"] = 'L'
middleStandardGeneticCode["YTG"] = 'L'
middleStandardGeneticCode["YTR"] = 'L'

middleStandardGeneticCode["CTR"] = 'L'
middleStandardGeneticCode["CTY"] = 'L'
middleStandardGeneticCode["CTS"] = 'L'
middleStandardGeneticCode["CTW"] = 'L'
middleStandardGeneticCode["CTK"] = 'L'
middleStandardGeneticCode["CTM"] = 'L'
middleStandardGeneticCode["CTB"] = 'L'
middleStandardGeneticCode["CTD"] = 'L'
middleStandardGeneticCode["CTH"] = 'L'
middleStandardGeneticCode["CTV"] = 'L'
middleStandardGeneticCode["CTN"] = 'L'

middleStandardGeneticCode["TCC"] = 'S'
middleStandardGeneticCode["TCT"] = 'S'
middleStandardGeneticCode["TCA"] = 'S'
middleStandardGeneticCode["TCG"] = 'S'

middleStandardGeneticCode["TCR"] = 'S'
middleStandardGeneticCode["TCY"] = 'S'
middleStandardGeneticCode["TCS"] = 'S'
middleStandardGeneticCode["TCW"] = 'S'
middleStandardGeneticCode["TCK"] = 'S'
middleStandardGeneticCode["TCM"] = 'S'
middleStandardGeneticCode["TCB"] = 'S'
middleStandardGeneticCode["TCD"] = 'S'
middleStandardGeneticCode["TCH"] = 'S'
middleStandardGeneticCode["TCV"] = 'S'
middleStandardGeneticCode["TCN"] = 'S'

middleStandardGeneticCode["AGT"] = 'S'
middleStandardGeneticCode["AGC"] = 'S'
middleStandardGeneticCode["AGY"] = 'S'

middleStandardGeneticCode["TAT"] = 'Y'
middleStandardGeneticCode["TAC"] = 'Y'
middleStandardGeneticCode["TAY"] = 'Y'


middleStandardGeneticCode["TGT"] = 'C'
middleStandardGeneticCode["TGC"] = 'C'
middleStandardGeneticCode["TGY"] = 'C'


middleStandardGeneticCode["TGG"] = 'W'


middleStandardGeneticCode["CCT"] = 'P'
middleStandardGeneticCode["CCC"] = 'P'
middleStandardGeneticCode["CCA"] = 'P'
middleStandardGeneticCode["CCG"] = 'P'

middleStandardGeneticCode["CCR"] = 'P'
middleStandardGeneticCode["CCY"] = 'P'
middleStandardGeneticCode["CCS"] = 'P'
middleStandardGeneticCode["CCW"] = 'P'
middleStandardGeneticCode["CCK"] = 'P'
middleStandardGeneticCode["CCM"] = 'P'
middleStandardGeneticCode["CCB"] = 'P'
middleStandardGeneticCode["CCD"] = 'P'
middleStandardGeneticCode["CCH"] = 'P'
middleStandardGeneticCode["CCV"] = 'P'
middleStandardGeneticCode["CCN"] = 'P'


middleStandardGeneticCode["CAT"] = 'H'
middleStandardGeneticCode["CAC"] = 'H'
middleStandardGeneticCode["CAY"] = 'H'


middleStandardGeneticCode["CAA"] = 'Q'
middleStandardGeneticCode["CAG"] = 'Q'
middleStandardGeneticCode["CAR"] = 'Q'


middleStandardGeneticCode["CGT"] = 'R'
middleStandardGeneticCode["CGC"] = 'R'
middleStandardGeneticCode["CGA"] = 'R'
middleStandardGeneticCode["CGG"] = 'R'

middleStandardGeneticCode["CGR"] = 'R'
middleStandardGeneticCode["CGY"] = 'R'
middleStandardGeneticCode["CGS"] = 'R'
middleStandardGeneticCode["CGW"] = 'R'
middleStandardGeneticCode["CGK"] = 'R'
middleStandardGeneticCode["CGM"] = 'R'
middleStandardGeneticCode["CGB"] = 'R'
middleStandardGeneticCode["CGD"] = 'R'
middleStandardGeneticCode["CGH"] = 'R'
middleStandardGeneticCode["CGV"] = 'R'
middleStandardGeneticCode["CGN"] = 'R'

middleStandardGeneticCode["AGA"] = 'R'
middleStandardGeneticCode["AGG"] = 'R'
middleStandardGeneticCode["AGR"] = 'R'

middleStandardGeneticCode["MGA"] = 'R'
middleStandardGeneticCode["MGG"] = 'R'
middleStandardGeneticCode["MGR"] = 'R'


middleStandardGeneticCode["ATT"] = 'I'
middleStandardGeneticCode["ATC"] = 'I'
middleStandardGeneticCode["ATA"] = 'I'

middleStandardGeneticCode["ATY"] = 'I'
middleStandardGeneticCode["ATW"] = 'I'
middleStandardGeneticCode["ATM"] = 'I'


middleStandardGeneticCode["ATG"] = 'M'


middleStandardGeneticCode["ACT"] = 'T'
middleStandardGeneticCode["ACC"] = 'T'
middleStandardGeneticCode["ACA"] = 'T'
middleStandardGeneticCode["ACG"] = 'T'

middleStandardGeneticCode["ACR"] = 'T'
middleStandardGeneticCode["ACY"] = 'T'
middleStandardGeneticCode["ACS"] = 'T'
middleStandardGeneticCode["ACW"] = 'T'
middleStandardGeneticCode["ACK"] = 'T'
middleStandardGeneticCode["ACM"] = 'T'
middleStandardGeneticCode["ACB"] = 'T'
middleStandardGeneticCode["ACD"] = 'T'
middleStandardGeneticCode["ACH"] = 'T'
middleStandardGeneticCode["ACV"] = 'T'
middleStandardGeneticCode["ACN"] = 'T'


middleStandardGeneticCode["AAT"] = 'N'
middleStandardGeneticCode["AAC"] = 'N'
middleStandardGeneticCode["AAY"] = 'N'


middleStandardGeneticCode["AAA"] = 'K'
middleStandardGeneticCode["AAG"] = 'K'

middleStandardGeneticCode["AAR"] = 'K'


middleStandardGeneticCode["GTT"] = 'V'
middleStandardGeneticCode["GTC"] = 'V'
middleStandardGeneticCode["GTA"] = 'V'
middleStandardGeneticCode["GTG"] = 'V'

middleStandardGeneticCode["GTR"] = 'V'
middleStandardGeneticCode["GTY"] = 'V'
middleStandardGeneticCode["GTS"] = 'V'
middleStandardGeneticCode["GTW"] = 'V'
middleStandardGeneticCode["GTK"] = 'V'
middleStandardGeneticCode["GTM"] = 'V'
middleStandardGeneticCode["GTB"] = 'V'
middleStandardGeneticCode["GTD"] = 'V'
middleStandardGeneticCode["GTH"] = 'V'
middleStandardGeneticCode["GTV"] = 'V'
middleStandardGeneticCode["GTN"] = 'V'


middleStandardGeneticCode["GCT"] = 'A'
middleStandardGeneticCode["GCC"] = 'A'
middleStandardGeneticCode["GCA"] = 'A'
middleStandardGeneticCode["GCG"] = 'A'

middleStandardGeneticCode["GCR"] = 'A'
middleStandardGeneticCode["GCY"] = 'A'
middleStandardGeneticCode["GCS"] = 'A'
middleStandardGeneticCode["GCW"] = 'A'
middleStandardGeneticCode["GCK"] = 'A'
middleStandardGeneticCode["GCM"] = 'A'
middleStandardGeneticCode["GCB"] = 'A'
middleStandardGeneticCode["GCD"] = 'A'
middleStandardGeneticCode["GCH"] = 'A'
middleStandardGeneticCode["GCV"] = 'A'
middleStandardGeneticCode["GCN"] = 'A'


middleStandardGeneticCode["GAT"] = 'D'
middleStandardGeneticCode["GAC"] = 'D'
middleStandardGeneticCode["GAY"] = 'D'


middleStandardGeneticCode["GAA"] = 'E'
middleStandardGeneticCode["GAG"] = 'E'
middleStandardGeneticCode["GAR"] = 'E'


middleStandardGeneticCode["GGT"] = 'G'
middleStandardGeneticCode["GGC"] = 'G'
middleStandardGeneticCode["GGA"] = 'G'
middleStandardGeneticCode["GGG"] = 'G'

middleStandardGeneticCode["GGR"] = 'G'
middleStandardGeneticCode["GGY"] = 'G'
middleStandardGeneticCode["GGS"] = 'G'
middleStandardGeneticCode["GGW"] = 'G'
middleStandardGeneticCode["GGK"] = 'G'
middleStandardGeneticCode["GGM"] = 'G'
middleStandardGeneticCode["GGB"] = 'G'
middleStandardGeneticCode["GGD"] = 'G'
middleStandardGeneticCode["GGH"] = 'G'
middleStandardGeneticCode["GGV"] = 'G'
middleStandardGeneticCode["GGN"] = 'G'

middleStandardGeneticCode["---"] = '-'



# def sequence_to_matrix(sequence):
dna_to_number_dictionry = {
         #A     c     G   T
    "A": [1.00, 0.00, 0.00, 0.00],
    "C": [0.00, 1.00, 0.00, 0.00],
    "G": [0.00, 0.00, 1.00, 0.00],
    "T": [0.00, 0.00, 0.00, 1.00],
    "U": [0.00, 0.00, 0.00, 1.00],
    "R": [0.50, 0.00, 0.50, 0.00],
    "Y": [0.00, 0.50, 0.00, 0.50],
    "S": [0.00, 0.00, 0.50, 0.50],
    "W": [0.50, 0.00, 0.00, 0.50],
    "K": [0.00, 0.00, 0.50, 0.50],
    "M": [0.05, 0.50, 0.00, 0.00],
    "B": [0.00, 0.33, 0.33, 0.33],
    "D": [0.33, 0.00, 0.33, 0.33],
    "H": [0.33, 0.33, 0.00, 0.33],
    "V": [0.33, 0.33, 0.33, 0.00],
    "N": [0.25, 0.25, 0.25, 0.25],
}

def dna_to_matix(sequence):
    matrix = np.zeros([ len(sequence),4], float)
    for i in range(0,  len(sequence)):
        matrix[i] = dna_to_number_dictionry.get(sequence[i])
    return matrix

dna_to_number_dictionry2 = {
    # A     C     G   T   ## 5 means missing value here
    "A": 1,
    "C": 2,
    "G": 3,
    "T": 4,
    "U": 4,
    "R": 5,
    "Y": 5,
    "S": 5,
    "W": 5,
    "K": 5,
    "M": 5,
    "B": 5,
    "D": 5,
    "H": 5,
    "V": 5,
    "N": 5,
}


def dna_to_matix2(sequence):
    matrix = []
    for i in range(0,  len(sequence)):
        matrix.append(dna_to_number_dictionry2.get(sequence[i]))
    return matrix

def dna_to_matix3(sequence):
    matrix = []
    for i in range(0,  len(sequence)):
        matrix.append(sequence[i])
    return matrix

def getAllPossibleWithIupac(seq):
    allPossibleForEachChar = []
    for thisChar in seq:
        allPossibleForThisChar = []
        for key in dnaIupacCode:
            if thisChar in dnaIupacCode[key]:
                allPossibleForThisChar.append(key)

        allPossibleForEachChar.append(allPossibleForThisChar)

    currentPossibleCombinations = []
    for allPossibleForThisChar in allPossibleForEachChar:
        if len(currentPossibleCombinations) > 0 :
            lastPossibleCombinations = currentPossibleCombinations
            currentPossibleCombinations = []
            for thisChar in allPossibleForThisChar:
                for i in lastPossibleCombinations:
                    thisCombination = i+thisChar
                    currentPossibleCombinations.append(thisCombination)
        else:
            for thisChar in allPossibleForThisChar:
                currentPossibleCombinations.append(thisChar)

    return currentPossibleCombinations


def check_splice_sites(s1t, s2t):
    for index in range(0, len(dornors)):
        dornor = dornors[index]
        acceptor = acceptors[index]
        allAcceptors = getAllPossibleWithIupac(acceptor)
        allDornors = getAllPossibleWithIupac(dornor)
        if (s1t in allDornors) & (s2t in allAcceptors):
            return True
    return False

def ifSpliceSitesOk(targetTranscript, targetGenome):
    for i in range(1, len(targetTranscript.Cds)):
        let = 0
        tst = 0
        s1t = ""
        s2t = ""
        if (targetTranscript.strand == "+"):
            let = targetTranscript.Cds[i-1][1]
            tst = targetTranscript.Cds[i][0]
            s1t = FastaFile.getSubSequence(targetGenome, targetTranscript.chromosome_name, let + 1, let + 2, "+")
            s2t = FastaFile.getSubSequence(targetGenome, targetTranscript.chromosome_name, tst - 2, tst - 1, "+")
        else:
            let = targetTranscript.Cds[i - 1][0]
            tst = targetTranscript.Cds[i][1]
            s1t = FastaFile.getSubSequence(targetGenome, targetTranscript.chromosome_name, let - 2, let - 1, "-")
            s2t = FastaFile.getSubSequence(targetGenome, targetTranscript.chromosome_name, tst + 1, tst + 2, "-")

        if not check_splice_sites(s1t, s2t):
            return False
    return True

def ifLengthDivisibleByThree(sequence):
    return 0 == len(sequence) % 3

def ifNewStopCodon(cdsSequence):
    j = 0
    while j < len(cdsSequence)-3:
        threeNaInFrame = cdsSequence[j : j+3];
        if (threeNaInFrame in mustStopCodons):
            return True
        j = j + 3

    return False


def ifEndWithStopCodon(cdsSequence):
    threeNaInFrame = cdsSequence[len(cdsSequence)- 3:len(cdsSequence)]
    if threeNaInFrame in possibleStopCodons:
        return True
    else:
        return False

def ifStartWithStartCodon(cdsSequence):
    threeNaInFrame = cdsSequence[0: 3];
    if threeNaInFrame in possibleStartCodons:
        return True
    else:
        return False

def checkOrfState( targetTranscript, targetGenome):
    metaInformation = []
    if_orf_conserved = True
    if ifSpliceSitesOk(targetTranscript,targetGenome):
        metaInformation.append("_spliceSitesConserved")
    else:
        metaInformation.append("_spliceSitesDestroyed")
        if_orf_conserved = False

    cdsSequenceString = targetTranscript.cds_sequence
    if (len(cdsSequenceString) < 3):
        metaInformation.append("_exonLengthLessThan3")
        if_orf_conserved = False
    else:
        metaInformation.append( "_exonLengthMoreThan3")
        if (ifLengthDivisibleByThree(cdsSequenceString)):
            metaInformation.append("_exonLengthIsDivisibleBy3")
        else:
            metaInformation.append("_exonLengthIsNotMultipleOf3")
            if_orf_conserved = False

        if ifNewStopCodon(cdsSequenceString):
            metaInformation.append( "_prematureStopCodon")
            if_orf_conserved = False
        else:
            metaInformation.append("_noPrematureStopCodon")

        if ifEndWithStopCodon(cdsSequenceString):
            metaInformation.append("_endWithStopCodon")
        else:
            metaInformation.append( "_notEndWithStopCodon")
            if_orf_conserved = False

        if ifStartWithStartCodon(cdsSequenceString):
            metaInformation.append( "_startWithStartCodon")
        else:
            metaInformation.append("_notStartWithStartCodon")
            if_orf_conserved = False
    if( if_orf_conserved ):
        metaInformation.append("_ConservedFunction")

    metaInformation.append("_local_")
    metaInformation.append(str(targetTranscript.start))
    metaInformation.append("-")
    metaInformation.append(str(targetTranscript.end))
    if( "+" ==  targetTranscript.strand):
        metaInformation.append("_positive")
    else:
        metaInformation.append("_negative")

    meta_informaiton = ''.join(metaInformation)
    targetTranscript.meta_informaiton = meta_informaiton
    targetTranscript.if_orf_conserved = if_orf_conserved


def getReverseComplementary(sequence):
    reversecomplementary = []
    for j in range(0, len(sequence)):
        i = len(sequence) - 1 - j
        c = sequence[i]
        if ('A' == c):
            c = 'T'
        elif ('T' == c):
            c = 'A'
        elif ('U' == c):
            c = 'A'
        elif ('C' == c):
            c = 'G'
        elif ('G' == c):
            c = 'C'
        elif ('R' == c):
            c = 'Y'
        elif ('Y' == c):
            c = 'R'
        elif ('K' == c):
            c = 'M'
        elif ('M' == c):
            c = 'K'
        elif ('B' == c):
            c = 'V'
        elif ('V' == c):
            c = 'B'
        elif ('D' == c):
            c = 'H'
        elif ('H' == c):
            c = 'D'
        reversecomplementary.append(c)
    return ''.join(reversecomplementary)
