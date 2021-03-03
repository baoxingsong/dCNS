"""
Created: 25/12/2017  00:40:39
Author: Baoxing Song
Email: songbaoxing168@163.com
This source code is partially adopted from https://github.com/bvilhjal/mixmogam/
And the R source code of emma is referred

modified: 31 Oct 2019
By Baoxing Song try to take the CNS present/absent as genotypic variant for GWAS analysis using a fixed model
about the f test implemented here https://towardsdatascience.com/fisher-test-for-regression-analysis-1e1687867259
"""

import numpy as np
from numpy import linalg
from scipy import stats
import warnings
import re
import statistics
def parse_cns_genotype_file(filename, absent_threshold=0.2, present_threshold=0.8, maf_count=5, min_total_count = 17):  # 25 peers and 3 pcs
    """
        read the cns length and coverage file and return a inter genotype table
    """
    individs = []
    with open(filename) as f:
        line = f.readline()
        l = list(map(str.strip, line.split()))
        for i in range(2, len(l)):
            individs.append(l[i])
    num_individs = len(individs)
    all_snps = {'chrs': [], 'positions': [],  'ends':[], 'snps': []}
    with open(filename) as f:
        for line_i, line in enumerate(f):
            if line_i > 0:
                l = list(map(str.strip, line.split()))
                m = re.search('^(\d+):(\d+)\-(\d+)', l[0])
                if (m != None):
                    chrom = int(m.group(1))
                    snp = np.zeros(num_individs, dtype='float')
                    cns_length = float(l[1])
                    an = 0 # absent count
                    pn = 0 # present count
                    for i in range(2, num_individs+2, 1):
                        nt = float(l[i])/float(cns_length)
                        if nt <= absent_threshold:
                            snp[i-2] = 0.0
                            an = an + 1
                        elif nt >=present_threshold:
                            snp[i-2] = 1.0
                            pn = pn + 1
                        else:
                            snp[i-2] = 3.0 #missing
                    all_snps['chrs'].append(chrom)
                    all_snps['positions'].append(int(m.group(2)))
                    all_snps['ends'].append(int(m.group(3)))
                    all_snps['snps'].append(snp)
    return all_snps, individs



def parse_cns_genotype_file_proportion(filename):
    """
        read the cns length and coverage file and return a inter genotype table
    """
    individs = []
    with open(filename) as f:
        line = f.readline()
        l = list(map(str.strip, line.split()))
        for i in range(2, len(l)):
            individs.append(l[i])
    num_individs = len(individs)
    all_snps = {'chrs': [], 'positions': [],  'ends':[], 'snps': []}
    with open(filename) as f:
        for line_i, line in enumerate(f):
            if line_i > 0:
                l = list(map(str.strip, line.split()))
                m = re.search('^(\d+):(\d+)\-(\d+)', l[0])
                if (m != None):
                    chrom = int(m.group(1))
                    snp = np.zeros(num_individs, dtype='float')
                    cns_length = float(l[1])
                    for i in range(2, num_individs+2, 1):
                        nt = float(l[i])/float(cns_length)
                        snp[i-2] = nt
                    all_snps['chrs'].append(chrom)
                    all_snps['positions'].append(int(m.group(2)))
                    all_snps['ends'].append(int(m.group(3)))
                    all_snps['snps'].append(snp)
    return all_snps, individs


def parse_gene_present_absent_file(filename, absent_threshold=0.4):
    """
        read the gene length and coverage file and return a inter gene present absent table
    """
    individs = {} # taxa name to index
    with open(filename) as f:
        line = f.readline()
        l = list(map(str.strip, line.split()))
        for i in range(3, len(l)):
            individs[l[i]] = i-3
    num_individs = len(individs)
    genes = {}
    all_genes = {'chrs': [], 'positions': [], 'ends': [], 'genes': []}
    with open(filename) as f:
        for line_i, line in enumerate(f):
            if line_i > 0:
                l = list(map(str.strip, line.split()))
                m = re.search('^(\S+):(\d+)\-(\d+)', l[1])
                genes[l[0]] = line_i-1
                if (m != None):
                    chrom = m.group(1)
                    gene = np.zeros(num_individs, dtype='float')
                    gene_length = float(l[2])
                    for i in range(3, num_individs+3, 1):
                        nt = float(l[i])/float(gene_length)
                        if nt <= absent_threshold:
                            gene[i-3] = 0.0
                        else:
                            gene[i-3] = 1.0
                    all_genes['chrs'].append(chrom)
                    all_genes['positions'].append(int(m.group(2)))
                    all_genes['ends'].append(int(m.group(3)))
                    all_genes['genes'].append(gene)
    return all_genes, individs, genes

#read genotype data
genopype, individs = parse_cns_genotype_file("../CNS") #the document for generating this file has been released at https://github.com/baoxingsong/CNSpublication/blob/master/CNS_analysis/howToUseCNSABVtoperformGWAS.html
genopype_proportion, individs_proportion = parse_cns_genotype_file_proportion("../CNS") #the document for generating this file has been released at https://github.com/baoxingsong/CNSpublication/blob/master/CNS_analysis/howToUseCNSABVtoperformGWAS.html


# get the genotypic variants taxa id to line number map. line number could be queried from the genotypic variants
print("cns reading done")
gene_present_absent, gene_present_absent_individs, gene_present_absent_genes = parse_gene_present_absent_file("../CNS")# here could be modified to exclude absent gene, here we are not doing that
print("gene reading done")

import sys
maf_count = 15
min_total_count = 35
f = open(sys.argv[1], 'rU') # this input is the gene expression level matrix obtained from https://www.nature.com/articles/nature25966
line_index = 0
phen_individs = {} # key is the taxa name and value is the the index
phen_dict = {}  #key is the phenotype id, value is a vector of phenotype
for line in f:
    if 0 == line_index:
        l = line.split()
        for i in range(0, len(l)):
            phen_individs[l[i]] = i
    else:
        d = []
        l = line.split()
        for i in range(1, len(l)):
            d.append(float(l[i]))
        pid = l[0]
        phen_dict[pid] = d
    line_index = line_index + 1

f.close()
#read phenotype data end
print("genotype reading done")
# get the PCS taxa id to line number map. line number could be queried from the genotypic variants
pcs = np.mat(np.loadtxt("../282set_PCs.txt", skiprows=1, usecols = (2,3,4))).astype("float") # this is the pcs matrix obtained from https://www.nature.com/articles/nature25966
pcs_id = {} # key is the taxa name and value is the id of taxa in the matrix
f1 = open("../282set_PCs.txt", 'rU')
i = 0
for line in f1:
    if i >0:
        l = list(map(str.strip, line.split()))
        pcs_id[l[0]] = i-1 # taxa id and the index of taxa in the matrix
    i = i+1

f1.close()
#pca reading end

peers = np.mat(np.loadtxt(sys.argv[2], skiprows=1, usecols = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25))).astype("float")   # here is the PEERs factor matrix obtained from https://www.nature.com/articles/nature25966
peers_id = {}  # key is the taxa name and value is the id of taxa in the matrix
f1 = open(sys.argv[2], 'rU')
i = 0
for line in f1:
    if i >0:
        l = list(map(str.strip, line.split()))
        peers_id[l[0].replace("282set_", "")] = i - 1 # taxa id and the index of taxa in the matrix
    i = i+1
f1.close()
print ("begein association analysis")
lin_depend_thres=1e-4


used_pav = set()

f = open("cns_gwas_forplot_small", 'rU') # this input is the gene expression level matrix obtained from https://www.nature.com/articles/nature25966
for line in f:
    l = list(map(str.strip, line.split()))
    gene = l[4]
    for snp_index in range(len(genopype['snps'])):
        if snp_index not in used_pav:
            if str(genopype['chrs'][snp_index]) == l[0] and str(genopype['positions'][snp_index]) == l[1] and str(genopype['ends'][snp_index]) == l[2] and gene in phen_dict:
                used_pav.add(snp_index)
                p = phen_dict[gene]
                an_phenotype = []
                pn_phenotype = []
                for i in range(len(genopype['snps'][snp_index])):
                    if (genopype['snps'][snp_index][i] != 3.0) and (individs[i] in pcs_id) and (individs[i] in peers_id) \
                            and (individs[i] in phen_individs) and (not (genopype['snps'][snp_index][i] == 0.0 \
                            and (individs[i] in gene_present_absent_individs) and gene in gene_present_absent_genes \
                            and gene_present_absent['genes'][gene_present_absent_genes[gene]][gene_present_absent_individs[individs[i]]] == 0.0)):
                        if genopype['snps'][snp_index][i] ==0:
                            an_phenotype.append(p[phen_individs[individs[i]]])
                        else:
                            pn_phenotype.append(p[phen_individs[individs[i]]])

                an_median = statistics.median(an_phenotype)
                pn_median = statistics.median(pn_phenotype)

                for i in range(len(genopype_proportion['snps'][snp_index])):
                    if individs[i] in phen_individs and abs(p[phen_individs[individs[i]]] - an_median) < abs(p[phen_individs[individs[i]]] - pn_median) :
                        print ("absence\t" + gene + "\t" + str(genopype['chrs'][snp_index]) + ":" + str(genopype['positions'][snp_index]) + "-" + str(genopype['ends'][snp_index])  + "\t" + str(genopype_proportion['snps'][snp_index][i]))
                    elif individs[i] in phen_individs:
                        print ("presence\t" + gene + "\t" + str(genopype['chrs'][snp_index]) + ":" + str(genopype['positions'][snp_index]) + "-" + str(genopype['ends'][snp_index])  + "\t" + str(genopype_proportion['snps'][snp_index][i]))
