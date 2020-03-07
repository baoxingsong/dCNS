#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 16:06:01 2016

@author: mm2842
"""
import sys
import csv
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import string
sys.path.append('/home/bs674/.local/lib/python3.6/site-packages/bhtsne/bh_tsne')
sys.path.append('/home/bs674/.local/lib/python3.6/site-packages/bhtsne/')
import bhtsne

from itertools import product
from collections import OrderedDict
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.decomposition import TruncatedSVD
from sklearn.utils import shuffle

'''To split a sequence into truncate into kmers'''
# seq = "cgagtctcacccaagagggttcaac"
# kmerlength = 5
# split_len(seq, kmerlength)
# ['cgagt', 'ctcac', 'ccaag', 'agggt', 'tcaac']
# TODO
def split_len(seq, kmerlength):
    return [''.join(x) for x in zip(*[list(seq[z::kmerlength]) for z in range(kmerlength)])]

# this only did complementary but did not do reverse
# kmer = "cgagtctcacccaagagggttcaac"
# complementarykmer(kmer)
# 'gctcagagtgggttctcccaagttg'
'''Reverse complementary'''
def complementarykmer(kmer):
    trans = str.maketrans('tagc', 'atcg')
    return kmer.translate(trans)

'''To slide throught a sequence and collect ALL the kmers'''
# kmer = "cgagtctcacccaagagggttcaac"
# kmersize = 4
#slide_string(kmer, kmersize)
# ['gtctcacccaagagggttcaac', 'agtctcacccaagagggttcaac', 'gagtctcacccaagagggttcaac']
# TODO
def slide_string(text,kmersize):
    n = kmersize-1
    sliced= []
    while n >= 1:
        sliced.append(text[n:])
        n=n-1
    return sliced

'''Map of all valid newtokens (no Ns)'''
# kmersize = 2
# createNewTokenIndexMap(kmersize)
# OrderedDict([('aa', 'aantt'), ('ac', 'acngt'), ('ag', 'agnct'), ('at', 'atnat'), ('ca', 'cantg'), ('cc', 'ccngg'), ('cg', 'cgncg'), ('ct', 'agnct'), ('ga', 'gantc'), ('gc', 'gcngc'), ('gg', 'ccngg'), ('gt', 'acngt'), ('ta', 'tanta'), ('tc', 'gantc'), ('tg', 'cantg'), ('tt', 'aantt')])
def createNewTokenIndexMap(kmersize):
    kmerIndexMap = OrderedDict()
    nucleotides = ["a", "c", "g", "t"]    
    kmerall = product(nucleotides, repeat=kmersize) #give all the possiable k-mers
    for i in kmerall:
        kmercomp = []
        kmer =''.join(i)
        kmercomp.append(kmer)
        kmercomprev = complementarykmer(kmer)[::-1]
        kmercomp.append(kmercomprev)
        kmerIndexMap[kmer] = "n".join(sorted(kmercomp)) # by sort aaa and ttt have the some result which is aaanttt
    return kmerIndexMap
    
'''Map of all valid kmers (no Ns)'''
# kmersize = 2
# createKmerIndexMap(kmersize)
# OrderedDict([('aa', 'aa'), ('ac', 'ac'), ('ag', 'ag'), ('at', 'at'), ('ca', 'ca'), ('cc', 'cc'), ('cg', 'cg'), ('ct', 'ct'), ('ga', 'ga'), ('gc', 'gc'), ('gg', 'gg'), ('gt', 'gt'), ('ta', 'ta'), ('tc', 'tc'), ('tg', 'tg'), ('tt', 'tt')])
def createKmerIndexMap(kmersize):
    kmerIndexMap = OrderedDict()
    nucleotides = ["a", "c", "g", "t"]    
    kmerall = product(nucleotides, repeat=kmersize)
    for i in kmerall:
        kmercomp = []
        kmer =''.join(i)
        kmercomp.append(kmer)
        kmerIndexMap[kmer] = kmer
    return kmerIndexMap


'''Collect from dictionaries seqs and tokenize them'''
# line = "cgagtctcacccaagagggttcaac"
# kmersize = 4
# kmerIndexMap = createNewTokenIndexMap(kmersize)
# seqToNewtokens(line, kmerIndexMap, kmersize)
# return 'cgagnctcg gagantctc acccngggt aagantctt acccngggt tcaanttga actcngagt ctcantgag cccantggg agagnctct aaccnggtt caacngttg agtcngact gtgantcac ccaanttgg cctcngagg gaacngttc agacngtct caccnggtg caagncttg agggnccct tgaanttca'
def seqToNewtokens(line, kmerIndexMap, kmersize):
    line = str.lower(line)
    if len(line) > 0:
        sentences = [ ]
        kmerlist = split_len(line, kmersize)
        for kmer in kmerlist:
            if "n" not in kmer:
                sentences.append(kmerIndexMap.get(kmer))
        for n in range(1,kmersize):
            kmerlist = split_len(line[n:], kmersize)
            for kmer in kmerlist:
                if "n" not in kmer:
                    sentences.append(kmerIndexMap.get(kmer))
    if not sentences:
        return "" #empty string
    else:
        return " ".join(sentences) #newtokens separated with spaces    

'''Collect from dictionaries seqs and tokenize them'''
def collectTokens(dictseqs, IndexMap, kmersize):
    listnewtokens = []
    for each in dictseqs:
        line = each['seq'].strip()
        line = str.lower(line)
        if len(line) > 0:
            sentences = [ ]
            kmerlist = split_len(line, kmersize)
            for kmer in kmerlist:
                if "n" not in kmer:
                    sentences.append(IndexMap.get(kmer))
            for n in range(1,kmersize):
                kmerlist = split_len(line[n:], kmersize)
                for kmer in kmerlist:
                    if "n" not in kmer:
                        sentences.append(IndexMap.get(kmer))
        if not sentences:
            continue
        else:
            listnewtokens.append(" ".join(sentences)) #newtokens separated with spaces
    return listnewtokens


'''Try-Catch block to test integer cast'''
def representsInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False


'''Process a file and return a list with dictionaries'''
def returnListDict(path, color):
    print ("processing: " + path)
    listresults = []
    with open(path) as bedfile:
        file_reader = csv.reader(bedfile, delimiter="\t")
        for line in file_reader:
            seq = str(line[1])
            label = str(line[0])
            seq_length = len(seq)
            # print("seq_length:" + str(seq_length))
            listresults.append({'seq':seq,'label':label,'color':color, 'seq_length':seq_length})
    return listresults


def plotTSNE(embedded_coords, embedded_colors, outfile, zz_target):
    fig = plt.figure(figsize=(10, 10))
    ax = plt.axes(frameon=False)
    plt.setp(ax, xticks=(), yticks=())
    plt.subplots_adjust(left=0.0, bottom=0.0, right=1.0, top=0.9,wspace=0.0, hspace=0.0)
    plt.scatter(embedded_coords[:, 0],embedded_coords[:, 1],c=embedded_colors,marker=".",alpha=0.1)
    plt.legend(handles=patches_array, loc='lower right')
    fig.savefig(outfile)
    plt.close(fig)
    f= open(outfile + ".txt","w+")
    for i in range(len(embedded_colors)):
        f.write( str(embedded_coords[i, 0]) + "\t" + str(embedded_coords[i, 1]) + "\t" + str(embedded_colors[i]) + "\t" + str(zz_target[i]) + "\n")
    f.close()



'''Collection define a container object for datasets: Dictionary-like object that exposes its keys as attributes'''
# This is the Bunch class from scikit-learn
# class renamed to avoid conflict with the scikit-learn code
class Collection(dict):
    def __init__(self, **kwargs):
        super(Collection, self).__init__(kwargs)
    def __setattr__(self, key, value):
        self[key] = value
    def __dir__(self):
        return self.keys()
    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError:
            raise AttributeError(key)
    def __setstate__(self, state):
        pass


print ("Omnia per numeros...")


file = ["./sequences_with_tag1",
        "./sequences_with_tag2",
        "./sequences_with_tag3",
        "./sequences_with_tag4",
        "./sequences_with_tag5",
        "./sequences_with_tag6",
        "./sequences_with_tag7"]




file = ["./chr8_sequences_with_tag_1",
        "./chr8_sequences_with_tag_3",
        "./chr8_sequences_with_tag_5",
        "./chr8_sequences_with_tag_7",
        "./chr8_sequences_with_tag_17",
        "./chr8_sequences_with_tag_19",
        "./chr8_sequences_with_tag_21",
        "./chr8_sequences_with_tag_23",
        "./chr8_sequences_with_tag_TE"]

regions = {}
color_vector = ['r','r','r','r', 'g','r','r','r', 'k'] #use this vector for non-genic regions - repeats black

i = 0
for filepath in file:
    regions[i]=returnListDict(filepath,color_vector[i])
    i = i + 1

kmersize_range = [5,2,8,3,4,6,7]

for kmersize in kmersize_range:
    print ( "kmersize: %d" % kmersize )
    if kmersize < 4:
        kmers = True
        print ( "Tokenizing data using kmers..." )
        indexMap = createKmerIndexMap(kmersize) ## without reverse complementary
    else:
        kmers = False
        print ( "Tokenizing data using newtokens..." )
        indexMap = createNewTokenIndexMap(kmersize) # with reverse complementary
    samples = len(regions)    
    print ("Printing labels and colors for %d samples ..." % samples )
    x = []; color_array = []; patches_array = []; target_array = []; seq_length_array = []
    for i,region in regions.items():
        patches_array.append(mpatches.Patch(color=region[i].get('color'), label=region[i].get('label')))
        tokenized_region = collectTokens(region,indexMap,kmersize)
        x += tokenized_region
        target_array.append(np.full((len(tokenized_region)), i, dtype=np.float))
        color_array.append(np.full((len(tokenized_region)), region[i].get('color'), dtype=np.string_))
        for re in regions[i]:
            seq_length_array.append(re['seq_length'])
        print ( '%d regions with %s color and %s label' % (len(tokenized_region),region[i].get('color'),region[i].get('label')) )
    y = np.concatenate(color_array)
    z = np.concatenate(target_array)
    zz = np.array(seq_length_array)

    X_data, y_colours, y_target, zz_target = shuffle(x, y, z, zz, random_state=42)
    print ( "Extracting features from the training data using TfidfVectorizer")
    vectorizer = TfidfVectorizer(min_df = 1, max_df = 1.0, sublinear_tf=True,use_idf=True) #vectorizer for kmer frequencies
    X_TFIDF = vectorizer.fit_transform(X_data) #byebye to the sequences - welcome to the index.
    kmer_names = vectorizer.get_feature_names()
    print ( '%d tokens in the vocabulary with k-size %d' % (len(kmer_names), kmersize) )
    print ( 'processing matrix with n_sampes %d and n_features %d' %  X_TFIDF.shape )
    if kmersize <= 4:
        maxima = True
        reducto = False
    else:
        maxima = False
        reducto = True
    if reducto:
        n_cmpnnts = 150
        svd = TruncatedSVD(n_components=n_cmpnnts, random_state=42)
        Xred_to_bhtsne = svd.fit_transform(X_TFIDF)
        np.savetxt(str(kmersize)+"svd.csv", Xred_to_bhtsne, delimiter="\t")
        print ( "Reducto! TruncatedSVD to bhtsne using n_components %d" %  n_cmpnnts )
        print ( svd.explained_variance_ratio_.sum() )
    if maxima:
        print ( "Maxima! All the matrix to bhtsne..." )
        Xmax_to_bhtsne=X_TFIDF.toarray()
    ndms = 2
    tht=0.5
    rndsd=42
    vrbs=True
    prplxt_range = [30,50]
    for prplxt in prplxt_range:
        print ('bhtsne params: perplexity %d and n_dims %d' %  (prplxt, ndms) )
        if reducto:
            X_embedded = np.asarray(list(bhtsne.run_bh_tsne(Xred_to_bhtsne, no_dims=ndms, perplexity=prplxt, theta=tht, randseed=rndsd, verbose=vrbs)))
            print ("Plotting Reducto...")
            if kmers:
                file_fig = str(kmersize)+"-kmersize_kmers_Reducto_Reg_True_length_"+str(prplxt)+"-per.png"
            else:
                file_fig = str(kmersize)+"-kmersize_newtoken_Reducto_Reg_True_length_"+str(prplxt)+"-per.png"
            file_fig = os.path.join("./",file_fig)
            y_c=[]
            for x in y_colours:
                y_c.append(str(x)[2])
            plotTSNE(X_embedded,y_c,file_fig, zz_target)
        if maxima:
            X_embedded = np.asarray(list(bhtsne.run_bh_tsne(Xmax_to_bhtsne, no_dims=ndms, perplexity=prplxt, theta=tht, randseed=rndsd, verbose=vrbs)))
            print ("Plotting Maxima...")
            if kmers:
                file_fig = str(kmersize)+"-kmersize_kmers_Maxima_Reg_True_length_"+str(prplxt)+"-per.png"
            else:
                file_fig = str(kmersize)+"-kmersize_newtoken_Maxima_Reg_True_length_"+str(prplxt)+"-per.png"
            file_fig = os.path.join("./",file_fig)
            y_c=[]
            for x in y_colours:
                y_c.append(str(x)[2])
            plotTSNE(X_embedded,y_c,file_fig, zz_target)
