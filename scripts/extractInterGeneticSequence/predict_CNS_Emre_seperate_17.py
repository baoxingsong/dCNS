from __future__ import print_function
import pandas as pd
import numpy as np
import random

import keras

from keras.datasets import mnist
from keras import backend as K
from keras.utils import np_utils
from keras.layers import *
from keras.optimizers import Adam
from keras.models import Sequential,Model
from keras.callbacks import EarlyStopping
import os
from keras.models import Model


from hyperopt import tpe, hp, fmin
from random import shuffle
import gc
from sklearn.metrics import accuracy_score, classification_report


import tensorflow as tf
config = tf.ConfigProto(allow_soft_placement=True)
config.gpu_options.per_process_gpu_memory_fraction = 0.6
config.gpu_options.allow_growth = True
sess = tf.Session(config=config)
K.set_session(sess)



# We need some global variables. Because couldn't pass them to Hyperopt.

global epochs
global Y
global ValIndexes
global trainIn
global epochs
global mxLen
global basePairs
global allIn
global mySize


space = {'hidden_layer1': hp.choice('hidden_layer1', [4, 8, 16, 32, 64, 128]),
         'hidden_layer2': hp.choice('hidden_layer2', [4, 8, 16, 32, 64, 128]),
 'drop_out1': hp.choice('drop_out1', [0, 0.05, 0.1, 0.2, 0.3]),
 'drop_out2': hp.choice('drop_out2', [0, 0.05, 0.1, 0.2, 0.3]),
 'filter_size1': hp.choice('filter_size1', [4, 8, 16, 32, 64, 128]),
 'kernel_size1': hp.choice('kernel_size1', [1, 3, 5, 8, 16, 32]),
 'pool_size1': hp.choice('pool_size1', [3, 5, 8, 16, 32]),
 'strides1': hp.choice('strides1', [1, 3, 5, 8]),
 'batch_size': hp.choice('batch_size', [16, 32, 64, 128]),
 'beta1': hp.uniform('beta1', 0.9, 0.9999),
 'beta2': hp.uniform('beta2', 0.9, 0.9999),
 'learning_rate': hp.choice('learning_rate', [0.1, 0.01, 0.001, 0.0001])}





myspace = {'hidden_layer1': [4, 8, 16, 32, 64, 128],
           'hidden_layer2': [4, 8, 16, 32, 64, 128],
           'drop_out1':  [0, 0.05, 0.1, 0.2, 0.3],
           'filter_size1': [4, 8, 16, 32, 64, 128],
           'kernel_size1': [1, 3, 5, 8, 16, 32],
           'pool_size1': [3, 5, 8, 16, 32],
           'strides1':  [1, 3, 5, 8],
           'drop_out2': [0, 0.05, 0.1, 0.2, 0.3],
           'filter_size2':  [4, 8, 16, 32, 64, 128],
           'kernel_size2':  [1, 3, 5, 8, 16, 32],
           'pool_size2':  [3, 5, 8, 16, 32],
           'strides2':  [1, 3, 5, 8],
           'hidden_layer3':  [4, 8, 16, 32, 64, 128],
           'drop_out3':  [0, 0.05, 0.1, 0.2, 0.3],
           'filter_size3':  [4, 8, 16, 32, 64, 128],
           'kernel_size3':  [1, 3, 5, 8, 16, 32],
           'pool_size3':  [3, 5, 8, 16, 32],
           'strides3':  [1, 3, 5, 8],
           'hidden_layer4':  [4, 8, 16, 32, 64, 128],
           'drop_out4': [0, 0.05, 0.1, 0.2, 0.3],
           'batch_size':  [16, 32, 64, 128],
           'learning_rate':  [0.1, 0.01, 0.001, 0.0001]
           }


def getTensor(Data, whichColumn):
    # Takes data as pandas data frame
    # Returns selected column as tensor. Selected column should have RNA sequences. Also returns max lenght of the RNA
    SEQs = getColumn(Data, Col=whichColumn)
    mxlen = len(max(SEQs.values[:,0],key=len))
    # Convert SEQs data frame to CNN compatible one-hot encoded.
    SEQTensor = getOneHotEncodedSet(SEQs, mxlen=mxlen)
    return SEQTensor, mxlen

def getrainTestValidationTSplits(size, valRate, numberOfChunks):
    myInVals = np.arange(size)
    random.shuffle(myInVals)
    numberOfValSamples = int(valRate*size)
    ValSamples = myInVals[:numberOfValSamples]
    TrTstSamples=myInVals[numberOfValSamples:]
    TrTstSamples=np.array_split(TrTstSamples, numberOfChunks)
    return TrTstSamples, ValSamples

def createModel(mxLen, args):
    hidden_layer1=args['hidden_layer1']
    hidden_layer2=args['hidden_layer2']
    drop_out1=args['drop_out1']
    drop_out2=args['drop_out2']
    filter_size1=args['filter_size1']
    kernel_size1=args['kernel_size1']
    pool_size1=args['pool_size1']
    strides1=args['strides1']
    modelIn, modelOut = buildModel(mxLen, width=4, filterSize=filter_size1, kernelSize=kernel_size1, poolSize=pool_size1, denseSize1=hidden_layer1, denseSize2=hidden_layer2, dropoutRate1=drop_out1, dropoutRate2=drop_out2, strides=strides1)
    finalModel = Model(inputs=modelIn, outputs=modelOut)
    return finalModel

def readData(which):
    path0 = "/media/bs674/1_8t/AndCns/kmerClusterCns/c1"
    rawData=pd.read_csv(path0)
    indexes = range(rawData.shape[0])
    cnt=(np.mod(indexes,2)==0)
    X0= rawData.values[cnt]
    path1 = "/media/bs674/1_8t/AndCns/kmerClusterCns/c7"
    rawData=pd.read_csv(path1)
    indexes = range(rawData.shape[0])
    cnt=(np.mod(indexes,2)==0)
    X1= rawData.values[cnt]
    Y = np.zeros((len(X0)+len(X1),2))
    Y[0:len(X0),0]=1
    Y[len(X0):,1] = 1
    X = pd.concat([pd.DataFrame(X0), pd.DataFrame(X1)])
    return pd.DataFrame(X.values), Y


def concatDf(dfList):
    # Concats list of data frames into one data frame
    df = dfList[0]
    for l in range(1, len(dfList)):
        df = pd.concat([df, dfList[l]], axis=0)
    return pd.DataFrame(df.values)

def filterDistMatrix(myUniqueSpecies, DM):
    #May be deleted! Not in use yet.
    specListinDistM= DM.columns.values
    res=np.zeros(DM.shape[0])
    for s in myUniqueSpecies:
        if s in specListinDistM:
            res[np.where(s==specListinDistM)]=1
    res=np.array(res, dtype=bool)
    dumDM = DM.values[:,res]
    dumDM2 = dumDM[res,:]
    return specListinDistM[res], dumDM2

def getColumn(df, Col):
    # Takes a pandas data frame and returns a column
    res = df[df.columns[Col]]
    return pd.DataFrame(res.values)

def  getTarget(df,targetColumn):
    # Takes data frame and target column number
    # Returns target as array
    ogt = getColumn(df, Col=targetColumn).values
    ogt=np.array(ogt)
    OGT=[]
    # Not optimal. A better way will be added.
    for a in ogt:
        OGT.append(a[0])
    ogt= np.array(OGT)
    return ogt

def getOneHotEncodedSet(SEQs, mxlen, depth=4):
    # Takes sqquence list, and max lenght of sequences
    # Returns one hot encoded sequence list. Shorter sequences are 0 padded.
    one_hot_seqs = []
    for seq in SEQs.values[:, 0]:
        one_hot_seqs.append(one_hot_encoding(seq, mxlen))
    one_hot_seqs = np.array(one_hot_seqs)
    return  one_hot_seqs.reshape(one_hot_seqs.shape[0], depth, mxlen, 1)

def one_hot_encoding(seq, mx):
    # Takes a sequence and max lenght
    # Returns one hot encoded 2-dimensional array.
    dict = {'A': [1, 0, 0, 0], 'T': [0, 1, 0, 0], 'C': [0, 0, 1, 0], 'G': [0, 0, 0, 1]}    #maybe chage this dict and change the CNN stripe of vertial as 2
    one_hot_encoded = np.zeros(shape=(4, mx))
    for i, nt in enumerate(seq):
        if nt.upper() == "A" or nt.upper() == "C" or nt.upper() == "G" or nt.upper() == "T" :
            one_hot_encoded[:,i] = dict[nt.upper()]
        else:
            continue
    return one_hot_encoded


def buildModel(window_size, width=4, filterSize=0, kernelSize=0, poolSize=0, denseSize1=0, denseSize2=0, dropoutRate1=0,  dropoutRate2=0, strides=0):
    # Builds a convolutional branch and returns input and output node
    #output node is not sized 1! It's an array.
    mInput = Input(shape=(width, window_size, 1))
    #Convolution
    print("kernelSize:" + str(kernelSize))
    model = Conv2D(filterSize, kernel_size=(1, kernelSize),  padding='valid', activation='relu')(mInput)
    model = MaxPooling2D(pool_size=(1, poolSize),strides=(1, strides), padding='same')(model)
    # model = Conv2D(filterSize, kernel_size=(1, kernelSize), padding='same', activation='relu')(model)
    # model = MaxPooling2D(pool_size=(1, poolSize), strides=(1, strides),padding='same')(model)
    # model = Conv2D(filterSize, kernel_size=(1, kernelSize), padding='same', activation='relu')(model)
    # model = MaxPooling2D(pool_size=(1, poolSize), strides=(1, strides),padding='same')(model)
    # Dense layers
    model = Flatten()(model)
    # model = Dense(denseSize, activation='relu')(model)
    # model = Dropout(dropoutRate)(model)
    # model = Dense(denseSize, activation='relu')(model)
    # model = Dropout(dropoutRate)(model)
    # model = Dense(denseSize, activation='relu')(model)
    # model = Dropout(dropoutRate)(model)
    model = Dense(denseSize1, activation='relu')(model)
    model = Dropout(dropoutRate1)(model)
    model = Dense(denseSize2, activation='relu')(model)
    model = Dropout(dropoutRate2)(model)
    mOutput = Dense(2, activation='softmax')(model)
    return mInput, mOutput

def hyperParameterOptimization( max_eval):
    best_params = fmin(objective_func, space,
                           algo=tpe.suggest, max_evals=max_eval)
    return best_params

def objective_func(args):
    # objective function is called by hyperopt package. It takes args dictionary which contains CNN parameters
    # Returns a loss. In our case it is defined as 1 - r squared of validation set
    # Because this function can get ony one parameter we used some global variables.
    global Y
    global ValIndexes
    global trainIn
    global epochs
    global mxLen
    global basePairs
    global mySize
    # Get variables from dictionary
    batch_size = args['batch_size']
    beta1 = args['beta1']
    beta2 = args['beta2']
    learning_rate = args['learning_rate']
    try: # sometimes hyperopt parameter combinations are not competable with our CNN model. To avoid crash, we use try-except
        #create model
        mymodel = createModel(mxLen, args)
        adam = keras.optimizers.Adam(lr=learning_rate, beta_1=beta1, beta_2=beta2)
        mymodel.compile(loss="binary_crossentropy",optimizer=adam, metrics=['accuracy'])
        valTr= trainIn # set train species
        valTs= ValIndexes # set validation species
        # True-False array that shows training set
        allIn=np.arange(mySize)
        valCont = (allIn == valTr[0])
        for unqs in valTr:
            valCont = np.logical_or(valCont, (allIn == unqs))
        PR=[]
        REAL=[]
        mymodel.fit(basePairs[valCont],  Y[valCont,:],  batch_size=batch_size, epochs=epochs, verbose=1)
        # Get predictions on validation set.
        for v in valTs:
            specContTest = (allIn == v)  # True-False array that shows validation species
            tey = Y[specContTest]
            prediction = mymodel.predict(basePairs[specContTest])
            PR.append(np.argmax(prediction)) # Species OGT prediction is obtained as median of RNA related predictions
            REAL.append(np.argmax(tey))
        PR = np.array(PR)
        REAL = np.array(REAL)
        ac = accuracy_score(REAL, PR)
    except: # If training did not work, provides a huge loss
        ac=-100
    del mymodel
    gc.collect()
    K.clear_session()
    return 1-ac

def performCNN(args):
    # objective function is called by hyperopt package. It takes args dictionary which contains CNN parameters
    # Returns a loss. In our case it is defined as 1 - r squared of validation set
    # Because this function can get ony one parameter we used some global variables.
    global mxLen
    global Y
    global trainIn
    global epochs
    global mySize
    global basePairs
    global allIn
    global testIn
    # Get variables from dictionary
    batch_size = args['batch_size']
    beta1 = args['beta1']
    beta2 = args['beta2']
    learning_rate = args['learning_rate']
    #create model
    mymodel = createModel(mxLen, args)
    adam = keras.optimizers.Adam(lr=learning_rate, beta_1=beta1, beta_2=beta2)
    mymodel.compile(loss="binary_crossentropy",optimizer=adam, metrics=['accuracy'])
    valTr= trainIn # set train species
    valTs= testIn # set validation species
    # True-False array that shows training set
    allIn=np.arange(mySize)
    valCont = (allIn == valTr[0])
    for unqs in valTr:
        valCont = np.logical_or(valCont, (allIn == unqs))
    PR=[]
    REAL=[]
    mymodel.fit(basePairs[valCont],  Y[valCont,:],  batch_size=batch_size, epochs=epochs, verbose=1)
    # Get predictions on validation set.
    for v in valTs:
        specContTest = (allIn == v)  # True-False array that shows validation species
        tey = Y[specContTest]
        prediction = mymodel.predict(basePairs[specContTest])
        PR.append(np.argmax(prediction))  # Species OGT prediction is obtained as median of RNA related predictions
        REAL.append(np.argmax(tey))
    PR = np.array(PR)
    REAL = np.array(REAL)
    del mymodel
    gc.collect()
    K.clear_session()
    return PR, REAL

def convertHyperopt(best_parameters):
    # Hyperopt retuns index of the selected hyper parameters.
    # This function takes indexes and returns best hyper parameters as dictionary.
    convertedDict={}
    for i in range(len(best_parameters)):
        itm = best_parameters.popitem()
        if itm[0]== 'beta1' or itm[0]== 'beta2':
            ky = itm[0]
            convertedDict[ky] = itm[1]
        else:
            ky=itm[0]
            vlind=itm[1]
            myarray = myspace[ky]
            myval=myarray[vlind]
            convertedDict[ky]=myval
    return convertedDict

def writefile(P, R, S, pth):
    # P is predictions, R is real values, S is species array. Pth is path.
    for i in range(len(R)):
        with open(pth , "a") as myfile:
            mystr = S[i] +","+ str(R[i]) + "," + str(P[i]) + "\n"
            myfile.write(mystr)
    return

####################################   START     ##############################################################
####################################   INPUTS    ###############################################################
whichTF='ZmTF2'

valRate = 0.05  # 0 to 1. Percentage of validation set
numberofChunks = 5 # To divide data set into groups after validation set is excluded.
#(Then numberofChunks -1 groups are used for training, 1 group is used for  testing. This process is repeated numberofChunks times.)

epochs=30 # CNN epochs
max_eval=50 # number of trials in hyper parameter optimization






# Read data from files
Data, Y = readData(whichTF)


basePairs, mxLen = getTensor(Data, 0)
mySize=basePairs.shape[0]

TrTstIndxes, ValIndexes = getrainTestValidationTSplits(mySize, valRate, numberofChunks)




first = 1 # to check if it is first training. Needed not to use hyperopt in each group. do not change
AllPR=[] # to keep all predictions
AllREAL=[] # to keep all real values related to predictions


for i in range(numberofChunks): # for each chunk
    trainIn=[]
    testIn = TrTstIndxes[i] # test species = ith group in the TrTstSpecies

    for j in range(numberofChunks): # This loop adds species to training set except ith group in the TrTstSpecies
        if j !=i:
            for asp in TrTstIndxes[j]:
                trainIn.append(asp)

    if first: # if it is first run then best hyper parameters needed to be found. Call hyperoptimization
        best_parameters = hyperParameterOptimization(max_eval)
        convertedParams=convertHyperopt(best_parameters) # best_parameters are indices. Convert appropiate format.
        first=0 # First run is finished. We get best parameters. And, do not want to do heyperopt in the next iteration.

    # Train the model and get predictions.
    myprediction, myrealValues = performCNN(convertedParams)

    # Store predictions' their real values and species names. This loop is not efficient but does not effect the performance
    # that much. Can be better.
    for l in range(len(myprediction)):
        AllPR.append(myprediction[l])
        AllREAL.append(myrealValues[l])
    # free keras memory
    gc.collect()
    K.clear_session()

#Calculate performance
AllPR=np.array(AllPR).ravel()
AllREAL=np.array(AllREAL).ravel()

print(classification_report(AllREAL, AllPR))
print(convertedParams)
