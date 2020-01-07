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



import tensorflow as tf
config = tf.ConfigProto(allow_soft_placement=True)
config.gpu_options.per_process_gpu_memory_fraction = 0.2
# config.gpu_options.allow_growth = True
sess = tf.Session(config=config)
K.set_session(sess)



########################################################################################################
# read DNA sequences
########################################################################################################
LEN_TAR=200

NTS=dict([["A",[1,0,0,0]],["C",[0,1,0,0]],["G",[0,0,1,0]],["T",[0,0,0,1]],["X",[0,0,0,0]]]) 
def seq2num(seq,LEN_TAR=LEN_TAR):
   if len(seq)<LEN_TAR:
      seq=seq+"X"*(LEN_TAR-len(seq))
   else:
      seq=seq[:LEN_TAR]
   nums=np.array([NTS[x] for x in seq])
   return nums
   
pos_seqs=[]
neg_seqs=[]

handle=open('toML','r')

for line in handle:
   x=seq2num(line.split('\t')[0])
   y=seq2num(line.split('\t')[1])
   pos_seqs.append(x)
   neg_seqs.append(y)

pos_seqs=np.array(pos_seqs)
neg_seqs=np.array(neg_seqs)

#########################################################################################
# train/test splitting
#########################################################################################

datalen=len(pos_seqs)
train_indices=random.sample(range(datalen),int(datalen*0.8))
test_indices=list(set(range(datalen))-set(train_indices))

pos_train = pos_seqs[train_indices]
pos_test = pos_seqs[test_indices]
neg_train = neg_seqs[train_indices]
neg_test = neg_seqs[test_indices]

predictors_train=np.concatenate([pos_train,neg_train])
predictors_test=np.concatenate([pos_test,neg_test])

targets_train=np.concatenate((np.repeat(1,len(pos_train)),np.repeat(0,len(neg_train))))
targets_test=np.concatenate((np.repeat(1,len(pos_test)),np.repeat(0,len(neg_test))))
targets_train=np_utils.to_categorical(targets_train,2)
targets_test=np_utils.to_categorical(targets_test,2)

#################################################################################################################
# model training
#################################################################################################################

def build():
   model=Sequential()

   model.add(Conv1D(128,kernel_size=6,padding='valid',activation='relu',input_shape=(LEN_TAR,4)))
   model.add(Conv1D(128,kernel_size=6,padding='valid',activation='relu'))
   model.add(MaxPooling1D(pool_size=3,padding='valid'))
   model.add(Dropout(0.1))

   model.add(Conv1D(128,kernel_size=6,padding='same',activation='relu'))
   model.add(Conv1D(128,kernel_size=6,padding='same',activation='relu'))
   model.add(MaxPooling1D(pool_size=3,padding='valid'))
   model.add(Dropout(0.1))

   model.add(Flatten())
   model.add(Dense(64,activation='relu'))
   model.add(Dropout(0.1))
   model.add(Dense(32,activation='relu'))
   model.add(Dropout(0.1))
   model.add(Dense(2,activation='softmax'))

   return model

callbacks=[EarlyStopping(monitor='val_loss',patience=3)]

model=build()
model.compile(loss='binary_crossentropy',optimizer='adam',metrics=['accuracy'])

model.fit(x = predictors_train,
          y = targets_train,
          validation_data = (predictors_test,targets_test),
          batch_size=512,
          epochs=100,
          shuffle=True,
          callbacks=callbacks)

accuracy=model.evaluate(x=predictors_test,y=targets_test)[1]
print (accuracy)

