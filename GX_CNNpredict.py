# -*- coding: utf-8 -*-
"""
Created on Jan.01.16

@author: ZiyiW
"""
from keras import models
from keras.models import load_model
import numpy as np
import h5py ##might need this
from matplotlib import pyplot as plt
import pdb

def CNNpredict(ute):
    # use SegNet model to make segmentation
    mymodel = 'myModel_utegrow_180201.h5'

    n_labels = 2
    # pdb.set_trace()
    img_h,img_w,img_d = np.shape(ute)

    mask = np.zeros(np.shape(ute))

    if(img_h!=128 or img_w!=128):
        raise Exception('Segmentation Image size should be 128 x 128 x n')

    # a deep matrix: in 3rd dimension, the label will be set to 1
    def de_label_map(myPred):
        myPred = np.reshape(myPred,[1,img_w,img_h,n_labels])
        return_map = np.zeros([img_w,img_h])
        for i in range(0,np.shape(myPred)[1]):
            for j in range(0,np.shape(myPred)[2]):
                return_map[i][j]=np.argmax(myPred[0,i,j,:])
        return return_map

    autoencoder = models.Sequential()

    autoencoder=load_model(mymodel)

    for i in range(0,img_d):
        ute_slice = np.fliplr(np.flipud(ute[:,:,i]))

        ute_thre = np.percentile(ute_slice, 99)

        ute_slice = np.divide(ute_slice,ute_thre)

        ute_slice[ute_slice > 1] = 1

        ute_slice = np.multiply(ute_slice,255)

        # pdb.set_trace()
        mask_slice = autoencoder.predict(np.reshape(ute_slice,(1,img_w,img_h,1)))
        mask_slice = de_label_map(mask_slice)

        mask_slice = np.fliplr(np.flipud(mask_slice))

        mask[:,:,i] = mask_slice

        # plt.subplot(1,2,1)
        # plt.imshow(mask_slice, origin="upper",cmap='grays')
        # plt.subplot(1,2,2)
        # plt.imshow(ute_slice, origin="upper",cmap='gray')

        # saveStr='Segment_'+str(i)+'.png'
        # plt.savefig(saveStr)

    mask = mask.astype(bool)

    return mask
