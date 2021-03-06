# -*- coding: utf-8 -*-
"""
Created on Jan.01.16

@author: ZiyiW
"""
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

from keras import models
from keras.models import load_model
import numpy as np
import h5py ##might need this

import pdb, os

def CNNpredict(ute):
    # use SegNet model to make segmentation
    mymodel = os.path.dirname(__file__)+'/myModel_utegrow_180201.h5'
    # mymodel = os.path.dirname(__file__)+'/myModel_utegrow_180403.h5'

    n_labels = 2
    # pdb.set_trace()
    img_h,img_w,img_d = np.shape(ute)

    mask = np.zeros(np.shape(ute))

    if(img_h!=128 or img_w!=128):
        raise Exception('Segmentation Image size should be 128 x 128 x n')

    # a deep matrix: in 3rd dimension, the label will be set to 1
    def de_label_map(myPred):
        myPred = np.reshape(myPred,[img_w,img_h,n_labels])
        return np.argmax(myPred,axis=2) # collapse the 3rd dimension

    autoencoder = models.Sequential()

    autoencoder=load_model(mymodel)

    for i in range(0,img_d):
        ute_slice = np.fliplr(np.flipud(ute[:,:,i]))

        ute_thre = np.percentile(ute_slice, 99)

        ute_slice = np.divide(ute_slice,ute_thre)

        ute_slice[ute_slice > 1] = 1

        ute_slice = np.multiply(ute_slice,255)

        # ute_slice = (ute_slice - 90.3463)/60.002483 # this is used for the newer model

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
