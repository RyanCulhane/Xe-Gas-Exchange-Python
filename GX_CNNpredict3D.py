# -*- coding: utf-8 -*-
"""
Created on Apirl.05.18

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
from GX_Map_utils import fullMontage
from GX_Benchmark_utils import timerfunc
import scipy.io as sio

fdata = '/home/ziyiw/Patients/005027/005027.mat'
mat_input = sio.loadmat(fdata)
ute = mat_input['ute']

img_mean = 79.43
img_std = 66.01

@timerfunc
def CNNpredict_3D(ute):
    # use SegNet model to make segmentation
    mymodel = os.path.dirname(__file__)+'/myModel_utegrow_180403.h5'

    n_labels = 2
    # pdb.set_trace()
    img_h,img_w,img_d = np.shape(ute)

    mask = np.zeros(np.shape(ute))
    mask1 = np.zeros(np.shape(ute))
    mask2 = np.zeros(np.shape(ute))

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
        ute_slice1 = np.fliplr(np.flipud(ute[i,:,:]))
        ute_slice2 = np.fliplr(np.flipud(ute[:,i,:]))
        # pdb.set_trace()

        ute_thre = np.percentile(ute_slice, 99)
        ute_slice = np.divide(ute_slice,ute_thre)
        ute_slice[ute_slice > 1] = 1
        ute_slice = np.multiply(ute_slice,255)

        ute_thre = np.percentile(ute_slice1, 99)
        ute_slice1 = np.divide(ute_slice1,ute_thre)
        ute_slice1[ute_slice1 > 1] = 1
        ute_slice1 = np.multiply(ute_slice1,255)

        ute_thre = np.percentile(ute_slice2, 99)
        ute_slice2 = np.divide(ute_slice2,ute_thre)
        ute_slice2[ute_slice2 > 1] = 1
        ute_slice2 = np.multiply(ute_slice2,255)

        ute_slice = (ute_slice - img_mean)/img_std
        ute_slice1 = (ute_slice1 - img_mean)/img_std
        ute_slice2 = (ute_slice2 - img_mean)/img_std
        # pdb.set_trace()
        mask_slice = autoencoder.predict(np.reshape(ute_slice,(1,img_w,img_h,1)))
        mask_slice = de_label_map(mask_slice)

        mask_slice1 = autoencoder.predict(np.reshape(ute_slice1,(1,img_w,img_h,1)))
        mask_slice1 = de_label_map(mask_slice1)

        mask_slice2 = autoencoder.predict(np.reshape(ute_slice2,(1,img_w,img_h,1)))
        mask_slice2 = de_label_map(mask_slice2)

        mask[:,:,i] = np.fliplr(np.flipud(mask_slice))
        mask1[i,:,:] = np.fliplr(np.flipud(mask_slice1))
        mask2[:,i,:] = np.fliplr(np.flipud(mask_slice2))

        print('finished '+str(i))

        # plt.subplot(1,2,1)
        # plt.imshow(mask_slice, origin="upper",cmap='grays')
        # plt.subplot(1,2,2)
        # plt.imshow(ute_slice, origin="upper",cmap='gray')

        # saveStr='Segment_'+str(i)+'.png'
        # plt.savefig(saveStr)

    mask3d = mask+mask1+mask2
    mask3d[mask3d<2] = 0
    mask3d[mask3d>0] = 1
    fullMontage(mask3d)
    pdb.set_trace()
    mask3d = mask3d.astype(bool)

    return mask3d

mask = CNNpredict_3D(ute)
pdb.set_trace()
