#  Ziyi, Jan 31. 2018
#
# This file decomp the dissolved-phase image using 1-point Dixon based method
#
# angle: in radians

import numpy as np
import nibabel as nib
import scipy.io as sio
import SimpleITK as sitk
import warnings
from matplotlib import pyplot as plt

import pdb

## Input **************************
fmask = 'BHUTE_Sub002102_FID49886_mask_grow_reg_vent.nii'
fdata = 'Sub002102_data.mat'

mat_input = sio.loadmat(fdata)
gas = mat_input['gasVol_highSNR']
dissolved = mat_input['dissolvedVol']

meanRbc2barrier = 0.457
mask = np.array(nib.load(fmask).get_data(),dtype='bool')

# check if the directions of images are the same
# pdb.set_trace()
# plt.subplot(321)
# plt.imshow(mask[:,:,55])
# plt.subplot(323)
# plt.imshow(mask[:,:,64])
# plt.subplot(325)
# plt.imshow(mask[:,:,80])
# plt.subplot(322)
# plt.imshow(abs(dissolved[:,:,55]))
# plt.subplot(324)
# plt.imshow(abs(dissolved[:,:,64]))
# plt.subplot(326)
# plt.imshow(abs(dissolved[:,:,80]))
# plt.show()

## apply delta angfe from RBC:barrier ******************
desired_angle = np.arctan2(meanRbc2barrier,1)

netVet_lung = np.sum(dissolved[mask])
current_angle = np.arctan2(np.imag(netVet_lung), np.real(netVet_lung))
delta_angle = desired_angle - current_angle

rotVol = np.multiply(dissolved,np.exp(1j*delta_angle))

## Correct B0 inhomogeneities
iterCount = 0
meanphase = 1 # arbitrary number

# shift the gas to have zero mean phase
while(abs(meanphase) > 1E-7):
    if(iterCount > 20): warnings.warn('can not converge in 20 iterations')
    diffphase = np.angle(gas)
    meanphase = np.mean(diffphase[mask])
    gas = np.multiply(gas,np.exp(-1j*meanphase))

diffphase = - np.angle(gas)

# apply difference phase
rotVol_B = np.multiply(rotVol,np.exp(1j*diffphase))

rbc = np.imag(rotVol_B)
barrier = np.real(rotVol_B)

sitk_rbc = sitk.GetImageFromArray(rbc)
sitk.WriteImage(sitk_rbc,"rbc.nii")
sitk_barrier = sitk.GetImageFromArray(barrier)
sitk.WriteImage(sitk_barrier,"barrier.nii")
# pdb.set_trace()
