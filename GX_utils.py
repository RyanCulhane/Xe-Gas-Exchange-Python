import numpy as np
import nibabel as nib
import scipy.io as sio
import SimpleITK as sitk
import warnings
from matplotlib import pyplot as plt

import pdb


def binning(volume,thresholds):

    # volume: mask_vented, thresholded 3D volume
    # thresholds: just the middle thresholds

    bvolume = np.ones(np.shape(volume))

    bvolume[(volume > 0) & (volume <= thresholds[0])] = 2

    for k in range(len(thresholds)-2):
        bvolume[(volume > thresholds[k+1]) & (volume <= thresholds[k+2])] = k+3

    bvolume[volume > thresholds[-1]] = len(thresholds)+2

    return bvolume

def dixonDecomp(gas_highSNR,dissolved,mask_vent,meanRbc2barrier):

    ## apply delta angfe from RBC:barrier ******************
    desired_angle = np.arctan2(meanRbc2barrier,1)

    netVet_lung = np.sum(dissolved[mask_vent])
    current_angle = np.arctan2(np.imag(netVet_lung), np.real(netVet_lung))
    delta_angle = desired_angle - current_angle

    rotVol = np.multiply(dissolved,np.exp(1j*delta_angle))

    ## Correct B0 inhomogeneities
    iterCount = 0
    meanphase = 1 # arbitrary number

    # shift the gas_highSNR to have zero mean phase
    while(abs(meanphase) > 1E-7):
        if(iterCount > 20): warnings.warn('can not converge in 20 iterations')
        diffphase = np.angle(gas_highSNR)
        meanphase = np.mean(diffphase[mask_vent])
        gas_highSNR = np.multiply(gas_highSNR,np.exp(-1j*meanphase))

    diffphase = - np.angle(gas_highSNR)

    # apply difference phase
    rotVol_B = np.multiply(rotVol,np.exp(1j*diffphase))

    rbc = np.imag(rotVol_B)
    barrier = np.real(rotVol_B)

    # sitk_rbc = sitk.GetImageFromArray(rbc)
    # sitk.WriteImage(sitk_rbc,"rbc.nii")
    # sitk_barrier = sitk.GetImageFromArray(barrier)
    # sitk.WriteImage(sitk_barrier,"barrier.nii")
    # pdb.set_trace()
    return rbc, barrier

def gasBinning(gas_highreso,bin_threshold,mask,percentile):
    ## binning for gas
    from GX_utils import binning

    gas_highreso = abs(gas_highreso)

    gas_thre = np.percentile(gas_highreso[mask], percentile)

    gas_highreso_m = np.divide(np.multiply(gas_highreso,mask),gas_thre)

    gas_highreso_m[gas_highreso_m > 1] = 1

    gas_binning = binning(gas_highreso_m, bin_threshold)

    # create ventilation mask
    mask_vent = gas_binning
    mask_vent[mask_vent<3] = 0
    mask_vent[mask_vent>0] = 1

    return gas_binning, mask_vent

def disBinning(discomp,gas_highSNR,bin_threshold,mask,cor=1):

    ## binning for rbc or barrier
    ground = 1e-5

    from GX_utils import binning

    comp2gas = np.zeros(np.shape(gas_highreso))

    comp2gas[mask] = np.divide(discomp[mask],abs(gas_highreso[mask]))
    comp2gas = np.multiply(comp2gas,cor)

    comp2gas[comp2gas < ground] = ground

    comp2gas_binning= binning(comp2gas, bin_threshold)

    return comp2gas_binning
