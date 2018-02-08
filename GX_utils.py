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

    gas_thre = np.percentile(gas_highreso[mask], percentile)

    gas_highreso_m = np.divide(np.multiply(gas_highreso,mask),gas_thre)

    gas_highreso_m[gas_highreso_m > 1] = 1

    gas_binning = binning(gas_highreso_m, bin_threshold)

    # create ventilation mask
    mask_vent = gas_binning
    mask_vent[mask_vent<3] = 0
    mask_vent[mask_vent>0] = 1
    mask_vent = mask_vent.astype(bool)

    return gas_binning, mask_vent

def disBinning(discomp,gas_highSNR,bin_threshold,mask,cor=1):

    ## binning for rbc or barrier
    ground = 1e-5

    from GX_utils import binning

    comp2gas = np.zeros(np.shape(gas_highSNR))

    comp2gas[mask] = np.divide(discomp[mask],gas_highSNR[mask])
    comp2gas = np.multiply(comp2gas,cor)

    comp2gas[comp2gas < ground] = ground

    comp2gas_binning= binning(comp2gas, bin_threshold)

    return comp2gas_binning

def register(gas_highreso,ute, mask):

    # register mask to gas_highreso, and apply the transform to UTE

    sitk_mask = sitk.GetImageFromArray(mask)
    sitk_gas = sitk.GetImageFromArray(gas_highreso)
    sitk_ute = sitk.GetImageFromArray(ute)

    ## registrate mask to gas_highreso
    # set up elastix filter for the warpping
    elastixImageFilter = sitk.ElastixImageFilter()
    elastixImageFilter.SetFixedImage(sitk_gas) # gas highreso
    elastixImageFilter.SetMovingImage(sitk_mask) # mask

    # set up parameters for the warpping, we use affine first and then use bspline interpolation for non-rigid warpping
    parameterMapVector = sitk.VectorOfParameterMap()
    parameterMapVector.append(sitk.GetDefaultParameterMap("affine"))
    elastixImageFilter.SetParameterMap(parameterMapVector)

    elastixImageFilter.Execute()

    # save warpped image
    sitk_mask_reg = elastixImageFilter.GetResultImage()

    np_mask_reg = sitk.GetArrayFromImage(sitk_mask_reg)

    np_mask_reg[np_mask_reg < 0.5] = 0
    np_mask_reg[np_mask_reg > 0] = 1
    np_mask_reg = np_mask_reg.astype(bool)

    # sitk_mask_reg = sitk.GetImageFromArray(np_mask_reg)
    # sitk_mask_reg = sitk.PermuteAxes(sitk_mask_reg,[2,1,0])
    # fix the direction change due to numpy array trans
    # sitk.WriteImage(sitk_mask_reg,p_mask_grow+"_reg.nii")

    # apply the same transform to the label image
    transformParameterMap = elastixImageFilter.GetTransformParameterMap()

    transformixImageFilter = sitk.TransformixImageFilter()
    transformixImageFilter.SetTransformParameterMap(transformParameterMap)
    transformixImageFilter.SetMovingImage(sitk_ute) # fixed ute

    transformixImageFilter.Execute()

    sitk_ute_reg = transformixImageFilter.GetResultImage()
    sitk_ute_reg = sitk.PermuteAxes(sitk_ute_reg,[2,1,0])

    np_ute_reg = sitk.GetArrayFromImage(sitk_ute_reg)
    # sitk.WriteImage(sitk_ute_reg,p_ute_recon+"_reg.nii")
    return np_ute_reg, np_mask_reg
