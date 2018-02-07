#  Ziyi, Feb.07.2018
#
# This file complete the binning for Gas high resolution
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
fmask = 'BHUTE_Sub002102_FID49886_mask_grow_reg.nii'
fdata = 'Sub002102_data.mat'

mat_input = sio.loadmat(fdata)
gas_highreso = mat_input['gasVol_highreso']

mask = np.array(nib.load(fmask).get_data(),dtype='bool')

from GX_defineColormaps import thre_vent as bin_threshold

num_bin = len(bin_threshold)+1

percentile = 0.99

stats = {}
## *********************************
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

pdb.set_trace()
