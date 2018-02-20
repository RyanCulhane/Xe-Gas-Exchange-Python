#  Ziyi, Jan 26. 2018
#
# This file register mask to gas_highreso, and apply the transform to ute_recon
#

import numpy as np
import nibabel as nib
import SimpleITK as sitk

import pdb

p_mask_grow = "BHUTE_mask_grow.nii"
p_ute_recon = "BHUTE_recon.nii"
p_gas_highreso = "Dixon_gas_recon_highreso.nii"

img_mask_grow = nib.load(p_mask_grow)
img_ute_recon = nib.load(p_ute_recon)
img_gas_highreso = nib.load(p_gas_highreso)

#convert from np to sitk image
np_mask = np.array(img_mask_grow.get_data())
np_gas = np.array(img_gas_highreso.get_data())
np_gas = np.divide(np_gas,np.max(np_gas)) # can use top percentile
np_ute = np.array(img_ute_recon.get_data())

sitk_mask = sitk.GetImageFromArray(np_mask)
sitk_gas = sitk.GetImageFromArray(np_gas)
sitk_ute = sitk.GetImageFromArray(np_ute)
## registrate mask to gas_highreso
# pdb.set_trace()
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

sitk_mask_reg = sitk.GetImageFromArray(np_mask_reg)
sitk_mask_reg = sitk.PermuteAxes(sitk_mask_reg,[2,1,0])
# fix the direction change due to numpy array trans
sitk.WriteImage(sitk_mask_reg,p_mask_grow+"_reg.nii")

# apply the same transform to the label image
transformParameterMap = elastixImageFilter.GetTransformParameterMap()

transformixImageFilter = sitk.TransformixImageFilter()
transformixImageFilter.SetTransformParameterMap(transformParameterMap)
transformixImageFilter.SetMovingImage(sitk_ute) # fixed ute

transformixImageFilter.Execute()

sitk_ute_reg = transformixImageFilter.GetResultImage()
sitk_ute_reg = sitk.PermuteAxes(sitk_ute_reg,[2,1,0])
sitk.WriteImage(sitk_ute_reg,p_ute_recon+"_reg.nii")
