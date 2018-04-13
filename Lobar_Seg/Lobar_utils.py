import numpy as np
import nibabel as nib
import scipy.io as sio
import SimpleITK as sitk
import nibabel as nib
import warnings
import math, pdb, os
import time
from GX_Map_utils import fullMontage
from scipy import io

#### input
# move_filename = '/home/ziyi/Lobe_Seg/ProtonLobeMasks/0004Lobes.nii'
# fix_filename = '/home/ziyi/Lobe_Seg/Sub005013.nii'
#
# img_fix = nib.load(fix_filename).get_data()

def gen_lobar_mask(move_filename, img_fix):
    ## need to align these 2 images
    img_move = nib.load(move_filename).get_data()

    pad0 = np.shape(img_fix)[0] - np.shape(img_move)[0]
    pad1 = np.shape(img_fix)[1] - np.shape(img_move)[1]
    pad2 = np.shape(img_fix)[2] - np.shape(img_move)[2]
    # padding the moving image
    img_move = np.pad(img_move,((pad0/2,pad0/2),(pad1/2,pad1/2),(pad2/2,pad2/2)), 'constant',constant_values=0)
    # img_move = np.rot90(np.flip(img_move,axis=2),k=-1) # use for old nii orientation
    img_move = np.rot90(img_move,k=1) # use for new orientation
    img_move_b = np.copy(img_move)
    img_move_b[img_move_b>0] = 1

    # img_move_b[img_move_b==0] = np.nan
    # img_fix[img_fix==0] = np.nan
    ##################
    # register img_move_b to img_fix, and apply the transform to img_follow

    sitk_move = sitk.GetImageFromArray(img_move_b)
    sitk_fix = sitk.GetImageFromArray(img_fix)
    sitk_follow = sitk.GetImageFromArray(img_move)

    # set up elastix filter for the warpping
    elastixImageFilter = sitk.ElastixImageFilter()
    elastixImageFilter.SetFixedImage(sitk_fix) # gas highreso
    elastixImageFilter.SetMovingImage(sitk_move) # img_move

    # set up parameters for the warpping, we use affine first
    parameterMapVector = sitk.VectorOfParameterMap()
    parameterMapVector.append(sitk.GetDefaultParameterMap("affine")) #"rigid"
    parameterMapVector.append(sitk.GetDefaultParameterMap("bspline")) #"nonrigid"
    elastixImageFilter.SetParameterMap(parameterMapVector)

    elastixImageFilter.Execute()

    # save warpped image
    sitk_move_reg = elastixImageFilter.GetResultImage()
    #
    # np_img_fix = sitk.GetArrayFromImage(sitk_fix)
    np_img_lobe_b = sitk.GetArrayFromImage(sitk_move_reg)
    # # pdb.set_trace()
    # np_img_move_b[np_img_move_b < 0.5] = 0
    # np_img_move_b[np_img_move_b > 0] = 1
    # np_img_move_b = np_img_move_b.astype(bool)

    # # sitk_move_reg = sitk.GetImageFromArray(np_img_move)
    # # sitk_move_reg = sitk.Permimg_followAxes(sitk_move_reg,[2,1,0])
    # # fix the direction change due to numpy array trans
    # # sitk.WriteImage(sitk_move_reg,p_img_move_grow+"_reg.nii")
    #
    # apply the same transform to the label image
    transformParameterMap = elastixImageFilter.GetTransformParameterMap()

    transformixImageFilter = sitk.TransformixImageFilter()
    transformixImageFilter.SetTransformParameterMap(transformParameterMap)
    transformixImageFilter.SetMovingImage(sitk_follow) # fixed img_follow

    transformixImageFilter.Execute()

    sitk_follow_reg = transformixImageFilter.GetResultImage()
    # sitk_follow_reg = sitk.Permimg_followAxes(sitk_follow_reg,[2,1,0])

    np_img_lobe = sitk.GetArrayFromImage(sitk_follow_reg)

    # converge into 5 discreet numbers
    np_img_lobe = np.multiply(np_img_lobe,img_fix)
    np_img_lobe[np_img_lobe<0.0] = 0.0
    np_img_lobe[(np_img_lobe>0.0) & (np_img_lobe<=1.5)] = 1.0
    np_img_lobe[(np_img_lobe>1.5) & (np_img_lobe<=2.5)] = 2.0
    np_img_lobe[(np_img_lobe>2.5) & (np_img_lobe<=3.5)] = 3.0
    np_img_lobe[(np_img_lobe>3.5) & (np_img_lobe<=4.5)] = 4.0
    np_img_lobe[np_img_lobe>4.5] = 5.0

    return np_img_lobe
