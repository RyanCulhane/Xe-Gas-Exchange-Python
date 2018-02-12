import numpy as np
import nibabel as nib
import scipy.io as sio
import SimpleITK as sitk
from skimage import color
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

    return gas_highreso_m, gas_binning, mask_vent

def disBinning(discomp,gas_highSNR,bin_threshold,mask,cor=1):

    ## binning for rbc or barrier
    ground = 1e-5

    from GX_utils import binning

    comp2gas = np.zeros(np.shape(gas_highSNR))

    comp2gas[mask] = np.divide(discomp[mask],gas_highSNR[mask])
    comp2gas = np.multiply(comp2gas,cor)

    comp2gas[comp2gas < ground] = ground

    comp2gas_binning= binning(comp2gas, bin_threshold)

    return comp2gas, comp2gas_binning

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

def mergeRGBandGray(ute_slice,binning_slice):
    # function combine the gray scale UTE with the RGB binning via HSV

    # construct RGB version of gray-level ute
    ute_slice_color = np.dstack((ute_slice,ute_slice,ute_slice))

    # Convert the input image and color mask to HSV
    ute_slice_hsv = color.rgb2hsv(ute_slice_color)
    binning_slice_hsv = color.rgb2hsv(binning_slice)

    # Replace the hue and saturation of the original image
    # with that of the color mask
    ute_slice_hsv[..., 0] = binning_slice_hsv[..., 0]
    ute_slice_hsv[..., 1] = binning_slice_hsv[..., 1]

    mask = ((binning_slice[:,:,0]==0) & (binning_slice[:,:,1]==0) & (binning_slice[:,:,2]==0))
    mask = ~mask

    ute_slice_hsv[mask,:] = binning_slice_hsv[mask,:]

    colormap = color.hsv2rgb(ute_slice_hsv)

    return colormap

def montage(Img):
    ## plot montage(2*8) of Img
    ## Img has to have 16 slices
    img_w = 128
    img_h = 128
    count = 16
    n_row = 2
    n_col = 8

    img_montage = np.zeros((n_row * img_h, n_col * img_w,3))

    image_id = 0
    for j in range(n_row):
        for k in range(n_col):
            if image_id >= count:
                break
            sliceN, sliceM = j * img_h, k * img_w
            img_montage[sliceN:sliceN + img_h, sliceM:sliceM + img_w,:] = Img[:, :, image_id,:]
            image_id += 1

    return img_montage

def makeMontage(bin_index, ute_reg, index2color, ind_start, ind_inter):
    ## make montage (2*8) from binning map and ute image
    ## the montage will pick the image from ind_start
    img_w,img_h,img_d = np.shape(bin_index)

    num_slice = 16
    ind_end = ind_start + ind_inter*num_slice

    colormap = np.zeros((img_w,img_h,num_slice,3))
    ind_slice = 0

    def convert_index2color(ind):
        return index2color[ind]

    ## convert each slice from index to RGB, then combine ute_reg and bin_index_RGB to HSV
    for k in range(ind_start, ind_end, ind_inter):

        ## convert bin_index to bin_rgb
        bin_rgb = map(convert_index2color, bin_index[:,:,k].flatten())
        bin_rgb = np.asarray(bin_rgb)
        bin_rgb = np.reshape(bin_rgb,(img_w,img_h,3))

        ## merge bin_rgb with ute_reg through hsv colorspace
        colorslice = mergeRGBandGray(ute_reg[:,:,k],bin_rgb)

        colormap[:,:,ind_slice,:] = colorslice
        ind_slice = ind_slice + 1

    ## make montage from the image stack
    img_montage = montage(colormap)
    return img_montage

def makeHistogram(data, color, x_lim, y_lim, num_bins, refer_fit, hist_name):
    ## plot histogram for the gas exchange ratio maps
    # make a thick frame
    rc('axes', linewidth=2)

    fig, ax = plt.subplots()

    # the histogram of the data
    data = data.flatten()
    data[data<0] = 0
    data = np.append(data, x_lim)

    weights = np.ones_like(data)/float(len(data))

    num, bins, patches = ax.hist(data, num_bins,color=color,weights=weights)

    # define and plot healthy reference line
    normal = refer_fit[0]*np.exp(-((bins-refer_fit[1])/refer_fit[2])**2);

    ax.plot(bins, normal, '--',color='k',linewidth=4)

    # ax.set_xlabel('Pixel Intensity', fontsize=26)
    ax.set_ylabel('Fraction of total pixels',fontsize=26)
    xlim((0,x_lim))
    ylim((0,y_lim))

    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(18)
        # tick.label1.set_fontweight('bold')
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(18)
        # tick.label1.set_fontweight('bold')

    # Tweak spacing to prevent clipping of ylabel
    fig.tight_layout()
    plt.savefig(hist_name)
