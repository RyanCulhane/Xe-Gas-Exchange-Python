import os
import numpy as np
import nibabel as nib
import scipy.io as sio
import SimpleITK as sitk
from skimage import color
import warnings
import math
from scipy.ndimage.morphology import binary_dilation
import pdfkit

from matplotlib import pyplot as plt

import pdb

def getIndexMaxOnes(arr):
    ## returns the starting index and ending index of the max consecutive ones
    # intitialize count
    cur_count = 0
    cur_sta = 0

    max_count = 0
    pre_state = 0

    index_sta = 0
    index_end = 0

    for i in range(0, np.size(arr)):

        if (arr[i] == 0):
            cur_count = 0
            if((pre_state == 1)&(cur_sta == index_sta)):
                index_end = i-1
            pre_state = 0

        else:
            if(pre_state == 0):
                cur_sta = i
                pre_state = 1
            cur_count+= 1
            if(cur_count>max_count):
                max_count = cur_count
                index_sta = cur_sta

    return index_sta,index_end

def decideStartInterval(mask):

    # determine the starting slice and the interval for the montage
    num_slice = 16.0

    sum_line = np.sum(np.sum(mask,axis=0),axis=0)

    binary_arr = sum_line>300

    ind_start, ind_end =  getIndexMaxOnes(binary_arr)

    ind_inter = np.ceil((ind_end-ind_start)/num_slice).astype(int)

    return ind_start, ind_inter

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
        iterCount = iterCount+1
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
    return rbc, barrier

def binning(volume,thresholds):

    # volume: mask_vented, thresholded 3D volume
    # thresholds: just the middle thresholds

    bvolume = np.ones(np.shape(volume))

    bvolume[(volume > 0) & (volume <= thresholds[0])] = 2

    for k in range(len(thresholds)-1):
        bvolume[(volume > thresholds[k]) & (volume <= thresholds[k+1])] = k+3

    bvolume[volume > thresholds[-1]] = len(thresholds)+2

    return bvolume

def gasBinning(gas_highreso,bin_threshold,mask,percentile):
    ## binning for gas
    gas_thre = np.percentile(gas_highreso[mask], percentile)

    gas_highreso_n = np.divide(np.multiply(gas_highreso,mask),gas_thre)

    gas_highreso_n[gas_highreso_n > 1] = 1

    gas_binning = binning(gas_highreso_n, bin_threshold)

    # create ventilation mask
    mask_vent = np.copy(gas_binning)
    mask_vent[mask_vent<3] = 0 # exclude VDP in the ventilation map
    # mask_vent[mask_vent>0] = 1
    mask_vent = mask_vent.astype(bool)

    return gas_highreso_n, gas_binning, mask_vent

def disBinning(discomp,gas_highSNR,bin_threshold, mask,cor=1):

    ## binning for rbc or barrier
    ground = 1e-5
    ground_gas = 1e-15
    # prevent divided by 0
    # gas_highSNR[(mask>0) & (gas_highSNR < ground_gas)] = ground_gas

    comp2gas = np.zeros(np.shape(gas_highSNR))

    comp2gas[mask] = np.divide(discomp[mask],gas_highSNR[mask])
    comp2gas = np.multiply(comp2gas,cor)

    #set negative values
    comp2gas[ comp2gas < 0] = ground

    comp2gas_binning= binning(comp2gas, bin_threshold)

    return comp2gas, comp2gas_binning

def register(gas_highreso, ute, mask):

    # register mask to gas_highreso, and apply the transform to UTE

    sitk_mask = sitk.GetImageFromArray(mask)
    sitk_gas = sitk.GetImageFromArray(gas_highreso)
    sitk_ute = sitk.GetImageFromArray(ute)

    ## registrate mask to gas_highreso
    # set up elastix filter for the warpping
    elastixImageFilter = sitk.ElastixImageFilter()
    elastixImageFilter.SetFixedImage(sitk_gas) # gas highreso
    elastixImageFilter.SetMovingImage(sitk_mask) # mask

    # set up parameters for the warpping, we use affine first
    parameterMapVector = sitk.VectorOfParameterMap()
    parameterMapVector.append(sitk.GetDefaultParameterMap("affine")) #"rigid"
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
    # sitk_ute_reg = sitk.PermuteAxes(sitk_ute_reg,[2,1,0])

    np_ute_reg = sitk.GetArrayFromImage(sitk_ute_reg)
    # sitk.WriteImage(sitk_ute_reg,p_ute_recon+"_reg.nii")
    return np_ute_reg.astype('float64'), np_mask_reg

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

def biasFieldCor(image, mask):

    pathInput = 'image.nii'
    pathMask = 'mask.nii'
    pathOutput = 'image_cor.nii'
    pathBiasField = 'biasfield.nii'

    pathN4 = './N4BiasFieldCorrection'

    # save the inputs into nii files so the execute N4 can read in
    nii_imge = nib.Nifti1Image(abs(image), np.eye(4))
    nii_mask = nib.Nifti1Image(mask.astype(float), np.eye(4))
    nib.save(nii_imge, pathInput)
    nib.save(nii_mask, pathMask)

    # cmd = pathN4+' -d 3 -i '+pathInput+' -s 2 -x '+pathMask+' -o ['+pathOutput+', '+pathBiasField+']'
    cmd = pathN4+' -d 3 -i '+pathInput+' -x '+pathMask+' -o ['+pathOutput+', '+pathBiasField+']'

    os.system(cmd)

    image_cor = np.array(nib.load(pathOutput).get_data())
    image_biasfield = np.array(nib.load(pathBiasField).get_data())

    # remove the generated nii files
    os.remove(pathInput)
    os.remove(pathMask)
    os.remove(pathOutput)
    os.remove(pathBiasField)

    return image_cor.astype('float64'), image_biasfield.astype('float64')

def fullMontage(X, colormap='gray'):
    # used for debug, plot the entire montage
    m, n, count = np.shape(X)
    mm = int(np.ceil(np.sqrt(count)))
    nn = mm
    M = np.zeros((mm * m, nn * n))

    image_id = 0
    for j in range(mm):
        for k in range(nn):
            if image_id >= count:
                break
            sliceN, sliceM = j * m, k * n
            M[sliceN:sliceN + n, sliceM:sliceM + m] = X[:, :, image_id]
            image_id += 1

    plt.figure()
    plt.imshow(M, cmap=colormap)
    plt.axis('off')
    plt.show(block=False)

    return M

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

def makeMontage(bin_index, ute_reg, index2color, ind_start, ind_inter, mon_name):
    ## make montage (2*8) from binning map and ute image
    ## the montage will pick the image from ind_start
    # normalize ute
    ute_thre = np.percentile(ute_reg, 99)
    ute_reg_m = np.divide(ute_reg,ute_thre)
    ute_reg_m[ute_reg_m>1] = 1

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
        colorslice = mergeRGBandGray(ute_reg_m[:,:,k],bin_rgb)

        colormap[:,:,ind_slice,:] = colorslice
        ind_slice = ind_slice + 1

    ## make montage from the image stack
    img_montage = montage(colormap)

    ## plot and save the montage
    plt.imshow(img_montage,interpolation='none')
    plt.axis('off')
    plt.savefig(mon_name,transparent = True,bbox_inches='tight',pad_inches=-0.1)

    return img_montage

def makeHistogram(data, color, x_lim, y_lim, num_bins, refer_fit, hist_name):
    ## plot histogram for the gas exchange ratio maps
    # make a thick frame
    from matplotlib.pyplot import rc, xlim, ylim
    # pdb.set_trace()
    rc('axes', linewidth=4)

    fig, ax = plt.subplots(figsize=(9, 6))

    # the histogram of the data
    # limit the range of data
    data = data.flatten()
    data[data<0] = 0
    data[data>x_lim] = x_lim
    data = np.append(data, x_lim)

    weights = np.ones_like(data)/float(len(data))

    # plot histogram
    num, bins, patches = ax.hist(data, num_bins,color=color,weights=weights,edgecolor='black')

    # define and plot healthy reference line
    normal = refer_fit[0]*np.exp(-((bins-refer_fit[1])/refer_fit[2])**2);
    ax.plot(bins, normal, '--',color='k',linewidth=4)

    # ax.set_xlabel('Pixel Intensity', fontsize=26)
    ax.set_ylabel('Fraction of Total Pixels',fontsize=35)
    xlim((0,x_lim))
    ylim((0,y_lim))

    plt.locator_params(axis='x', nbins=4, nticks=4)

    # barrier/RBC add a tick on the end
    from GX_defineColormaps import histogram_ticks
    if('bar_hist.png' in hist_name):
        xticks = histogram_ticks['bar_xticks']
        yticks = histogram_ticks['bar_yticks']

    if('rbc_hist.png' in hist_name):
        xticks = histogram_ticks['rbc_xticks']
        yticks = histogram_ticks['rbc_yticks']

    if('ven_hist.png' in hist_name):
        xticks = histogram_ticks['ven_xticks']
        yticks = histogram_ticks['ven_yticks']

    # plt.locator_params(axis='y', nbins=4, nticks=4)
    xticklabels = ['{:.1f}'.format(x) for x in xticks]
    yticklabels = ['{:.2f}'.format(x) for x in yticks]

    plt.xticks(xticks, xticklabels,fontsize=40)
    plt.yticks(yticks, yticklabels,fontsize=40)


    # Tweak spacing to prevent clipping of ylabel
    fig.tight_layout()
    # pdb.set_trace()
    plt.savefig(hist_name)

def fSNR_3T(image,mask):

    # x-y border exclusion
    xybe = 0
    mini_cube_dim = [8,8,8]

    my_dim = np.shape(image)

    # dilate the mask to analyze noise area away from the signal
    def util(x):
        return int((math.ceil(x*0.025)*2+1))

    dilate_struct = np.ones((util(my_dim[0]), util(my_dim[1]), util(my_dim[2])))
    noise_mask = binary_dilation(mask,dilate_struct).astype(bool)

    # exclue border too
    if(xybe>0):
        noise_mask[0:xybe,:,:] = True
        noise_mask[-xybe:,:,:] = True
        noise_mask[:,0:xybe,:] = True
        noise_mask[:,-xybe:,:] = True

    noise_temp = np.copy(image)
    noise_temp[noise_mask] = np.nan

    # set up for using mini noise cubes to go through the image and calculate std for noise
    (mini_x,mini_y,mini_z) = mini_cube_dim

    n_noise_vox = mini_x*mini_y*mini_z

    mini_vox_std = 0.75*n_noise_vox # minimul number of voxels to calculate std

    stepper = 0
    total = 0
    std_dev_mini_noise_vol = []

    for ii in range(0, my_dim[0]/mini_x):
        for jj in range(0, my_dim[1]/mini_y):
            for kk in range(0, my_dim[2]/mini_z):

                mini_cube_noise_dist = noise_temp[ii*mini_x:(ii+1)*mini_x, jj*mini_y:(jj+1)*mini_y, kk*mini_z:(kk+1)*mini_z]

                mini_cube_noise_dist = mini_cube_noise_dist[~np.isnan(mini_cube_noise_dist)]

                # only calculate std for the noise when it is long enough
                if(len(mini_cube_noise_dist)>mini_vox_std):
                    std_dev_mini_noise_vol.append(np.std(mini_cube_noise_dist,ddof=1))
                    stepper = stepper+1

                total = total+1

    image_noise = np.median(std_dev_mini_noise_vol)
    image_signal = np.average(image[mask])

    SNR = image_signal/image_noise
    SNR_Rayleigh = SNR*0.66

    return SNR, SNR_Rayleigh, image_signal, image_noise

def binStats(rawdata, bindata, mask, mask_all, key, recondata=None):

    statsbox = {}
    maskall = np.sum(mask_all).astype('float')

    statsbox[key+'_defect'] = np.divide(np.sum((bindata == 2)),maskall)
    statsbox[key+'_low'] = np.divide(np.sum((bindata == 3)),maskall)
    statsbox[key+'_mean'] = np.average(abs(rawdata[mask]))

    if ((key == 'rbc')|(key == 'ven')):
        statsbox[key+'_high'] = np.divide(np.sum((bindata == 6)|(bindata == 7)),maskall)
    else:
        statsbox[key+'_high'] = np.divide(np.sum((bindata == 8)|(bindata == 9)),maskall)

    if ((key == 'rbc')|(key == 'bar')):
        statsbox[key+'_negative'] = np.divide(np.sum((recondata < 0)&(mask > 0)),maskall)
        statsbox[key+'_SNR'],_,_,_ = fSNR_3T(recondata, mask_all)
    else:
        _, SNR_Rayleigh, _, _ = fSNR_3T(recondata, mask_all)
        statsbox[key+'_SNR'] = SNR_Rayleigh

    return statsbox

def pdfScaleMerge(input_file1, input_file2, output_file, scale=1):

    from pyPdf import PdfFileWriter, PdfFileReader

    with open(input_file1, "rb") as in_f1:
        with open(input_file2,'rb') as in_f2:
            input1 = PdfFileReader(in_f1)
            input2 = PdfFileReader(in_f2)
            output = PdfFileWriter()

            p = input1.getPage(0)
            p.scale(4,4) # scale it up by a factor of 2
            output.addPage(p)

            p = input2.getPage(0)
            p.scale(4,4) # scale it up by a factor of 2
            output.addPage(p)

            with open(output_file, "wb") as out_f:
                output.write(out_f)

def genHtmlPdf(Subject_ID, data_dir, RBC2barrier, stats_box):
    # reder html using the templates and stats
    temp_clinical = "html_tmp/temp_clinical.html"
    temp_technical = "html_tmp/temp_technical.html"

    report_clinical = "report_clinical.html"
    report_technical = "report_technical.html"

    def adj_format1(x):
        num = np.around((x)*100,decimals=0).astype(int)
        if((num==0)):
            return '<1'
        else:
            return num.astype(str)

    def adj_format2(x):
        return np.around(x,decimals=2).astype(str)


    html_parameters = {
        'Subject_ID': Subject_ID[:3]+'-'+Subject_ID[3:],
        'inflation': np.around(stats_box['inflation'],decimals=2),
        'RBC2barrier': np.around(RBC2barrier,decimals=3).astype(str),
        'ven_defect': adj_format1(stats_box['ven_defect']),
        'ven_low': adj_format1(stats_box['ven_low']),
        'ven_high': adj_format1(stats_box['ven_high']),
        'ven_mean': adj_format2(stats_box['ven_mean']),
        'ven_SNR': adj_format2(stats_box['ven_SNR']),
        'bar_defect': adj_format1(stats_box['bar_defect']),
        'bar_low': adj_format1(stats_box['bar_low']),
        'bar_high': adj_format1(stats_box['bar_high']),
        'bar_mean': adj_format2(stats_box['bar_mean']),
        'bar_SNR': adj_format2(stats_box['bar_SNR']),
        'bar_negative': adj_format1(stats_box['bar_negative']),
        'rbc_defect': adj_format1(stats_box['rbc_defect']),
        'rbc_low': adj_format1(stats_box['rbc_low']),
        'rbc_high': adj_format1(stats_box['rbc_high']),
        'rbc_mean': adj_format2(stats_box['rbc_mean']),
        'rbc_SNR': adj_format2(stats_box['rbc_SNR']),
        'rbc_negative': adj_format1(stats_box['rbc_negative']),

        'ven_montage':data_dir+'/ven_montage.png',
        'bar_montage':data_dir+'/bar_montage.png',
        'rbc_montage':data_dir+'/rbc_montage.png',

        'ven_hist':data_dir+'/ven_hist.png',
        'bar_hist':data_dir+'/bar_hist.png',
        'rbc_hist':data_dir+'/rbc_hist.png'
    }

    from GX_defineColormaps import referece_stats
    html_parameters.update(referece_stats)

    # generate pdf from html
    # reder html second
    with open(temp_clinical, 'r') as f:
        data = f.read()
        rendered = data.format(**html_parameters)

    with open(report_clinical, 'w') as o:
        o.write(rendered)

    with open(temp_technical, 'r') as f:
        data = f.read()
        rendered = data.format(**html_parameters)

    with open(report_technical, 'w') as o:
        o.write(rendered)

    options = {
        'page-width':300,#320,
        'page-height':150,#160,
        'margin-top': 1,
        'margin-right': 0.1,
        'margin-bottom': 0.1,
        'margin-left': 0.1,
        'dpi':300,
        # 'zoom':2,
        'disable-smart-shrinking':'',
        'encoding': "UTF-8",
        }

    # generate and scale PDF
    pdf_clinical = data_dir+'/report_clinical_'+Subject_ID+'.pdf'
    pdf_technical = data_dir+'/report_technical_'+Subject_ID+'.pdf'
    pdf_merge = data_dir+'/report_'+Subject_ID+'.pdf'

    # pdb.set_trace()

    pdfkit.from_file(report_clinical, pdf_clinical, options=options)
    pdfkit.from_file(report_technical,pdf_technical,options=options)

    # pdb.set_trace()
    os.remove(report_technical)
    os.remove(report_clinical)

    # scale and merge pdf to convert to ppt
    pdfScaleMerge(input_file1=pdf_clinical, input_file2=pdf_technical, output_file=pdf_merge, scale=1)

    # generate ppt
    genPPT(pdf_file=pdf_merge)

    # remove files generated when making ppt
    os.remove(data_dir+'/report_'+Subject_ID+'_1.pdf')
    os.remove(data_dir+'/report_'+Subject_ID+'_1.jpg')
    os.remove(data_dir+'/report_'+Subject_ID+'_2.pdf')
    os.remove(data_dir+'/report_'+Subject_ID+'_2.jpg')

def genPPT(pdf_file):
    # generate ppt from pdf
    from ppt_pdf.cli_pdf_to_ppt import PdfToPpt

    PdfToPpt(pdf_file=pdf_file).execute()
