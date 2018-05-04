import numpy as np
import nibabel as nib
import scipy.io as sio
import SimpleITK as sitk
import nibabel as nib
import math, pdb, os, glob, time, warnings
# from GX_Map_utils import fullMontage
from scipy import io
import pandas as pd

from Lobar_utils import gen_lobar_mask

move_filename = '/home/ziyiw/Gas_Exchange/Lobar_Seg/ProtonLobeMasks/0004Lobes.nii'
data_dir = '/home/ziyiw/Patients/'

Subject_IDs = os.listdir(data_dir)

# Subject_IDs_no = ['005018', '004003A', '005010A', '004003A_Rep', '005006']
# create an empty pandas table
columns = ['left_low','left_up', 'right_low', 'right_mid','right_up','nleft_low','nleft_up', 'nright_low', 'nright_mid','nright_up']
ven_sum_table = pd.DataFrame(columns = columns)
bar_sum_table = pd.DataFrame(columns = columns)
rbc_sum_table = pd.DataFrame(columns = columns)

for Subject_ID in Subject_IDs:

    fdata = data_dir+Subject_ID+'/'+Subject_ID+'.mat'
    mat_input = sio.loadmat(fdata)

    rbc2gas = mat_input['rbc2gas']
    bar2gas = mat_input['bar2gas']
    ventilation = mat_input['ventilation']

    # if(Subject_ID in Subject_IDs_no):
    #     mask = mat_input['mask_reg']
    #     mask_lobe = gen_lobar_mask(move_filename = move_filename, img_fix = mask).astype(int)
    #
    #     # save the mask to the .mat file for each subject
    #     mat_input['mask_lobe'] = mask_lobe
    #     sio.savemat(fdata, mat_input)
    # else:
    try:
        mask_lobe = mat_input['mask_lobe']
    except:
        mask = mat_input['mask_reg']
        mask_lobe = gen_lobar_mask(move_filename = move_filename, img_fix = mask).astype(int)

        # save the mask to the .mat file for each subject
        mat_input['mask_lobe'] = mask_lobe
        sio.savemat(fdata, mat_input)

    num_pix = np.size(mask_lobe[mask_lobe>0])

    for lobe_index in range(1,6):

        sum_ven = np.sum(ventilation[mask_lobe==lobe_index])
        sum_rbc = np.sum(rbc2gas[mask_lobe==lobe_index])
        sum_bar = np.sum(bar2gas[mask_lobe==lobe_index])

        per_pix = np.size(mask_lobe[mask_lobe==lobe_index])/(num_pix*1.0)

        ven_sum_table.loc[Subject_ID,columns[lobe_index-1]] = sum_ven
        bar_sum_table.loc[Subject_ID,columns[lobe_index-1]] = sum_bar
        rbc_sum_table.loc[Subject_ID,columns[lobe_index-1]] = sum_rbc

        ven_sum_table.loc[Subject_ID,columns[lobe_index-1+5]] = per_pix
        bar_sum_table.loc[Subject_ID,columns[lobe_index-1+5]] = per_pix
        rbc_sum_table.loc[Subject_ID,columns[lobe_index-1+5]] = per_pix

    print('Completed Subject '+Subject_ID)

ven_sum_table.to_csv('ven_sum_new')
bar_sum_table.to_csv('bar_sum_new')
rbc_sum_table.to_csv('rbc_sum_new')
pdb.set_trace()
