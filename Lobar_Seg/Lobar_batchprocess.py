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
data_dir = '/home/ziyiw/Patients_old/'

Subject_IDs = os.listdir(data_dir)

Subject_IDs_no = ['005018', '004003A', '005010A', '004003A_Rep', '005006']
# create an empty pandas table
columns = ['llow','lup', 'rlow', 'rmid','rup']
ven_table = pd.DataFrame(columns = columns)
bar_table = pd.DataFrame(columns = columns)
rbc_table = pd.DataFrame(columns = columns)

for Subject_ID in Subject_IDs:

    if(Subject_ID in Subject_IDs_no):
        continue

    fdata = data_dir+Subject_ID+'/'+Subject_ID+'.mat'
    mat_input = sio.loadmat(fdata)

    rbc2gas = mat_input['rbc2gas']
    bar2gas = mat_input['bar2gas']
    ventilation = mat_input['ventilation']
    mask = mat_input['mask_reg']

    mask_lobe = gen_lobar_mask(move_filename = move_filename, img_fix = mask).astype(int)

    # save the mask to the .mat file for each subject
    mat_input['mask_lobe'] = mask_lobe
    sio.savemat(fdata, mat_input)

    for lobe_index in range(1,6):

        mean_ven = np.average(ventilation[mask_lobe==lobe_index])
        mean_rbc = np.average(rbc2gas[mask_lobe==lobe_index])
        mean_bar = np.average(bar2gas[mask_lobe==lobe_index])

        ven_table.loc[Subject_ID,columns[lobe_index-1]] = mean_ven
        bar_table.loc[Subject_ID,columns[lobe_index-1]] = mean_bar
        rbc_table.loc[Subject_ID,columns[lobe_index-1]] = mean_rbc

    print('Completed Subject '+Subject_ID)

ven_table.to_csv('ven_table_old')
bar_table.to_csv('bar_table_old')
rbc_table.to_csv('rbc_table_old')
pdb.set_trace()
