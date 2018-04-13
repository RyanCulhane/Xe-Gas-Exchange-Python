import numpy as np
import nibabel as nib
import scipy.io as sio
import SimpleITK as sitk
import nibabel as nib
import math, pdb, os, glob, time, warnings
from GX_Map_utils import fullMontage
from scipy import io
import pandas as pd

from Lobar_utils import gen_lobar_mask

move_filename = '/home/ziyi/Lobe_Seg/ProtonLobeMasks/0004Lobes.nii'
data_dir = '/home/ziyi/Patients/'

Subject_IDs = os.listdir(data_dir)

for Subject_ID in Subject_IDs:

    fdata = data_dir+Subject_ID+'/'+Subject_ID+'.mat'
    mat_input = sio.loadmat(fdata)

    rbc2gas = mat_input['rbc2gas']
    bar2gas = mat_input['bar2gas']
    ventilation = mat_input['ventilation']
    mask = mat_input['mask_reg']

    # create an empty pandas table
    columns = ['llow','lup', 'rlow', 'rmid','rup']
    ven_table = pd.DataFrame(columns = columns)
    bar_table = pd.DataFrame(columns = columns)
    rbc_table = pd.DataFrame(columns = columns)

    mask_lobe = gen_lobar_mask(move_filename = move_filename, img_fix = mask).astype(int)

    for lobe_index in range(1,6):

        mean_ven = np.average(ventilation[mask_lobe==lobe_index])
        mean_rbc = np.average(rbc2gas[mask_lobe==lobe_index])
        mean_bar = np.average(bar2gas[mask_lobe==lobe_index])

        ven_table.loc[Subject_ID,columns[lobe_index-1]] = mean_ven
        bar_table.loc[Subject_ID,columns[lobe_index-1]] = mean_bar
        rbc_table.loc[Subject_ID,columns[lobe_index-1]] = mean_rbc

    ven_table.to_csv('ven_table')
    bar_table.to_csv('bar_table')
    rbc_table.to_csv('rbc_table')
    pdb.set_trace()
