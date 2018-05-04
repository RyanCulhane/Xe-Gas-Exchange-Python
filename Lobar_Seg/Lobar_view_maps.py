# plot montage of all the lobar maps and save to a director
# ziyi Apirl 15. 2018

import os, sys, pdb
from GX_Map_utils import fullMontage, decideStartInterval
import scipy.io as sio
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

data_dir = '/home/ziyiw/Patients_old/'
plot_dir = '/home/ziyiw/Gas_Exchange/Lobar_Seg/Plots/'

Subjects_ID = os.listdir(data_dir)

for Subject_ID in Subjects_ID:
  Subject_ID = '005015'
  fdata = data_dir+Subject_ID+'/'+Subject_ID+'.mat'

  mat_input = sio.loadmat(fdata)
  try:
      lobe_mask = mat_input['mask_lobe']
      # lobe_mask = mat_input['mask_reg']
      pdb.set_trace()
      import nibabel as nib
      nii_mask = nib.Nifti1Image(lobe_mask.astype(float), np.eye(4))
      nib.save(nii_mask, 'lobe_mask_005015')
      pdb.set_trace()

      ind_start, ind_inter = decideStartInterval(lobe_mask.astype(bool))
      lobe_mask_plot = lobe_mask[:,:,ind_start::ind_inter]
      lobe_mask_plot = lobe_mask_plot[:,:,:16]
      fullMontage(lobe_mask_plot)
      plt.title(Subject_ID)
      plt.savefig(plot_dir+Subject_ID+'_mask')
      plt.close()
  except:
      continue
