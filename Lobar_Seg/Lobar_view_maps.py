# plot montage of all the lobar maps and save to a director
# ziyi Apirl 15. 2018

import os, sys, pdb
from GX_Map_utils import fullMontage, decideStartInterval
import scipy.io as sio
import matplotlib.pyplot as plt
import numpy as np

data_dir = '/home/ziyiw/Patients/'
plot_dir = '/home/ziyiw/Gas_Exchange/Lobar_Seg/Plots/'

Subjects_ID = os.listdir(data_dir)

for Subject_ID in Subjects_ID:
  fdata = data_dir+Subject_ID+'/'+Subject_ID+'.mat'
  
  mat_input = sio.loadmat(fdata)
  lobe_mask = mat_input['mask_lobe']
  
  ind_start, ind_inter = decideStartInterval(lobe_mask.astype(bool))
  lobe_mask_plot = lobe_mask[:,:,ind_start::ind_inter]
  lobe_mask_plot = lobe_mask_plot[:,:,:16]
  pdb.set_trace()
  fullMontage(lobe_mask_plot)
  plt.savefig(plot_dir+Subject_ID)
  

  # yingxuan
