import numpy as np
import scipy.io as sio
import warnings
from matplotlib import pyplot as plt
import pdb

## Input **************************
fdata = 'bin_data.mat'

mat_input = sio.loadmat(fdata)
bin_color = mat_input['binning']
ute = mat_input['ute']
ute = np.divide(ute, np.amax(ute))

ute = np.flip(np.flip(np.flip(ute,1),0),2).astype('float64')
bin_color = np.flip(np.flip(np.flip(bin_color,1),0),2).astype('float64')

from GX_utils import makeMontage
from GX_defineColormaps import long_index2color

img_montage = makeMontage(bin_index = bin_color,
                          ute_reg = ute,
                          index2color = long_index2color,
                          ind_start = 22,
                          ind_inter = 5)

# my_dpi = 182
# plt.figure(figsize=(128*8/my_dpi, 128*2/my_dpi), dpi=my_dpi)
plt.imshow(img_montage)
plt.axis('off')
plt.savefig("testMontage.png",transparent = True,bbox_inches='tight',pad_inches=-0.1)
# pdb.set_trace()
