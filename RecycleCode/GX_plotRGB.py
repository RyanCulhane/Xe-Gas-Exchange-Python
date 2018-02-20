import numpy as np
import pdb
import nibabel as nib
from matplotlib import pyplot as plt

indexmap = np.zeros((128,128))

indexmap[:3,:] = 1
indexmap[3:6,:] = 2
indexmap[6:9,:] = 3
indexmap[9:12,:] = 4
indexmap[12:15,:] = 5
indexmap[15:18,:] = 6
indexmap[18:21,:] = 7
indexmap[21:24,:] = 8
indexmap[24:128,:] = 0

from GX_defineColormaps import long_index2color

def index2color(ind):
    return long_index2color[ind]

colormap = map(index2color, indexmap.flatten())

colormap = np.asarray(colormap)
colormap = np.reshape(colormap,(128,128,3))

fute = 'BHUTE_Sub002102_FID49886_recon.nii'

ute = np.array(nib.load(fute).get_data())

ute = ute[:,:,60]

plt.imshow(colormap)
plt.imshow(ute,cmap='gray')
plt.show()
pdb.set_trace()
