from GX_utils import fullMontage
from scipy.ndimage.morphology import binary_dilation
import scipy.io as sio
import nibabel as nib
import numpy as np
import math
import pdb

##### input
fdata = 'Sub002102_data.mat'
mat_input = sio.loadmat(fdata)
image = abs(mat_input['gasVol_highreso'])

fmask = 'BHUTE_Sub002102_FID49886_mask_grow.nii'
mask = np.array(nib.load(fmask).get_data()).astype(bool)

#####
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

pdb.set_trace()
