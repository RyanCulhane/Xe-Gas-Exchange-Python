import numpy as np
import pdb
import time
import re

from GX_Recon_utils import generate_traj, read_twix_hdr, recon
from twix_parser import readTwix
## input
twix_file = '005012_BHUTE.dat'
###################################################################
## read in data
scans, evps = readTwix(twix_file)

data = np.asarray([x.data for x in scans])

# parse hdr "MEAS" for more information
meas_dict = read_twix_hdr(evps[2][1])

npts = np.shape(data)[1]
nFrames = np.shape(data)[0]

gen_traj_dict = {
    'npts': npts,
    'nFrames':nFrames,
    'traj_type':3, #halton Spiral
    'dwell_time': float(meas_dict['alDwellTime'].split()[0])/1000.0,
    'oversampling': 2,
    'ramp_time': float(meas_dict['RORampTime']),
    'plat_time': 2500,
    'decay_time': 60,
    'del_x': 0,
    'del_y': 0,
    'del_z': 0,
}

x, y, z = generate_traj(**gen_traj_dict)

def vectorize(x):
    return np.reshape(x,(npts*nFrames,1))

traj = np.squeeze(0.5*np.stack((vectorize(x), vectorize(y), vectorize(z)), axis=-1))

kernel_sharpness = 0.15
kernel_extent = 7*kernel_sharpness

recon_dict = {
    'traj': traj,
    'data': np.reshape(data,(npts*nFrames,1)),
    'kernel_sharpness': kernel_sharpness,
    'kernel_extent': kernel_extent,
    'overgrid_factor': 3,
    'n_pipe_iter': 20,
    'image_size': (npts, npts, npts),
    'verbosity': 1,
}

reconVol = recon(**recon_dict)
