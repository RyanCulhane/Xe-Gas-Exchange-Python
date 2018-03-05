import numpy as np
import pdb
import time
import re

from GX_Recon_classmap import Gaussian, L2Proximity, MatrixSystemModel, IterativeDCF, LSQgridded
from GX_Recon_utils import gen_traj_c, generate_radial_1D_traj, read_twix_hdr
from twix_parser import readTwix
## input
twix_file = '005012_BHUTE.dat'
###################################################################
## read in data
scans, evps = readTwix(twix_file)

data = np.asarray([x.data for x in scans])

# parse hdr "MEAS" for more information
meas_dict = read_twix_hdr(evps[2][1])

dwell_time = float(meas_dict['alDwellTime'].split()[0])/1000.0
ramp_time = float(meas_dict['RORampTime'])
npts = np.shape(data)[1]
nFrames = np.shape(data)[0]

traj_para = {
    'npts': npts,
    'dwell_time': dwell_time,
    'oversampling': 2,
    'ramp_time': ramp_time,
    'plat_time': 2500,
    'decay_time': 60,
}

verbosity = 1
kernel_sharpness = 0.15
kernel_extent = 7*kernel_sharpness
overgrid_factor = 3
n_pipe_iter = 20
traj_type = 3 # 3 for halton_spiral
image_size = (npts, npts, npts)

del_x = 0
del_y = 0
del_z = 0

# Generate and vectorize traj and data
if(del_x == del_y == del_z):
    traj_para.update({'grad_delay_time':del_x})
    radial_distance_x = generate_radial_1D_traj(**traj_para)
    radial_distance_y = radial_distance_x
    radial_distance_z = radial_distance_x
else:
    traj_para.update({'grad_delay_time':del_x})
    radial_distance_x = generate_radial_1D_traj(**traj_para)
    traj_para.update({'grad_delay_time':del_y})
    radial_distance_y = generate_radial_1D_traj(**traj_para)
    traj_para.update({'grad_delay_time':del_z})
    radial_distance_z = generate_radial_1D_traj(**traj_para)

x, y, z  = gen_traj_c(nFrames, traj_type)

x = np.array([radial_distance_x]).transpose().dot(np.array([x])).transpose()
y = np.array([radial_distance_y]).transpose().dot(np.array([y])).transpose()
z = np.array([radial_distance_z]).transpose().dot(np.array([z])).transpose()

def vectorize(x):
    return np.reshape(x,(npts*nFrames,1))

traj = np.squeeze(0.5*np.stack((vectorize(x), vectorize(y), vectorize(z)), axis=-1))

data = np.reshape(data,(npts*nFrames,1))

## Starting reconstruction
kernel_obj = Gaussian(kernel_extent = kernel_extent, kernel_sigma = kernel_sharpness, verbosity = verbosity)

prox_obj = L2Proximity(kernel_obj= kernel_obj, verbosity = verbosity)

system_obj = MatrixSystemModel(proximity_obj=prox_obj, overgrid_factor=overgrid_factor, image_size=image_size, traj=traj, verbosity=verbosity)

dcf_obj = IterativeDCF(system_obj = system_obj, dcf_iterations = n_pipe_iter, verbosity = verbosity)

recon_obj = LSQgridded(system_obj = system_obj, dcf_obj = dcf_obj, verbosity = verbosity)

reconVol = recon_obj.reconstruct(data = data, traj = traj)
