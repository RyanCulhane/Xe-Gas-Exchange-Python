import ctypes as ct
import pdb
import scipy.io as sio
import numpy as np
from numpy.ctypeslib import ndpointer
import time
import os

# fdata = 'gridding_inputoutput.mat'
#
# mat_input = sio.loadmat(fdata)
#
# input1 = mat_input['traj']
# input2 = mat_input['input2_kernelExtent_overgridF']
# input3 = mat_input['input3_uint32_matrixS']
# input4 = mat_input['input4']
# pdb.set_trace()
#wrap up c function:
# void sparse_gridding_distance(double *coords, double kernel_width,
# 		unsigned int npts, unsigned int ndims,
# 		unsigned int *output_dims,
# 		double *nonsparse_sample_indices,
# 		double  *nonsparse_voxel_indices,
# 		double *nonsparse_distances,
# 		unsigned int *n_nonsparse_entries,
# 		unsigned int max_size,
# 		int force_dim)

_traj = ct.CDLL(os.getcwd()+'/libtraj.so')
_traj.gen_traj.argtypes = (ct.c_long,ct.c_long)

# output 3 coordinates of the trajectory points
num_projs = 5
traj_type = 5

output_size = 3*num_projs

_traj.gen_traj.restype = ndpointer(dtype=ct.c_double, shape=(output_size,))

result = _traj.gen_traj(ct.c_long(num_projs), ct.c_long(traj_type))

x = result[:num_projs]
y = result[num_projs:2*num_projs]
z = result[2*num_projs:3*num_projs]

print x
print y
print z
pdb.set_trace()
