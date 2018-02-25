import ctypes as ct
import pdb
import scipy.io as sio
import numpy as np
from numpy.ctypeslib import ndpointer
import time

fdata = 'gridding_inputoutput.mat'

mat_input = sio.loadmat(fdata)

input1 = mat_input['traj']
input2 = mat_input['input2_kernelExtent_overgridF']
input3 = mat_input['input3_uint32_matrixS']
input4 = mat_input['input4']
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

_sparse = ct.CDLL('/home/ziyi/Gas_Exchange/libsparse.so')

_sparse.sparse_gridding_distance.argtypes = (\
 ct.POINTER(ct.c_double),ct.c_double, ct.c_uint, ct.c_uint, ct.POINTER(ct.c_uint),\
 ct.POINTER(ct.c_double),ct.POINTER(ct.c_double),ct.POINTER(ct.c_double),\
 ct.POINTER(ct.c_uint),ct.c_uint,ct.c_int)

def sparse_gridding(traj, kernel_para, matrix_size, force_dim):
    global _sparse

    npts,ndim = np.shape(traj)
    # flatten traj to a list for input
    traj = traj.flatten().tolist()
    kernel_para = kernel_para.flatten()[0]
    matrix_size = matrix_size.flatten().tolist()
    force_dim = force_dim.flatten()[0]

    num_coord = len(traj)
    num_matrixsize = len(matrix_size)

    # calculate max size of the output indices
    max_nNeighbors = 1
    for dim in range(0,ndim):
        max_nNeighbors = int(max_nNeighbors*(kernel_para+1))

    max_size = npts*max_nNeighbors

    # create empty output
    sample_indices = [0.0] * max_size
    voxel_indices = [0.0] * max_size
    distances = [0.0] * max_size
    nSparsePoints = [0] * 1

    # define argument types
    coord_type = ct.c_double * num_coord
    outputsize_type = ct.c_uint * num_matrixsize
    outindices_type = ct.c_double * max_size
    n_nonsparse_entries_type = ct.c_uint * 1

    # set_result
    _sparse.sparse_gridding_distance.restype = ndpointer(dtype=ct.c_double, shape=(max_size*3,))

    result = _sparse.sparse_gridding_distance(\
    coord_type(*traj), ct.c_double(kernel_para), ct.c_uint(npts),\
    ct.c_uint(ndim), outputsize_type(*matrix_size),\
    outindices_type(*sample_indices), outindices_type(*voxel_indices),\
    outindices_type(*distances), \
    n_nonsparse_entries_type(*nSparsePoints),\
    ct.c_uint(max_size),ct.c_int(force_dim)
    )

    sample_indices = result[:max_size]
    voxel_indices = result[max_size:2*max_size]
    distances = result[2*max_size:3*max_size]
    return sample_indices, voxel_indices, distances

time_start = time.time()
sample_indices, voxel_indices, distances = sparse_gridding(input1, input2, input3, input4)
time_end = time.time()
print(time_end-time_start)
pdb.set_trace()
