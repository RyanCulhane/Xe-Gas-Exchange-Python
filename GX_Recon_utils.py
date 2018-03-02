from scipy.stats import norm
import scipy.sparse as sps
import numpy as np
import pdb
import ctypes as ct
import os
from numpy.ctypeslib import ndpointer

def sparse_gridding_c(traj, kernel_para, matrix_size, force_dim):

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

    _sparse = ct.CDLL(os.getcwd()+'/libsparse.so')

    _sparse.sparse_gridding_distance.argtypes = (\
     ct.POINTER(ct.c_double),ct.c_double, ct.c_uint, ct.c_uint, ct.POINTER(ct.c_uint),\
     ct.POINTER(ct.c_double),ct.POINTER(ct.c_double),ct.POINTER(ct.c_double),\
     ct.POINTER(ct.c_uint),ct.c_uint,ct.c_int)

     # construct a wrapper function for the c code function
    def sparse_gridding(traj, kernel_para, matrix_size, force_dim):

        # global _sparse

        npts,ndim = np.shape(traj)
        # flatten traj to a list for input
        traj = traj.flatten().tolist()
        kernel_para = kernel_para.flatten()[0]
        matrix_size = matrix_size.astype(int).flatten().tolist()

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

        # set_result to return an array of numbers
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

    sample_indices, voxel_indices, distances = sparse_gridding(traj, kernel_para, matrix_size, force_dim)

    return sample_indices, voxel_indices, distances

def gen_traj_c(num_projs, traj_type):
    # generate xyz coordinates for the trajectory samples based on the number of projs and traj type
    # traj_type:
    # 1: Spiral
    # 2. Halton
    # 3. Haltonized Spiral
    # 4. ArchimedianSeq
    # 5. Double Golden Mean

    # output 3 coordinates of the trajectory points
    output_size = 3*num_projs

    _traj = ct.CDLL(os.getcwd()+'/libtraj.so')
    _traj.gen_traj.argtypes = (ct.c_long,ct.c_long)
    _traj.gen_traj.restype = ndpointer(dtype=ct.c_double, shape=(output_size,))

    result = _traj.gen_traj(ct.c_long(num_projs), ct.c_long(traj_type))

    x = result[:num_projs]
    y = result[num_projs:2*num_projs]
    z = result[2*num_projs:3*num_projs]

return x, y, z
