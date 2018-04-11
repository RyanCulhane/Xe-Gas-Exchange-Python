from scipy.stats import norm
import scipy.sparse as sps
import numpy as np
import pdb
import ctypes as ct
import os
import re
from numpy.ctypeslib import ndpointer
from GX_Twix_parser import readTwix
from matplotlib import pyplot as plt

from GX_Benchmark_utils import timerfunc # decorator for timing a function
from memory_profiler import profile # decorator to monitor memory usage

def generate_radial_1D_traj(dwell_time, grad_delay_time, ramp_time, plat_time, decay_time, npts,oversampling):
    # generate 1d radial distance array based on the timing and the amplitude and the gradient
    grad_delay_npts = grad_delay_time/dwell_time
    ramp_npts = ramp_time/dwell_time
    plat_npts = plat_time/dwell_time
    decay_npts = decay_time/dwell_time

    pts_vec = np.array(range(0, npts))

    # assume a isotropic recon
    def calcRadDist(out_size):
        # calculate sample number of each region boundary
        ramp_start_pt = grad_delay_npts
        plat_start_pt = ramp_start_pt+ramp_npts
        decay_start_pt = plat_start_pt+plat_npts
        decay_end_pt = decay_start_pt+decay_npts

        # calculate binary mask for each region
        in_delay = pts_vec<ramp_start_pt
        in_ramp = (pts_vec>=ramp_start_pt) & (pts_vec<plat_start_pt)
        in_plat = (pts_vec>=plat_start_pt) & (pts_vec<decay_start_pt)
        in_decay = (pts_vec>=decay_start_pt) & (pts_vec<decay_end_pt)

        # calculate times in each region
        ramp_pts_vec = np.multiply((pts_vec-ramp_start_pt),in_ramp)
        plat_pts_vec = np.multiply((pts_vec-plat_start_pt),in_plat)
        decay_pts_vec = np.multiply((pts_vec-decay_start_pt),in_decay)

        # calculate the gradient amplitude  over time(assume plateau is 1)
        ramp_g = ramp_pts_vec/ramp_npts
        plat_g = in_plat
        decay_g = np.multiply((1.0-decay_pts_vec/decay_npts),in_decay)

        # calculate radial position (0.5)
        ramp_dist = 0.5*np.multiply(ramp_pts_vec,ramp_g)
        plat_dist = 0.5*ramp_npts*in_plat + np.multiply(plat_pts_vec,plat_g)
        decay_dist = (0.5*ramp_npts+plat_npts)*in_decay + np.multiply(in_decay,np.multiply(decay_pts_vec*0.5,(1.0+decay_g)))

        radial_distance = (ramp_dist+plat_dist+decay_dist)/out_size

        return radial_distance

    return calcRadDist(npts)

# @profile
def sparse_gridding_c(traj, kernel_para, matrix_size, force_dim):

    #wrap up c function:
    # void sparse_gridding_distance(double *coords, double kernel_width,
    # 		unsigned int npts, unsigned int ndims,
    # 		unsigned int *output_dims,
    # 		unsigned int *n_nonsparse_entries,
    # 		unsigned int max_size,
    # 		int force_dim)

    _sparse = ct.CDLL(os.getcwd()+'/libsparse.so')

    _sparse.sparse_gridding_distance.argtypes = (\
     ct.POINTER(ct.c_double),ct.c_double, ct.c_uint, ct.c_uint, ct.POINTER(ct.c_uint),\
     ct.POINTER(ct.c_uint),ct.c_uint,ct.c_int)


    npts,ndim = np.shape(traj)
    # flatten traj to a list for input
    traj = traj.flatten().tolist()
    kernel_para = kernel_para
    matrix_size = matrix_size.astype(int).flatten().tolist()

    num_coord = len(traj)
    num_matrixsize = len(matrix_size)

    # calculate max size of the output indices
    max_nNeighbors = 1
    for dim in range(0,ndim):
        max_nNeighbors = int(max_nNeighbors*(kernel_para+1))

    max_size = npts*max_nNeighbors

    # create empty output
    nSparsePoints = [0] * 1

    # define argument types
    coord_type = ct.c_double * num_coord
    outputsize_type = ct.c_uint * num_matrixsize
    n_nonsparse_entries_type = ct.c_uint * 1

    # set_result to return an array of numbers
    _sparse.sparse_gridding_distance.restype = ndpointer(dtype=ct.c_double, shape=(max_size*3,))

    result = _sparse.sparse_gridding_distance(\
    coord_type(*traj), ct.c_double(kernel_para), ct.c_uint(npts),\
    ct.c_uint(ndim), outputsize_type(*matrix_size),\
    n_nonsparse_entries_type(*nSparsePoints),\
    ct.c_uint(max_size),ct.c_int(force_dim)
    )

    sample_indices = result[:max_size]
    voxel_indices = result[max_size:2*max_size]
    distances = result[2*max_size:3*max_size]

    return sample_indices, voxel_indices, distances

# @profile
@timerfunc
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

def read_twix_hdr(bufferr):
    # parse the buffer string into a dictionary
    # remove empty lines
    p = re.compile('\n\s*\n')
    bufferr = p.sub('',bufferr)

    # split ascconv and xprotocco
    p = re.compile('### ASCCONV BEGIN[^\n]*\n(.*)\s### ASCCONV END ###')
    split_list = p.split(bufferr,1)

    # just take xprot and forget about ascconv
    if(len(split_list) == 1):
        xprot = split_list[0]
    elif(len(split_list)==2):
        # ascconv is not parsed at this moment
        ascconv = split_list[0]
        xprot = split_list[1]
    else:
        raise Exception('Twix file has multiple Ascconv')

    # parse xprot
    p = re.compile('<Param(?:Bool|Long|String)\."(\w+)">\s*{([^}]*)')
    token_list = p.findall(xprot)

    p = re.compile('<ParamDouble\."(\w+)">\s*{\s*(<Precision>\s*[0-9]*)?\s*([^}]*)')
    token_list = token_list + p.findall(xprot)

    # parse into name + value and return
    twix_dict = {}

    for token in token_list:
        name = token[0]

        p = re.compile('("*)|( *<\w*> *[^\n]*)')
        value = p.sub(token[-1],'').strip()

        p = re.compile('\s*')
        value = p.sub(value,' ')

        twix_dict[name] = value

    return(twix_dict)

def generate_traj(dwell_time,ramp_time, plat_time, decay_time, npts, oversampling, del_x, del_y, del_z, nFrames, traj_type):
    # Generate and vectorize traj and data
    traj_para = {
        'npts': npts,
        'dwell_time': dwell_time,
        'oversampling': oversampling,
        'ramp_time': ramp_time,
        'plat_time': plat_time,
        'decay_time': decay_time,
    }

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

    return x, y, z

def remove_noise_rays(data, x, y, z, thre_snr):
    # remove noisy FID rays in Dixon image
    nFrames = np.shape(data)[0]
    thre_dis = thre_snr*np.average(abs(data[:,:5]))
    max_tail = np.amax(abs(data[:,10:]),axis=1)
    good_index = max_tail<thre_dis

    n_Frames_good = np.sum(good_index)
    data = data[good_index]
    x = x[good_index]
    y = y[good_index]
    z = z[good_index]

    return data, x, y, z, n_Frames_good

def complex_align(x):
    return(np.flip(np.flip(np.flip(np.transpose(x,(2,1,0)),2),1),0))

def recon(data, traj, kernel_sharpness, kernel_extent, overgrid_factor, image_size, n_pipe_iter, verbosity):

    from GX_Recon_classmap import Gaussian, L2Proximity, MatrixSystemModel, IterativeDCF, LSQgridded

    ## Starting reconstruction
    kernel_obj = Gaussian(kernel_extent = kernel_extent, kernel_sigma = kernel_sharpness, verbosity = verbosity)

    prox_obj = L2Proximity(kernel_obj= kernel_obj, verbosity = verbosity)

    system_obj = MatrixSystemModel(proximity_obj=prox_obj, overgrid_factor=overgrid_factor, image_size=image_size, traj=traj, verbosity=verbosity)

    dcf_obj = IterativeDCF(system_obj = system_obj, dcf_iterations = n_pipe_iter, verbosity = verbosity)

    recon_obj = LSQgridded(system_obj = system_obj, dcf_obj = dcf_obj, verbosity = verbosity)

    reconVol = recon_obj.reconstruct(data = data, traj = traj)

    return(reconVol)

# @profile
@timerfunc
def recon_ute(twix_file):
    # recon the ute file, input the path of the Siemens twix file
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
        'oversampling': 3,
        'ramp_time': float(meas_dict['RORampTime']),
        'plat_time': 2500,
        'decay_time': 60,
        'del_x': 0,
        'del_y': 0,
        'del_z': 0,
    }

    x, y, z = generate_traj(**gen_traj_dict)

    def vectorize(x):
        return np.reshape(x,(np.prod(np.shape(x)),1))

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
        'verbosity': 0,
    }
    print('Starting recon UTE')
    uteVol = recon(**recon_dict)

    uteVol = np.transpose(uteVol,(2,1,0))
    # uteVol = complex_align(uteVol)

    return(uteVol)

@timerfunc
def recon_dixon(twix_file):
    ## recon_dixon images
    ###################################################################
    ## read in data
    scans, evps = readTwix(twix_file)

    data_dixon = np.asarray([x.data for x in scans[:-2]]) # the last 2 are spectrums

    data_spect = np.asarray(scans[-1].data)
    data_dis = data_dixon[3::2,:]
    data_gas = data_dixon[2::2,:] # the first gas is contaminated, so we throw it away

    # parse hdr "MEAS" for more information
    meas_dict = read_twix_hdr(evps[2][1])

    npts = np.shape(data_dixon)[1]
    nFrames = np.shape(data_dixon)[0]+2 # the last 2 are spectrums
    TE90 = float(meas_dict['alTE'].split()[0])

    gen_traj_dict = {
        'npts': npts,
        'nFrames':np.floor(nFrames/2).astype(int),
        'traj_type':3, #halton Spiral
        'dwell_time': float(meas_dict['alDwellTime'].split()[0])/1000.0,
        'oversampling': 3,
        'ramp_time': float(meas_dict['RORampTime']),
        'plat_time': 2500,
        'decay_time': 60,
        'del_x': 24-13,
        'del_y': 24-14,
        'del_z': 24-9,
    }

    x, y, z = generate_traj(**gen_traj_dict)

    # the last 2 are used for spectroscopy
    x = x[1:-1:,:]
    y = y[1:-1:,:]
    z = z[1:-1:,:]

    def vectorize(x):
        return np.reshape(x,(np.prod(np.shape(x)),1))

    data_dis, x_dis, y_dis, z_dis, nFrames_dis = remove_noise_rays(data = data_dis, x = x, y = y, z = z, thre_snr = 0.5)
    print("Threw away bad dissolved FID: "+str(nFrames/2-2-nFrames_dis))
    data_gas, x_gas, y_gas, z_gas, nFrames_gas = remove_noise_rays(data = data_gas, x = x, y = y, z = z, thre_snr = 0.5)
    print("Threw away bad gas FID: "+str(nFrames/2-2-nFrames_gas))

    ## recon gas for high SNR and high resolution
    traj = np.squeeze(0.5*np.stack((vectorize(x_gas), vectorize(y_gas), vectorize(z_gas)), axis=-1))
    recon_dict_highSNR = {
        'traj': traj,
        'data': np.reshape(data_gas,(npts*nFrames_gas,1)),
        'kernel_sharpness': 0.14,
        'kernel_extent': 9*0.14,
        'overgrid_factor': 3,
        'n_pipe_iter': 20,
        'image_size': (npts, npts, npts),
        'verbosity': 0,
    }
    print('Starting recon gas high SNR')
    gasVol_highSNR = recon(**recon_dict_highSNR)

    recon_dict_highreso = {
        'traj': traj,
        'data': np.reshape(data_gas,(npts*nFrames_gas,1)),
        'kernel_sharpness': 0.32,
        'kernel_extent': 9*0.32,
        'overgrid_factor': 3,
        'n_pipe_iter': 20,
        'image_size': (npts, npts, npts),
        'verbosity': 0,
    }
    print('Starting recon gas high resolution')
    gasVol_highreso = recon(**recon_dict_highreso)

    ## recon dissolved for high SNR
    traj = np.squeeze(0.5*np.stack((vectorize(x_dis), vectorize(y_dis), vectorize(z_dis)), axis=-1))
    recon_dict_highSNR['traj'] = traj
    recon_dict_highSNR['data'] = np.reshape(data_dis,(npts*nFrames_dis,1))

    print('Starting recon dissolved')
    dissolvedVol = recon(**recon_dict_highSNR)

    gasVol_highreso = np.transpose(gasVol_highreso,(2,1,0))
    gasVol_highSNR = np.transpose(gasVol_highSNR,(2,1,0))
    dissolvedVol = np.transpose(dissolvedVol,(2,1,0))

    return gasVol_highSNR, gasVol_highreso, dissolvedVol, TE90
