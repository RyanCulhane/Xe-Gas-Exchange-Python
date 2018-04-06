import pyfftw
import scipy.signal
import numpy as np
from timeit import Timer
import pdb

# a = pyfftw.empty_aligned((128, 64), dtype='complex128')
# b = pyfftw.empty_aligned((128, 64), dtype='complex128')
#
# a[:] = numpy.random.randn(128, 64) + 1j*numpy.random.randn(128, 64)
# b[:] = numpy.random.randn(128, 64) + 1j*numpy.random.randn(128, 64)
#
# t = Timer(lambda: scipy.signal.fftconvolve(a, b))
#
# print('Time with scipy.fftpack: %1.3f seconds' % t.timeit(number=100))
#
# # Monkey patch fftpack with pyfftw.interfaces.scipy_fftpack
# scipy.fftpack = pyfftw.interfaces.scipy_fftpack
# scipy.signal.fftconvolve(a, b) # We cheat a bit by doing the planning first
#
# # Turn on the cache for optimum performance
# pyfftw.interfaces.cache.enable()
#
# print('Time with monkey patched scipy_fftpack: %1.3f seconds' %
#        t.timeit(number=100))
recont = np.load('recon_t.npy')
reconf = np.load('recon_f.npy')
#
#
# a = pyfftw.empty_aligned(np.shape(recont),dtype='complex128')
# a[:] = recont
# b = pyfftw.interfaces.numpy_fft.ifft(a)


import pyfftw
import time

time_start = time.time()
a = pyfftw.empty_aligned(np.shape(recont), dtype='complex128')
b = pyfftw.empty_aligned(np.shape(recont), dtype='complex128')
fft_object_c = pyfftw.FFTW(a, b, axes=(0,1,2))
recon_ff = fft_object_c(recont)/np.prod(np.shape(recont))
time_end = time.time()
print('FFT using FFTW: '+str(time_end-time_start))

time_start = time.time()
recon_npf = np.fft.ifftn(recont)
time_end = time.time()
print('FFT using nump: '+str(time_end-time_start))

pdb.set_trace()
