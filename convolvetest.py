from scipy import signal
import matplotlib.pyplot as plt
from scipy import misc
import numpy as np
import pyfits as pf
plt.ion()

constant = 2#2*np.sqrt(2*np.log(2))   # constant used to change FWHM to sigma

# assuming the units in arcsec for resolution and pixelsize 
pixelsize = 1

FWHM1 = 40/pixelsize 	# for initial
sigma1 = FWHM1/constant
FWHM2 = 55/pixelsize 	# for expected
sigma2 = FWHM2/constant

inital = np.outer(signal.general_gaussian(100,1, sigma1), signal.general_gaussian(100,1, sigma1))
Expected = np.outer(signal.general_gaussian(100,1, sigma2), signal.general_gaussian(100,1, sigma2))

# FWHM and sigma for kernel
FWHM = np.sqrt(FWHM2**2 - FWHM1**2)
sigma = FWHM/(constant)

# making the kernel map
kernel = np.outer(signal.general_gaussian(100,1, sigma), signal.general_gaussian(100,1, sigma))/((sigma)*np.sqrt(2*np.pi))

# convolved image
convolved = signal.fftconvolve(inital, kernel, mode='same')

# plotting
fig, (ax_orig, ax_kernel, ax_convolved, expected) = plt.subplots(1, 4)
ax_orig.imshow(inital, cmap='gray')
ax_orig.set_title('Original')
ax_orig.set_axis_off()
ax_kernel.imshow(kernel, cmap='gray')
ax_kernel.set_title('Gaussian kernel')
ax_kernel.set_axis_off()
ax_convolved.imshow(convolved, cmap='gray')
ax_convolved.set_title('convolved')
ax_convolved.set_axis_off()
expected.imshow(Expected, cmap='gray')
expected.set_title('Expected')
expected.set_axis_off()
fig.show()

plt.figure()
plt.plot(convolved[50]/np.max(convolved[50]), label = 'convolved')
plt.plot(Expected[50]/np.max(Expected[50]), label = 'expected')
plt.legend()
