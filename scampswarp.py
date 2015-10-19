import os
from numpy import *
import subprocess
import docopt

"""
scampswarp - creates config files for scamp and swarp and contains a 
             a function to perform warp on image to improve astrometry
Requires the following packages: os, numpy, subprocess, docopt
Requires the following files: check_apassoverlap.py

usage: scampswarp [-h]

Options:
    -h, --help              Show this screen
    -v, --verbose           If True, print extra text

"""

catalogchecks = {}
catalogchecks['APASS'] = 'check_apassoverlap.py'

arguments = docopt.docopt(__doc__)

VERBOSE = arguments['--verbose']

sextractor_config = '''
    ANALYSIS_THRESH 1.5
	BACK_FILTERSIZE 3
	BACKPHOTO_TYPE GLOBAL
	BACK_SIZE 128
	CATALOG_NAME test.cat
	CATALOG_TYPE FITS_LDAC
	CLEAN Y
	CLEAN_PARAM 1.
	DEBLEND_MINCONT 0.005
	DEBLEND_NTHRESH 32
	DETECT_MINAREA 5
	DETECT_THRESH {detect_thresh}
	DETECT_TYPE CCD
	FILTER Y
	FILTER_NAME {filter_name}
	FLAG_IMAGE flag.fits
	GAIN 0.0
	MAG_GAMMA 4.
	MAG_ZEROPOINT 0.0
	MASK_TYPE CORRECT
	MEMORY_BUFSIZE 1024
	MEMORY_OBJSTACK 3000
	MEMORY_PIXSTACK 300000
	PARAMETERS_NAME {parameters_name}
	PHOT_APERTURES 5
	PHOT_AUTOPARAMS 2.5, 3.5
	PIXEL_SCALE 2.0
	SATUR_LEVEL 50000.
	SEEING_FWHM 1.5
	STARNNW_NAME {starnnw_name}
	VERBOSE_TYPE {verbose_type}
'''

scamp_params = """XWIN_IMAGE
YWIN_IMAGE

ERRAWIN_IMAGE
ERRBWIN_IMAGE
ERRTHETAWIN_IMAGE

XWIN_IMAGE
YWIN_IMAGE

FLUX_AUTO

FLUXERR_AUTO

FLAGS
#FLAGS_WEIGHT
#IMAFLAGS_ISO
FLUX_RADIUS
"""

default_conv = """CONV NORM
# 3x3 ``all-ground'' convolution mask with FWHM = 2 pixels.
1 2 1
2 4 2
1 2 1
"""

default_nnw = """NNW
# Neural Network Weights for the SExtractor star/galaxy classifier (V1.3)
# inputs:	9 for profile parameters + 1 for seeing.
# outputs:	``Stellarity index'' (0.0 to 1.0)
# Seeing FWHM range: from 0.025 to 5.5'' (images must have 1.5 < FWHM < 5 pixels)
# Optimized for Moffat profiles with 2<= beta <= 4.

 3 10 10  1

-1.56604e+00 -2.48265e+00 -1.44564e+00 -1.24675e+00 -9.44913e-01 -5.22453e-01  4.61342e-02  8.31957e-01  2.15505e+00  2.64769e-01
 3.03477e+00  2.69561e+00  3.16188e+00  3.34497e+00  3.51885e+00  3.65570e+00  3.74856e+00  3.84541e+00  4.22811e+00  3.27734e+00

-3.22480e-01 -2.12804e+00  6.50750e-01 -1.11242e+00 -1.40683e+00 -1.55944e+00 -1.84558e+00 -1.18946e-01  5.52395e-01 -4.36564e-01 -5.30052e+00
 4.62594e-01 -3.29127e+00  1.10950e+00 -6.01857e-01  1.29492e-01  1.42290e+00  2.90741e+00  2.44058e+00 -9.19118e-01  8.42851e-01 -4.69824e+00
-2.57424e+00  8.96469e-01  8.34775e-01  2.18845e+00  2.46526e+00  8.60878e-02 -6.88080e-01 -1.33623e-02  9.30403e-02  1.64942e+00 -1.01231e+00
 4.81041e+00  1.53747e+00 -1.12216e+00 -3.16008e+00 -1.67404e+00 -1.75767e+00 -1.29310e+00  5.59549e-01  8.08468e-01 -1.01592e-02 -7.54052e+00
 1.01933e+01 -2.09484e+01 -1.07426e+00  9.87912e-01  6.05210e-01 -6.04535e-02 -5.87826e-01 -7.94117e-01 -4.89190e-01 -8.12710e-02 -2.07067e+01
-5.31793e+00  7.94240e+00 -4.64165e+00 -4.37436e+00 -1.55417e+00  7.54368e-01  1.09608e+00  1.45967e+00  1.62946e+00 -1.01301e+00  1.13514e-01
 2.20336e-01  1.70056e+00 -5.20105e-01 -4.28330e-01  1.57258e-03 -3.36502e-01 -8.18568e-02 -7.16163e+00  8.23195e+00 -1.71561e-02 -1.13749e+01
 3.75075e+00  7.25399e+00 -1.75325e+00 -2.68814e+00 -3.71128e+00 -4.62933e+00 -2.13747e+00 -1.89186e-01  1.29122e+00 -7.49380e-01  6.71712e-01
-8.41923e-01  4.64997e+00  5.65808e-01 -3.08277e-01 -1.01687e+00  1.73127e-01 -8.92130e-01  1.89044e+00 -2.75543e-01 -7.72828e-01  5.36745e-01
-3.65598e+00  7.56997e+00 -3.76373e+00 -1.74542e+00 -1.37540e-01 -5.55400e-01 -1.59195e-01  1.27910e-01  1.91906e+00  1.42119e+00 -4.35502e+00

-1.70059e+00 -3.65695e+00  1.22367e+00 -5.74367e-01 -3.29571e+00  2.46316e+00  5.22353e+00  2.42038e+00  1.22919e+00 -9.22250e-01 -2.32028e+00


 0.00000e+00 
 1.00000e+00 
"""   

def scampswarp(image,interptype = 'LANCZOS3',cattype = 'APASS'):
    """
    image:       name of image to scamped and swarped
    intertype:   interpolation method

    Updates image and resuffixes it with _ss.fits - saves
    in original image location.
    """
    # Extract directory
    bits = image.split('/')
    name = bits.pop(-1)
    imdir = '/'.join(bits)+'/'
    # Choose output image name
    imageout = name.split('.fits')[0]+'_ss.fits'
    # Choose catalog name
    catalog = image.split('.fits')[0]+'.cat'
    # Check that the chosen catalog overlaps with data
    os.system('python '+catalogchecks[cattype]+' -v -s /opt/sextractor/2.8.6/bin/sex -r /mnt/scratch-lustre/njones/SURP2015/nscripts/APASS save '+image)
    # Run sextractor
    os.system('sex -c default.sex '+image)
    # rename catalog
    os.system('mv test.cat '+catalog)
    # Run scamp
    os.system('scamp -c scamp.default '+catalog)
    # Run swarp
    os.system('swarp -c swarp.default '+image)
    # Rename output image
    os.system('mv coadd.fits '+imdir+imageout)
    # Check overlap with catalog again
    os.system('python '+catalogchecks[cattype]+' -v -s /opt/sextractor/2.8.6/bin/sex -r /mnt/scratch-lustre/njones/SURP2015/nscripts/APASS save '+imdir+imageout)



if __name__ == '__main__': 
    # Choose background sigma detection threshold
    detect_thresh = 50
    # Choose config file names. These are the default names, changing them
    # means updating the config files with new file names
    sextractor_config_name = "scamp.sex"
    params_name = "scamp.param"
    nnw_name = "default.nnw"
    conv_name = "default.conv"
    scamp_config_name = "default.scamp"
    swarp_config_name = "default.swarp"

    if VERBOSE:
        verbose_type = "NORMAL"

    # Write config files for sextractor
    fp = open(sextractor_config_name, "w+")
    fp.write(sextractor_config.format(detect_thresh=detect_thresh, filter_name=conv_name,
        parameters_name=params_name, starnnw_name=nnw_name, verbose_type=verbose_type))
    fp.close()
    fp = open(params_name, "w+")
    fp.write(scamp_params)
    fp.close()
    fp = open(conv_name, "w+")
    fp.write(default_conv)
    fp.close()
    fp = open(nnw_name, "w+")
    fp.write(default_nnw)
    fp.close()
    # Create a default SCAMP config file
    subprocess.call("scamp -d > {config}".format(config=scamp_config_name), shell=True)

    # Create a default SWARP config file
    subprocess.call("swarp -d > {config}".format(config=swarp_config_name), shell=True)
