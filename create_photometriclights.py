#!/usr/bin/env python

"""create_photometriclights.py -- based on APASS catalog, flatten light frames. 

Usage:
    zeropoint2 [-h] [-v] [-c] [-u] [-p] [-x NUMBER] [-X NUMBER] [-y NUMBER] [-Y NUMBER] [-m NUMBER] [-t NUMBER] [-e ARCSEC] [-d NUMBER] [-f FILTER] [-r DIRECTORY] [-s LOCATION] [-o DIRECTORY] [-k] <catalog>

Options:
    -h, --help                                  Show this screen
    -v, --verbose                               Show extra information [default: False]      
    -e ARCSEC, --error=ARCSEC                   Maximum allowable separation in arcsec [default: 2.2]
    -d NUMBER, --dropsig NUMBER                 Sources with zeropoints outside this number of sigmas of mean is excluded [default: 2.3] 
    -p --plot                                   Plot the data and save in same directory as input catalog or image.
    -c, --cache                                 Create/use APASS catalog cache [default: False]
    -u, --update                                Update image headers with information from the catalog
    -f FILTER, --filter FILTER                  Name of the filter ('g', 'r'); By default, FILTNAM is read from fits file if supplied. 
    -x NUMBER, --xmin NUMBER
    -X NUMBER, --xmax NUMBER
    -y NUMBER, --ymin NUMBER
    -Y NUMBER, --ymax NUMBER
    -m NUMBER, --magmin NUMBER
    -t NUMBER, --threshold NUMBER               Number of sigma to threshold [default: 10]
    -r DIRECTORY, --refcat DIRECTORY            The directory containing the reference APASS catalogs [default: /mnt/gfsproject/naiad/jlzhang/APASS/]
    -s LOCATION, --sex LOCATION                 Location of SExtractor executable [default: /opt/sextractor/2.8.6/bin/sex]
    -o DIRECTORY, --outputdir DIRECTORY         Output directory name  [default: .]
    -k, --kfit                                  Correct for the colour term
    

Examples:
    python zeropoint2.py -v -s C:\\Users\\abrah_000\\Desktop\\test.cat
    python zeropoint2.py -v --gband C:\\Users\\abrah_000\\Desktop\\test.cat
    python zeropoint2.py -v -c --airmasscorrect /Volumes/savehere/ image.fits

"""
import os.path
import math
import sys
import itertools
import subprocess
import re

import docopt
import numpy
import pandas as pd
import scipy.stats
import matplotlib.pyplot as plt
import pylab as P

from scipy import spatial
from astropy.io import ascii
from astropy.io import fits
from astropy.stats import sigma_clip

from datetime import datetime
import astropy.wcs as wcs

import warnings
from astropy.modeling import models, fitting

from scipy.optimize import curve_fit

def print_verbose_string(printme):
    print >> sys.stderr, "VERBOSE: %s" % printme


# Exception functions

class ZeroPointError(Exception):
    """Base class for exceptions in this module."""
    pass

class FilterError(ZeroPointError):
    """Raised when an unknown filter is specified

    Attributes:
        msg -- explanation of the error
    """
    def __init__(self, msg):
        self.msg = msg


# Returns the full path to the APASS file corresponding to a given declination.

def SelectAPASSFile(declination,apass_dir):

    # We only have files for a limited range of declinations for DR8
    #if declination < 20:
    #    return(None)

    # Declinations greater than 90 make no sense
    #if declination > 90:
    #    return(None)

    # Return the filename stub by rounding down to the
    # nearest 5 degrees in declination. So 52.1 => 50,
    # 58.1 => 55, etc.
    file_stub = int(5 * int(float(declination)/5))
    # if dec is less than 10, e.g. 0 or 5, file name is 00 or 05 and not 0 or 5.
    if file_stub <10:
        file_stub = '0'+str(file_stub)

    # We only have files for a limited range of declinations for DR8; dec>90 doesn't make sense
    if declination > 20 and declination < 90:
        apass_file = apass_dir + "zp%s_8.sum" % file_stub

    # We remainder declination range is in DR7; 
    if declination > 0 and declination < 20:
        apass_file = apass_dir + "zp%s_7.sum" % file_stub

    if os.path.isfile(apass_file):
        return(apass_file)
    else:
        return(None)


# Loads an APASS catalog in a specific RA range and return the data as a set of narrays.

def LoadAPASSFile(filename, min_ra, max_ra, verbose=False):

    colspecs = [(0,10),(11,21),(22,28),(29,39),(40,46),(47,51),(52,56), (57,63),(64,70),
        (71,77),(78,84),(85,91),(92,98),(99,105),(106,112),(113,119),(120,126),(127,133),(134,140)]
    names = ['name', 'ra', 'raerr', 'dec', 'decerr', 'nobs', 'mobs', 'V', 'BV', 'B', 'g', 'r', 'i',
        'Verr',  'BVerr', 'Berr',  'gerr',  'rerr','ierr']

    if verbose:
        print_verbose_string("Reading in %s" % filename)

    aacat = pd.read_fwf(filename, colspecs=colspecs, header=None, names=names, comment='#' )
    if verbose:
        print_verbose_string("Number of objects read in: %d" % len(aacat.ra))

    # Cull sample in RA - note that the Pandas Series object has been refactored so it subclasses NDFrame rather than ndarray
    # which is why we need to call .values below. If you are using an old version of Pandas you may need to remove .values. 
    # See: http://stackoverflow.com/questions/21822988/what-series-method-replaced-searchsorted
    indices = []
    if verbose:
        print_verbose_string("Culling data in right ascension...")
    if ( (max_ra - min_ra) > 180.0 ):
        if verbose:
            print "RA range of data 'wraps around' the 24h mark."
        start1 = numpy.searchsorted( aacat.ra.values, max_ra, 'left' )
        end1 = len( aacat.ra ) - 1
        indices1 = numpy.arange( start1, end1 )
        start2 = 0
        end2 = numpy.searchsorted( aacat.ra.values, min_ra, 'right' )
        indices2 = numpy.arange( start2, end2 )
        indices = numpy.concatenate( (indices1,indices2) )
    else:
        start = numpy.searchsorted( aacat.ra.values, min_ra, 'left' )
        end = numpy.searchsorted( aacat.ra.values, max_ra, 'right' )
        indices = numpy.arange( start, end )
    if verbose:
        print "Number of objects remaining: %d" % len(indices)  

    # Convert from NDFrame to ndarray data type
    name = aacat.name[indices].values
    ra = aacat.ra[indices].values
    dec = aacat.dec[indices].values
    g = aacat.g[indices].values
    gerr = aacat.gerr[indices].values
    r = aacat.r[indices].values
    rerr = aacat.rerr[indices].values

    return(name, ra, dec, g, r, gerr, rerr)

# SExtractor configuration files

#CHECKIMAGE_TYPE MINIBACKGROUND FILTERED
sextractor_config = """
    ANALYSIS_THRESH 1.5
    BACK_FILTERSIZE 3
    BACKPHOTO_TYPE GLOBAL
    BACK_SIZE 128
    CATALOG_NAME test.cat
    CATALOG_TYPE ASCII_HEAD
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
    GAIN 1.0
    MAG_GAMMA 4.
    MAG_ZEROPOINT 0.0
    MASK_TYPE CORRECT
    MEMORY_BUFSIZE 4096
    MEMORY_OBJSTACK 30000
    MEMORY_PIXSTACK 3000000
    PARAMETERS_NAME {parameters_name}
    PHOT_APERTURES 5
    PHOT_AUTOPARAMS 2.5, 3.5
    PIXEL_SCALE 2.85
    SATUR_LEVEL 50000.
    SEEING_FWHM 2.5
    STARNNW_NAME {starnnw_name}
    VERBOSE_TYPE {verbose_type}
"""

sextractor_params = """NUMBER
FLUX_AUTO
FLUXERR_AUTO
FLUX_APER
FLUXERR_APER
X_IMAGE
Y_IMAGE
X_WORLD
Y_WORLD
FLUX_RADIUS
FLAGS
CLASS_STAR
MAG_AUTO
MAGERR_AUTO
MAG_ISO
MAGERR_ISO
BACKGROUND
A_IMAGE
B_IMAGE
THETA_IMAGE
THETA_SKY
ISOAREA_IMAGE
FWHM_IMAGE
ISOAREAF_IMAGE
ELLIPTICITY
"""

default_conv = """CONV NORM
# 3x3 ``all-ground'' convolution mask with FWHM = 2 pixels.
1 2 1
2 4 2
1 2 1
"""

default_nnw = """NNW
# Neural Network Weights for the SExtractor star/galaxy classifier (V1.3)
# inputs:   9 for profile parameters + 1 for seeing.
# outputs:  ``Stellarity index'' (0.0 to 1.0)
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


def create_catalog(image_name, detect_thresh=10):

    # Create a config, param, conv, nnw file for Sextractor
    sextractor_config_name = "/tmp/scamp.sex"
    params_name = "/tmp/scamp.param"
    nnw_name = "/tmp/default.nnw"
    conv_name = "/tmp/default.conv"
    catalog_name = "/tmp/cz.cat"
    if verbose:
        verbose_type = "NORMAL"
    else:
        verbose_type = "QUIET"
    fp = open(sextractor_config_name, "w")
    fp.write(sextractor_config.format(detect_thresh=detect_thresh, filter_name=conv_name,
        parameters_name=params_name, starnnw_name=nnw_name, verbose_type=verbose_type))
    fp.close()
    fp = open(params_name, "w")
    fp.write(sextractor_params)
    fp.close()
    fp = open(conv_name, "w")
    fp.write(default_conv)
    fp.close()
    fp = open(nnw_name, "w")
    fp.write(default_nnw)
    fp.close()
    
    backname         = image_name.split('.')[0]+'_miniback.fits'
    filteredbackname = image_name.split('.')[0]+'_filteredback.fits'
    subprocess.call(sex_loc+" -c {config} -CATALOG_NAME {catalog} -CHECKIMAGE_NAME '{back_name} {filteredback_name}' {image}".format(config=sextractor_config_name, catalog=catalog_name, back_name=backname, filteredback_name=filteredbackname, image=image_name), shell=True)

    return catalog_name

def getimsize(image_name):
    h = fits.getheader(image_name)
    naxis1=h['NAXIS1']
    naxis2=h['NAXIS2']
    return naxis1,naxis2

def save_image(image, factor_map, newimage_dir, suffix='pc', savemap=False):
    
    # Divide image by factor map
    old_data    = fits.getdata(image)
    header  = fits.getheader(image)
    new_data    = old_data/factor_map

    # Include all header entries in old image into the new image
    hduP            = fits.PrimaryHDU(new_data)
    hdulist_final   =fits.HDUList(hduP)
    hdulist_final[0].header.update(header)

    # Save photometry corrected frame
    image_basename  = image.split('/')[-1]
    image_basename  = image_basename.split('.')[-2]
    new_image_name  = newimage_dir + '/' + image_basename + '_'+suffix +'.fits' 
    if os.path.exists(new_image_name):
        os.remove(new_image_name)
    hdulist_final.writeto(new_image_name)

    if savemap:
        basename    = image.split('/')[-1]
        # Save factor map
        factor_map_name  = newimage_dir + '/' + basename.split('.')[-2] + '_'+suffix+'map_'+'.fits'
        hduP        = fits.PrimaryHDU(factor_map)
        hdulist_final   =fits.HDUList(hduP)
        hdulist_final[0].header.update(header)
        if os.path.exists(factor_map_name):
            os.remove(factor_map_name)
        if verbose:
            print 'factor map is saved to ' + factor_map_name
        hdulist_final.writeto(factor_map_name)
        
    return new_image_name

def curve_fit_wrapper(m0_mADU,C,y_sig):
    # This application of linear least squares doesn't using weights
    # Fit data, calculate parameters
    # popt, pconv = curve_fit(magfunc, C, m0_mADU)
    popt2, pconv2 = curve_fit(magfunc, C, m0_mADU,absolute_sigma=True,sigma=y_sig)
    z1 = popt2[0]
    m  = popt2[1]
    # Calculate errors
    N           = len(m0_mADU)
    #z1_var      = pconv[0][0]
    #m_var       = pconv[1][1]
    z1_var2     = pconv2[0][0]
    m_var2       = pconv2[1][1]
    y_fit       = magfunc(C,z1,m)
    res1        = y_fit - m0_mADU
    res1_var     = numpy.sum( res1**2 ) / (N-1)
    # Goodness of Fit measure R2
    y_mean = numpy.mean(m0_mADU)
    ESS = numpy.sum( (y_fit-y_mean)**2 )
    TSS = numpy.sum( (m0_mADU-y_mean)**2 )
    Rsq = ESS/TSS
    return z1, z1_var2, m, m_var2, res1, res1_var, Rsq

def magfunc(x, z0, m):
    return z0 + m*x

def surffit(x,y,var,deg=1):
    # Fit a tilted plane to var maps
    p_init = models.Polynomial2D(degree=deg)
    fit_p  = fitting.LevMarLSQFitter()
    
    with warnings.catch_warnings():
        # Ignore model linearity warning from the fitter
        warnings.simplefilter('ignore')
        p = fit_p(p_init, x, y, var)

    return p

####################### BODY OF PROGRAM STARTS HERE ########################

if __name__ == "__main__":

    arguments = docopt.docopt(__doc__)

    # Mandatory argument
    catalog = arguments['<catalog>']

    # Non-mandatory options without arguments
    verbose = arguments['--verbose']
    cache = arguments['--cache']
    update = arguments['--update']
    show_plot = arguments['--plot']
    correctcolour = arguments['--kfit']

    # Non-mandatory options with arguments
    band = arguments['--filter']
    saveim_dir = arguments['--outputdir']

    maximum_error = arguments['--error']
    if maximum_error == None:
        maximum_error = 1.5;
    maximum_error = float(maximum_error)/3600.0   # Convert to degrees

    xmin = arguments['--xmin']
    if xmin == None:
        xmin = 0
    else:
        xmin = float(xmin)

    xmax = arguments['--xmax']
    if xmax == None:
        xmax = 1000000
    else:
        xmax = float(xmax)

    ymin = arguments['--ymin']
    if ymin == None:
        ymin = 0
    else:
        ymin = float(ymin)

    ymax = arguments['--ymax']
    if ymax == None:
        ymax = 1000000
    else:
        ymax = float(ymax)

    magmin = arguments['--magmin']
    if magmin == None:
        magmin = 0 
    else:
        magmin = float(magmin)

    detect_sigma = arguments['--threshold']
    if detect_sigma == None:
        detect_sigma = 10
    detect_sigma    = float(detect_sigma)

    sex_loc         = arguments['--sex']
    clip_sigma      = arguments['--dropsig']

    apass_dir = arguments['--refcat']

    if verbose:
        print arguments
        print_verbose_string( "Maximum permissible positon error is: %f deg" % maximum_error )

# Create SExtractor catalog.

    if re.search('.fits$', catalog):
        # Need to run SExtractor to generate the catalog
        if band == None:
            hdulist = fits.open(catalog)
            hdu = hdulist[0]
            hdulist.close()

            if 'FILTNAM' not in hdu.header.keys():
                if hdu.header['SERIALNO'] in ['83F010687','83F010783','83F010820','83F010827']:
                    band = 'SloanG'
                elif hdu.header['SERIALNO'] in ['83F010692','83F010730','83F010784','83F010826']:
                    band = 'SloanR'
            else:
                band = hdu.header['FILTNAM']

            if re.search('SloanG',band):
                band = 'g-band'
            elif re.search('SloanR',band):
                band = 'r-band'
            else:
                raise FilterError("Filter unknown. Specify a filter using the -f keyword.")
            if verbose:
                print_verbose_string( "Filter determined from header: %s" % band )
        elif re.search('SloanG',band):
            band = 'g-band'
        elif re.search('^g',band):
            band = 'g-band'
        elif re.search('SloanR',band):
            band = 'r-band'
        elif re.search('^r',band):
            band = 'r-band'
        else:
            raise FilterError("Filter unknown. Specify a filter using the -f keyword. Known filters: g, r.")

        cat = create_catalog(catalog,detect_thresh=detect_sigma)

    #Load catalog data 
    #(number, flux_auto, fluxerr_auto, fluxerr_aper, x_image, y_image, flux_radius, flags, class_star)
    sexdata = ascii.read(cat) 
    number = sexdata['NUMBER']
    flux_auto = sexdata['FLUX_AUTO']
    fluxerr_auto = sexdata['FLUXERR_AUTO']
    flux_aper = sexdata['FLUX_APER']
    fluxerr_aper = sexdata['FLUXERR_APER']
    x_image = sexdata['X_IMAGE']
    y_image = sexdata['Y_IMAGE']
    x_world = sexdata['X_WORLD']
    y_world = sexdata['Y_WORLD']
    flux_radius = sexdata['FLUX_RADIUS']
    flags = sexdata['FLAGS']
    class_star = sexdata['CLASS_STAR']
    fwhm_image = sexdata['FWHM_IMAGE']
    ellip      = sexdata['ELLIPTICITY']
    theta_image= sexdata['THETA_IMAGE']

# Load AAVSO APASS catalog (name,ra,dec,g,r)
    
    # Determine the minimum and maximum image positions
    min_x_image = min(x_image)
    max_x_image = max(x_image)
    min_y_image = min(y_image)
    max_y_image = max(y_image)
    x_centre = min_x_image + (max_x_image - min_x_image)/2.0
    y_centre = min_y_image + (max_y_image - min_y_image)/2.0

    # Determine the declination range spanned by the data
    # Plus a little: the same cached file is used for dithered images
    min_ra = min(x_world)-1.1
    max_ra = max(x_world)+1.1
    min_dec = min(y_world)-1.1
    max_dec = max(y_world)+1.1
    x_center = min_x_image + (max_x_image - min_x_image)/2.0
    y_center = min_y_image + (max_y_image - min_y_image)/2.0
    if verbose:
        print_verbose_string( "Min right ascension: %f" % min_ra ) 
        print_verbose_string( "Max right ascension: %f" % max_ra ) 
        print_verbose_string( "Min declination: %f" % min_dec )
        print_verbose_string( "Max declination: %f" % max_dec )

    # Load the AAVSO APASS catalog, either from the cache (fast) or from the original files (slow). 
    if cache and os.path.isfile('apass_cache.dat'):
        (name,ra,dec,g,r,gerr,rerr) = numpy.loadtxt('apass_cache.dat')
        if verbose:
            print_verbose_string( "Loading cached catalog data" )
    else:
        # Load the APASS file which corresponds to the minimum declination in the field of view
        apass_file = SelectAPASSFile( min_dec, apass_dir)
        (name,ra,dec,g,r,gerr,rerr) = LoadAPASSFile( apass_file, min_ra, max_ra, verbose=verbose )

        # If the file for the maximum declination is different from the file for the minimum declination then
        # load that one too and augment the arrays accordingly
        apass_file_extra = SelectAPASSFile( max_dec , apass_dir)

        if not (apass_file == apass_file_extra):
            (name_extra,ra_extra,dec_extra,g_extra,r_extra,gerr_extra,rerr_extra) = LoadAPASSFile(apass_file_extra, min_ra, max_ra, verbose=verbose )
            name = numpy.concatenate( (name, name_extra) )
            ra = numpy.concatenate( (ra, ra_extra) )
            dec = numpy.concatenate( (dec, dec_extra) )
            g = numpy.concatenate( (g, g_extra) )
            r = numpy.concatenate( (r, r_extra) )
            gerr = numpy.concatenate( (gerr, gerr_extra) )
            rerr = numpy.concatenate( (rerr, rerr_extra) )
        if cache:
            numpy.savetxt('apass_cache.dat',(name,ra,dec,g,r,gerr,rerr))
    
    # Select the magnitude corresponding to the known filter. If a sum of g and r
    # is being used then guesstimate the summed magnitude.
    try:
        if band == "g-band":
            catalog_mag = g
            catalog_mag_err = gerr
        elif band == "r-band":
            catalog_mag = r
            catalog_mag_err = rerr
        elif band == "sum":
            catalog_mag = -2.5*numpy.log10((10.0**(-0.4*g) + 10.0**(-0.4*r))/2.0)
            catalog_mag_err = -2.5*numpy.log10((10.0**(-0.4*gerr) + 10.0**(-0.4*rerr))/2.0)
        else:
            raise ZeroPointError 

    except ValueError:
        if verbose:
            print_verbose_string( "Arithmetic error encountered." )

# Match Dragonfly SExtractor catalog sources with those in APASS.

    # Build the KDTree
    DragonflyRADec = [x_world,y_world]
    DragonflyRADec = numpy.array(numpy.transpose(DragonflyRADec))
    AAVSORADec = [ra,dec]
    AAVSORADec = numpy.array(numpy.transpose(AAVSORADec))

    AAVSOTree = spatial.KDTree(AAVSORADec)
    nn = AAVSOTree.query(DragonflyRADec)
    dist = numpy.array(nn[0])
    index = numpy.array(nn[1])

    # Extract catalog entries in APASS catalog that correspond to Dragonfly catalog sources
    ra_matched   = numpy.array(ra)[index]
    dec_matched  = numpy.array(dec)[index]
    catalog_mag_matched = numpy.array(catalog_mag)[index] 
    catalog_mag_err_matched = numpy.array(catalog_mag_err)[index]  
    catalog_gr_matched  = numpy.array(g - r)[index]
    catalog_g_matched       = numpy.array(g)[index]
    catalog_r_matched       = numpy.array(r)[index]

# Remove certain sources from analysis to make trimmed catalog list

    # Condition (use bitwise "&" and not boolean "and")
    # Retain only a subset of closely matching pairs
    # If the catalog magnitude is unknown it gets a value of 99, weed those out
    # Don't keep Dragonfly stars that have negative flux
    condition = numpy.array(    (dist < maximum_error) &
                            (x_image > xmin) & (x_image < xmax) & (y_image > ymin) & (y_image < ymax) &
                            (catalog_mag_matched > magmin) & (catalog_mag_matched < 50) & (flux_auto > 0) &
                            (catalog_gr_matched > -5) & (catalog_gr_matched < 5) &
                            (fwhm_image < 5)
                        )
    # Extract Dragonfly and APASS catalog entries that fit condition
    # Dragonfly: 
    # (number, flux_auto, fluxerr_auto, fluxerr_aper, x_image, y_image, flux_radius, flags, class_star)
    # APASS: (name,ra,dec; catalog_mag)
    ii = numpy.array(numpy.where(condition)[0]) 
    # _c means matched & trimmed: APASS and Dragonfly matched AND retains only subset of good stars
    flux_auto_c    = numpy.array(flux_auto[ii])
    fluxerr_auto_c = numpy.array(fluxerr_auto[ii])
    flux_aper_c     = numpy.array(flux_aper[ii])
    x_image_c       = numpy.array(x_image[ii])
    y_image_c       = numpy.array(y_image[ii])
    fwhm_image_c    = numpy.array(fwhm_image[ii])
    ellip_c         = numpy.array(ellip[ii])
    theta_image_c   = numpy.array(theta_image[ii])

    # Turn Dragonfly catalog flux into magnitude
    #mag_aper_c  = -2.5*numpy.log10(flux_aper_c)
    mag_auto_c  = -2.5*numpy.log10(flux_auto_c)   

    # Calculate radius of source from centre of image
    rad_c       = numpy.sqrt( (x_image_c - x_centre)**2 + (y_image_c - y_centre)**2 )

    # Grab the right APASS catalog information
    catalog_mag_c   = numpy.array(catalog_mag_matched[ii])
    catalog_gr_c    = numpy.array(catalog_gr_matched[ii])
    catalog_mag_err_c = numpy.array(catalog_mag_err_matched[ii])
    catalog_g_c    = numpy.array(catalog_g_matched[ii])
    catalog_r_c    = numpy.array(catalog_r_matched[ii])

# Sigma clip using simple zeropoints before sending into fitting for various affects

    #zp_c        = catalog_mag_c - mag_aper_c
    zp_c        = catalog_mag_c - mag_auto_c
    clip_sigma  = float(clip_sigma)
    clipped_zp  = sigma_clip(zp_c, sig=clip_sigma, iters = 5)
    outliers    = clipped_zp.mask
    clipped_zp1  = clipped_zp[~outliers].data

    clipped_flux_auto      = numpy.ma.masked_array(flux_auto_c,mask=outliers)[~outliers].data
    clipped_fluxerr_auto   = numpy.ma.masked_array(fluxerr_auto_c,mask=outliers)[~outliers].data
    clipped_flux_aper       = numpy.ma.masked_array(flux_aper_c,mask=outliers)[~outliers].data
    clipped_x_image         = numpy.ma.masked_array(x_image_c,mask=outliers)[~outliers].data
    clipped_y_image         = numpy.ma.masked_array(y_image_c,mask=outliers)[~outliers].data
    clipped_fwhm_image      = numpy.ma.masked_array(fwhm_image_c,mask=outliers)[~outliers].data
    clipped_ellip           = numpy.ma.masked_array(ellip_c,mask=outliers)[~outliers].data
    clipped_theta_image     = numpy.ma.masked_array(theta_image_c,mask=outliers)[~outliers].data
    clipped_mag_auto      = numpy.ma.masked_array(mag_auto_c,mask=outliers)[~outliers].data

    clipped_rad           = numpy.ma.masked_array(rad_c,mask=outliers)[~outliers].data

    clipped_catalog_mag     = numpy.ma.masked_array(catalog_mag_c,mask=outliers)[~outliers].data
    clipped_catalog_gr      = numpy.ma.masked_array(catalog_gr_c,mask=outliers)[~outliers].data
    clipped_catalog_mag_err = numpy.ma.masked_array(catalog_mag_err_c,mask=outliers)[~outliers].data
    clipped_catalog_g      = numpy.ma.masked_array(catalog_g_c,mask=outliers)[~outliers].data
    clipped_catalog_r      = numpy.ma.masked_array(catalog_r_c,mask=outliers)[~outliers].data


# Zeropoints and residues before doing any fitting.

    m0_mADU = clipped_catalog_mag-clipped_mag_auto #same as clipped_zp1    
    z01          = clipped_zp1.mean()
    z0_err1      = clipped_zp1.std()
    z0_var1      = z0_err1**2
    clipped_res1    = clipped_zp1 - z01

    # Measure dispersion in values
    z0MAD1       = numpy.median(   numpy.absolute(clipped_zp1 - numpy.median(clipped_zp1))    )
    q75, q25 = numpy.percentile(clipped_zp1, [75 ,25])
    z0IQR1        = q75-q25

# Fit colour term and remove 

    if correctcolour:
        clipped_catalog_flux_err = 10**(clipped_catalog_mag_err/-2.5)
        clipped_fluxerr_auto = clipped_fluxerr_auto
        clipped_fluxerr_total = numpy.sqrt(clipped_catalog_flux_err**2 + clipped_fluxerr_auto**2)
        y_sig = numpy.abs(-2.5*numpy.log10(clipped_fluxerr_total))
        (z02,z0_var2,k2,k_var2,clipped_res2,res2_var,Rsq2)=curve_fit_wrapper(m0_mADU,clipped_catalog_gr,y_sig)
        m0_mADU_kC = m0_mADU - k2*clipped_catalog_gr
        clipped_zp2 = z02 + clipped_res2

        # Measure dispersion in values
        z0MAD2       = numpy.median(   numpy.absolute(clipped_zp2 - numpy.median(clipped_zp2))    )
        q75, q25 = numpy.percentile(clipped_zp2, [75 ,25])
        z0IQR2        = q75-q25

    else: #skip colour term calculations
        clipped_res2 = clipped_res1
        clipped_zp2  = clipped_zp1


# Fit residues with surface and divide it out

    p_res = surffit(clipped_x_image,clipped_y_image,clipped_res2, deg=1)
    (naxis1,naxis2) = getimsize(catalog)
    x = numpy.array(map(float,range(naxis1)))
    y = numpy.array(map(float,range(naxis2)))
    xv, yv = numpy.meshgrid(x, y)
    model_res_mag = p_res(xv,yv)
    model_res_factor = 10**(model_res_mag/(-2.5))

    # save new image
    new_image=save_image(catalog,model_res_factor,saveim_dir,suffix='pcapass',savemap=True)

# Calculate new zeropoint and scatter
    
    clipped_zp3 = clipped_zp2 - p_res(clipped_x_image,clipped_y_image)
    z03         = clipped_zp3.mean()
    z0_err3      = clipped_zp3.std()
    z0_var3      = z0_err3**2
    clipped_res3    = clipped_zp3 - z03

    # Measure dispersion in values
    z0MAD3       = numpy.median(   numpy.absolute(clipped_zp3 - numpy.median(clipped_zp3))    )
    q75, q25 = numpy.percentile(clipped_zp3, [75 ,25])
    z0IQR3        = q75-q25

# Fit radial surface to photometry residues and divide it out

# If "--update" (and a .fits was provided), add the information as keywords to the image

    if update and re.search('.fits$', new_image):
        print new_image + ' header will be updated with addition parameters'
        # If update and fits file provided, do the updating!
        data, header = fits.getdata(new_image, header=True)
        header['ZP1'] = (z01, 'Dragonfly Pipeline')
        header['ZPRMS1'] = (z0_var1, 'Dragonfly Pipeline')
        header['ZPMAD1'] = (z0MAD1, 'Dragonfly Pipeline')
        header['ZPIQR1'] = (z0IQR1, 'Dragonfly Pipeline')
        if correctcolour:
            header['ZP2'] = (z02, 'Dragonfly Pipeline')
            header['Z0FITRMS']=(z0_var2, 'z0_var2, Dragonfly Pipeline')
            header['ZPRMS2'] = (res2_var, 'res2_var Dragonfly Pipeline')
            header['ZPMAD2'] = (z0MAD2, 'Dragonfly Pipeline')
            header['ZPIQR2'] = (z0IQR2, 'Dragonfly Pipeline')
            header['KTERM']=(k2, 'colour term k*(g-r), Dragonfly Pipeline')
            header['KTERMVAR']=(k_var2, 'Dragonfly Pipeline')
            header['RSQ2']=(Rsq2, 'Dragonfly Pipeline')
        header['ZP3'] = (z03, 'Dragonfly Pipeline')
        header['ZPRMS3'] = (z0_var3, 'sqrt(z0 var) Dragonfly Pipeline')
        header['ZPMAD3'] = (z0MAD3, 'Dragonfly Pipeline')
        header['ZPIQR3'] = (z0IQR3, 'Dragonfly Pipeline')
        header['SIGCLIP'] = (clip_sigma, 'Sources with simple zeropoints within this number of sigmas is included in the final zeropoint calculation, Dragonfly Pipeline')
        header['ZPNGOOD'] = len(clipped_zp)
        header['ZPNREJ'] = ( len(zp_c) - len(clipped_zp))
        header['FWHM'] = sexdata['FWHM_IMAGE'].mean()
        header['FWHMRMS'] = sexdata['FWHM_IMAGE'].std()
        header['AXRATIO'] = (sexdata['B_IMAGE']/sexdata['A_IMAGE']).mean()
        header['AXRRMS'] = (sexdata['B_IMAGE']/sexdata['A_IMAGE']).std()
        header['THMEAN'] = (sexdata['THETA_IMAGE']).mean()
        header['THRMS'] = (sexdata['THETA_IMAGE']).std()
        header['SKYLVL'] = (sexdata['BACKGROUND']).mean()
        header['SKYRMS'] = (sexdata['BACKGROUND']).std()    
        header['SWARPSCL'] = (10**(z03/2.5)) / (10**(27/2.5))
        if ((sexdata['BACKGROUND']).mean() > 0.0):
            header['SKYSB'] = z03 - 2.5*math.log10((sexdata['BACKGROUND']).mean()) + 5.0*math.log10(2.85)
        fits.writeto(new_image,data,header,clobber=True) 

# Optionally plot the results

    if show_plot and correctcolour:
        
        print ' '
        print 'plots will be made showing zeropoint distribution'

    # Figure out histogram limits
        final_minz = 1000;
        final_maxz = -1000;
        for z in [clipped_zp1, clipped_zp2, clipped_zp3]:
            minz = min(z)
            maxz = max(z)
            if minz < final_minz:
                final_minz = minz
            if maxz > final_maxz:
                final_maxz = maxz
        final_minz = final_minz-0.05
        final_maxz = final_maxz+0.05

    # plot zeropoints (clipped) 

        zp2plot = clipped_zp3
        cm = plt.cm.get_cmap('RdYlBu')
        plt.figure(figsize=(22,24))

        # m0-mADU or z0+kC+f(x,y)+g(f,y)+residue distribution across FOV
        plt.subplot(4,3,1)
        yy = clipped_res1
        zmin = yy.mean() - 3*yy.std()
        zmax = yy.mean() + 3*yy.std()
        sc = plt.scatter(clipped_x_image, clipped_y_image, c=yy, vmin=zmin, vmax=zmax, s=30, cmap=cm)
        plt.axis( [min(x_image),max(x_image),min(y_image),max(y_image)] )
        plt.colorbar(sc,label='kC+f(x,y)+g(x,y)+residue')
        plt.xlabel('X_IMAGE')
        plt.ylabel('Y_IMAGE')

        # histogram of z01 
        plt.subplot(4,3,2)
        n, bins, patches = P.hist(clipped_zp1, 20)
        P.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
        P.xlim([final_minz, final_maxz])
        r = max(clipped_zp1) - min(clipped_zp1)
        r = ("%.3f" % r)
        rms = clipped_zp1.std()
        rms = ("%.5f" % rms)
        plt.xlabel('z0+kC+f(x,y)+g(x,y)+residue; R: '+r+' std: '+rms)
        plt.ylabel('N')

        # Add title to this subplot as title of whole plot
        filename = str(catalog)
        filename = filename.split('/')[-1]
        plt.title(  filename+ '; clipped zeropoints at sigma = ' +str(clip_sigma) + ' model: m0=m_ADU+z0+kC+f(x,y)+g(f,y)+res'  )

        # magnitude vs zeropoint1 scatter plot
        plt.subplot(4,3,3)
        sc = plt.scatter(clipped_catalog_mag, clipped_zp1, c=clipped_rad, vmin=numpy.min(clipped_rad), vmax=numpy.max(clipped_rad), s=20, cmap=cm,alpha=0.4)
        plt.colorbar(sc,label='Radius')
        ymax = max(clipped_zp1) +0.1
        ymin = min(clipped_zp1) -0.1
        plt.axis( [min(clipped_catalog_mag),max(clipped_catalog_mag),ymin,ymax] )
        plt.xlabel('mag')
        plt.ylabel('z0+kC+f(x,y)+g(x,y)+residue')


            


        # m0-mADU-kC or z0+f(x,y)+g(f,y)+residue distribution across FOV
        plt.subplot(4,3,4)
        yy = clipped_res2
        zmin = yy.mean() - 3*yy.std()
        zmax = yy.mean() + 3*yy.std()
        sc = plt.scatter(clipped_x_image, clipped_y_image, c=yy, vmin=zmin, vmax=zmax, s=30, cmap=cm)
        plt.axis( [min(x_image),max(x_image),min(y_image),max(y_image)] )
        plt.colorbar(sc,label='f(x,y)+g(x,y)+residue')
        plt.xlabel('X_IMAGE')
        plt.ylabel('Y_IMAGE')

        # histogram of z02
        plt.subplot(4,3,5)
        n, bins, patches = P.hist(clipped_zp2, 20)
        P.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
        P.xlim([final_minz, final_maxz])
        r = max(clipped_zp2) - min(clipped_zp2)
        r = ("%.3f" % r)
        rms = clipped_zp2.std()
        rms = ("%.5f" % rms)
        plt.xlabel('z0+f(x,y)+g(x,y)+residue; R: '+r+' std: '+rms)
        plt.ylabel('N')

        # g-r colour vs zeropoint2 scatter plot AND fitted colour term line
        plt.subplot(4,3,6)
        xx = [min(clipped_catalog_gr),max(clipped_catalog_gr)]
        yy = z02 + numpy.array(xx)*k2
        plt.plot(xx,yy)
        sc = plt.scatter(clipped_catalog_gr, clipped_zp1, c=y_sig, vmin=numpy.min(y_sig), vmax=numpy.max(y_sig), s=20, cmap=cm,alpha=0.4)
        plt.colorbar(sc,label='error in m0-mAUD')
        ymax = max(clipped_zp1) +0.1
        ymin = min(clipped_zp1) -0.1
        plt.axis( [min(clipped_catalog_gr),max(clipped_catalog_gr),ymin,ymax] )
        plt.xlabel('g-r catalog mag')
        plt.ylabel('z0+kC+f(x,y)+g(x,y)+residue')





        # m0-mADU-kC-f(x,y) or z0+g(f,y)+residue distribution across FOV
        plt.subplot(4,3,7)
        yy = clipped_res3
        zmin = yy.mean() - 3*yy.std()
        zmax = yy.mean() + 3*yy.std()
        sc = plt.scatter(clipped_x_image, clipped_y_image, c=yy, vmin=zmin, vmax=zmax, s=30, cmap=cm)
        plt.axis( [min(x_image),max(x_image),min(y_image),max(y_image)] )
        plt.colorbar(sc,label='g(x,y)+residue')
        plt.xlabel('X_IMAGE')
        plt.ylabel('Y_IMAGE')

        # histogram of z03
        plt.subplot(4,3,8)
        n, bins, patches = P.hist(clipped_zp3, 20)
        P.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
        P.xlim([final_minz, final_maxz])
        r = max(clipped_zp3) - min(clipped_zp3)
        r = ("%.3f" % r)
        rms = clipped_zp3.std()
        rms = ("%.5f" % rms)
        plt.xlabel('z0+g(x,y)+residue; R: '+r+' std: '+rms)
        plt.ylabel('N')

        # Fit to residue map
        plt.subplot(4,3,9)
        zmin = clipped_res3.mean() - 3*clipped_res3.std()
        zmax = clipped_res3.mean() + 3*clipped_res3.std()
        sc = plt.scatter(clipped_x_image, clipped_y_image, c=p_res(clipped_x_image,clipped_y_image), vmin=zmin, vmax=zmax, s=30, cmap=cm)
        plt.axis( [min(x_image),max(x_image),min(y_image),max(y_image)] )
        plt.colorbar(sc,label='Fit to m0-mADU-kC')
        plt.xlabel('X_IMAGE')
        plt.ylabel('Y_IMAGE')





        # if subplot(5,3) #10-12 inclusive is res4 map, z4 distribution, radius vs z3


  
        z0=z03
        res = clipped_res3
        # radius vs z0+residue 
        plt.subplot(4,3,10)
        xx = [min(clipped_rad),max(clipped_rad)]
        yy = [z0,z0] 
        plt.plot(xx,yy)
        sc=plt.scatter(clipped_rad,z0+res,c=clipped_catalog_gr, vmin=numpy.min(clipped_catalog_gr), vmax=numpy.max(clipped_catalog_gr), s=20, cmap=cm,alpha=0.3)
        ymax = max(z0+res) +0.06
        ymin = min(z0+res) -0.06
        plt.axis( [min(clipped_rad),max(clipped_rad),ymin,ymax] )
        plt.colorbar(sc,label='Colour')
        plt.xlabel('Radius')
        plt.ylabel('z0+residue')

        # fwhm vs z0+residue
        plt.subplot(4,3,11)
        xx = [min(clipped_fwhm_image),max(clipped_fwhm_image)]
        yy = [z0,z0] 
        plt.plot(xx,yy)
        sc=plt.scatter(clipped_fwhm_image,z0+res,c=clipped_rad, vmin=numpy.min(clipped_rad), vmax=numpy.max(clipped_rad), s=20, cmap=cm,alpha=0.3)
        ymax = max(z0+res) +0.06
        ymin = min(z0+res) -0.06
        plt.axis( [min(clipped_fwhm_image),max(clipped_fwhm_image),ymin,ymax] )
        plt.colorbar(sc,label='Radius')
        plt.xlabel('FWHM')
        plt.ylabel('z0+residue')

        # ellipticity vs z0+residue
        plt.subplot(4,3,12)
        xx = [min(clipped_ellip),max(clipped_ellip)]
        yy = [z0,z0] 
        plt.plot(xx,yy)
        sc=plt.scatter(clipped_ellip,z0+res,c=clipped_rad, vmin=numpy.min(clipped_rad), vmax=numpy.max(clipped_rad), s=20, cmap=cm,alpha=0.3)
        ymax = max(z0+res) +0.06
        ymin = min(z0+res) -0.06
        plt.axis( [min(clipped_ellip),max(clipped_ellip),ymin,ymax] )
        plt.colorbar(sc,label='Radius')
        plt.xlabel('Ellipticity')
        plt.ylabel('z0+residue')

        ## theta vs z0+residue
        #plt.subplot(3,4,12)
        #xx = [min(clipped_theta_image),max(clipped_theta_image)]
        #yy = [z02,z02] 
        #plt.plot(xx,yy)
        #sc=plt.scatter(clipped_theta_image,z02+clipped_z0_residues2,c=clipped_rad, vmin=numpy.min(clipped_rad), vmax=numpy.max(clipped_rad), s=20, cmap=cm,alpha=0.3)
        #ymax = max(z02+clipped_z0_residues2) +0.06
        #ymin = min(z02+clipped_z0_residues2) -0.06
        #plt.axis( [min(clipped_theta_image),max(clipped_theta_image),ymin,ymax] )
        #plt.colorbar(sc,label='Radius')
        #plt.xlabel('Theta')
        #plt.ylabel('z0+residue')

        # save plots
        saveloc = str(catalog)
        saveloc = saveloc.split('.')[0]+'_clipped_apass.png'
        print 'saving second plot in: '+saveloc
        print ' '
        plt.savefig(saveloc)


        # plot colour versus zeropoint with colourbar showing catalog brightness of star in g and r bands
        plt.figure(figsize=(20,18))
        cm = plt.cm.get_cmap('RdYlBu')

        plt.subplot(2,1,1)
        xx = [min(clipped_catalog_gr),max(clipped_catalog_gr)]
        yy = z02 + numpy.array(xx)*k2
        plt.plot(xx,yy)
        sc = plt.scatter(clipped_catalog_gr, clipped_zp1, c=clipped_catalog_g, vmin=numpy.min(clipped_catalog_g), vmax=numpy.max(clipped_catalog_g), s=30, cmap=cm,alpha=0.5)
        plt.colorbar(sc,label='APASS g magnitude')
        plt.xlabel('g-r catalog mag')
        plt.ylabel('z0+kC+f(x,y)+g(x,y)+residue')

        plt.subplot(2,1,2)
        xx = [min(clipped_catalog_gr),max(clipped_catalog_gr)]
        yy = z02 + numpy.array(xx)*k2
        plt.plot(xx,yy)
        sc = plt.scatter(clipped_catalog_gr, clipped_zp1, c=clipped_catalog_r, vmin=numpy.min(clipped_catalog_r), vmax=numpy.max(clipped_catalog_r), s=30, cmap=cm,alpha=0.5)
        plt.colorbar(sc,label='APASS r magnitude')
        plt.xlabel('g-r catalog mag')
        plt.ylabel('z0+kC+f(x,y)+g(x,y)+residue')

        # save plots
        saveloc = str(catalog)
        saveloc = saveloc.split('.')[0]+'_colourapass.png'
        print 'saving second plot in: '+saveloc
        print ' '
        plt.savefig(saveloc)

