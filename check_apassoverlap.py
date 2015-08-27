#!/usr/bin/env python

"""check_apassoverlap -- Plot the overlapping stars in the APASS catalog and an input image with wcs; specific catalog columns are needed.  :

Usage:  
    check_apassoverlap [-h] [-v] [-c] [-t THREASHOLD] [-e ARCSEC] [-s SEXLOCATION] [-r DIRECTORY] (save | show) <image>

Options:
    -h, --help                              Show this screen
    -v, --verbose                           Show extra information [default: False]
    -c, --cache                             Create/use APASS catalog cache [default: False]
    -t NUMBER, --threshold NUMBER           Number of sigma to threshold [default: 15] 
    -e ARCSEC, --error=ARCSEC               Maximum allowable separation in arcsec [default: 1.5]
    -s SEXLOCATION, --sex SEXLOCATION       Location of SExtractor executable [default: /opt/local/bin/sex]
    -r DIRECTORY, --refcat DIRECTORY        The directory containing the reference APASS catalogs [default: /Volumes/dragonfly/data/dragonflysurvey/APASS]

Examples:
    python check_apassoverlap.py -v save input.fits
"""

import docopt
import sys
import re
from astropy.io import ascii
import matplotlib.pyplot as plt
import subprocess
import os
import pandas as pd
import numpy
from scipy import spatial

def print_verbose_string(printme):
    # Print information to standard error
    print >> sys.stderr, "VERBOSE: %s" % printme

def SelectAPASSFile(declination,apass_dir):

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
        print 'This apass catalog file has been selected: '+apass_file
        return(apass_file)
    else:
        print 'This apass catalog file does not exist: '+apass_file
        return(None)

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
    r = aacat.r[indices].values
    gerr = aacat.gerr[indices].values
    rerr = aacat.rerr[indices].values

    return(name, ra, dec, g, r,gerr,rerr)

# SExtractor configuration files

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
    MEMORY_BUFSIZE 1024
    MEMORY_OBJSTACK 3000
    MEMORY_PIXSTACK 300000
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
    if verbose:
        print 'Run the following SExtractor command: ' 
        print (sex_loc+" -c {config} -CATALOG_NAME {catalog} {image}".
                format(config=sextractor_config_name, catalog=catalog_name, image=image_name)
                )
    # Run SExtractor command
    subprocess.call(sex_loc+" -c {config} -CATALOG_NAME {catalog} {image}".format(config=sextractor_config_name, catalog=catalog_name, image=image_name), shell=True)
    return catalog_name

def remove_tmpfiles():
    sextractor_config_name = "/tmp/scamp.sex"
    params_name = "/tmp/scamp.param"
    nnw_name = "/tmp/default.nnw"
    conv_name = "/tmp/default.conv"
    catalog_name = "/tmp/cz.cat"
    os.remove(sextractor_config_name)
    os.remove(params_name)
    os.remove(nnw_name)
    os.remove(conv_name)
    os.remove(catalog_name)






####################### BODY OF PROGRAM STARTS HERE ########################
if __name__ == "__main__":

# Collect Arguments

    arguments = docopt.docopt(__doc__)

    # Mandatory
    image  = arguments['<image>']
    saveit      = arguments['save']
    showit      = arguments['show']

    # Optional
    verbose     = arguments['--verbose']
    cache = arguments['--cache']
    detect_sigma    = arguments['--threshold']
    maximum_error   = arguments['--error']
    maximum_error = float(maximum_error)/3600.0   # Convert to degrees
    sex_loc         = arguments['--sex']
    apass_dir = arguments['--refcat']
    apass_dir = apass_dir + '/'

    if verbose:
        print arguments

# Load catalog file; generate first if input is fits file

    if verbose and re.search('.cat$', image):
        print 'The following catalog columns are needed for input image catalog: '
        print sextractor_params
        print ''

    if re.search('.cat$', image) == None:
        # Need to run SExtractor to generate the catalog
        cat = create_catalog(image,detect_thresh=detect_sigma)
    else:
        # The catalog is just the input image file
        cat = image

    # Load catalog file
    sexdata = ascii.read(cat) 
    number = sexdata['NUMBER']
    flux_auto = sexdata['FLUX_AUTO']
    fluxerr_auto = sexdata['FLUXERR_AUTO']
    flux_aper = sexdata['FLUX_APER']
    fluxerr_aper = sexdata['FLUXERR_APER'], 
    x_image = sexdata['X_IMAGE']
    y_image = sexdata['Y_IMAGE']
    x_world = sexdata['X_WORLD']
    y_world = sexdata['Y_WORLD']
    flux_radius = sexdata['FLUX_RADIUS']
    flags = sexdata['FLAGS']
    class_star = sexdata['CLASS_STAR']

# Load APASS catalog file(s)

    # Determine the declination range spanned by the data
    min_ra = min(x_world)-1.2
    max_ra = max(x_world)+1.2
    min_dec = min(y_world)-1.2
    max_dec = max(y_world)+1.2
    if verbose:
        print_verbose_string( "Min right ascension: %f" % min_ra ) 
        print_verbose_string( "Max right ascension: %f" % max_ra ) 
        print_verbose_string( "Min declination: %f" % min_dec )
        print_verbose_string( "Max declination: %f" % max_dec ) 

    if cache and os.path.isfile('apass_cache.dat'):
        (name,ra,dec,g,r,gerr,rerr) = numpy.loadtxt('apass_cache.dat')
        if verbose:
            print_verbose_string( "Loading cached catalog data" )
    else:
        # Load the APASS file which corresponds to the minimum declination in the field of view
        apass_file = SelectAPASSFile( min_dec, apass_dir)
        (name,ra,dec,g,r,gerr,rerr) = LoadAPASSFile( apass_file, min_ra-1, max_ra+1, verbose=verbose )
        # If the file for the maximum declination is different from the file for the minimum declination then
        # load that one too and augment the arrays accordingly
        apass_file_extra = SelectAPASSFile( max_dec , apass_dir)
        if not (apass_file == apass_file_extra):
            (name_extra,ra_extra,dec_extra,g_extra,r_extra,gerr_extra,rerr_extra) = LoadAPASSFile( apass_file_extra, min_ra-1, max_ra+1, verbose=verbose )
            name = numpy.concatenate( (name, name_extra) )
            ra = numpy.concatenate( (ra, ra_extra) )
            dec = numpy.concatenate( (dec, dec_extra) )
            g = numpy.concatenate( (g, g_extra) )
            r = numpy.concatenate( (r, r_extra) )
            gerr = numpy.concatenate( (gerr,gerr_extra) )
            rerr = numpy.concatenate( (rerr, rerr_extra) )
        if cache:
            numpy.savetxt('apass_cache.dat',(name,ra,dec,g,r,gerr,rerr))

    # Extract RA and DECS from image cat and APASS cat
    DragonflyRADec = [x_world,y_world]
    DragonflyRADec = numpy.array(numpy.transpose(DragonflyRADec))
    AAVSORADec = [ra,dec]
    AAVSORADec = numpy.array(numpy.transpose(AAVSORADec))  

# Match APASS sources to image catalog sources

    # Build KDTree in order to match image SExtractor catalog and APASS catalog targets by position  
    AAVSOTree = spatial.KDTree(AAVSORADec)
    nn = AAVSOTree.query(DragonflyRADec)    # Finding which AAVSO RA DECs are closest to each DragonflyRADec
    dist = numpy.array(nn[0])               # Dist between closest match
    index = numpy.array(nn[1])              # AAVSO cat index that matches each input image cat entry in order

    if verbose:
        print 'Total number of stars extracted from APASS catalog: ' + str(len(ra))
        print 'Total number of stars extracted from SExtracted image catalog: ' + str(len(dist))

    # Extract APASS sources that correspond to image catalog sources
    aara   =numpy.array(ra[index])
    aadec  =numpy.array(dec[index])

# Plot positions, selected sources and deviations
    plt.figure(figsize=(16,12))
    plt.suptitle('Plots of positions of image and APASS catalog sources; Blue box shows image FOV')
    cm = plt.cm.get_cmap('RdYlBu')

    # Plot surrounding AAVSO (APASS) sources' RA, DEC and input image (Dragonfly) sources' RA DEC
    cond = numpy.array((g<14))
    jj=numpy.array(numpy.where(cond)[0])
    ra_bright    = ra[jj]
    dec_bright   = dec[jj]
    if verbose:
        print 'Total number of stars extracted from APASS catalog brighter than g mag 14: ' + str(len(ra_bright))
    plt.subplot(2,2,1)
    plt.scatter(ra_bright,dec_bright,marker='o',color='m',label='apass mag < 14')
    plt.scatter(x_world,y_world,marker='x',label='dragonfly')
    plt.legend(bbox_to_anchor=(0., 1.01, 1., .101), loc=3, ncol=2, mode="expand", borderaxespad=0.)
    plt.xlabel('dec')
    plt.ylabel('ra')
    plt.plot([min(x_world),min(x_world)],[min(y_world),max(y_world)],color='b')
    plt.plot([max(x_world),max(x_world)],[min(y_world),max(y_world)],color='b')
    plt.plot([min(x_world),max(x_world)],[min(y_world),min(y_world)],color='b')
    plt.plot([min(x_world),max(x_world)],[max(y_world),max(y_world)],color='b')
    plt.axis( [min(x_world)-1,max(x_world)+1,min(y_world)-1,max(y_world)+1] )
    plt.axis('tight')

    # Plot all matched APASS sources and dragonfly sources
    plt.subplot(2,2,2)
    plt.scatter(aara,aadec,marker='o',color='m',label='apass mag < 14')
    plt.scatter(x_world,y_world,marker='x',label='dragonfly')
    plt.legend(bbox_to_anchor=(0., 1.01, 1., .101), loc=3, ncol=2, mode="expand", borderaxespad=0.)
    plt.xlabel('dec')
    plt.ylabel('ra')
    plt.plot([min(x_world),min(x_world)],[min(y_world),max(y_world)],color='b')
    plt.plot([max(x_world),max(x_world)],[min(y_world),max(y_world)],color='b')
    plt.plot([min(x_world),max(x_world)],[min(y_world),min(y_world)],color='b')
    plt.plot([min(x_world),max(x_world)],[max(y_world),max(y_world)],color='b')
    plt.axis( [min(x_world)-0.1,max(x_world)+0.1,min(y_world)-0.1,max(y_world)+0.1] )
    plt.axis('tight')

    # Plot again after removing outliers
    cond=numpy.array((dist<maximum_error) & (flux_auto>0))
    ii=numpy.array(numpy.where(cond)[0])
    aara_out    = aara[ii]
    aadec_out   = aadec[ii]
    x_world_out = x_world[ii]
    y_world_out = y_world[ii]
    if verbose:
        print 'Total number of stars extracted from APASS catalog that mathced well with image SExtracted catalog: ' + str(len(aara_out))
    if len(aara_out) < 20:
        print ''
        print '*********Problem: there are not many stars matching between image and APASS catalog!*********'
        print ''
    plt.subplot(2,2,3)
    plt.scatter(aara_out,aadec_out,marker='o',color='m',label='apass with dist<'+ str(maximum_error*3600)+'"')
    plt.scatter(x_world_out,y_world_out,marker='x',label='dragonfly')
    plt.legend(bbox_to_anchor=(0., 1.01, 1., .101), loc=3, ncol=2, mode="expand", borderaxespad=0.)
    plt.xlabel('dec')
    plt.ylabel('ra')
    plt.plot([min(x_world),min(x_world)],[min(y_world),max(y_world)],color='b')
    plt.plot([max(x_world),max(x_world)],[min(y_world),max(y_world)],color='b')
    plt.plot([min(x_world),max(x_world)],[min(y_world),min(y_world)],color='b')
    plt.plot([min(x_world),max(x_world)],[max(y_world),max(y_world)],color='b')
    plt.axis( [min(x_world)-0.1,max(x_world)+0.1,min(y_world)-0.1,max(y_world)+0.1] )
    plt.axis('tight')

    # Plot Dragonfly sources and colour by distance error
    plt.subplot(2,2,4)
    dist_arcsec = dist*3600.0
    #zmin = 0 - 0.8*dist_arcsec.std()
    #zmax = 0 + 0.8*dist_arcsec.std()
    zmin = 0
    zmax = 0 + 6.0
    sc=plt.scatter(x_world,y_world,marker='x', c=dist_arcsec, vmin=zmin, vmax=zmax, s=35, cmap=cm,label='dragonfly')
    plt.legend(bbox_to_anchor=(0., 1.01, 1., .101), loc=3, ncol=1, mode="expand", borderaxespad=0.)
    plt.xlabel('dec')
    plt.ylabel('ra')
    plt.plot([min(x_world),min(x_world)],[min(y_world),max(y_world)],color='b')
    plt.plot([max(x_world),max(x_world)],[min(y_world),max(y_world)],color='b')
    plt.plot([min(x_world),max(x_world)],[min(y_world),min(y_world)],color='b')
    plt.plot([min(x_world),max(x_world)],[max(y_world),max(y_world)],color='b')
    plt.axis( [min(x_world)-0.1,max(x_world)+0.1,min(y_world)-0.1,max(y_world)+0.1] )
    plt.colorbar(sc,label='Offset from APASS location (arcsec)')
    plt.axis('tight')

    if saveit:
        saveloc = image.split('.')[0]+'_apassoverlap.png'
        plt.savefig(saveloc)
    if showit:
        plt.show()

    

# Clean up by removing all tmp files
    remove_tmpfiles()

