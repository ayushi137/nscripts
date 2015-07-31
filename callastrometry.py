"""

callastrometry - contains functions to handle WCS header information

Requires the following modules: os, docopt, astropy

Contains the following functions: callastrometry, scrubwcsheader

Usage:
callastrometry [-h]

Options:
    -h, --help          Show this screen

"""

########################## IMPORT PACKAGES ###########################
import os
import docopt
from astropy import wcs
from astropy.io import fits

docopt.docopt(__doc__)

########################## DATA LISTS ###########################

# List of file extensions produced by astrometry
fileextensions = ['.axy','.corr','-indx.xyls','.match','.new','-objs.png',
                  '.rdls','.solved']
# List of dictionary keys corresponding to WCS solutions
keys = ['CTYPE1','CTYPE2','CRPIX1','CRVAL1','CRPIX2','CRVAL2','CD1_1',
        'CD1_2','CD2_1','CD2_2','EQUINOX']

########################## FUNCTIONS ###########################
def callastrometry(fname,generate=False,filekeep=False):
    """
    Find WCS solution for an image: the transformation between pixels and R.A./dec

    fname:      name of image file on which to perform astrometry
    generate:   Boolean that specifies whether to force astrometry solution or not
                (kwarg, default = False)
    filekeep:   Boolean that specifies whether to keep auxiliary files from
                astrometry solution (kwarg, default = False)

    Returns nothing explicitly, implicitly saves update file
    """
    # Read in fits header
    header = fits.getheader(fname)
    # If all WCS entries already in header and generation not forced, done
    if set(keys) < set(header.keys()) and generate == False:
        print 'Astrometry already done'
    # Otherwise, do astrometry
    else:
        # Run astrometry.net
        command = 'solve-field --no-fits2fits --use-sextractor --cpulimit 20 {0}'.format(fname)
        os.system(command)
        trimfname = fname.split('.fits')[0]
        if filekeep == False:
            for ext in fileextensions:
                os.system('rm -f '+trimfname+ext)
        # Read in WCS information
        lightwcs = fits.getheader(trimfname+'.wcs')
        w = wcs.WCS(lightwcs)
        # Get x,y pixel scales
        scales = wcs.utils.proj_plane_pixel_scales(w)*3600
        # Read in file data
        data,light = fits.getdata(fname,header=True)
        # Update header
        light['PSCALY']=scales[0]
        light['PSCALX']=scales[1]
        for key in keys:
            light[key] = lightwcs[key]
        # Save updated file
        fits.writeto(fname,data,light,clobber=True)
        if filekeep == False:
            os.system('rm -f '+trimfname+'.wcs')

def scrubwcsheader(fname,ind=0):
    """
    Removes all WCS information from header

    fname:      name of file from which to remove WCS solution
    ind:        index of image in fits file

    Returns nothing explicitly, implicitly saves update file
    """
    # Open file and read header
    light = fits.open(fname)
    header = light[ind].header
    # Remove WCS keys
    for key in keys:
        try:
            del header[key]
        except (ValueError,KeyError):
            continue
    data = light[ind].data
    light.close()
    # Save updated file
    fits.writeto(fname,data,header,clobber = True)
