"""
callAPASS - contains a function to create APASS catalogues

Requires the following modules: urllib, docopt, numpy, os

Contains the following functions: callAPASS

Usage:
callAPASS [-h]

Options:
    -h, --help          Show this screen

"""

########################## IMPORT PACKAGES ###########################

import urllib as url
import docopt
from numpy import *
import os

docopt.docopt(__doc__)

########################## FUNCTIONS ###########################

def callAPASS(ra,dec,fov,APASSdir):
    """
    Finds an appropriate APASS catalogue for an image and saves it in APASSdir

    ra:         right ascension in degrees
    dec:        declination in degrees
    fov:        field of view in degrees (max 15)
    APASSdir:
    
    Saves an array (data), where each column contains a column of APASS data
    data[:,0] = right ascension in degrees
    data[:,1] = error in right ascension in arcseconds
    data[:,2] = declination in degrees
    data[:,3] = error in declination in arcseconds
    data[:,4] = number of observations
    data[:,5] = Johnson V magnitude
    data[:,6] = error in Johnson V magnitude
    data[:,7] = Johnson B magnitude
    data[:,8] = error in Johnson B magnitude
    data[:,9] = Sloan g' magnitude
    data[:,10] = error in Sloan g' magnitude
    data[:,11] = Sloan r' magnitude
    data[:,12] = error in Sloan r' magnitude
    data[:,13] = Sloan i' magnitude
    data[:,14] = error in Sloan i' magnitude

    """
    # Create URL by which to access catalogue
    address='http://www.aavso.org/cgi-bin/apass_download.pl?ra={0}&dec={1}&radius={2}&outtype=0'.format(ra,dec,fov)
    catalog = url.urlopen(address)
    # Do a lot of terrible data parsing
    data = catalog.read()
    catalog.close()
    data = data.replace('<Td><font size=-1>','')
    data = data.replace('</td>','')
    data = data[577:len(data)-256]
    data = data.split('\n<tr>')
    data = array(data)[where(array(data) != '')]
    data = array([i.replace(' ','').replace('NA','-1').split('\n\t') for i in data])
    data = delete(data,5,axis=1)
    data = data.astype('float')
    # Save catalogue to file
    savetxt(APASSdir+'ra{0}dec{1}fov{2}.apass'.format(ra,dec,fov),data)
