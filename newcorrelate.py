#!/usr/bin/env python

"""
correlate - correlate Dragonfly images against Herschel images

Expects a particular file system set up for Dragonfly images: 
<parent directory>/<object name>/<subdirectory>/
This may be followed by further subdirectories by date, or contain the files
themselves.

Requires the following software: sextractor, astrometry.net
Requires the following packages: numpy, scipy, astropy, docopt, matplotlib, os
Requires the following files:    photometry.py, resconvolve.py, maskdata.py
                                 regrid.py, backgroundplane.py

Usage:
correlate [-hvlgpcrmb] [-d DIRECTORY] [-u DIRECTORIES] [-o OBJECTNAMES] 
	      [-s DIRECTORY] [-f FILEPATHS] [-x DIRECTORY] [-a DIRECTORY]

Options:
    -h, --help                        Show this screen
    -v, --verbose                     Print processing steps and extra 
                                      information
    -l, --lowmemory                   If True, use a version of regridding code
                                      that does not require high memory usage

    -d DIRECTORY, --IOdir DIRECTORY   Location of input files (parent directory)
                                      [default: /mnt/gfsproject/naiad/njones/moddragonflydata/]
    -u LIST, --unmod LIST             A list of subdirectories of --IOdir 
    								  containing raw files and calibration frames
    								  [default: /raw_lights/, /calframes/]
    -o OBJECTS, --objects OBJECTS     Provide a list of object subdirectory 
                                      names as a string
                                      [default: PGM_1_2, spi1_1]
    -s DIRECTORY, --sub DIRECTORY     Subdirectory containing flat-fielded and
                                      dark subtracted images
                                      [default: /NCalibrated/] 
    -f FILES, --filelist FILES        Provide a list of file paths to process
                                      as a string. If empty, processes all file
                                      in input directory
                                      [default: ]
    -x DIRECTORY, --cross DIRECTORY   Location of files to correlate with
    								  [default: ../herschel/]
    -a DIRECTORY, --apass DIRECTORY   Location of APASS catalogues
                                      [default: APASS/]
Testing Options:
    -g, --generate                    If False, do not generate data from
                                      any of the following substeps unless that
                                      data is missing
                                      If True, force all data to be generated
    -n, --darksub					  If False, do not create master darks and 
    								  subtract them unless files are missing
    -q, --flatfield 				  If False, do not create master flats and 
    								  divide by them unless files are missing
    -p, --photometry                  If False, do not perform photometry on an
                                      image unless photometry data is missing
                                      If True, force photometry to be done
    -c, --convolve                    If False, do not convolve an image unless
                                      convolved image missing
                                      If True, force photometry to be done
    -r, --regrid                      If False, do not regrid an image unless
                                      regridded image is missing
                                      If True, force regridding to be done
    -m, --mask                        If False, do not mask an image unless
                                      masked image is missing
                                      If True, force masking to be done
    -b, --backsub                     If False, do not background plane subtract
                                      an image unless subtracted is missing
                                      If True, force background plane 
                                      subtraction to be done


"""

########################### IMPORT BASE PACKAGES ############################

import os
import docopt
import numpy as np
nmax = np.max
nmin = np.min
from numpy import *
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib.colors import LogNorm

########################### COMMAND LINE ARGUMENTS #############################

arguments = docopt.docopt(__doc__)

# Non-mandatory options without arguments

VERBOSE = arguments['--verbose']
LOWMEMORY = arguments['--lowmemory']

# Non-mandatory options with arguments

directory = arguments['--IOdir']
rawdirs = arguments['--unmod']
rawdir = rawdirs.split(', ')[0]
caldir = rawdirs.split(', ')[1]
subdir = arguments['--sub']
APASSdir = arguments['--apass']
filelist = arguments['--filelist']
objects = arguments['--objects']
corrdir = arguments['--cross']

# Testing options

GENERATE = arguments['--generate']
DARKSUB = arguments['--darksub']
FLATFIELD = arguments['--flatfield']
PHOTOMETRY = arguments['--photometry']
CONVOLVE = arguments['--convolve']
REGRID = arguments['--regrid']
MASK = arguments['--mask']
BACKSUB = arguments['--backsub']

########################### IMPORT FUNCTIONS ############################

from photometry import *
from resconvolve import resconvolve
from maskdata import maskdata
from regrid import reshape
from backgroundplane import subBGplane

if LOWMEMORY:
	from regrid import regrid_lowmemory as regrid
elif not LOWMEMORY:
	from regrid import regrid

'''
#################################### FUNCTIONS #################################

def sexcalls(f,di,odi,cdi,bdi):
	"""
	Run source extractor on f, producing a catalogue, and object and background maps

	f:		file name to run source extractor on
	di:		directory path to f
	odi:	directory to save object map in
	cdi:	directory to save catalogue in
	bdi:	directory to save background map in

	Returns nothing explicitly, but runs source extractor on the file

	"""
	# Split to make file name for catalogue, object map and background map filenames
    fname = f.split('.fits')[0]
    # Construct source extractor call
    objsexcall = 'sex -PARAMETERS_NAME photo.param -CATALOG_NAME '+cdi+fname+'.cat'
    objsexcall += ' -CHECKIMAGE_TYPE OBJECTS, BACKGROUND -CHECKIMAGE_NAME '+odi+fname
    objsexcall += '_objects.fits '+di+f+', '+bdi+fname+'_background.fits '+di+f
    os.system(objsexcall)

def hist2d(x,y,nbins,maskval = 0,saveloc = '',labels=[]):
	"""
	Creates a 2D histogram from data given by numpy's histogram

	x,y:		two 2D arrays to correlate, the second masked
	nbins:		number of bins
	maskval:	value that indicates masked areas
				(kwarg, default = 0)
	saveloc:	place to save histogram plot - if unspecified, do not
				save plot (kwarg, default = '')
	labels:		labels for histogram plot, with the following format
				[title,xlabel,ylabel,zlabel] - if unspecified, do 
				not label plot (kwarg, default = [])

	Returns the edges of the histogram bins and the 2D histogram

	"""
	# Remove NAN values
    a = where(isnan(x) == False)
    x = x[a]
    y = y[a]
    b = where(isnan(y) == False)
    x = x[b]
    y = y[b]
    # Remove masked areas
    c = where(y != maskval)
    x = x[c]
    y = y[c]
    # Create histogram
    H,xedges,yedges = histogram2d(x,y,bins=nbins)
    # Reorient appropriately
    H = rot90(H)
    H = flipud(H)
    # Mask zero value bins
    Hmasked = ma.masked_where(H==0,H)

    # Begin creating figure
    plt.figure(figsize=(12,10))
    # Use logscale
    plt.pcolormesh(xedges,yedges,Hmasked,
                   norm = LogNorm(vmin = Hmasked.min(),vmax = Hmasked.max()))
    cbar = plt.colorbar()
    # Add labels
    if labels != []:
        title,xlabel,ylabel,zlabel = labels
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        cbar.ax.set_ylabel(zlabel)
    # Save plot
    if saveloc != '':
    	plt.savefig(saveloc)
    plt.close()
    # Return histogram
    return xedges,yedges,Hmasked

############################# DIRECTORIES ##################################

# location of images to correlate against
herdir = '../herschel/'

# specify herschel resolution at each wavelength
SPIRE = {'PSW':17.6,'PMW':23.9,'PLW':35.2}
spirekeys = SPIRE.keys()

# if necessary, specify observation date subdirectory names for each object 
# in a dictionary
objectdates = {}
objectdates['PGM_1_2'] = ['2013-09-23','2013-09-24','2013-09-25','2013-09-26',
                          '2013-09-27']
objectdates['spi1_1'] = ['2014-10-31','2014-11-01','2014-11-19',
                         '2014-11-22','2014-11-23','2014-11-24']
keys = objectdates.keys()

# indicate object names in herschel file name formatting
objectnames = {}
objectnames['spi1_1'] = 'spider_SPIRE_'
objectnames['PGM_1_2'] = 'draco_'

# source extractor object maps
objdir = '/objects/'
# source extractor background maps
bakdir = '/background/'
# source extractor catalogue
catdir = '/catalogue/'

# OUTPUT DIRECTORIES
# location of helper plots from zero point magnitude calculations
magdir = '/magnitudecalcs/'
# location of photometered fits files
outdir = '/photometered/'
# location of regridded images
regdir = '/regrid/'
# location of background plane subtracted images
bsudir = '/backgroundsub/'
# location of correlation plots
cordir = '/correlations/'

# if any of the output directories are missing, create them now
dirlist = [outdir,regdir,cordir,magdir,bsudir]
for d in dirlist:
    for key in keys:
        if os.path.isdir(directory+key+d) == False:
            os.system('mkdir '+directory+key+d)

################################ STEPS #####################################

# DATE FORMATTING
DATES = True
# if set to true, directory structure is assumed to include sudirectories for
# each date of observations



for key in keys:
    print 'OBJECT '+key
    dates = objectdates[key]
#    dates = ['2013-09-23']
    for date in dates:
        print 'DATE '+date
        di = directory+key+subdir+date+'/'
        odi = directory+key+objdir+date+'/'
        cdi = directory+key+catdir+date+'/'
        bdi = directory+key+bakdir+date+'/'
        pdi = directory+key+outdir+date+'/'
        rdi = directory+key+regdir+date+'/'
        sdi = directory+key+cordir+date+'/'
        mdi = directory+key+magdir+date+'/'
        gdi = directory+key+bsudir+date+'/'
        hdi = herdir+key+'/'
        dis = [pdi,rdi,sdi,mdi,gdi]
        for d in dis:
            if os.path.isdir(d) == False:
                os.system('mkdir '+d)
        files = os.listdir(di)
        files.sort()
 #       files = ['83F010826_17_light_ds_ff.fits']
        for f in files:
            try:
                print 'FILE '+di+f
                spl = f.split('.fits')[0]
                catdata = loadtxt(cdi+spl+'.cat')
                dflyimage = fits.open(di+f)
                dflydata = dflyimage[0].data
                dflyheader = dflyimage[0].header
                dflyimage.close()
                pscalx = dflyheader['PSCALX']
                pscaly = dflyheader['PSCALY']
                pixscale = (pscalx+pscaly)/2.
                catdata = catdata[where(catdata[:,8] == 0)]
                pfwhm = mean(catdata[:,6])
                afwhm = mean(catdata[:,7])*3600.
                #print pfwhm,pfwhm*2.85,afwhm
                dflybeam = pfwhm*pixscale
                dflyheader['FWHM'] = (dflybeam, 'arcseconds')
                if dflyheader['FILTNAM'] == 'Pol':
                    print 'polarization data'
                    continue              
                pname = spl+'_photo.fits'
                outplots = mdi+spl
                print 'Do photometry'
                if os.path.isfile(pdi+pname) == True and generate == False:
                    pdata,dflyheader = fits.getdata(pdi+pname,header=True)
                    zp = dflyheader['M0']
                elif os.path.isfile(pdi+pname) != True or generate == True:
                    pdata,dflyheader,zp = photometry(di+f,cdi+spl+'.cat',APASSdir,pdi+pname,
                                                     create = True,plot = outplots)
                if os.path.isfile(odi+spl+'_photo_objects.fits') != True or generate == True:
                    sexcalls(pname,pdi,odi,cdi,bdi)
                dflyheader['kJpADU'] = (float(tokjypersr(1,pixscale,zp)),'kJy/sr per ADU/pixel')
                dflyheader['FWHM'] = (dflybeam, 'arcseconds')
                n = where(isnan(pdata) == True)
                pdata[n] = 0
                pstar = fits.getdata(odi+spl+'_photo_objects.fits')
                for skey in spirekeys:
                    if 'draco' in objectnames[key]:
                        hername = hdi+'reindex_'+objectnames[key]
                        hername += skey+'_sanepic.fits'
                    elif 'spider' in objectnames[key]:
                        hername = hdi+'reindex_'+objectnames[key]
                        hername += skey+'_map_sanepic.fits'
                    pspl = pname.split('.fits')[0]
                    herbeam = SPIRE[skey]
                    cname = pspl+'_convto'+str(herbeam)+'.fits'
                    ocname = pspl+'_convto'+str(herbeam)+'_objects.fits'
                    print 'Convolve data'
                    if os.path.isfile(pdi+cname) == True and generate == False:
                        cdata,dflyheader = fits.getdata(pdi+cname,header=True)
                    elif os.path.isfile(pdi+cname) != True or generate == True:
                        cdata,dflyheader = resconvolve(pdata,dflybeam,
                                                       herbeam,
                                                       outfile = pdi+cname,
                                                       header = dflyheader)
                    if os.path.isfile(odi+ocname) == True and generate == False:
                        ocdata,oheader = fits.getdata(odi+ocname,header=True)
                    elif os.path.isfile(odi+ocname) != True or generate == True:
                        ocdata,oheader = resconvolve(pstar,dflybeam,herbeam,
                                                     outfile = odi+ocname,
                                                     header = dflyheader)
                    if dflyheader == 0:
                        raise ValueError('Convolution failed')
                    cutoff = 10*median(ocdata)
                    print 'cutoff = ',cutoff
                    mspl = cname.split('.fits')[0]
                    mname = mspl+'_mask.fits'
                    ds9 = '/home/njones/ds9 '+pdi+mname+' -histequ -zoom to fit &'
                    print 'Regrid and mask'
                    if os.path.isfile(pdi+mname) == True and generate == False:
                        mdata,dflyheader = fits.getdata(pdi+mname,header=True)
                        target = fits.getdata(hername)
                    if os.path.isfile(pdi+mname) != True or generate == True or dflyheader['MASKCUT'] - cutoff > 1e-3:
                        cdata,target = regrid(pdi+cname,hername)
                        ocdata,target = regrid(odi+ocname,hername)
                        mdata,dflyheader = maskdata(cdata,ocdata,cutoff,
                                                    outfile = pdi+mname,
                                                    header = dflyheader)
                    #os.system(ds9)
                    print 'Reshaping'
                    mspl = mname.split('.fits')[0]
                    ress = sqrt((dflybeam/s2f)**2+(herbeam/s2f)**2)
                    r = reshape(mdata,mdata,2*ress)
                    t = reshape(target,mdata,2*ress)
                    dflyheader['kJpADU'] = (float(tokjypersr(1,pixscale,zp)),'kJy/sr per ADU/pixel')
                    dflyheader['FWHM'] = (dflybeam, 'arcseconds')
                    rname = mspl+'_regrid.fits'
                    hname = hername.split('.fits')[0]
                    hname = hname+'_reshaped.fits'
                    fits.writeto(rdi+rname,r,dflyheader,clobber=True)
                    rspl = rname.split('.fits')[0]
                    saveloc = sdi+mspl+'_'+skey+'.png'
                    title = 'Correlation between Dragonfly and Herschel'
                    xlabel = 'Herschel [MJy/sr]'
                    ylabel = 'Dragonfly [kJy/sr]'
                    zlabel = 'Pixels'
                    labels = [title,xlabel,ylabel,zlabel]
                    H = hist2d(t,r,50,saveloc,labels = labels)
                    bname = rspl+'_backsub.fits'
                    p0 = [2,1,1,1]
                    newr,bg = subBGplane(r,t,p0)
                    if isinstance(newr,float) != True:
                        dflyheader['BACKSUB'] = 'TRUE'
                        fits.writeto(gdi+bname,newr,dflyheader,clobber=True)
                        fits.writeto(bdi+spl+'_bgplane.fits',bg,dflyheader,clobber=True)
                        print 'Plotting'
                        bspl = bname.split('.fits')[0]
                        saveloc = sdi+bspl+'_'+skey+'.png'
                        title = 'Correlation between Dragonfly and Herschel'
                        xlabel = 'Herschel [MJy/sr]'
                        ylabel = 'Dragonfly [kJy/sr]'
                        zlabel = 'Pixels'
                        labels = [title,xlabel,ylabel,zlabel]
                        negs = where(r > 0)
                        H = hist2d(t[negs],newr[negs],50,saveloc,labels = labels)
                    elif isinstance(newr,float) == True:
                        print 'Failed background subtraction'
                    #fits.writeto(hname,t,clobber=True)
            except (ValueError,TypeError) as e:
                print e
                print di+f+' failed to correlate'
'''
