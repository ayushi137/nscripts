"""
photometry - a set of functions to perform photometry on an image

Requires the following packages: numpy, scipy, astropy, matplotlib, os
                                 docopt
Requires the following files: callAPASS.py

Contains the following functions: magnitude, residuals, tokjypersr,
                                  APASScheck, photocheck, photometry
                                  catheader
"""

############################## IMPORT PACKAGES ##########################

import urllib as url
import numpy as np
nmin = np.min
nmax = np.max
from astropy.io import fits
from astropy import units as u
from astropy import wcs
from numpy import *
import os
from scipy.optimize import leastsq
import matplotlib.pyplot as plt
import matplotlib

############################## FORMAT PLOTS ##########################


font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 20}

matplotlib.rc('font', **font)


############################## FUNCTIONS ##########################

def callAPASS(ra,dec,APASSdir):
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
    data = np.array(data)[np.where(np.array(data) != '')]
    data = np.array([i.replace(' ','').replace('NA','-1').split('\n\t') for i in data])
    data = np.delete(data,5,axis=1)
    data = data.astype('float')
    # Save catalogue to file
    savetxt(APASSdir+'ra{0}dec{1}fov{2}.apass'.format(ra,dec,fov),data)


def magnitude(m0,counts):
    """
    Converts counts [ADU] to apparent magnitude
    
    m0:         zero point magnitude
    counts:     value(s) in ADU to be converted

    Returns an array of apparent magnitudes in the same shape as counts.

    """
    return m0 - 2.5*log10(counts)

def residuals(m0,counts,realms):
    """
    Finds the absolute difference between catalogue magnitudes and those
        calculated with magnitude(m0,counts)

    m0:         zero point magnitude
    counts:     value(s) in ADU to be converted
    realms:     catalogue magnitude(s) in the same shape as counts

    Returns an array in the same shape as counts and realms

    """
    return abs(realms - magnitude(m0,counts))

def tokjypersr(image,pixscale,m0):
    """
    Converts values from units of ADU/pixel to units of kJy/sr

    image:      value(s) in ADU to be converted
    pixscale:   arcseconds per pixel
    m0:         zero point magnitude

    Returns an array with units kJy/sr in the same shape as image

    """
    if isinstance(pixscale,float):
        tojy = (3631*u.Jy)*(image/(pixscale*u.arcsec)**2)*10**(-0.4*m0)
    if isinstance(pixscale,(list,ndarray)):
        tojy = (3631*u.Jy)*(image/(pixscale[0]*u.arcsec*pixscale[1]*u.arcsec))*10**(-0.4*m0)
    return array(tojy.to(u.kJy/u.sr))

def catheader(catname):
    """
    Reads the header of the catalogue file and creates a dictionary specifying
        Python indices of each row

    catname:    Path to catalogue file

    Reuturns dictionary of header information

    """
    # Create dictionary
    header = {}
    # Cycle through lines
    with open(catname) as f:
        for line in f:
            # Read only commented lines, then quit
            if '#' in line:
                L = line.split()
                header[L[2]] = int(L[1])-1
            else:
                break
    return header

def APASScheck(APASSdir,scanl,scanr,scand,scanu,racf,deccf,newfov=7,fulloutput = False):
    """
    Check if the relevant part of the APASS catalogue is already in 
        APASSdir. If not, download it.

    APASSdir:      directory containing APASS catalogue files
    scanl:         minimum R.A. value for sources in source 
                   extractor catalogue
    scanr:         maximum R.A. value for sources in source 
                   extractor catalogue
    scand:         minimum declination value for sources in source 
                   extractor catalogue
    scanu:         maximum declination value for sources in source 
                   extractor catalogue
    racf:          central R.A. value from .fits header
    deccf:         central declination value from .fits header
    newfov:        determine fov of new APASS catalogue 
                   (kwarg, default = 7)

    Returns array representing table of APASS catalogue.

    """
    # List all files in APASSdir
    fs = os.listdir(APASSdir)
    # Boolean switch indicating whether a file should be used
    use = False
    i = 0
    # While not to be used, check each file
    while not use:
        try:
            # Slice up the file name and extract R.A., dec and fov
            s = fs[i].replace('ra','').replace('dec','-')
            s = s.replace('fov','-').replace('.apass','')
            s = s.split('-')
            ra = float(s[0])
            dec = float(s[1])
            fov = int(s[2])
            # Find the distance between the centre of the APASS catalogue
            # and each of the extreme R.A. and dec values
            rad1 = sqrt((scanl-ra)**2 + (scand-dec)**2)
            rad2 = sqrt((scanl-ra)**2 + (scanu-dec)**2)
            rad3 = sqrt((scanr-ra)**2 + (scand-dec)**2)
            rad4 = sqrt((scanr-ra)**2 + (scanu-dec)**2)
            # If all distances less than the fov, use current file
            if int(rad1) <= fov and int(rad2) <= fov and int(rad3) <= fov and int(rad4) <= fov:
                print 'Successfully found extant APASS catalogue'
                # Load APASS catalogue
                if fulloutput == False:
                	apass = loadtxt(APASSdir+fs[i],usecols = (0,2,9,11))
                # Remove any row with a -1 entry
                # -1 entry implies original information missing
               		apass = apass[where((apass[:,0] != -1) & (apass[:,1] != -1) & (apass[:,2] != -1) & (apass[:,3] != -1))]
               	if fulloutput == True:
               		apass = loadtxt(APASSdir+fs[i])
               		apass = apass[where((apass[:,0] != -1) & (apass[:,1] != -1) & (apass[:,2] != -1) & (apass[:,3] != -1))]
                # End while loop
                use = True
            # If this file does not work, move to the next
            else:
                i+=1
        # If the end of available files is reached create catalogue
        except IndexError:
            print 'Creating APASS catalogue'
            # Create catalogue
            callAPASS(racf,deccf,APASSdir,fov = newfov)
            print 'Made file ra{0}dec{1}fov{2}.apass'.format(racf,deccf,newfov)
            # Read in created catalogue
            if fulloutput == False:
            	apass = loadtxt(APASSdir+'ra{0}dec{1}fov{2}.apass'.format(racf,deccf,newfov),usecols = (0,2,9,11))
            # Remove any row with a -1 entry
            # -1 entry implies original information missing
           		apass = apass[where((apass[:,0] != -1) & (apass[:,1] != -1) & (apass[:,2] != -1) & (apass[:,3] != -1))] 
           	if fulloutput == True:
           		apass = loadtxt(APASSdir+'ra{0}dec{1}fov{2}.apass'.format(racf,deccf,newfov))
           		apass = apass[where((apass[:,0] != -1) & (apass[:,1] != -1) & (apass[:,2] != -1) & (apass[:,3] != -1))] 
    return apass

def photocheck(outname):
    """
    Checks whether photometry has already been done on outname

    outname:    path to where photometered file is to be saved

    Returns False if outname exists, and has the correct units
    otherwise return True.
    """
    # By default, return True unless changed
    create = True
    # If file exists, read header
    if os.path.isfile(outname) == True:
        o = fits.open(outname)
        oheader = o[0].header
        try:
            # If units keyword is set, check
            # whether it is kJy/sr - if not, set True
            if oheader['UNITS'] == 'kJy/sr':
                create = False
            elif oheader['UNITS'] == 'ADU':
                create = True
        # If units keyword not set, set True
        except KeyError:
            create = True
    # If file does not exist, set True
    elif os.path.isfile(outname) == False:
        create = True
    return create

def photometry(fname,catname,APASSdir,outname,threshold=30./3600.,stdcutoff=2,
               plot=False, create = 0, m0 = 27):
    """
    Performs photometry on an image

    fname:      path to file to perform photometry on
    catname:    path to source catalogue created by source extractor
    APASSdir:   directory where APASS catalogues are
    outname:    path to file to save photometered data
    threshold:  specifies how close catalogue and images need to be 
                in arcseconds (kwarg, default = 30./3600.)
    stdcutoff:  multiple of the standard deviation in the residual source
                magnitude distribution beyond which to cutout sources
                (kwarg, default = 2)
    m0:         initial guess for zero point magnitude fitting
                (kwarg, default = 27)
      


    """
    # Import data and header
    f = fits.open(fname)
    data = f[0].data
    header = f[0].header
    f.close()
    # If not force creation, check if creation needs to be forced
    if create == 0:
        create = photocheck(outname)
    # If creation forced, process photometry
    if create == True:
        # Read in relevant header data
        # Find filter type
        filt = header['FILTNAM']
        # Use WCS solution to find coordinates of center of field
        fwcs = wcs.WCS(header)
        cf = fwcs.wcs_pix2world([[data.shape[0]/2,data.shape[1]/2]],0)
        ra = cf[0,0]
        dec = cf[0,1]
        # Find pixel scales
        pscalex = header['PSCALX']
        pscaley = header['PSCALY']
        pixscale = [pscalex,pscaley]
 
        # Read in catalogue information
        sexcat = loadtxt(catname)
        cheader = catheader(catname)
        # Remove all entries that sextractor flagged
        sexcat = sexcat[where(sexcat[:,cheader['FLAGS']] == 0)]
        # Grab relevant columns
        source_ras = sexcat[:,cheader['ALPHA_J2000']]
        source_decs = sexcat[:,cheader['DELTA_J2000']]
        source_fluxes = sexcat[:,cheader['FLUX_BEST']]
        
        # Find R.A./dec limits of detected sources
        scanl = min(source_ras)
        scanr = max(source_ras)
        scand = min(source_decs)
        scanu = max(source_decs)
        
        # Get APASS catalogue
        apass = APASScheck(APASSdir,scanl,scanr,scand,scanu,ra,dec)
        # Truncate catalogue to only sources within the image
        goodapass = where((apass[:,0] > scanl) & (apass[:,0] < scanr) & (apass[:,1] > scand) & (apass[:,1] < scanu))
        apass = apass[goodapass]
        # Extract columns
        ras = apass[:,0]
        decs = apass[:,1]
        sloangs = apass[:,2]
        sloanrs = apass[:,3]
        
        # initialize lists to match sources
        apassinds = []
        sexcainds = []
        probapass = []
        
        # for sources detected by source extractor
        for i in range(len(sexcat)):
            # find the distance between the source and every
            # catalogue star
            ra_diff = ras - source_ras[i]
            dec_diff = decs - source_decs[i]
            rad = sqrt(ra_diff**2 + dec_diff**2)
            # create a list of possible matches, catalogue stars
            # that are closer than the threshold
            possiblematch = [r for r in rad if r < threshold]
            # if any matches are found
            if possiblematch != []:
                # choose the closest catalogue star as the match
                apassind = where(rad==min(possiblematch))[0][0]
                # if that catalogue star not already matched, add it 
                # to the list
                if apassind not in apassinds:
                    apassinds.append(apassind)
                    sexcainds.append(i)
                elif apassind in apassinds:
                    probapass.append(apassind)

        # If a catalogue star was matched to two source extractor
        # stars, remove both of those stars
        for a in probapass:
            try:
                bad = apassinds.index(a)
                del apassinds[bad]
                del sexcainds[bad]
            except ValueError:
                continue

        # Choose appropriate fluxes
        if filt == 'SloanG':
            fluxes = sloangs
        if filt == 'SloanR':
            fluxes = sloanrs
    
        # Start preliminary fit
        # Retriever star fluxes from in ADU from sextractor and in
        # magnitudes from APASS
        imagefluxes = source_fluxes[sexcainds]
        apassmagnis = fluxes[apassinds]
        if plot != False:
            plt.figure(figsize = (12,10))
            plt.plot(apassmagnis,-2.5*log10(imagefluxes),'.')
            plt.xlabel('g')
            plt.ylabel('uncorrected magnitudes')
            plt.title('Uncorrected image source magnitude vs sloan magnitudes')
            plt.savefig(plot+'_rawsample.png')
            plt.close()

        # Make preliminary fit
        zpmag = leastsq(residuals, m0, args = (imagefluxes,apassmagnis),full_output=1)
        success = zpmag[-1]
        allowedsuccess = [1,2,3]
        if success in allowedsuccess:
        # Find difference between preliminary fit and APASS magnitudes, and remove
        # outliers
            res = residuals(zpmag[0],imagefluxes,apassmagnis)
            good = where(res <= stdcutoff*std(res))
            
            imagefluxes = imagefluxes[good]
            apassmagnis = apassmagnis[good]
            
        # Perform second iteration of fitting
            
            zpmagi2 = leastsq(residuals,zpmag[0],args = (imagefluxes,apassmagnis),full_output=1)
    
        # Plot some statistics of the fit
            if plot != False:
                plt.figure(figsize =(12,10))
                plt.plot(apassmagnis,magnitude(zpmagi2[0],imagefluxes),'.')
                plt.xlabel('g')
                plt.ylabel('$m_0$ - 2.5$log_{10}(ADU)$')
                plt.title('Final sample of sources used to fit for the zero point magnitude')
                plt.savefig(plot+'_sample.png')
                plt.close()
                plt.figure(figsize=(12,10))
                plt.plot(apassmagnis,apassmagnis+2.5*log10(imagefluxes),'.')
                plt.axhline(zpmagi2[0],color = 'red',linewidth = 4)
                plt.title('Final $m_0$ = {0}'.format(zpmagi2[0]))
                plt.ylabel('($2.5log_{10}$(ADU)) + g')
                plt.xlabel('g')
                plt.savefig(plot+'_m0.png')
                plt.close()
        # Determine whether fit was successful
            success = zpmagi2[-1]
            
            if success in allowedsuccess:
                # Convert image units
                kJperADU = float(tokjypersr(1,pixscale,zpmagi2[0]))
                data *= kJperADU
                header['UNITS'] = 'kJy/sr'
                header['M0'] = float(zpmagi2[0])
                header['kJpADU'] = (kJperADU,'kJy/sr per ADU/pixel')
            elif not success in allowedsuccess:
                print 'No zero point magnitude found'
                header['UNITS'] = 'ADU'
                header['M0'] = 'N/A'
        # Write data to file
            fits.writeto(outname,data,header,clobber=True)
            return data,header,zpmagi2[0]
        elif success not in allowedsuccess:
            print 'No zero point magnitude found'
            header['UNITS'] = 'ADU'
            header['M0'] = 'N/A'
            fits.writeto(outname,data,header,clobber=True)
            return data,header,zpmag[0]
    if create == False:
        return fits.getdata(outname),header,0
