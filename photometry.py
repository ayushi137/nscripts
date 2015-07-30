import numpy as np
nmin = np.min
nmax = np.max
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from numpy import *
import os
from scipy.optimize import leastsq
from callAPASS import callAPASS
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import gridspec

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 20}

matplotlib.rc('font', **font)

gain = 0.37

def magnitude(m0,counts):
    return m0 - 2.5*log10(counts)

def residuals(m0,counts,realms):
    return abs(realms - magnitude(m0,counts))

def tokjypersr(image,pixscale,m0):
    tojy = (3631*u.Jy)*(image/(pixscale*u.arcsec)**2)*10**(-0.4*m0)
    return array(tojy.to(u.kJy/u.sr))

def APASScheck(APASSdir,scanl,scanr,scand,scanu,racf,deccf):
    fs = os.listdir(APASSdir)
    use = 0
    i = 0
    while use == 0:
        try:
            s = fs[i].replace('ra','').replace('dec','-').replace('fov','-').replace('.apass','')
            s = s.split('-')
            ra = float(s[0])
            dec = float(s[1])
            fov = int(s[2])
            rad1 = sqrt((scanl-ra)**2 + (scand-dec)**2)
            rad2 = sqrt((scanl-ra)**2 + (scanu-dec)**2)
            rad3 = sqrt((scanr-ra)**2 + (scand-dec)**2)
            rad4 = sqrt((scanr-ra)**2 + (scanu-dec)**2)
            if int(rad1) <= fov and int(rad2) <= fov and int(rad3) <= fov and int(rad4) <= fov:
                print 'Successfully found extant APASS catalogue'
                apass = loadtxt(APASSdir+fs[i],usecols = (0,2,9,11))
                apass = apass[where((apass[:,0] != -1) & (apass[:,1] != -1) & (apass[:,2] != -1) & (apass[:,3] != -1))]
                use = 1
            if abs(ra-racf) < 1e-5 and abs(dec-deccf) < 1e-5:
                print 'Successfully found extant APASS catalogue'
                apass = loadtxt(APASSdir+fs[i],usecols = (0,2,9,11))
                apass = apass[where((apass[:,0] != -1) & (apass[:,1] != -1) & (apass[:,2] != -1) & (apass[:,3] != -1))]
                use = 1  
            else:
                i+=1
        except IndexError:
            print 'Creating APASS catalogue'
            callAPASS(racf,deccf,5,APASSdir)
            apass = loadtxt(APASSdir+'ra{0}dec{1}fov5.apass'.format(racf,deccf),usecols = (0,2,9,11))
            apass = apass[where((apass[:,0] != -1) & (apass[:,1] != -1) & (apass[:,2] != -1) & (apass[:,3] != -1))]  
    return apass

def photocheck(outname):
    create = True
    if os.path.isfile(outname) == True:
        o = fits.open(outname)
        oheader = o[0].header
        try:
            if oheader['UNITS'] == 'kJy/sr':
                create = False
            elif oheader['UNITS'] == 'ADU':
                create = True
        except KeyError:
            create = True
    elif os.path.isfile(outname) == False:
        create = True
    return create

def photometry(fname,catname,APASSdir,outname,threshold=30./3600.,stdcutoff=2,
               pixscale=2.85,plot=False,create = 0):
    print 'Unload input'
    f = fits.open(fname)
    data = f[0].data
    header = f[0].header
    f.close()
    if create == 0:
        create = photocheck(outname)
    print create
    if create == True:
        #data /= gain
        filt = header['FILTNAM']
        ra = header['RA']
        dec = header['DEC'].replace('d',':')
        coords = SkyCoord(ra,dec,unit = (u.hourangle, u.deg))
        ra = coords.ra.deg
        dec = coords.dec.deg
 
        print 'Read in source catalogue'
        sexcat = loadtxt(catname)
        sexcat = sexcat[where(sexcat[:,8] == 0)]
        source_ras = sexcat[:,4]
        source_decs = sexcat[:,5]
        source_fluxes = sexcat[:,1]
        
        scanl = min(source_ras)
        scanr = max(source_ras)
        scand = min(source_decs)
        scanu = max(source_decs)
        
        print 'Check for APASS data'
        apass = APASScheck(APASSdir,scanl,scanr,scand,scanu,ra,dec)
        goodapass = where((apass[:,0] > scanl) & (apass[:,0] < scanr) & (apass[:,1] > scand) & (apass[:,1] < scanu))
        apass = apass[goodapass]
        ras = apass[:,0]
        decs = apass[:,1]
        sloangs = apass[:,2]
        sloanrs = apass[:,3]
        print len(apass),len(sexcat)
        
        print 'Matching source catalogue with APASS'
        apassinds = []
        sexcainds = []
        discrepancy = []
        probsexca = []
        probapass = []
        probdiscrepancy = []
        
        for i in range(len(sexcat)):
            ra_diff = ras - source_ras[i]
            dec_diff = decs - source_decs[i]
            rad = sqrt(ra_diff**2 + dec_diff**2)
            possiblematch = [r for r in rad if r < threshold]
            if possiblematch != []:
                apassind = where(rad==min(possiblematch))[0][0]
                if apassind not in apassinds:
                    apassinds.append(apassind)
                    sexcainds.append(i)
                    discrepancy.append(possiblematch)
                elif apassind in apassinds:
                    probapass.append(apassind)
                    probsexca.append(i)
                    probdiscrepancy.append(possiblematch)

        for a in probapass:
            try:
                bad = apassinds.index(a)
                del apassinds[bad]
                del sexcainds[bad]
            except ValueError:
                continue
        
        print filt
        if filt == 'SloanG':
            fluxes = sloangs
        if filt == 'SloanR':
            fluxes = sloanrs
    
        print 'Making preliminary fit'
        imagefluxes = source_fluxes[sexcainds]
        apassmagnis = fluxes[apassinds]
        if plot == True:
            plt.figure(figsize = (12,10))
            plt.plot(apassmagnis,-2.5*log10(imagefluxes),'.')
            plt.xlabel('g')
            plt.ylabel('uncorrected magnitudes')
            plt.title('Uncorrected image source magnitude vs sloan magnitudes')
        
        m0 = 27
        zpmag = leastsq(residuals, m0, args = (imagefluxes,apassmagnis),full_output=1)
    
        print 'Scrubbing outliers and making final fit'
        res = residuals(zpmag[0],imagefluxes,apassmagnis)
        good = where(res <= stdcutoff*std(res))
    
        imagefluxes = imagefluxes[good]
        apassmagnis = apassmagnis[good]
        
        zpmagi2 = leastsq(residuals,zpmag[0],args = (imagefluxes,apassmagnis),full_output=1)
    
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
        success = zpmagi2[-1]
        allowedsuccess = [1,2,3,4]
        if success in allowedsuccess:
            print 'Successful fit'
            data= tokjypersr(data,pixscale,zpmagi2[0])
            kJperADU = float(tokjypersr(1,pixscale,zpmagi2[0]))
            header['UNITS'] = 'kJy/sr'
            header['M0'] = float(zpmagi2[0])
            header['kJpADU'] = (kJperADU,'kJy/sr per ADU/pixel')
        if not success in allowedsuccess:
            print 'No zero point magnitude found'
            header['UNITS'] = 'ADU'
            header['M0'] = 'N/A'
        fits.writeto(outname,data,header,clobber=True)
        return data,header,zpmagi2[0]
    if create == False:
        return fits.getdata(outname),header,0
