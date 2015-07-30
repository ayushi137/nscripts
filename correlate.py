import numpy as np
nmax = np.max
nmin = np.min
from numpy import *
import matplotlib.pyplot as plt
from photometry import *
from resconvolve import resconvolve
from maskdata import maskdata
from regrid import regrid,reshape
import os
from astropy.io import fits
from allsexd import sexcalls
from matplotlib.colors import LogNorm
from backgroundplane import subBGplane
#plt.ion()
generate = False

def hist2d(x,y,nbins,saveloc,labels=[]):
    a = where(isnan(x) == False)
    x = x[a]
    y = y[a]
    b = where(isnan(y) == False)
    x = x[b]
    y = y[b]
    c = where(y != 0)
    x = x[c]
    y = y[c]
    H,xedges,yedges = histogram2d(x,y,bins=nbins)
    H = rot90(H)
    H = flipud(H)
    Hmasked = ma.masked_where(H==0,H)
    #print Hmasked
    plt.figure(figsize=(12,10))
    plt.pcolormesh(xedges,yedges,Hmasked,
                   norm = LogNorm(vmin = Hmasked.min(),vmax = Hmasked.max()))
    cbar = plt.colorbar()
    if labels != []:
        title,xlabel,ylabel,zlabel = labels
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        cbar.ax.set_ylabel(zlabel)
    plt.savefig(saveloc)
    plt.close()
    return xedges,yedges,Hmasked
    
directory = '/mnt/gfsproject/naiad/njones/moddragonflydata/'
subdir = '/NCalibrated/'
objdir = '/objects/'
catdir = '/catalogue/'
outdir = '/photometered/'
APASSdir = '/mnt/gfsproject/naiad/njones/APASS/'
herdir = '../herschel/'
bakdir = '/background/'
regdir = '/regrid/'
cordir = '/correlations/'
magdir = '/magnitudecalcs/'
bsudir = '/backgroundsub/'
masdir = '/masks/'

dirlist = [outdir,regdir,cordir,magdir,bsudir,masdir]
for d in dirlist:
    if os.path.isdir(directory+'PGM_1_2'+d) == False:
        os.system('mkdir '+directory+'PGM_1_2'+d)
    if os.path.isdir(directory+'spi1_1'+d) == False:
        os.system('mkdir '+directory+'spi1_1'+d)



SPIRE = {'PSW':17.6,'PMW':23.9,'PLW':35.2}
spirekeys = SPIRE.keys()

objectdates = {}
objectdates['PGM_1_2'] = ['2013-09-23','2013-09-24','2013-09-25','2013-09-26',
                          '2013-09-27']
objectdates['spi1_1'] = ['2014-10-31','2014-11-01','2014-11-19',
                         '2014-11-22','2014-11-23','2014-11-24']
objectnames = {}
objectnames['spi1_1'] = 'spider_SPIRE_'
objectnames['PGM_1_2'] = 'draco_'

s2f = 2*sqrt(2*log(2))

keys = objectdates.keys()
keys = ['PGM_1_2']
for key in keys:
    print 'OBJECT '+key
    dates = objectdates[key]
    dates = ['2013-09-23']
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
        kdi = directory+key+masdir+date+'/'
        hdi = herdir+key+'/'
        dis = [pdi,rdi,sdi,mdi,gdi,kdi]
        for d in dis:
            if os.path.isdir(d) == False:
                os.system('mkdir '+d)
        files = os.listdir(di)
        files.sort()
        files = ['83F010826_17_light_ds_ff.fits']
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
                dflyheader['kJpADU'] = (float(tokjypersr(1,pixscale,zp)),'kJy/sr per ADU/pixel')

                #if os.path.isfile(odi+spl+'_photo_objects.fits') != True or generate == True:
                    
#sexcalls(pname,pdi,odi,cdi,bdi)
                pstar = fits.getdata(odi+spl+'_objects.fits')
                pstar *= dflyheader['kJpADU']
                fits.writeto(odi+spl+'_photo_objects.fits',pstar,dflyheader,clobber=True)
                dflyheader['FWHM'] = (dflybeam, 'arcseconds')
                n = where(isnan(pdata) == True)
                pdata[n] = 0
                pstar = fits.getdata(odi+spl+'_photo_objects.fits')
                for skey in spirekeys:
                    print 'Herschel ', skey
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
                        cdata,dflyheader = resconvolve(pdata,dflybeam/s2f,
                                                       herbeam/s2f,
                                                       outfile = pdi+cname,
                                                       header = dflyheader)
                    if os.path.isfile(odi+ocname) == True and generate == False:
                        ocdata,oheader = fits.getdata(odi+ocname,header=True)
                    elif os.path.isfile(odi+ocname) != True or generate == True:
                        ocdata,oheader = resconvolve(pstar,dflybeam/s2f,herbeam/s2f,
                                                     outfile = odi+ocname,
                                                     header = dflyheader)
                    if dflyheader == 0 or oheader == 0:
                        raise ValueError('Convolution failed')
                    cutoff = 10*median(ocdata)
                    howcutoff = 'Ten times the median = 10^'
                    print 'cutoff = ',cutoff
                    nbins = 200
                    values,base = histogram(cdata.flatten(),bins=nbins)
                    cumulative = cumsum(values)
                    tlen = len(cdata.flatten())
                    plt.figure(figsize = (12,10))
                    plt.title('Cumulative distribution of Pixel Brightness')
                    plt.subplot(211)
                    plt.plot(base[:-1],tlen-cumulative,linewidth = 3)
                    plt.xlim(max(base[:-1]),min(base[:-1]))
                    plt.axvline(cutoff,color = 'red',
                                label = howcutoff+str(round(log10(cutoff),2)))
                    plt.legend(loc = 'best')
                    plt.ylabel('Number of Pixels')
                    plt.subplot(212)
                    plt.semilogy(base[:-1],tlen-cumulative,linewidth = 3)
                    plt.xlim(max(base[:-1]),min(base[:-1])) 
                    plt.axvline(cutoff,color = 'red')
                    plt.xlabel('kJy/sr')
                    plt.ylabel('Number of Pixels')
                    plt.savefig('cumulativetest_'+skey+'.png')
                    plt.close()
 
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
                        nbins = 200
                        values,base = histogram(newr[where(r!= 0)],bins=nbins)
                        cumulative = cumsum(values)
                        tlen = len(newr[where(r!=0)])
                        plt.figure(figsize = (12,10))
                        plt.title('Cumulative distribution of Pixel Brightness')
                        plt.subplot(211)
                        plt.plot(base[:-1],tlen-cumulative,linewidth = 3)
                        plt.xlim(max(base[:-1]),min(base[:-1]))
                        #plt.axvline(cutoff,color = 'red',
                        #            label = howcutoff+str(round(log10(cutoff),2)))
                        #plt.legend(loc = 'best')
                        plt.ylabel('Number of Pixels')
                        plt.subplot(212)
                        plt.semilogy(base[:-1],tlen-cumulative,linewidth = 3)
                        plt.xlim(max(base[:-1]),min(base[:-1])) 
                        #plt.axvline(cutoff,color = 'red')
                        plt.xlabel('kJy/sr')
                        plt.ylabel('Number of Pixels')
                        plt.savefig('cumulativetest_bsub_'+skey+'.png')
                        plt.close()
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
                    
