from astropy.io import fits
import os

ind = 'IMAGE'
herdir = '../herschel/PGM_1_2/'
files = os.listdir(herdir)

K_PtoE = {} # from http://herschel.esac.esa.int/Docs/SPIRE/html/spire_om.html
K_PtoE['PSW'] = 91.289
K_PtoE['PMW'] = 51.799
K_PtoE['PLW'] = 24.039

for f in files:
    if 'reindex_' not in f:
        h = fits.open(herdir+f)
        hdata = h[ind].data
        hhder = h[ind].header
        h.close()
        if 'spi1_1' in herdir:
            multiple = f.split('_')[-3]
        elif 'PGM_1_2' in herdir:
            multiple = f.split('_')[-2]
        newf = 'reindex_'+f
        hdata = K_PtoE[multiple]*hdata
        hhder['UNITS'] = 'MJy/sr'
        fits.writeto(herdir+newf, hdata, hhder,clobber = True)
    
 
