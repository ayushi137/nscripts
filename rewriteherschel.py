from astropy.io import fits
import docopt
import os

"""
To reindex Herschel files, making them more Python friendly.
This is currently unique to Draco and Spider clouds only.

Usage:
rewriteherschel [-h] [-d DIRECTORY]

Options:
    -h, --help                      Show this screen

    -d DIRECTORY, --dir DIRECTORY   Choose directory with Herschel files     

"""

arguments = docopt.docopt(__doc__)

herdir = arguments['--dir']

# Index to access image files in Herschel
ind = 'IMAGE'
# List all Herschel files
files = os.listdir(herdir)

# Create dictionary to convert raw Herschel units to MJy/sr
K_PtoE = {} # from http://herschel.esac.esa.int/Docs/SPIRE/html/spire_om.html
K_PtoE['PSW'] = 91.289
K_PtoE['PMW'] = 51.799
K_PtoE['PLW'] = 24.039

codewavelength = {}
codewavelength[250] = 'PSW'
codewavelength[350] = 'PMW'
codewavelength[500] = 'PLW'

# Cycle through all files
for f in files:
    # If file is not already reindexed
    if 'reindex_' not in f:
        # Read in data from file
        h = fits.open(herdir+f)
        hdata = h[ind].data
        hhder = h[ind].header
        h.close()
        # Find the wavelength being used
        multiple = K_PtoE[codewavelength[hhder['wavelength']]]
        newf = 'reindex_'+f
        hdata = K_PtoE[multiple]*hdata
        hhder['UNITS'] = 'MJy/sr'
        fits.writeto(herdir+newf, hdata, hhder,clobber = True)
    
 
