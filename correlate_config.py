from astropy import units as u

# specify telescope details
scope_lat = {'Dragonfly':32.902836*u.deg}
scope_long = {'Dragonfly':-105.528350*u.deg}
scope_alt = {'Dragonfly':2214*u.m}

# specify herschel resolution at each wavelength
SPIRE = {'PSW':17.6,'PMW':23.9,'PLW':35.2}
spirekeys = SPIRE.keys()

config_data = {
	
	'scope':'Dragonfly', # telescope being used
	'objects':['/PGM_1_2/','/spi1_1/'], # names of cloud directories
	'raws':'/raw_lights/', # name of raw file directory
	'cals':'/calframes/', # name of calibration file directory
	'dsff':'/darksub_flatfield/', # name of calibrated file directory (will be created)
	'rwcorrfiles':False, # rewrites files to be correlated against (Herschel specific)
	'badframesdir':'/removedframes', # place to put bad files
	'cutoff_type':'median', # rule to apply to data to calculate cutoff (median/mean)
	'cutoff_mult':10, # multiplicative factor to apply to cutoff
	'cutoff':0 # if zero, multiplies cutoff_mult by cutoff_type(data) to calculate cutoff - otherwise hard cutoff			
}

# OUTPUT DIRECTORIES

# master dark location
damdir = '/dark_masters/'
# master flat location
flmdir = '/flat_masters/'
# dark subtracted files
dsmdir = '/dark_subtracted/'

# source extractor object maps
objdir = '/objects/'
# source extractor background maps
bakdir = '/background/'
# source extractor catalogue
catdir = '/catalogue/'

# location of helper plots from zero point magnitude calculations
magdir = '/magnitudecalcs/'
# location of photometered fits files
phodir = '/photometered/'
# location of regridded images
regdir = '/regrid/'
# location of background plane subtracted images
bsudir = '/backgroundsub/'
# location of correlation plots
cordir = '/correlations/'