from numpy import *
import matplotlib.pyplot as plt
plt.ion()

outlierremove = True

stats = loadtxt('stats/stats.txt',dtype = 'str')

objects = stats[:,0]
times = stats[:,1].astype(float)
serials = stats[:,2]
colours = stats[:,3]
expns = stats[:,4]
altitudes = stats[:,5].astype(float)
airmasses = stats[:,6].astype(float)
m0s = stats[:,7].astype(float)
fwhms = stats[:,8].astype(float)
dates = stats[:,9]

spi = where(objects == 'spi1_1')
dra = where(objects == 'PGM_1_2')

inds = {}
inds['spi'] = spi
inds['dra'] = dra

names = {}
names['spi'] = 'Spider'
names['dra'] = 'Draco'

figinds = {}
figinds['spi'] = 0
figinds['dra'] = 5

reds = ['r','lightcoral','darkred','darkorange','indianred']
greens = ['g','lime','darkgreen','darkseagreen','greenyellow']

for key in names.keys():
	serial = serials[inds[key]]
	colour = colours[inds[key]]
	m0 = m0s[inds[key]]
	am = airmasses[inds[key]]
	fwhm = fwhms[inds[key]]
	time = times[inds[key]]
	date = dates[inds[key]]
	if outlierremove:
		good = where((fwhm < 30) & (m0 > 26.5))
		serial = serial[good]
		colour = colour[good]
		m0 = m0[good]
		am = am[good]
		fwhm = fwhm[good]
		time = time[good]
	cams = unique(serial)
	cr = 0
	cg = 0
	plt.figure(1+figinds[key],figsize = (10,8))
	plt.title(names[key])
	plt.xlabel('Time [Exposure Number]')
	plt.ylabel('$m_0$',fontsize = 20)
	plt.figure(2+figinds[key],figsize = (10,8))
	plt.title(names[key])
	plt.xlabel('Time [Exposure Number]')
	plt.ylabel('FWHM ["]')
	plt.figure(3+figinds[key],figsize = (10,8))
	plt.title(names[key])
	plt.xlabel('FWHM ["]')
	plt.ylabel('$m_0$',fontsize = 20)	
	plt.figure(4+figinds[key],figsize = (10,8))
	plt.title(names[key])	
	plt.xlabel('Airmass')
	plt.ylabel('FWHM ["]')	
	plt.figure(5+figinds[key],figsize = (10,8))
	plt.title(names[key])	
	plt.xlabel('Airmass')
	plt.ylabel('$m_0$',fontsize = 20)	
	for cam in cams:
		i = where(cam == serial)
		c = colour[i][0]
		if c == 'SloanG':
			color = greens[cg]
			cg += 1
		elif c == 'SloanR':
			color = reds[cr]
			cr += 1
		order = time[i].argsort()
		plt.figure(1+figinds[key])
		plt.plot(m0[i][order],'o',color = color,markersize = 10,label = cam)
		plt.figure(2+figinds[key])
		plt.plot(fwhm[i][order],'o',color = color,markersize = 10,label = cam)
		plt.figure(3+figinds[key])
		plt.plot(fwhm[i],m0[i],'o',color = color,markersize = 10,label = cam)
		plt.figure(4+figinds[key])
		plt.plot(am[i],fwhm[i],'o',color = color,markersize = 10,label = cam)
		plt.figure(5+figinds[key])
		plt.plot(am[i],m0[i],'o',color = color,markersize = 10,label = cam)
	plt.figure(1+figinds[key])
	plt.legend(loc = 'best')
	if not outlierremove:
		plt.savefig('stats/{0}_m0_time.png'.format(names[key]))
	if outlierremove:
		plt.savefig('stats/no_out_{0}_m0_time.png'.format(names[key]))
	plt.figure(2+figinds[key])
	plt.legend(loc = 'best')
	if not outlierremove:
		plt.savefig('stats/{0}_FWHM_time.png'.format(names[key]))
	if outlierremove:
		plt.savefig('stats/no_out_{0}_FWHM_time.png'.format(names[key]))
	plt.figure(3+figinds[key])
	plt.legend(loc = 'best')
	if not outlierremove:
		plt.savefig('stats/{0}_FWHM_m0.png'.format(names[key]))
	if outlierremove:
		plt.savefig('stats/no_out_{0}_FWHM_m0.png'.format(names[key]))
	plt.figure(4+figinds[key])
	plt.legend(loc = 'best')
	if not outlierremove:
		plt.savefig('stats/{0}_airmass_FWHM.png'.format(names[key]))
	if outlierremove:
		plt.savefig('stats/no_out_{0}_airmass_FWHM.png'.format(names[key]))
	plt.figure(5+figinds[key])
	plt.legend(loc = 'best')
	if not outlierremove:
		plt.savefig('stats/{0}_airmass_m0.png'.format(names[key]))
	if outlierremove:
		plt.savefig('stats/no_out_{0}_airmass_m0.png'.format(names[key]))
