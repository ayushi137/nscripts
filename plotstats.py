from numpy import *
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import leastsq

def m0am(p,am):
	m,b = p
	return m*am + b

def residuals(p,am,m0):
	return abs(m0-m0am(p,am))

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
slopes = stats[:,10].astype(float)

order = times.argsort()
objects = objects[order]
times = times[order]
serials = serials[order]
colours = colours[order]
expns = expns[order]
altitudes = altitudes[order]
airmasses = airmasses[order]
m0s = m0s[order]
fwhms = fwhms[order]
dates = dates[order]
slopes = slopes[order]

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
figinds['dra'] = 6

reds = ['r','lightcoral','darkred','darkorange','indianred']
greens = ['g','lime','darkgreen','darkseagreen','greenyellow']

for key in names.keys():
	if key == 'dra':
		daydiv = 0.1
	if key == 'spi':
		daydiv = 0.05
	serial = serials[inds[key]]
	colour = colours[inds[key]]
	m0 = m0s[inds[key]]
	am = airmasses[inds[key]]
	fwhm = fwhms[inds[key]]
	time = times[inds[key]]
	time -= np.min(time)
	date = dates[inds[key]]
	slope = slopes[inds[key]]
	diff = abs(roll(time,1)-time)
	transition = where(diff > 0.1)
	transition = delete(transition, 0)
	daysplits = []
	for t in transition:
		indstochange = arange(t,len(time))
		for j in indstochange:
			time[j] -= diff[t]
			time[j] += daydiv
		try:
			daysplits.append(0.5*(time[t]-time[t-1])+time[t-1])
		except IndexError:
			continue
	if outlierremove:
		good = where((fwhm < 30) & (m0 > 26.5))# & (slope > 0.5) & (slope < 4))
		serial = serial[good]
		colour = colour[good]
		m0 = m0[good]
		am = am[good]
		fwhm = fwhm[good]
		time = time[good]
		slope = slope[good]
		date = date[good]
	cams = unique(serial)
	cr = 0
	cg = 0
	plt.figure(1+figinds[key],figsize = (10,8))
	plt.title(names[key])
	plt.xlabel('Time')
	plt.ylabel('$m_0$',fontsize = 20)
	plt.figure(2+figinds[key],figsize = (10,8))
	plt.title(names[key])
	plt.xlabel('Time')
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
	plt.figure(6+figinds[key],figsize = (10,8))
	plt.title(names[key])	
	plt.xlabel('Time')
	plt.ylabel('Slope')
	for cam in cams:
		i = where(cam == serial)
		c = colour[i][0]
		if c == 'SloanG':
			color = greens[cg]
			cg += 1
		elif c == 'SloanR':
			color = reds[cr]
			cr += 1
		plt.figure(1+figinds[key])
		plt.plot(time[i],m0[i],'o',color = color,markersize = 10,label = cam)
		plt.figure(2+figinds[key])
		plt.plot(time[i],fwhm[i],'o',color = color,markersize = 10,label = cam)
		plt.figure(3+figinds[key])
		plt.plot(fwhm[i],m0[i],'o',color = color,markersize = 10,label = cam)
		plt.figure(4+figinds[key])
		plt.plot(am[i],fwhm[i],'o',color = color,markersize = 10,label = cam)
		plt.figure(5+figinds[key])
		plt.plot(am[i],m0[i],'o',color = color,markersize = 10,label = cam)
		plt.figure(6+figinds[key])
		plt.plot(time[i],slope[i],'o',color=color,markersize = 10,label = cam)
	for day in daysplits:
		plt.figure(1+figinds[key])
		plt.axvline(day,color='k',linewidth = 4)
		plt.figure(2+figinds[key])
		plt.axvline(day,color='k',linewidth = 4)
		plt.figure(6+figinds[key])
		plt.axvline(day,color='k',linewidth = 4)
	plt.figure(1+figinds[key])
	plt.legend(loc = 'best',fontsize = 10)
	if not outlierremove:
		plt.savefig('stats/{0}_m0_time.png'.format(names[key]))
	if outlierremove:
		plt.savefig('stats/no_out_{0}_m0_time.png'.format(names[key]))
	plt.figure(2+figinds[key])
	plt.legend(loc = 'best',fontsize = 10)
	if not outlierremove:
		plt.savefig('stats/{0}_FWHM_time.png'.format(names[key]))
	if outlierremove:
		plt.savefig('stats/no_out_{0}_FWHM_time.png'.format(names[key]))
	plt.figure(3+figinds[key])
	plt.legend(loc = 'best',fontsize = 10)
	if not outlierremove:
		plt.savefig('stats/{0}_FWHM_m0.png'.format(names[key]))
	if outlierremove:
		plt.savefig('stats/no_out_{0}_FWHM_m0.png'.format(names[key]))
	plt.figure(4+figinds[key])
	plt.legend(loc = 'best',fontsize = 10)
	if not outlierremove:
		plt.savefig('stats/{0}_airmass_FWHM.png'.format(names[key]))
	if outlierremove:
		plt.savefig('stats/no_out_{0}_airmass_FWHM.png'.format(names[key]))
	plt.figure(5+figinds[key])
	plt.legend(loc = 'best',fontsize = 10)
	if not outlierremove:
		plt.savefig('stats/{0}_airmass_m0.png'.format(names[key]))
	if outlierremove:
		plt.savefig('stats/no_out_{0}_airmass_m0.png'.format(names[key]))
	plt.figure(6+figinds[key])
	plt.legend(loc = 'best',fontsize = 10)
	if not outlierremove:
		plt.savefig('stats/{0}_slope_time.png'.format(names[key]))
	if outlierremove:
		plt.savefig('stats/no_out_{0}_slope_time.png'.format(names[key]))
	plt.close('all')
	for d in date:
		greenm0s = m0[where((colour == 'SloanG') & (date == d))]
		greenams = am[where((colour == 'SloanG') & (date == d))]
		redm0s = m0[where((colour == 'SloanR') & (date == d))]
		redams = am[where((colour == 'SloanR') & (date == d))]
		p0green = [-0.14,27.5]
		p0red = [-0.11,27.3]
		if len(greenm0s) != 1 and len(redm0s) != 1:
			pgreen = leastsq(residuals,p0green,args = (greenams,greenm0s))
			pred = leastsq(residuals,p0red,args = (redams,redm0s))
			amls = arange(np.min(am),np.max(am),0.01)
			plt.figure()
			plt.plot(greenams,greenm0s,'go',markersize = 10,label = 'k = {0}'.format(pgreen[0][0]))
			plt.plot(amls,m0am(pgreen[0],amls),'k',linewidth = 3)
			plt.plot(redams,redm0s,'ro',markersize = 10,label = 'k = {0}'.format(pred[0][0]))
			plt.plot(amls,m0am(pred[0],amls),'k',linewidth = 3)
			plt.xlabel('Airmass')
			plt.ylabel('$m_0$',fontsize = 20)
			plt.title(names[key])
			plt.legend(loc = 'best')
			plt.savefig('stats/{0}_{1}_extinction.png'.format(names[key],d))
			plt.close()
			print names[key]
			print 'Green ',pgreen[0]
			print 'Red ',pred[0]
