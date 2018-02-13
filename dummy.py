'''
way to test different things
'''

import argparse
from astropy.io import fits
import numpy as np
import fitsio
from matplotlib import pyplot as plt

''' 
test writing to fits file


l1 = np.ones((10))
l2 = np.zeros((10))
col = fits.Column(name='ones',format='D', array=l1)
col2 = fits.Column(name='zeros',format='D', array=l2)
hdulist = fits.BinTableHDU.from_columns([col,col2])
header = hdulist.header
hdulist.writeto('test.fits', overwrite=True)
'''

''' 
test close pair code
'''
from cattools import *
f = fitsio.read('/Users/ashleyross/fitsfiles/ELG.v5_10_7.latest.fits')
chunkl = ['eboss21','eboss22']
w = ((f['chunk'] == chunkl[0] ) | (f['chunk'] == chunkl[1] )) &   (f['isdupl']==False)
print f[w].size,f.size
# l = np.zeros((f.size))
a = f[w]
cc = fits.Column(name='CHUNK',format='A10',array=a['chunk'])
cols = [cc]
hdulist = fits.BinTableHDU.from_columns(cols)
header = hdulist.header
hdulist.writeto('test.fits', overwrite=True)
f = fitsio.read('test.fits')
print np.unique(f['CHUNK'])
# ws = (a['isdupl']==True)# & (a['hasfiber'] == 1)
# print a[ws].size
# #print sum(a['hasfiber']), sum(a['Z_reliable']), a[a['sector_TSR']==0].size
# #plt.plot(a[a['sector_TSR']==0]['ra'],a[a['sector_TSR']==0]['dec'],'k,')
# #plt.show()
# cpw,cp,cp_match,compreal,compboss = fiber_collision(a)
# 
# unique_sectors = np.unique(a['sector'])
# mc = []
# ac = []
# emc = []
# eac = []
# sl = []
# w2 = (a['sector_TSR'] != 0)
# for sector in unique_sectors:
# 	wi = (a['sector']==sector) & w2
# 	mcs = np.mean(a['sector_TSR'][wi])
# 	emcs = np.std(a['sector_TSR'][wi])#/sqrt(a['sector_TSR'][wi].size)
# 	acs = np.mean(compreal[wi])
# 	eacs = np.std(compreal[wi])#/sqrt(compreal[wi].size)
# 	mc.append(mcs)
# 	ac.append(acs)
# 	emc.append(emcs)
# 	eac.append(eacs)
# 	if emcs > .01:
# 		print sector,emcs
# print min(emc),max(emc),min(eac),max(eac)	
#plt.errorbar(ac,mc,emcs,fmt='ko')
#xl = [0,1]
#plt.plot(xl,xl,'r-')
#plt.show()

#ws = (compreal == 0)
#print np.unique(a['sector'][ws])
#plt.plot(a['ra'][ws],a['dec'][ws],'ko')
#plt.show()
#ws = (a['sector_TSR']==0) & (compreal != 0)
#print a[ws].size
#print np.unique(a['sector'])
# for i in range(0,a.size):
# 	if a[i]['sector_TSR'] == 0:
# 		sect = a[i]['sector']
# 		print sect
# 		wi = (a['sector']==sect)
# 		print np.mean(a['sector_TSR'][wi])

#plt.plot(compreal[w2],a['sector_TSR'][w2],'ko')
#xl = [0,1]
#plt.plot(xl,xl,'r-')
#plt.show()
# ncpw = 0
# ncp = 0
# ncpwdup = 0
# ndup = 0
# for i in range(0,cpw.size):
# 	if cp[i] == 3:
# 		ncp += 1.
# 	if cp[i] == 1:
# 		ncpw += cpw[i]	
# 	if a[i]['isdupl'] == True:
# 		ncpwdup += 	cpw[i]
# 		ndup += 1.
# 	if cpw[i] > 5:
# 		plt.plot(a['ra'][cp_match==i],a['dec'][cp_match==i],'ro')
# 		plt.plot(a['ra'][cp==1],a['dec'][cp==1],'ko')
# 		plt.xlim(a[i]['ra']-1/30.,	a[i]['ra']+1/30.)
# 		plt.ylim(	a[i]['dec']-1/30.,	a[i]['dec']+1/30.)
# 		plt.show()
		
#print min(cpw),max(cpw),min(cp),max(cp),f.size,a.size,ncpw,ncp/float(a.size),ncpw/float(a.size),ncpwdup,ndup
