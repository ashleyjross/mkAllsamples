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
l = np.zeros((f.size))
a = f[f['isdupl']==False]
print a.size
print sum(a['hasfiber']), sum(a['Z_reliable'])
cpw,cp = fiber_collision(a)
ncpw = 0
ncp = 0
ncpwdup = 0
ndup = 0
for i in range(0,cpw.size):
	if cp[i] == 3:
		ncp += 1.
	if cp[i] == 1:
		ncpw += cpw[i]	
	if a[i]['isdupl'] == True:
		ncpwdup += 	cpw[i]
		ndup += 1.
	if cpw[i] > 5:
		plt.plot(a['ra'][cp==3],a['dec'][cp==3],'ro')
		plt.plot(a['ra'][cp==1],a['dec'][cp==1],'ko')
		plt.xlim(a[i]['ra']-1/30.,	a[i]['ra']+1/30.)
		plt.ylim(	a[i]['dec']-1/30.,	a[i]['dec']+1/30.)
		plt.show()
		
# print min(cpw),max(cpw),min(cp),max(cp),f.size,a.size,ncpw,ncp/float(a.size),ncpw/float(a.size),ncpwdup,ndup
