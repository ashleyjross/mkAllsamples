'''
Take Anand's files and make simplified ELG catalogs
'''

dir = '/Users/ashleyross/fitsfiles/'
from astropy.io import fits
import numpy as np
import fitsio #reading in entire file with fitsio is much faster than with astropy


def mkgalELGsimp(reg='SGC',v='v5_10_7',compl=0.5,zmin=0,zmax=1.5,vo='test'):
#zmin=.6,zmax=1.1,samp='21',v='v5_10_7',c='sci',app='.fits',compl=0,compls=0,fkp='fkp',wm='wstar',zm=''):
	#from healpix import healpix, radec2thphi
	#if c == 'sci': #AJR uses this define directory for machine he uses
	#	dir = dirsci
	if reg == 'SGC':
		chunkl = ['eboss21','eboss22']			
	if reg == 'NGC':
		chunkl = ['eboss23','eboss25']			

	ffkp = np.loadtxt('nbarELG'+reg+v+'.dat').transpose()
	ral = []
	decl = []
	zl = []
	fkpweightl = []
	sysweightl = []
	#for j in range(0,len(chunkl)):
		#f = fits.open(dir+'eboss'+chunkl[j]+'.'+v+'.latest.fits')[1].data
		#f = fitsio.read(dir+'eboss'+chunkl[j]+'.'+v+'.latest.fits')
	f = fitsio.read(dir+'ELG.'+v+'.latest.fits')
	for i in range(0,len(f)):
		#z = f.Z[i]
		z = f[i]['Z']
		w =1.
		m = 0
	#if f[i]['dec'] < .5 and f[i]['ra'] > 350 and f[i]['ra'] < 355:
	#	m = 1
		if f[i]['sector_TSR'] < compl:
			m = 1
		kz = 0
		if f[i]['Z_reliable'] == True:
			kz = 1
		kc = 0
		if f[i]['chunk'] == chunkl[0] or f[i]['chunk'] == chunkl[1]:
			kc = 1		
		if z > zmin and z < zmax and m == 0 and kz == 1 and f[i]['isdupl'] == False and kc ==1 and f[i]['hasfiber'] == True:# and f[i]['depth_ivar_r'] > 100 and f[i]['depth_ivar_g'] > 300 and f[i]['depth_ivar_z'] > 50:
			ra,dec = f[i]['ra'],f[i]['dec']
			ral.append(ra)
			decl.append(dec)
			zl.append(z)
			zind = int(z/.01)
			fkpw = ffkp[-1][zind]
			fkpweightl.append(fkpw)
			sysweightl.append(1.) #none of these yet

	print len(ral)
	ral = np.array(ral)
	decl = np.array(decl)
	zl = np.array(zl)
	fkpweightl = np.array(fkpweightl)
	sysweightl = np.array(sysweightl)

	cols = []
	RAc = fits.Column(name='RA',format='D', array=ral)
	cols.append(RAc)
	DECc = fits.Column(name='DEC',format='D', array=decl)
	cols.append(DECc)
	Zc = fits.Column(name='Z',format='D', array=zl)
	cols.append(Zc)
	fkpc = fits.Column(name='WEIGHT_FKP',format='D', array=fkpweightl)
	cols.append(fkpc)
	sysc = fits.Column(name='WEIGHT_SYSTOT',format='D', array=sysweightl)
	cols.append(sysc)
	
	hdulist = fits.BinTableHDU.from_columns(cols)
	header = hdulist.header
	hdulist.writeto(dir+'ELG'+reg+vo+'.dat.fits', overwrite=True)
	
	return True		

def mkranELGsimp(reg='SGC',v='v5_10_7',compl=0.5,vo='test'):
	from random import random
	if reg == 'SGC':
		chunkl = ['eboss21','eboss22']			
	if reg == 'NGC':
		chunkl = ['eboss23','eboss25']			
	ffkp = np.loadtxt('nbarELG'+reg+v+'.dat').transpose()
	dat = fitsio.read(dir+'ELG'+reg+vo+'.dat.fits')
	print len(dat['Z'])
	ral = []
	decl = []
	zl = []
	fkpweightl = []
	sysweightl = []
	#for j in range(0,len(chunkl)):
	#f = fits.open(dir+'eboss'+chunkl[j]+'.'+v+'.latest.fits')[1].data
	#f = fitsio.read(dir+'eboss'+chunkl[j]+'.'+v+'.latest.rands.fits')
	f = fitsio.read(dir+'ELG.'+v+'.latest.rands.fits')
	print len(f), compl
	for i in range(0,len(f)):
		#z = f.Z[i]
		kc = 0
		if f[i]['chunk'] == chunkl[0] or f[i]['chunk'] == chunkl[1]:
			kc = 1		

		if f[i]['sector_TSR'] >= compl and kc == 1:
			indz = int(random()*len(dat['Z']))
			z = dat[indz]['Z']
			fkpw = dat[indz]['WEIGHT_FKP']
			wsys = dat[indz]['WEIGHT_SYSTOT']*f[i]['sector_TSR']#*f[i]['plate_SSR']
			ra,dec = f[i]['ra'],f[i]['dec']
			ral.append(ra)
			decl.append(dec)
			zl.append(z)
			fkpweightl.append(fkpw)
			sysweightl.append(wsys) 
	print sum(sysweightl)/10000.
	print len(ral)
	ral = np.array(ral)
	decl = np.array(decl)
	zl = np.array(zl)
	fkpweightl = np.array(fkpweightl)
	sysweightl = np.array(sysweightl)

	cols = []
	RAc = fits.Column(name='RA',format='D', array=ral)
	cols.append(RAc)
	DECc = fits.Column(name='DEC',format='D', array=decl)
	cols.append(DECc)
	Zc = fits.Column(name='Z',format='D', array=zl)
	cols.append(Zc)
	fkpc = fits.Column(name='WEIGHT_FKP',format='D', array=fkpweightl)
	cols.append(fkpc)
	sysc = fits.Column(name='WEIGHT_SYSTOT',format='D', array=sysweightl)
	cols.append(sysc)
	
	hdulist = fits.BinTableHDU.from_columns(cols)
	header = hdulist.header
	hdulist.writeto(dir+'ELG'+reg+vo+'.ran.fits', overwrite=True)

	return True
