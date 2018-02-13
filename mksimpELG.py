'''
Take Anand's files and make simplified ELG catalogs
'''

dir = '/Users/ashleyross/fitsfiles/'
from astropy.io import fits
import numpy as np
import fitsio #reading in entire file with fitsio is much faster than with astropy
from cattools import fiber_collision

def mkgalELG_specfull(reg='SGC',v='v5_10_7',vo='test'):
	'''
	Split into regions
	add collider information
	remove duplicates
	classify objects
	'''

	if reg == 'SGC':
		chunkl = ['eboss21','eboss22']			
	if reg == 'NGC':
		chunkl = ['eboss23','eboss25']			
	f = fitsio.read(dir+'ELG.'+v+'.latest.fits')
	#a = f[f['isdupl']==False] #remove duplicates
	w = (f['chunk'] == chunkl[0] ) | (f['chunk'] == chunkl[1] ) &   (f['isdupl']==False) #select wanted chunks and remove duplicates
	freg = f[w] #apply selection
	#columns for files
	ral = freg['ra'] #RA
	decl = freg['dec'] #DEC
	zl = freg['Z'] #Z, redshift
	sectorl = freg['sector'] #SECTOR
	platel = freg['PLATE'] #PLATE
	idl = freg['decals_objid'] #DECALS_OBJID
	fl = np.zeros(freg.size, dtype=int) #OBS_FLAG, flag for observation: 0 missed, 1 successful eBOSS redshift, 3 close pair, 4 star, 7 failure

	cpw,cp,cpind,comptot,compboss = fiber_collision(freg) #get collision information, sector stats
	#cpw contains traditional close pair weight
	#cp contains flag for successful redshift (1), missed (0), close pair (3)
	#cpind is id of colliding galaxy
	#comptot should be n_fib/n_targ and should match Anand's calculation
	#compboss should match
	
	for i in range(0,freg.size):
		of = cp[i]
		if freg[i]['Z_reliable'] == False:
			of = 7
		if freg[i]['CLASS'] == 'STAR': #I think this order is the one we want
			of = 4
		fl[i] = of
	

	cols = []
	RAc = fits.Column(name='RA',format='D', array=ral)
	cols.append(RAc)
	DECc = fits.Column(name='DEC',format='D', array=decl)
	cols.append(DECc)
	Zc = fits.Column(name='Z',format='D', array=zl)
	cols.append(Zc)
	wcpc = fits.Column(name='WEIGHT_CP',format='D', array=cpw)
	cols.append(wcpc)	
	wtsr = fits.Column(name='COMP_TSR',format='D', array=comptot)
	cols.append(wtsr)	
	wcb = fits.Column(name='COMP_BOSS',format='D', array=compboss)
	cols.append(wcb)	
	cc = fits.Column(name='CHUNK',format='A10',array=freg['chunk'])
	cols.append(cc)

	sc = fits.Column(name='SECTOR',format='K',array=sectorl)
	cols.append(sc)
	pc = fits.Column(name='PLATE',format='K',array=platel)
	cols.append(pc)
	idc = fits.Column(name='DECALS_OBJID',format='K',array=idl)
	cols.append(idc)
	cidc = fits.Column(name='COLLIDER_DECALS_OBJID',format='K',array=cpind)
	cols.append(cidc)

	fc = fits.Column(name='TYPE_FLAG',format='K',array=fl)
	cols.append(fc)
	#don't change the name of these columns
	morecols = ['plate_rSN2','XFOCAL','YFOCAL','FIBERID','galdepth_g','galdepth_r','galdepth_z','ebv','psfsize_g','psfsize_r','psfsize_z']
	for name in morecols:
		form = 'D'
		if name == 'FIBERID':
			form = 'K'
		c = fits.Column(name=name,format=form,array=freg[name])
		cols.append(c)
	hdulist = fits.BinTableHDU.from_columns(cols)
	header = hdulist.header
	hdulist.writeto(dir+'ELG'+reg+vo+'full.dat.fits', overwrite=True)
	
	return True		


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
	gdepthl = [] #to maybe be used for weights 
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
			gdepthl.append(f[i]['galdepth_g'])

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
	depthc = fits.Column(name='galdepth_g',format='D', array=gdepthl)
	cols.append(depthc)
	
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

def mkranELG_zgdepth(reg='SGC',v='v5_10_7',compl=0.5,vo='test',sub=1):
	#set sub to greater integer values if you want to test things
	from random import random
	if reg == 'SGC':
		chunkl = ['eboss21','eboss22']			
	if reg == 'NGC':
		chunkl = ['eboss23','eboss25']			
	ffkp = np.loadtxt('nbarELG'+reg+v+'.dat').transpose()
	dat = fitsio.read(dir+'ELG'+reg+vo+'.dat.fits')
	dat.sort(order='galdepth_g')
	depthl  = [] #subsample by factor of 100 to make smaller list to get percentiles
	for i in range(0,dat.size,100):
		depthl.append(dat[i]['galdepth_g'])
	depthl = np.array(depthl)
	#depthl = dat['galdepth_g']#.sort()
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
	ndat = dat.size
	for i in range(0,len(f),sub):
		#z = f.Z[i]
		kc = 0
		if f[i]['chunk'] == chunkl[0] or f[i]['chunk'] == chunkl[1]:
			kc = 1		

		if f[i]['sector_TSR'] >= compl and kc == 1:
			depth = f[i]['galdepth_g']
			#per = dat[dat['galdepth_g']>depth].size/float(ndat)
			per = depthl[depthl>depth].size/float(ndat)
			#perd = 90
			#peru = 100
			#if per > 0.05 and per < 0.95:
			#	perd = (per-.05)*100
			#	peru = (per+.05)*100
			#elif per <= 0.05:
			#	perd = 0
			#	peru = 10
			#dd = np.percentile(depthl,perd)
			#du = np.percentile(depthl,peru)
			#w = (dat['galdepth_g']>dd) & (dat['galdepth_g']<du)				
			#zdat = dat[w] 
			if per <= 0.05:
				per = 0.05
			if per >= 0.95:
				per = .95 #to 
			perd = (per-.05)
			peru = (per+.05)

			indu = int(peru*ndat)
			indd = int(perd*ndat)
			nd = indu-indd					
			indz = int(random()*nd)
			
			#zdat = dat[indd:indu]
			#z = zdat[indz]['Z']
			#fkpw = zdat[indz]['WEIGHT_FKP']
			#wsys = zdat[indz]['WEIGHT_SYSTOT']*f[i]['sector_TSR']#*f[i]['plate_SSR']
			while indd+indz >= ndat:
				print indu,indd,indz,per,perd,peru,i
				indz = indz -1
			z = dat[indd+indz]['Z']
			fkpw = dat[indd+indz]['WEIGHT_FKP']
			wsys = dat[indd+indz]['WEIGHT_SYSTOT']*f[i]['sector_TSR']#*f[i]['plate_SSR']
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
	hdulist.writeto(dir+'ELG'+reg+vo+'gdepth.ran.fits', overwrite=True)

	return True
