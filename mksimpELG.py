'''
Take Anand's files and make simplified ELG catalogs
'''

dir = '/Users/ashleyross/fitsfiles/'
from astropy.io import fits
import numpy as np
import fitsio #reading in entire file with fitsio is much faster than with astropy
from cattools import fiber_collision
from matplotlib import pyplot as plt
from math import *

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
	gdepthl = []
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
			gdepthl.append(f[i]['galdepth_g'])			
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
	dc = fits.Column(name='galdepth_g',format='D', array=gdepthl)
	cols.append(dc)
	
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
	ndep = 0
	for i in range(0,dat.size,100):
		depthl.append(dat[i]['galdepth_g'])
		ndep += 1.
	depthl = np.array(depthl)
	gdepthl = [] #list to write back out
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
			per = 1.-depthl[depthl>depth].size/ndep
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
			gdepthl.append(f[i]['galdepth_g'])	
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
	dc = fits.Column(name='galdepth_g',format='D', array=gdepthl)
	cols.append(dc)
	
	hdulist = fits.BinTableHDU.from_columns(cols)
	header = hdulist.header
	hdulist.writeto(dir+'ELG'+reg+vo+'gdepth.ran.fits', overwrite=True)

	return True

def testnz(reg,ver='test',zmin=0,zmax=1.5):
	f = fitsio.read(dir+'ELG'+reg+ver+'gdepth.ran.fits')
	fd = fitsio.read(dir+'ELG'+reg+ver+'.dat.fits')
	print np.percentile(f['galdepth_g'],50),np.percentile(fd['galdepth_g'],50)
	spl = np.percentile(fd['galdepth_g'],50)
	nzl = []
	nzlr = []
	nzlh = []
	nzlhr = []

	zl = []
	ni = int(100*(zmax-zmin))
	for i in range(0,ni):
		zl.append(zmin+.005+.01*i)
		nzl.append(0)
		nzlr.append(0)
		nzlh.append(0)
		nzlhr.append(0)

	for i in range(0,f.size):
		zind = int((f[i]['Z']-zmin)*100)
		if f[i]['galdepth_g'] > spl:
			nzlhr[zind] += 1.
		else:
			nzlr[zind] += 1.
	for i in range(0,fd.size):
		zind = int((fd[i]['Z']-zmin)*100)
		if fd[i]['galdepth_g'] > spl:
			nzlh[zind] += 1.
		else:
			nzl[zind] += 1.
	plt.plot(zl,np.array(nzl)/sum(nzl),	zl,np.array(nzlr)/sum(nzlr))
	plt.show()
	plt.plot(zl,np.array(nzlh)/sum(nzlh),zl,	np.array(nzlhr)/sum(nzlhr))
	plt.show()
	plt.plot(zl,np.array(nzlr)/sum(nzlr),zl,	np.array(nzlhr)/sum(nzlhr))
	plt.show()

	return True
	
			
def ngvdepth(reg,ver,zmin,zmax,sysmin=0,sysmax=10000,sys='galdepth_g'):
	#sample is the sample being used, e.g. 'lrg'
	stl = []
	wstl = []
	errl = []
	binng = []
	binnr = []
	binnrd = []
	nsysbin = 10 #10 bins are being used
	for i in range(0,nsysbin):
		binng.append(0)
		binnr.append(0)
		binnrd.append(0)

	
	nr = 0
	nrd = 0
	sysm = float(nsysbin)/(sysmax-sysmin)
	bs = 0
	bsr = 0
	bsrd = 0

	f = fitsio.read(dir+'ELG'+reg+ver+'.ran.fits') #read rands, no depth adjustment
	for i in range (0,len(f)):
		if f[i]['Z'] > zmin and f[i]['Z'] < zmax:
			w = f[i]['WEIGHT_SYSTOT']*f[i]['WEIGHT_FKP']
			sysv = f[i]['galdepth_g']
			bins = int((sysv-sysmin)*sysm)
			if bins >= 0 and bins < nsysbin:
				binnr[bins] += w
			else:
				bsr += w

			nr += w
	print min(f[sys]),max(f[sys])
	print nr

	f = fitsio.read(dir+'ELG'+reg+ver+'gdepth.ran.fits') #read rands, with depth adjustment
	for i in range (0,len(f)):
		if f[i]['Z'] > zmin and f[i]['Z'] < zmax:
			w = f[i]['WEIGHT_SYSTOT']*f[i]['WEIGHT_FKP']
			sysv = f[i]['galdepth_g']
			bins = int((sysv-sysmin)*sysm)
			if bins >= 0 and bins < nsysbin:
				binnrd[bins] += w
			else:
				bsrd += w

			nrd += w
	print min(f[sys]),max(f[sys])
	print nrd


	f = fitsio.read(dir+'ELG'+reg+ver+'.dat.fits') #read galaxy file
	
	no = 0
	zm = 0
	nt = 0
	avebs = 0
	sysstr = sys.split('_')
	
	for i in range (0,len(f)):
		z = f[i]['Z']
		if z > zmin and z < zmax:
			no += 1
			w = f[i]['WEIGHT_SYSTOT']*f[i]['WEIGHT_FKP']
			sysv = f[i]['galdepth_g']
			bins = int((sysv-sysmin)*sysm)
			
					
			if bins >= 0 and bins < nsysbin:
				binng[bins] += 1.*w
			else:
				bs += w #count numbers outside of sysmin/sysmax
				avebs += w*sysv
			zm += w*z
			nt += w
			
	if bs > 0:
		avebs = avebs/bs
	print avebs,sysmin		
	print 'total number, weighted number'
	print no,nt
	print 'mean redshift'
	print zm/nt

	print 'total number of randoms/objects '+str(nr)+'/'+str(nt)
	print 'number of randoms/objects outside tested range '+str(bsr)+'/'+str(bs)			
	ave = nt/nr
	aved = nt/nrd
	print 'average number of objects per random is '+ str(ave)+' '+str(aved)
	fs = open('n'+'gELG'+reg+'_'+ver+'_mz'+str(zmin)+'xz'+str(zmax)+'v'+sys+'.dat','w')
	xl = []
	yl = []
	el = []
	yld = []
	eld = []
	for i in range(0,nsysbin):
		sysv = sysmin + 1./(2.*sysm) + i/sysm
		if binnr[i] > 0:
			ns = binng[i]/binnr[i]/ave
			nsd = binng[i]/binnrd[i]/aved
			nse = sqrt(binng[i]/(binnr[i])**2./(ave)**2.+(binng[i]/ave)**2./(binnr[i])**3.) #calculate poisson error
			nsed = sqrt(binng[i]/(binnrd[i])**2./(aved)**2.+(binng[i]/aved)**2./(binnrd[i])**3.) #calculate poisson error
		else:
			ns = 1. #write out 1.0 1.0 if no pixels at given value of sys
			nse = 1.	
			nsd = 1.
			nsed = 1.	
		fs.write(str(sysv)+' '+str(ns)+' '+str(nse)+' '+str(nsd)+' '+str(nsed)+'\n')
		xl.append(sysv)
		yl.append(nsd)
		el.append(nsed)
		yld.append(ns)
		eld.append(nse)

	fs.close()
	chin = sum((np.array(yl)-1.)**2./np.array(el)**2.)
	print chin
	plt.errorbar(xl,yl,el,fmt='ko')
	plt.errorbar(xl,yld,eld,fmt='^r')
	plt.show()	
	return True
