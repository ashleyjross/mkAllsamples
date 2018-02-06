import os
import astropy.io.fits as fits
import sys
import matplotlib.pyplot as plt
import numpy as np
from pydl.pydlutils.yanny import yanny
from matplotlib.patches import Ellipse
import matplotlib
import subprocess
from matplotlib import gridspec
import datetime
import colorsys
import re
sys.path.append('/uufs/chpc.utah.edu/common/home/u0992342/Scripts/elgredshiftflag/trunk/python/')
import HandleReducedELGPlate
from matplotlib.backends.backend_pdf import PdfPages
from astropy.convolution import convolve, Box1DKernel
import speclite
import astropy.units as u
from astropy.coordinates import SkyCoord
import mangle
import healpy as hp



STILTSCMD='java -jar -Xmx4096M /uufs/chpc.utah.edu/common/home/u0992342/software/stilts.jar '
pid      = str(os.getpid())
TMPDIR   = os.getenv('MYSAS_DIR')+'/Catalogs/tmpdir/'

# appending fits tables (with the same format!)
# http://docs.astropy.org/en/stable/io/fits/usage/table.html#appending-tables
def append_fits(INFITSLIST,OUTFITS,COLNAMES=None):
	nfits = len(INFITSLIST)
	nrowslist = np.zeros(nfits,dtype='int')
	hdulist   = np.zeros(nfits,dtype='object')
	for i in xrange(nfits):
		hdulist[i]   = fits.open(INFITSLIST[i])
		nrowslist[i] = hdulist[i][1].data.shape[0]
	nrows = np.sum(nrowslist) # total nb of obj.
	print 'nrows = '+str(nrows)
	# if no COLNAMES is provided (COLNAMES==None), we take the columns from INFITSLIST[0]
	if (COLNAMES==None):
		COLNAMES = hdulist[0][1].columns.names
	# building the table structure: we take the columns (name,format) from INFITSLIST[0]
	collist = []
	for name in COLNAMES:
		print name
		ind = np.where(np.array(hdulist[0][1].columns.names)==name)[0][0]
		collist.append(fits.Column(name=name,format=hdulist[0][1].columns[ind].format))
	cols  = fits.ColDefs(collist)
	hdu   = fits.BinTableHDU.from_columns(cols, nrows=nrows)
	for i in xrange(nfits):
		indstart = np.sum(nrowslist[:i])
		indend   = indstart + nrowslist[i]
		print i, INFITSLIST[i], indstart, indend
		for colname in COLNAMES:
			hdu.data[colname][indstart:indend] = hdulist[i][1].data[colname]
	hdu.writeto(OUTFITS,clobber=True)
	return True


# SDSS tiling radius [degrees]
def get_tile_rad():
	return 1.49


# tiling directory
def get_tile_dir(CHUNK):
	# see Hee-Jong email "Re: eboss25 tiling is ready" [19/01/2018 20:43]
	return '/uufs/chpc.utah.edu/common/home/sdss05/software/svn.sdss.org/repo/eboss/ebosstilelist/trunk/outputs/'+CHUNK+'/'


# palette of unique colors
# http://stackoverflow.com/questions/470690/how-to-automatically-generate-n-distinct-colors
def get_colors(num_colors):
	colors=[]
	for i in np.arange(0., 360., 360. / num_colors):
		hue = i/360.
		lightness = (50 + np.random.rand() * 10)/100.
 		saturation = (90 + np.random.rand() * 10)/100.
		colors.append(colorsys.hls_to_rgb(hue, lightness, saturation))
	return colors



# function to read the tile properties
def read_tiles(parfile):
	a    = yanny(parfile)
	tile = np.array(a['STRUCT1']['tile']).astype('str')
	ra   = np.array(a['STRUCT1']['racen'])
	dec  = np.array(a['STRUCT1']['deccen'])
	return tile,ra,dec


# function to grab the footprint, and xlim,ylim for plots
def get_footprint(chunk):
	if (chunk=='eboss21'):
		footprint_ra  = np.array([317,360,360,317,317])
		footprint_dec = np.array([-2.0,-2.0,2.0,2.0,-2.0])
		xlim = [315,362]
		ylim = [-3,3]
	if (chunk=='eboss22'):
		footprint_ra  = np.array([0,45,45,0,0])
		footprint_dec = np.array([-5,-5,5,5,-5])
		xlim = [-1,46]
		ylim = [-6,6]
	if (chunk=='eboss23'):
		footprint_ra  = np.array([126.0, 136.5, 136.5, 157.0, 157.0, 142.5, 142.5, 126.0, 126.0])
		footprint_dec = np.array([16.0,  16.0,  13.8,  13.8,  27.0,  27.0,  29.0,  29.0,  16.0])
		xlim = [125,158]
		ylim = [12,31]
	if (chunk=='eboss25'):
		footprint_ra = np.array([131.0,166.0,166.0,157.0,157.0,142.5,142.5,131.0,131.0])
		footprint_dec= np.array([32.5, 32.5, 23.0, 23.0, 27.0, 27.0, 29.0, 29.0, 32.5 ])
		xlim = [130,167]
		ylim = [22,34]
	return footprint_ra, footprint_dec, xlim, ylim


# function to remove duplicates, if any
def rmv_dups(specfits,zquant,verbose=False):
	# inputs:
	# - specfits: fits catalog, with columns PLUG_RA, PLUG_DEC, ZWARNING, RCHI2, PLATE, MJD, FIBERI
	# - zquant: 'Z' or 'Z_NOQSO'
	# - verbose: if verbose=True, print info
	# output:
	# - keepind: np.array of True/False
	tmpfits = TMPDIR+'tmp.fits_'+pid
	# identifying duplicates
	tmpstr = (STILTSCMD+' tmatch1 action=identify '+
				'matcher=sky values="PLUG_RA PLUG_DEC" params=0.1 '+
				'in='+specfits+' ifmt=fits '+
				'out='+tmpfits+' ofmt=fits')
	print tmpstr
	subprocess.call(tmpstr, shell=True)
	# if no duplicates, no STILTS output catalogue...
	if not os.path.isfile(tmpfits):
		hdu     = fits.open(specfits)
		keepind = np.ones(len(hdu[1].data['Z']), dtype='bool')
	else:
		# keeping non duplicates and, for duplicates, the one with the best redshift
		hdu     = fits.open(tmpfits)
		groupid = hdu[1].data['GroupID']
		groupsize=hdu[1].data['GroupSize']
		if (zquant=='Z'):
			zwarning= hdu[1].data['ZWARNING']
			rchi2   = hdu[1].data['RCHI2']
		else:
			zwarning= hdu[1].data['ZWARNING_NOQSO']
			rchi2   = hdu[1].data['RCHI2_NOQSO']
		if (verbose==True):
			plate   = hdu[1].data['PLATE']
			mjd     = hdu[1].data['MJD']
			fiberid = hdu[1].data['FIBERID']
		## non-duplicates
		keepind = (groupsize<0) # empty -> -2147483648
		## for duplicates, the one with the best redshift
		groupid_dups = np.unique(groupid[(groupsize>=2)])
		for gid in groupid_dups:
			dupsind = np.arange(len(rchi2))[(groupid==gid)]                 # indexes of the dups from this group
			if (verbose==True):
				for ind in dupsind:
					print gid, ind, plate[ind], mjd[ind], fiberid[ind], zwarning[ind], rchi2[ind]
			dupsind   = dupsind[np.argsort(zwarning[dupsind])]              # sorting by increasing zwarning
			# if one has zwarning==0, we put it first
			if (zwarning[dupsind][0]==0):
				dupsind   = dupsind[(zwarning[dupsind]==zwarning[dupsind][0])]  # selecting dups with smallest zwarning
			dupsind   = dupsind[np.argsort(np.abs(rchi2[dupsind]-1.))]      # sorting by abs(rchi2-1)
			keepind[dupsind[0]] = True                                      # adding the dups with lowest abs(rchi2-1)
			if (verbose==True):
				print 'we keep:'
				print gid, ind, plate[dupsind[0]], mjd[dupsind[0]], fiberid[dupsind[0]], zwarning[dupsind[0]], rchi2[dupsind[0]]
				print ''
		## cleaning
		subprocess.call('rm '+tmpfits, shell=True)
	return keepind



def analysis_indiv(infits, zmin, zmax, outpng):
	# inputs
	# output
	# is it an indiv. PLATE-MJD or the file with all plates?
	if (infits.split('/')[-1].split('.')[-2]=='all'):
		isall = 1
	else:
		isall = 0
	# reading the fits
	hdu          = fits.open(infits)
	data         = hdu[1].data
	fiberid      = data['FIBERID']
	SN_MEDIAN_ALL= data['SN_MEDIAN_ALL']
	Z            = data['Z']
	ZflagJ       = data['Z_reliable']
	Z_NOQSO      = data['Z_NOQSO']
	Z_NOQSOflagJ = data['Z_NOQSO_reliable']

	# plate infos
	plate      = data['plate']
	mjd        = data['mjd']
	run2D      = data['run2D']
	tmp, tmpind= np.unique(plate,return_index=True)
	indlist    = tmpind[np.argsort(mjd[tmpind])] # sorting by increasing MJD
	plate_uniq = plate[indlist]
	mjd_uniq   = mjd[indlist]   # only 1 mjd per plate here
	run2D_uniq = run2D[indlist]
	n_uniq     = len(plate_uniq)
	nexp_uniq  = np.zeros(n_uniq,dtype='int')
	for i in xrange(n_uniq):
		# spPlate fits file to get the nexp center of the plate 
		spPlate = (os.getenv('EBOSS_ROOT')+'/spectro/redux/'+run2D_uniq[i]+'/'+
					str(plate_uniq[i])+'/'+
					'spPlate-'+str(plate_uniq[i])+'-'+str(mjd_uniq[i])+'.fits')
		# plate nexp
		tmphdu       = fits.open(spPlate)
		nexp_uniq[i] = tmphdu[0].header['NEXP_B1']
	info_uniq = np.chararray(n_uniq,itemsize=17)
	for i  in xrange(n_uniq):
		info_uniq[i] = str(plate_uniq[i])+'-'+str(mjd_uniq[i])+'('+str(nexp_uniq[i]).zfill(2)+'exp)'

	# splitting by spectro
	spec0     = (fiberid>=1)   & (fiberid<=1000)
	spec1     = (fiberid>=1)   & (fiberid<=500)
	spec2     = (fiberid>=501) & (fiberid<=1000)
	fcol      = ['0.5','#6699ff','#ff99cc']
	fill      = [True,False,False]
	hatch     = ['-','\\','/']
	col       = ['0.1','b','r']
	alpha     = [0.3,0.6,0.6]
	label     = ['Spec1&2','Spec1','Spec2']

	# efficiency (1st line: Z; 2nd line: Z_NOQSO)
	nobj      = np.zeros(3,dtype='int')
	Eff       = np.zeros((2,3)) # no zflag cut
	EffJ      = np.zeros((2,3)) # Johan zflag cut
	for i in xrange(3):
		exec 'tmpspec = spec'+str(i)
		# N
		nobj[i]    = len(Z[tmpspec])
		# zmin<z<zmax
		for iz in xrange(2):
			zquant    = ['Z','Z_NOQSO'][iz]	
			# no ZflagJ cut
			exec 'tmpzok = (tmpspec) & ('+zquant+'>=zmin) & ('+zquant+'<=zmax)'
			Eff[iz,i] = 100. * float(len(Z[tmpzok])) / nobj[i]
			# with ZflagJ cut
			exec 'tmpzok = (tmpspec) & ('+zquant+'>=zmin) & ('+zquant+'<=zmax) & ('+zquant+'flagJ)'
			EffJ[iz,i] = 100. * float(len(Z[tmpzok])) / nobj[i]
	# starting plot
	fig    = plt.figure(figsize=(30,7))
	gs     = gridspec.GridSpec(1,3,wspace=0.2)

	## loop on plots
	for ip in xrange(3):
		## SN_MEDIAN_ALL
		if (ip==0):
			binwidth  = 0.1
			binmin    = 0.0
			binmax    = 7.0
			xquant    = SN_MEDIAN_ALL
			zquant    = 'Z'
			tmpZflagJ = ZflagJ
			xlabel    = 'SN_MEDIAN_ALL'
			xlim      = [0,3]
			ylim      = [0,2]
			lablist   = ['N','SN_MEDIAN_ALL']
		## Z
		if (ip==1):
			binwidth  = 0.05
			binmin    = 0.0
			binmax    = 7.0
			xquant    = Z
			zquant    = 'Z'
			tmpZflagJ = ZflagJ
			xlabel    = 'Z'
			xlim      = [0,2]
			ylim      = [0,4.5]
			lablist   = ['N',str(zmin)+'<Z<'+str(zmax)]
		## Z_NOQSO
		if (ip==2):
			binwidth  = 0.05
			binmin    = 0.0
			binmax    = 7.0
			xquant    = Z_NOQSO
			zquant    = 'Z_NOQSO'
			tmpZflagJ = Z_NOQSOflagJ
			xlabel    = 'Z_NOQSO'
			xlim      = [0,2]
			ylim      = [0,4.5]
			lablist   = ['N',str(zmin)+'<Z_NOQSO<'+str(zmax)]
		### binning
		bins      = np.arange(binmin,binmax+binwidth,binwidth)
		nbins     = len(bins)-1
		bincenters= bins[0:nbins]+binwidth/2.
		ax = fig.add_subplot(gs[ip])
		## loop on spectro
		for i in xrange(3):
			# no ZflagJ cut
			exec 'tmpspec = spec'+str(i)
			tmpX = xquant[tmpspec]
			hist, bla = np.histogram(tmpX,bins=bins,weights=tmpX*0+1./(nobj[i]*binwidth))
			ax.bar(bincenters,hist,width=binwidth,
					color=fcol[i],edgecolor=fcol[i],
					fill=fill[i],hatch=hatch[i],
					alpha=alpha[i],label=label[i])
			# ZflagJ cut
			exec 'tmpspec = (spec'+str(i)+') & (tmpZflagJ)'
			tmpX  = xquant[tmpspec]
			hist, bla = np.histogram(tmpX,bins=bins,weights=tmpX*0+1./(nobj[i]*binwidth))
			# closing the histogram
			tmpx = np.concatenate(( np.array([bincenters[0],bincenters[0]]),
									bincenters+binwidth,
									binwidth+np.array([bincenters[-1],bincenters[-1]])))
			tmpy = np.concatenate(( np.array([0,hist[0]]),
									hist,
									np.array([hist[-1],0])))
			ax.step(tmpx,tmpy,linewidth=2,color=col[i],label=label[i]+' + '+zquant+'flagJ')
		## infos
		tmpxl = 0.01
		tmpdx = 0.12
		tmpy  = 0.95
		tmpdy = -0.05
		fontsize = 15
		for lab in lablist:
			tmpx = 0.35
			ax.text(tmpxl,tmpy,lab,transform=ax.transAxes,color='k',fontsize=fontsize)
			## loop on spectros
			for i in xrange(3):
				if (lab=='N'):
					ax.text(tmpx,tmpy,'%.0f'%nobj[i],transform=ax.transAxes,color=fcol[i],fontsize=fontsize,ha='center')
					tmpx += tmpdx
					exec 'tmpn = len(Z[(spec'+str(i)+') & (tmpZflagJ)])'
					ax.text(tmpx,tmpy,'%.0f'%tmpn,transform=ax.transAxes,color=col[i],fontsize=fontsize,ha='center')
					tmpx += tmpdx
				if (lab=='SN_MEDIAN_ALL'):
					exec 'tmpspec = spec'+str(i)
					tmpSN = np.mean(SN_MEDIAN_ALL[tmpspec])
					ax.text(tmpx,tmpy,'%.2f'%tmpSN,transform=ax.transAxes,color=fcol[i],fontsize=fontsize,ha='center')
					tmpx += tmpdx
					exec 'tmpspec = (spec'+str(i)+') & (tmpZflagJ)'
					tmpSN = np.mean(SN_MEDIAN_ALL[tmpspec])
					ax.text(tmpx,tmpy,'%.2f'%tmpSN,transform=ax.transAxes,color=col[i],fontsize=fontsize,ha='center')
					tmpx += tmpdx
				if (lab==str(zmin)+'<'+zquant+'<'+str(zmax)):
					exec 'tmpzok = (spec'+str(i)+') & ('+zquant+'>=zmin) & ('+zquant+'<=zmax)'
					tmpeff = 100. * float(len(Z[tmpzok])) / nobj[i]
					ax.text(tmpx,tmpy,'%.1f'%tmpeff+'%',transform=ax.transAxes,color=fcol[i],fontsize=fontsize,ha='center')
					tmpx += tmpdx
					exec 'tmpzok = (spec'+str(i)+') & ('+zquant+'>=zmin) & ('+zquant+'<=zmax) & (tmpZflagJ)'
					tmpeff = 100. * float(len(Z[tmpzok])) / nobj[i]
					ax.text(tmpx,tmpy,'%.1f'%tmpeff+'%',transform=ax.transAxes,color=col[i],fontsize=fontsize,ha='center')
					tmpx += tmpdx
			tmpy = tmpy + tmpdy
		## stuff
		ax.tick_params(axis='both', which='major', labelsize=15)
		ax.set_xlabel(xlabel,fontsize=15)
		ax.set_ylabel('"Normalized" counts',fontsize=15)
		ax.set_xlim(xlim)
		ax.set_ylim(ylim)
		ax.grid(True)
		ax.legend(loc=1,fontsize=15,bbox_to_anchor=(1.00, 0.75))
		if ((xlabel=='Z') | (xlabel=='Z_NOQSO')):
			ax.plot(np.zeros(2)+zmin,[0,3],'-',c='y',linewidth=2)
			ax.plot(np.zeros(2)+zmax,[0,3],'-',c='y',linewidth=2)
		ax.set_title(infits.split('/')[-1])
	# saving the plot
	plt.savefig(outpng,bbox_inches='tight')
	plt.close()

	return True



def addZNOQSOinfo(infits, spZallfits, quantlist):
	# inputs:
	# - infits: input catalog with columns FIBERID, ZNUM_NOQSO
	# - spZallfits: corresponding spZall file
	# - quantlist: comma separated string, e.g. 'RCHI2,VDISPZ'
	#    function built to add 'RCHI2'... be sure that 
	#    the quant+'_NOQSO' doesn't already exist in the spZbest!
	# output: 
	# - reads infits
	# - for each FIBERID, grap the ZNUM_NOQSO
	# - appends to infits a new column RCHI2_NOQSO,
	#    with RCHI2 from spallfits for the considered FIBERID and ZNUM_NOQSO

	quantlist = quantlist.split(',')

	# reading infits
	hdu       = fits.open(infits)
	cols      = hdu[1].columns
	fiberid   = hdu[1].data['FIBERID']
	znum_noqso= hdu[1].data['ZNUM_NOQSO']
	nobj      = len(fiberid)

	# reading spZall
	spZall_hdu     = fits.open(spZallfits)
	spZall_fiberid = spZall_hdu[1].data['FIBERID']
	spZall_nobj    = len(spZall_fiberid)
	for quant in quantlist:
		# spZall quant
		exec 'spZall_'+quant+' = spZall_hdu[1].data["'+quant+'"]'
		# initializing (to have the correct type)
		exec quant+'_znoqso = np.copy(hdu[1].data["'+quant+'"])'

	for i in xrange(nobj):
		# spZall indexes for the considered FIBERID
		tmp       = (spZall_fiberid==fiberid[i])
		tmpind    = np.arange(spZall_nobj)[tmp]
		# index corresponding to ZNUM_NOQSO
		spZall_ind = tmpind[znum_noqso[i]-1]
		# attributing the spAll quantity measurement
		for quant in quantlist:
			exec quant+'_znoqso[i] = spZall_'+quant+'[spZall_ind]'

	# adding spZall columns to infits
	newcollist = []
	for quant in quantlist:
		exec 'newcollist.append(fits.Column(name="'+quant+'_NOQSO",format=hdu[1].columns["'+quant+'"].format,array='+quant+'_znoqso))'
	newcols = fits.ColDefs(newcollist)
	newhdu  = fits.BinTableHDU.from_columns(cols+newcols)
	newhdu.writeto(infits,clobber=True)
	return True


def addZALLinfo(infits, znumlist, spZallfits, quantlist, namelist):
	# inputs:
	# - infits: input fits
	# - numlist: np.array with "Best fit redshift/classification index" [1-index]
	# - spZallfits: corresponding spZall file
	# - znumlist: 
	# - quantlist: comma separated string, e.g. 'RCHI2,VDISPZ'
	#    function built to add 'RCHI2'...
	# - namelist: comma separated string with the column names for quantlist
	#   !  be sure that elements from namelist doesn't already exist in the spZbest !
	# output: 
	# - reads infits
	# - for each FIBERID, grap the requested znum
	# - appends to infits new columns (namelist) for the quantlist quantities from znumlist in spZallfits

	quantlist = quantlist.split(',')
	namelist  = namelist.split(',')

	# reading infits
	hdu       = fits.open(infits)
	cols      = hdu[1].columns
	fiberid   = hdu[1].data['FIBERID']
	nobj      = len(fiberid)
	# initializing (to have the correct type)
	for quant in quantlist:
		exec quant+' = np.copy(hdu[1].data["'+quant+'"])'

	# reading spZall
	spZall_hdu     = fits.open(spZallfits)
	spZall_fiberid = spZall_hdu[1].data['FIBERID']
	for quant in quantlist:
		# spZall quant
		exec 'spZall_'+quant+' = spZall_hdu[1].data["'+quant+'"]'

	for i in xrange(nobj):
		# spZall indexes for the considered FIBERID
		tmpind     = np.where(spZall_fiberid==fiberid[i])[0]
		# index corresponding to ZNUM
		spZall_ind = tmpind[znumlist[i]-1]
		# attributing the spAll quantity measurement
		for quant in quantlist:
			exec quant+'[i] = spZall_'+quant+'[spZall_ind]'

	# adding spZall columns to infits
	newcollist = []
	for quant,name in zip(quantlist,namelist):
		exec 'newcollist.append(fits.Column(name="'+name+'",format=hdu[1].columns["'+quant+'"].format,array='+quant+'))'
	newcols = fits.ColDefs(newcollist)
	newhdu  = fits.BinTableHDU.from_columns(cols+newcols)
	newhdu.writeto(infits,clobber=True)
	return True


def MergePlateInfo(PLATE,MJD,CHUNK,RUN1D,RUN2D):
	# inputs
	# output
	OUTROOT   = os.getenv('MYELG_DIR')+'/'+CHUNK+'/cats/'+CHUNK+'.'+RUN2D+'.'+PLATE+'-'+MJD
	PLATEDIR  = os.getenv('EBOSS_ROOT')+'/spectro/redux/'+RUN2D+'/'+PLATE+'/'
	# settings
	TSPHOTCAT = os.getenv('MYELG_DIR')+'/ELG_TS_master.fits' # TS photometric catalogue
	SPPLATE   = PLATEDIR+'/spPlate-'+PLATE+'-'+MJD+'.fits'
	SPZBEST   = PLATEDIR+'/'+RUN1D+'/spZbest-'+PLATE+'-'+MJD+'.fits'
	SPZALL    = PLATEDIR+'/'+RUN1D+'/spZall-'+PLATE+'-'+MJD+'.fits' 
	# EBOSS_TARGET1
	if ((CHUNK=='eboss21') | (CHUNK=='eboss22')):
		EBOSSBIT = 43
	if ((CHUNK=='eboss23') | (CHUNK=='eboss25')):
		EBOSSBIT = 44
	TILEFITS = get_tile_dir(CHUNK)+'/final-'+CHUNK+'.fits'

	# we compute the fits file with zQ,zCont
	SPZELGF_Z       = OUTROOT+'.ELGflag.Z.fits'
	SPZELGF_Z_NOQSO = OUTROOT+'.ELGflag.ZNOQSO.fits'
	print 'SPZELGF_Z=',SPZELGF_Z
	## computing SPZELGF_ Z and SPZELGF_ZNOQSO
	for zquant in ['Z','Z_NOQSO']:
		exec ('a = HandleReducedELGPlate.HandleReducedELGPlate('+
				'SPZELGF_'+zquant+',PLATEDIR,plate=PLATE,mjd=MJD,zquant=zquant,run1d=RUN1D)')
		a.load()
		a.fitLineFluxes()
		a.assign_ELG_redshift_line_flag()
		a.assign_ELG_redshift_continuum_flag()
		a.save_result()

	tmpfits  = TMPDIR+'tmp.fits_'+pid
	tmp2fits = TMPDIR+'tmp2.fits_'+pid
	tmpasc   = TMPDIR+'tmp.asc_'+pid
	# adding to spZbest: XFOCAL,YFOCAL,CARTID
	spBhdu  = fits.open(SPZBEST) # spZbest
	spPhdu  = fits.open(SPPLATE) # spPlate
	newhdu  = fits.BinTableHDU.from_columns(
				spBhdu[1].columns + 
				spPhdu['PLUGMAP'].columns['XFOCAL']+
				spPhdu['PLUGMAP'].columns['YFOCAL']+
				spPhdu['PLUGMAP'].columns['EBOSS_TARGET1']+
				spPhdu['PLUGMAP'].columns['EBOSS_TARGET_ID']+
				fits.Column(name='CARTID',format='K',array=[spPhdu[0].header['CARTID'] for x in spBhdu[1].data['PLUG_RA']]))
	newhdu.writeto(tmpfits,clobber=True)
	# selecting ELG targets only
	hdu      = fits.open(tmpfits)
	ebosstag = hdu[1].data['EBOSS_TARGET1']
	bitsel   = ((ebosstag & 2**EBOSSBIT)>0)
	fits.update(tmpfits,hdu[1].data[bitsel],1)	

	# adding Z_NOQSO RCHI2
	hdu        = fits.open(tmpfits)
	ZNUM_NOQSO = hdu[1].data['ZNUM_NOQSO']
	a  = addZALLinfo(tmpfits, ZNUM_NOQSO, SPZALL, 'RCHI2', 'RCHI2_NOQSO')
	# adding second best fit infos
	a  = addZALLinfo(tmpfits, ZNUM_NOQSO*0+2, SPZALL, 'Z,CLASS,RCHI2', 'Z_SCND,CLASS_SCND,RCHI2_SCND')

	# adding ELG zFlag if the file is there; flags from Z
	for zquant in ['Z','Z_NOQSO']:
		exec 'SPZELGF = SPZELGF_'+zquant
		## to rename columns and adding a column zquant+'_reliable'
		tmpstr = "echo 'addcol "+zquant+"_reliable \"zQ>=2 || (zQ>=1 && zCont>0) || (zQ>=0 && zCont>=2.5)\"' > "+tmpasc
		print tmpstr
		subprocess.call(tmpstr, shell=True)
		tmpstr = (STILTSCMD+' tpipe in='+SPZELGF+' ifmt=fits cmd=\'meta Name\' | '+
			'grep -v \\+ | grep -v Name | '+
			'awk \'{print "colmeta -name '+zquant+'_"$2 " "$2}\' >> '+tmpasc)
		print tmpstr
		subprocess.call(tmpstr, shell=True)
		## cross-matching
		tmpstr = (STILTSCMD+' tmatch2 matcher=exact join=all1 '+
			'in1='+tmpfits+' ifmt1=fits values1="FIBERID" '+
			'icmd2=@'+tmpasc+' '+
			'in2='+SPZELGF+' ifmt2=fits values2="'+zquant+'_FIBERID" '+
			'out='+tmp2fits+' ofmt=fits')
		print tmpstr
		subprocess.call(tmpstr, shell=True)
		subprocess.call('mv '+tmp2fits+' '+tmpfits, shell=True)
		subprocess.call('rm '+tmpasc, shell=True)

	# crossing with TSPHOTCAT
	tmpstr = (STILTSCMD+' tmatch2 matcher=exact join=1and2 '+	# no duplicates on a given plate
				'in1='+tmpfits+'   ifmt1=fits values1="EBOSS_TARGET_ID" '+
				"icmd1='delcols \"EBOSS_TARGET1\";' "+			# already in TSPHOTCAT
				'in2='+TSPHOTCAT+' ifmt2=fits values2="EBOSS_TARGET_ID" '+
				'out='+OUTROOT+'.fits ofmt=fits '+
				"ocmd='colmeta -name EBOSS_TARGET_ID EBOSS_TARGET_ID_1; delcols \"EBOSS_TARGET_ID_2\";'")
	print tmpstr
	subprocess.call(tmpstr, shell=True)
	subprocess.call('rm '+tmpfits, shell=True)

	# adding hasfiber, plate_ssr and systematic weights
	hdu     = fits.open(OUTROOT+'.fits')
	cols    = hdu[1].columns
	data    = hdu[1].data
	ra      = data['ra']
	dec     = data['dec']
	zwarn   = data['ZWARNING']
	zok     = data['Z_reliable']
	CLASS   = data['CLASS']
	hasfiber= np.ones(len(ra),dtype='int')
	tmp     =  (((zwarn & 2**1)>0) |	# not LITTLE_COVERAGE
				((zwarn & 2**7)>0) |	# not UNPLUGGED
				((zwarn & 2**8)>0) |	# not BAD_TARGET
				((zwarn & 2**9)>0))		# not NODATA
	hasfiber[tmp] = -1
	tmp    = (hasfiber==1) & (zok) & (CLASS!='STAR')
	ssr    = float(len(ra[tmp]))/float(len(ra[hasfiber==1]))
	systweight = get_hpsystweight(ra,dec,CHUNK)
	newcols = [cols[x] for x in cols.names]
	newcols.append(fits.Column(name='WEIGHT_SYSTOT',format='E',array=systweight))
	newcols.append(fits.Column(name='hasfiber',     format='I',array=hasfiber))
	newcols.append(fits.Column(name='plate_SSR',    format='E',array=ra*0.+ssr))
	newhdu = fits.BinTableHDU.from_columns(newcols)
	newhdu.writeto(OUTROOT+'.fits',clobber=True)

	return True



def plot_syst(CHUNK,RUN2D,SYSTFITS):

	# settings
	chunkdir= os.getenv('MYELG_DIR')+'/'+CHUNK+'/'
	infits  = chunkdir+'cats/'+CHUNK+'.'+RUN2D+'.meta.fits'
	outroot = chunkdir+'plots/'+CHUNK+'.'+RUN2D+'.radec'

	# title
	hdu     = fits.open(infits)
	tmp     = (hdu[1].data['ISLATEST']==1)
	title   = (CHUNK+' , '+RUN2D+' ['+
				str(np.min(hdu[1].data['MJD'][tmp]))+r'$\leq$MJD$\leq$'+str(np.max(hdu[1].data['MJD'][tmp]))+']')

	# catalog with systematics values
	hdu        = fits.open(SYSTFITS)
	data       = hdu[1].data
	systra     = data['hpra']
	systdec    = data['hpdec']
	systlist   = np.array(['nstar','ebv','psfsize_g','psfsize_r','psfsize_z','galdepth_g','galdepth_r','galdepth_z'])
	for quant in systlist:
		exec 'syst'+quant+' = data["hp'+quant+'"]'
	print str(len(systra))+' syst obj. before cutting'
	if (CHUNK=='eboss21'):
		tmp = (systra>317) & (systra<360) & (np.abs(systdec)<2.)
	if (CHUNK=='eboss22'):
		tmp = (systra>0) & (systra<45) & (np.abs(systdec)<5.)
	if (CHUNK=='eboss23'):
		tmp = (((systra>126.0) & (systra<157.0) & (systdec>13.8) & (systdec<29.0)) & 
				((systra>136.5) | (systdec>16.0)) & ((systra<142.5) | (systdec<27.0)))
	if (CHUNK=='eboss25'):
		tmp = (	((systra>131)   & (systra<166) & (systdec>29) & (systdec<32.5)) |
				((systra>142.5) & (systra<166) & (systdec>27) & (systdec<29)) |
				((systra>157)   & (systra<166) & (systdec>13) & (systdec<27)))
	for quant in np.append(np.array(['ra','dec']),systlist):
		exec 'syst'+quant+' = syst'+quant+'[tmp]'
	print str(len(systra))+' syst obj. after cutting'

    # reading metafits to get plate infos
	hdu      = fits.open(infits)
	tmpplate = hdu[1].data['PLATE']
	bla, tmpind = np.unique(tmpplate, return_index=True)
	platelist= tmpplate[tmpind]
	rac      = hdu[1].data['RADEG'][tmpind]
	decc     = hdu[1].data['DECDEG'][tmpind]
	nplate   = len(platelist)

	for quant in systlist:
		fig       = plt.figure(figsize=(20,7))
		ax        = fig.add_subplot(111)
		# cbar settings
		exec 'c = syst'+quant
		if (quant=='nstar'):
			cbarlabel = 'PS1 stellar density [/deg2]'
			clim      = [1500,2500]
		if (quant=='ebv'):
			cbarlabel = 'e(b-v)'
			clim      = [0,0.1]
		if (quant=='galdepth_g'):
			cbarlabel = quant[-1]+'-band 5 sigma extended-source depth [ABmag]'
			if ((CHUNK=='eboss21') | (CHUNK=='eboss22')):
				clim = [23.5,25.5]
			if ((CHUNK=='eboss23') | (CHUNK=='eboss25')):
				clim = [23.0,25.0]
		if (quant=='galdepth_r'):
			cbarlabel = quant[-1]+'-band 5 sigma extended-source depth [ABmag]'
			if ((CHUNK=='eboss21') | (CHUNK=='eboss22')):
				clim = [23.5,25.5]
			if ((CHUNK=='eboss23') | (CHUNK=='eboss25')):
				clim = [22.5,24.5]
		if (quant=='galdepth_z'):
			cbarlabel = quant[-1]+'-band 5 sigma extended-source depth [ABmag]'
			if ((CHUNK=='eboss21') | (CHUNK=='eboss22')):
				clim = [22.0,24.0]
			if ((CHUNK=='eboss23') | (CHUNK=='eboss25')):
				clim = [21.5,23.5]
		if (quant[0:3]=='psf'):
			cbarlabel = quant[-1]+'-band PSF FWHM [arcsec]'
			clim      = [0.8,2.0]
		SC = ax.scatter(systra,systdec,c=c,s=40,edgecolor='none',alpha=0.8,vmin=clim[0],vmax=clim[1],cmap=matplotlib.cm.jet)
		# plate number
		for i in xrange(nplate):
			ax.text(rac[i],decc[i],platelist[i],fontsize=10,fontweight='bold',rotation='vertical',va='center',ha='center')
			# plate contours
			e = Ellipse(xy=np.array([rac[i],decc[i]]),width=3.0/np.cos(decc[i]/180.*np.pi),height=3.0,angle=0)
			e.set_color('k')
			e.set_facecolor('none')
			e.set_linewidth(2)
			e.set_alpha(1.0)
			e.set_clip_box(ax.bbox)
			ax.add_artist(e)
		# ELG footprint + xlim,ylim
		footprint_ra,footprint_dec,xlim,ylim = get_footprint(CHUNK)
		ax.plot(footprint_ra,footprint_dec,c='k',linewidth=4,zorder=10)
		ax.set_xlim(xlim)
		ax.set_ylim(ylim)
		# colorbar
		cbarticks     = np.linspace(clim[0],clim[1],11)
		cbar_ylab     = ['%.2f' % x for x in cbarticks]
		cbar_ylab[0]  = '<'+cbar_ylab[0]
		cbar_ylab[-1] = '>'+cbar_ylab[-1]
		cbar = plt.colorbar(SC,ticks=cbarticks)
		cbar.ax.set_yticklabels(cbar_ylab)
		cbar.set_label(cbarlabel)
		cbar.set_clim(clim)
		ax.grid(True)
		ax.set_xlabel('R.A.')
		ax.set_ylabel('DEC.')
		ax.set_title(title)
		plt.savefig(outroot+'.'+quant+'.png',bbox_inches='tight')
		plt.close()
		print outroot+'.'+quant+'.png done'

	return True



def analysis_all(OUTROOT,FITSLIST,ZQUANT,METAFITS,RUN2D,XQUANTLIST=['RSN2']):
	# input:
	# output:

	# GMT date
	GMTdate = 'GMT '+datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")

	nfits   = len(FITSLIST)

	# grabbing plate infos, SN_MEDIAN_ALL, and efficiency
	hdu     = fits.open(METAFITS)
	ZMIN    = hdu[0].header['ZMIN']
	ZMAX    = hdu[0].header['ZMAX']
	data    = hdu[1].data
	# corresponding indexes in metafits
	indlist = np.zeros(nfits,dtype='int')
	for i in xrange(nfits):
		indlist[i] = np.arange(len(data['FITSFILE']))[(data['FITSFILE']==FITSLIST[i])][0]
	# meta quantities
	plate = data['PLATE']  [indlist]
	mjd   = data['MJD']    [indlist]
	nexp  = data['NEXP']   [indlist]
	nobs  = data['NOBSEXP'][indlist]
	exec 'eff   = data["EFFJ'+re.sub('Z','',ZQUANT)+'"][indlist]'
	for quant in XQUANTLIST:
		exec quant+' = data["'+quant+'"][indlist]'
	# one color for each plate
	plate_uniq    = np.unique(plate)
	print 'plate_uniq = ', plate_uniq
	platecol_uniq = get_colors(len(plate_uniq))
	
	# plotting Eff vs. X
	for xquant in XQUANTLIST:
		exec 'x = '+xquant
		if (xquant=='SN_MEDIAN_G'):
			xlim  = [0.3,1.2]
		if (xquant=='SN_MEDIAN_ALL'):
			xlim  = [0.3,1.6]
		if (xquant=='nexp'):
			xlim  = [0,20]
		if (xquant=='SN_FOII'):
			xlim  = [-1,20]
		if (xquant=='SN2EXTI'):
			xlim   = [0,100]
		if (xquant=='RSN2'):
			xlim   = [0,60]
		ylim  = [30,80]
		col   = ['0.5','b','r']
		alpha = [1.0,0.2,0.2]
		label = ['Spec1+2','Spec1','Spec2']
		title = str(ZMIN)+'<'+ZQUANT+'<'+str(ZMAX)+'+zflagJ'
		# figure with all plates
		figall, axall = plt.subplots()
		for i in xrange(len(plate_uniq)):
			ind_i = np.arange(nfits)[((plate==plate_uniq[i]) & (eff[:,0]>0))]
			# loop on spectros
			#for j in xrange(3):
			for j in xrange(1):
				# sorting by increasing x
				ind_i = ind_i[np.argsort(x[ind_i,j])]
				## plotting as stars when QUALITY=bad
				ind_good = []
				ind_bad  = []
				for tmpind in ind_i:
					tmpstr = ("awk '{if ($1=="+str(plate_uniq[i])+" && $2=="+str(mjd[tmpind])+") print $7;}' "+
								os.getenv('EBOSS_ROOT')+"/spectro/redux/"+RUN2D+"/platelist-mjdsort.txt")
					p1 = subprocess.Popen(tmpstr,stdout=subprocess.PIPE, shell=True)
					tmpqual = p1.communicate()[0].strip()
					if (tmpqual=='good'):
						ind_good.append(tmpind)
					if (tmpqual=='bad'):
						ind_bad.append(tmpind)
				## plot with all plates
				if (i==0):
					axall.scatter(x[ind_i,j],   100.*eff[ind_i,j],   marker='o',s=20, c=col[j],alpha=0.5,     label=label[j])
					axall.scatter(x[ind_good,j],100.*eff[ind_good,j],marker='s',s=50, c='b',alpha=alpha[j],label=label[j]+' [QUALITY="good"]')
					axall.scatter(x[ind_bad ,j],100.*eff[ind_bad, j],marker='*',s=100,c='r',alpha=alpha[j],label=label[j]+' [QUALITY="bad"]') 
				else:
					axall.scatter(x[ind_i,j],   100.*eff[ind_i,j],   marker='o',s=20, c=col[j],alpha=0.5)
					axall.scatter(x[ind_good,j],100.*eff[ind_good,j],marker='s',s=50, c='b',alpha=alpha[j])
					axall.scatter(x[ind_bad ,j],100.*eff[ind_bad, j],marker='*',s=100,c='r',alpha=alpha[j])
				axall.plot   (x[ind_i,j],100.*eff[ind_i,j],color=platecol_uniq[i],alpha=0.3)
		axall.set_title(title)
		if (xquant=='RSN2'):
			axall.set_xlabel('RSN2 [scale factor 2.58]')
		else:
			axall.set_xlabel(xquant)
		axall.set_ylabel('Efficiency [%]')
		axall.set_xlim(xlim)
		axall.set_ylim(ylim)
		axall.grid(True)
		axall.legend(loc=4)
		plt.savefig(OUTROOT+'.EFFvs'+xquant+'.png',bbox_inches='tight')
		plt.close()
	return True



def create_metafits(OUTFITS,FITSLIST,ZMIN,ZMAX,RUN2D,RUN1D):
	nfits   = len(FITSLIST)
	print nfits
	# grabbing plate infos, SN_MEDIAN_ALL, and efficiency
	SPZBESTLIST = np.array(['' for x in xrange(nfits)],dtype='object')
	EXPNUM      = np.array(['' for x in xrange(nfits)],dtype='object')
	for quant in ['RADEG','DECDEG','SEEING50','AIRMASS']:
		exec quant+' = np.zeros(nfits)'
	for quant in ['PLATE','MJD','NOBSEXP','NEXP','ISLATEST']:
		exec quant+' = np.zeros(nfits,dtype="int")'
	NOBJ   = np.zeros((nfits,3),dtype='int')
	for quant in [	'SN_MEDIAN_U','SN_MEDIAN_G','SN_MEDIAN_R','SN_MEDIAN_I','SN_MEDIAN_Z',
					'SN_MEDIAN_ALL','RSN2',
					'SNFOII','SNFOII_NOQSO',
					'EFFJ','EFFJ_NOQSO',
					'FRACZOK','FRACZSTAR',
					'FRACZOK_NOQSO','FRACZSTAR_NOQSO']:
		exec quant+' = np.zeros((nfits,3))'
	for i in xrange(nfits):
		# reading fits
		print FITSLIST[i]
		infits  = FITSLIST[i]
		hdu     = fits.open(infits)
		data    = hdu[1].data
		hasfiber= data['hasfiber']
		fiberid = data['FIBERID']
 		snmedall= data['SN_MEDIAN_ALL']
		snmedu  = data['SN_MEDIAN'][:,0]
		snmedg  = data['SN_MEDIAN'][:,1]
		snmedr  = data['SN_MEDIAN'][:,2]
		snmedi  = data['SN_MEDIAN'][:,3]
		snmedz  = data['SN_MEDIAN'][:,4]
		for ZQUANT in ['Z','Z_NOQSO']:
			exec ZQUANT+' = data["'+ZQUANT+'"]'
			exec 'CLASS'+re.sub('Z','',ZQUANT)+'    = data["CLASS'+re.sub('Z','',ZQUANT)+'"]'
			exec ZQUANT+'_reliable                  = data["'+ZQUANT+'_reliable"]'
			exec 'zCont'+re.sub('Z','',ZQUANT)+'    = data["'+ZQUANT+'_zCont"]'
			exec 'FOII' +re.sub('Z','',ZQUANT)+'    = data["'+ZQUANT+'_O2_3728_flux"]'
			exec 'FOIIerr' +re.sub('Z','',ZQUANT)+' = data["'+ZQUANT+'_O2_3728_fluxErr"]'
		# plate infos
		PLATE[i]       = int(FITSLIST[i].split('/')[-1].replace('-','.').split('.')[2])
		MJD[i]         = int(FITSLIST[i].split('/')[-1].replace('-','.').split('.')[3])
		platedir       = os.getenv('EBOSS_ROOT')+'/spectro/redux/'+RUN2D+'/'+str(PLATE[i])+'/'
		SPZBESTLIST[i] = platedir+RUN1D+'/spZbest-'+str(PLATE[i])+'-'+str(MJD[i])+'.fits'
		spPlate        = platedir+'spPlate-'+str(PLATE[i])+'-'+str(MJD[i])+'.fits'
		spPhdu         = fits.open(spPlate)
		spPhdr         = spPhdu[0].header
		# plate ra,dec, seeing, airmass
		for quant in ['RADEG','DECDEG','SEEING50','AIRMASS']:
			exec quant+'[i] = spPhdr["'+quant+'"]'
		# nobsexp
		spPlancomb = platedir+'spPlancomb-'+str(PLATE[i])+'-'+str(MJD[i])+'.par'
		tmpstr     = 'grep science '+spPlancomb+' | wc -l'
		p1         = subprocess.Popen(tmpstr, stdout=subprocess.PIPE, shell=True)
		NOBSEXP[i] = int(p1.communicate()[0].strip())
		# nexp
		NEXP[i]    = spPhdr['NEXP_B1']

		## SN2 values from Vivek: looping on indiv. exposures to sum the sn2
		## looping on indiv. exposures
		for j in xrange(spPhdr['NEXP_B1']):
			# exposure number
			tmpnum   = spPhdr['EXPID'+str(1+j).zfill(2)].split('-')[1]
			# updating the list
			EXPNUM[i]+= tmpnum+','
			# at what MJD was this exposure observed?
			tmpstr  = 'grep '+tmpnum+' '+platedir+'spPlan2d-'+str(PLATE[i])+'-?????.par | awk \'{print $3}\''
			p1      = subprocess.Popen(tmpstr, stdout=subprocess.PIPE, shell=True)
			tmpmjd  = p1.communicate()[0].strip('\n')
			## RSN2 re-scaling after MJD>=57712
			## see Vivek's email : [eboss-pipeline 5862] new scale factor for ELG plates]
			## "I have retuned the scale factor for ELG observations from 1.467 to 2.58."
			if (int(tmpmjd)<57712):
				tmpcoeff = 2.58 / 1.467
			else:
				tmpcoeff = 1.
			# SOS fits files
			r1fits   = os.getenv('BOSS_SOS')+'/'+tmpmjd+'/sci-'+str(PLATE[i])+'-r1-'+tmpnum.zfill(8)+'.fits'
			RSN2[i,1]+= tmpcoeff * fits.open(r1fits)[0].header['FRAMESN2']
			r2fits   = os.getenv('BOSS_SOS')+'/'+tmpmjd+'/sci-'+str(PLATE[i])+'-r2-'+tmpnum.zfill(8)+'.fits'
			RSN2[i,2]+= tmpcoeff * fits.open(r2fits)[0].header['FRAMESN2']
		## removing the last ','
		EXPNUM[i] = EXPNUM[i][0:-1]
		RSN2[i,0] = (RSN2[i,1]+RSN2[i,2])/2.

		# spec1 & spec2
		spec0 = (fiberid>=1)   & (fiberid<=1000)	# spectro 1+2
		spec1 = (fiberid>=1)   & (fiberid<=500)		# spectro 1
		spec2 = (fiberid>=501) & (fiberid<=1000)	# spectro 2
		#
		tmpfiberok = (hasfiber==1)					# valid fiber
		for j in xrange(3):
			exec 'NOBJ[i,j] = len(Z[spec'+str(j)+'])'
			# SN_MEDIAN
			exec 'SN_MEDIAN_U[i,j]   = np.median(snmedu[spec'+str(j)+'])'
			exec 'SN_MEDIAN_G[i,j]   = np.median(snmedg[spec'+str(j)+'])'
			exec 'SN_MEDIAN_R[i,j]   = np.median(snmedr[spec'+str(j)+'])'
			exec 'SN_MEDIAN_I[i,j]   = np.median(snmedi[spec'+str(j)+'])'
			exec 'SN_MEDIAN_Z[i,j]   = np.median(snmedz[spec'+str(j)+'])'
			# SN_MEDIAN_ALL
			exec 'SN_MEDIAN_ALL[i,j] = np.median(snmedall[spec'+str(j)+'])'
			for ZQUANT in ['Z','Z_NOQSO']:
				ZQUANTSTR = re.sub('Z','',ZQUANT) # 'Z'->'' ; 'Z_NOQSO'->'_NOQSO'
				exec 'tmpj     = (spec'+str(j)+')'
				exec 'tmpstar  = (CLASS'+ZQUANTSTR+'=="STAR")'
				exec 'tmpzok   = ('+ZQUANT+'_reliable)'
				exec 'tmpzvalok= ('+ZQUANT+'>=ZMIN) & ('+ZQUANT+'<=ZMAX)'
				exec 'tmpfoiiok= (FOII'+ZQUANTSTR+'>-90.)'# valid FOII meas.
				# SNFOII
				tmp = (tmpj) & (tmpfiberok) & (tmpfoiiok)
				exec 'SNFOII'+ZQUANTSTR+'[i,j] = np.median(FOII'+ZQUANTSTR+'[tmp]/FOIIerr'+ZQUANTSTR+'[tmp])'
				# reliable z (not star)
				tmpn  = float(len(Z[(tmpj) & (tmpfiberok) & (~tmpstar)]))
				tmpnok= float(len(Z[(tmpj) & (tmpfiberok) & (~tmpstar) & (tmpzok)]))
				exec 'FRACZOK'+ZQUANTSTR+'[i,j] = tmpnok/tmpn'
				# stars
				tmpn  = float(len(Z[(tmpj) & (tmpfiberok)]))
				tmpnok= float(len(Z[(tmpj) & (tmpfiberok) & (tmpstar)]))
				exec 'FRACZSTAR'+ZQUANTSTR+'[i,j] = tmpnok/tmpn'
				# 0.7<z_reliable<1.1
				tmpn  = float(len(Z[(tmpj) & (tmpfiberok) & (~tmpstar)]))
				tmpnok= float(len(Z[(tmpj) & (tmpfiberok) & (~tmpstar) & (tmpzok) & (tmpzvalok)]))
				exec 'EFFJ'+ZQUANTSTR+'[i,j] = tmpnok/tmpn'
	# marking the latest file, for each plate
	for p in np.unique(PLATE):
		tmpind = np.arange(nfits)[PLATE==p]
		tmpmjd = MJD[tmpind]
		ISLATEST[tmpind[np.argmax(tmpmjd)]] = 1
	# sorting by increasing plate
	indsort = np.argsort(PLATE)
	# creating the outfits
	cols_list = []
	cols_list.append(fits.Column(name='SPZBESTFILE',
								 format='A'+str(np.max(np.array([len(x) for x in SPZBESTLIST]))),
								 array=SPZBESTLIST[indsort]))
	cols_list.append(fits.Column(name='FITSFILE',
								 format='A'+str(np.max(np.array([len(x) for x in FITSLIST]))),
								 array=np.array(FITSLIST)[indsort]))
	cols_list.append(fits.Column(name='EXPNUM',
								 format='A'+str(np.max(np.array([len(x) for x in EXPNUM]))),
								 array=EXPNUM[indsort]))
	for quant in ['PLATE','MJD','NOBSEXP','NEXP','ISLATEST']:
		exec 'cols_list.append(fits.Column(name="'+quant+'",format="K",array='+quant+'[indsort]))'
	for quant in ['RADEG','DECDEG','SEEING50','AIRMASS']:
		exec 'cols_list.append(fits.Column(name="'+quant+'",format="E",array='+quant+'[indsort]))'
	for quant in [	'NOBJ',
					'SN_MEDIAN_U','SN_MEDIAN_G','SN_MEDIAN_R','SN_MEDIAN_I','SN_MEDIAN_Z',
					'SN_MEDIAN_ALL','RSN2',
					'SNFOII','SNFOII_NOQSO',
					'EFFJ','EFFJ_NOQSO',
					'FRACZOK','FRACZSTAR',
                    'FRACZOK_NOQSO','FRACZSTAR_NOQSO']:
		if (quant=='NOBJ'):
			tmpformat = '3K'
		else:
			tmpformat = '3E'
		exec 'cols_list.append(fits.Column(name="'+quant+'",format="'+tmpformat+'",array='+quant+'[indsort,:]))'
	cols = fits.ColDefs(cols_list)
	# header
	prihdr = fits.Header()
	prihdr['DATE'] = 'GMT '+datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S")
	prihdr['RUN2D']= RUN2D
	prihdr['RUN1D']= RUN1D
	prihdr['ZMIN'] = ZMIN
	prihdr['ZMAX'] = ZMAX
	prihdu = fits.PrimaryHDU(header=prihdr)
	hdu    = fits.BinTableHDU.from_columns(cols)
	hdulist= fits.HDUList([prihdu, hdu])
	hdulist.writeto(OUTFITS,clobber=True)

	return True
	


def fast_launch(INFITS,FASTDIR,FASTROOT):
	CURRDIR   = os.getcwd()
	# reading fits
	hdu   = fits.open(INFITS)
	data  = hdu[1].data
	id    = np.arange(len(data))
	zspec = data['Z']
	# AB mags
	m_g  = data['g']
	m_r  = m_g - data['gr']
	m_z  = m_r - data['rz']
	m_w1 = m_r - data['rw1']
	m_w2 = m_r - data['rw2']
	# fluxes in nanomaggies (ABmag = 22.5 - 2.5*log10(flux))
	for filt in ['g','r','z','w1','w2']:
		exec 'f_'+filt+' = 10.**(-0.4*(m_'+filt+'-22.5))'
	# flux errors
	for filt in ['g','r','z','w1','w2']:
		exec 'e_'+filt+' = f_'+filt+' / data["snr_'+filt+'"]'
	#
	for filt in ['g','r','z','w1','w2']:
		exec 'tmp = (np.isfinite(f_'+filt+')==False) | (f_'+filt+'<0)'
		exec 'f_'+filt+'[tmp] = -99'
		exec 'e_'+filt+'[tmp] = -99'
	# writing FAST input cat
	f = open(FASTDIR+FASTROOT+'.cat','w')
	f.write('# id z_spec F224 E224 F225 E225 F227 E227 F229 E229 F230 E230\n')
	for i in xrange(len(id)):
		f.write(str(id[i])+'\t'+str(zspec[i])+'\t'+
				str(f_g [i])+'\t'+str(e_g [i])+'\t'+
				str(f_r [i])+'\t'+str(e_r [i])+'\t'+
				str(f_z [i])+'\t'+str(e_z [i])+'\t'+
				str(f_w1[i])+'\t'+str(e_w1[i])+'\t'+
				str(f_w2[i])+'\t'+str(e_w2[i])+'\n')
	f.close()
	# launching fast
	os.chdir(FASTDIR)
	tmpstr = 'idl -e fast -args '+FASTDIR+'/'+FASTROOT+'.param'
	print tmpstr
	subprocess.call(tmpstr, shell=True)
	os.chdir(CURRDIR)



def fast_add(INFITS,FASTDIR,FASTROOT):
	# number of objects
	nobj = len(fits.open(INFITS)[1].data)
	# FAST ABmag corresponding to the best-fit spectrum
	## reading the filters, to compute fitted AB mag
	decam_g = speclite.filters.load_filters('decam2014-g')
	decam_r = speclite.filters.load_filters('decam2014-r')
	decam_z = speclite.filters.load_filters('decam2014-z')
	wise_w1 = speclite.filters.load_filters('wise2010-W1')
	wise_w2 = speclite.filters.load_filters('wise2010-W2')
	## array to store the mags
	instrlist= ['decam','decam','decam','wise','wise']
	filtlist = ['g','r','z','w1','w2']
	for filt in filtlist:
		exec 'fast_'+filt+'mag = np.zeros(nobj)'
	## loop on objects
	for i in xrange(nobj):
		### best-fit model (1e-19.erg.s-1.cm-2.A-1)
		ll, ff = np.loadtxt(FASTDIR+'/BEST_FITS/'+FASTROOT+'_'+str(i)+'.fit', unpack=True)
		### ABmag
		for instr,filt in zip(instrlist,filtlist):
			exec ('fast_'+filt+'mag[i]  = '+
					instr+'_'+filt+'.get_ab_magnitudes('+
					'ff*10**(-19.)*u.erg/u.s/u.cm**2/u.AA,'+
					'wavelength=ll*u.AA)._data[0][0]')
	## storing to a fits file
	collist = []
	for filt in filtlist:
		exec ('collist.append(fits.Column('+
				'name="fast_'+filt+'mag" ,format="E", '+
				'array=fast_'+filt+'mag))')
	cols = fits.ColDefs(collist)
	hdu  = fits.BinTableHDU.from_columns(cols)
	hdu.writeto('tmp.fits_'+pid,clobber=True)
	# adding FAST measurements to the fits
	tmpstr = (STILTSCMD+' tmatchn nin=3 matcher=exact '+
			'in1='+INFITS+" ifmt1=fits values1='$0' "+
			'in2='+FASTDIR+'/'+FASTROOT+".fout ifmt2=ascii values2='$0' "+
			'in3=tmp.fits_'+pid+" ifmt3=fits values3='$0' "+
			'out=tmp.match.fits_'+pid+' ofmt=fits '+
			"icmd2='colmeta  -name fast_id id; "+
					"colmeta -name fast_z z; "+
					"colmeta -name fast_ltau ltau; "+
					"colmeta -name fast_metal metal; "+
					"colmeta -name fast_lage lage; "+
					"colmeta -name fast_Av Av; "+
					"colmeta -name fast_lmass lmass; "+
					"colmeta -name fast_lsfr lsfr; "+
					"colmeta -name fast_lssfr lssfr; "+
					"colmeta -name fast_la2t la2t; "+
					"colmeta -name fast_chi2 chi2;'")
	print tmpstr
	subprocess.call(tmpstr, shell=True)
	subprocess.call('mv tmp.match.fits_'+pid+' '+INFITS, shell=True)
	subprocess.call('rm tmp.fits_'+pid, shell=True)
	return True



def fast_plot(LATESTFITS,CHUNK,FASTDIR,REDUXDIR,RUN2D,OUTPDF,funit='flam'):
	# reading fits
	hdu  = fits.open(LATESTFITS)
	data =  hdu[1].data
	# AB mags
	m_g  = data['g']
	m_r  = m_g - data['gr']
	m_z  = m_r - data['rz']
	m_w1 = m_r - data['rw1']
	m_w2 = m_r - data['rw2']
	# fluxes in nanomaggies (ABmag = 22.5 - 2.5*log10(flux))
	for filt in ['g','r','z','w1','w2']:
		exec 'f_'+filt+' = 10.**(-0.4*(m_'+filt+'-22.5))'
	# flux errors
	for filt in ['g','r','z','w1','w2']:
		exec 'e_'+filt+' = f_'+filt+' / data["snr_'+filt+'"]'
	# spectro
	zspec= data['Z']
	plate= data['plate'].astype('str')
	mjd  = data['mjd'].astype('str')
	fiber= data['fiberid'].astype('str')
	zQ   = data['Z_zQ']
	zCont= data['Z_zCont']
	lmass= data['fast_lmass']
	ltau = data['fast_ltau']
	lage = data['fast_lage']
	Av   = data['fast_Av']
	for filt in ['g','r','z','w1','w2']:
		exec 'fast_'+filt+' = data["fast_'+filt+'mag"]'
	### 0.7<zok<1.1
	tmp  = (	(zspec>0.7) & (zspec<1.1) &
			((zQ>=2) | ((zQ>=1) & (zCont>0)) | ((zQ>=0.0) & (zCont>=2.5))))
	sel  = np.arange(len(zspec))[tmp]
	### sorting by increasing lmass
	sel  = sel[np.argsort(lmass[sel])]
	### picking N gals
	tmpn = 1e2
	sel  = sel[np.round(np.linspace(0,len(sel)-1,tmpn)).astype('int')]
	print str(len(sel))+' objects plotted'
	## grabbing lambda_c (grzW1W2)
	lc, bla = np.loadtxt(FASTDIR+'/BEST_FITS/'+CHUNK+'.'+RUN2D+'.latest_0.input_res.fit', unpack=True)
	## f and ef are in nanomaggies (ABZP=22.5)
	if (funit=='flam'):
		## we convert then to 1e-19.erg.s-1.cm-2.A-1 (https://en.wikipedia.org/wiki/AB_magnitude)
		## f[1e-19.erg.s-1.cm-2.A-1] = f[nanommagies] * coeff
		coeff_decals = 1e19 * 10.**(-0.4*(22.5-8.90)) / (3.34*10.**4*lc**2.)
		ylabel = 'Flux [1e-19 erg/s/cm2/A]'
		ylim   = [1e0,1e3]
		legloc = 3
	if (funit=='fnu'):
		## mAB = -2.5*log10(f[microJansky])+23.9
		coeff_decals = 10.**(-0.4*(22.5-23.90))
		ylabel ='Flux [microJy]'
		ylim   = [1e0,1e2]
		legloc = 4
	with PdfPages(OUTPDF) as pdf:	
		for i  in xrange(len(sel)):
			seli = sel[i]
			fig,ax = plt.subplots()
			ll, ff = np.loadtxt(FASTDIR+'/BEST_FITS/'+CHUNK+'.'+RUN2D+'.latest_'+str(seli)+'.fit', unpack=True)
			if (funit=='flam'):
				coeff_fast   = 1.0
			if (funit=='fnu'):
				coeff_fast   = 10.**(-9.) * 3.34*ll**2.
			## best-fit model
			ax.plot(ll,ff*coeff_fast,label='FAST Best-fit model (logM='+'%.2f'%lmass[seli]+')')
			## grzw1w2 photometry
			ax.errorbar(lc,  np.array([f_g[seli],f_r[seli],f_z[seli],f_w1[seli],f_w2[seli]])*coeff_decals,
						yerr=np.array([e_g[seli],e_r[seli],e_z[seli],e_w1[seli],e_w2[seli]])*coeff_decals,
						fmt=None,c='r',ecolor='r',elinewidth=2,capthick=2,label='DECaLS photometry')
			## observed spectrum (stored in '1E-17 erg/cm^2/s/Ang')
			## cf. http://www.sdss.org/dr12/tutorials/quicklook/
			spPlate = REDUXDIR+'/'+RUN2D+'/'+plate[seli]+'/spPlate-'+plate[seli]+'-'+mjd[seli]+'.fits'
			sphdu   = fits.open(spPlate)
			c0      = sphdu[0].header['coeff0']
			c1      = sphdu[0].header['coeff1']
			spnpix  = sphdu[0].header['naxis1']
			spl     = 10.**(c0+c1*np.arange(spnpix))
			if (funit=='flam'):
				coeff_eboss = 1.e2
			if (funit=='fnu'):
				coeff_eboss = 10**(-7.) * 3.34*spl**2.
			spf     = sphdu[0].data[fiber[seli].astype('int')-1,:]
			spivar  = sphdu[1].data[fiber[seli].astype('int')-1,:]
			spf_smoothed = convolve(spf*(spivar>0),Box1DKernel(50))
			ax.plot(spl,convolve(spf*coeff_eboss*(spivar/coeff_eboss**2>0),Box1DKernel(50)),
						c='#e67300',linewidth=0.4,alpha=0.5,label='eBOSS Spectrum',zorder=2)
			## plot stuff
			ax.set_xlim(3e3,6e4)
			ax.set_xticks([3e3,4e3,5e3,6e3,7e3,8e3,9e3,1e4,2e4,3e4,4e4,5e4])
			ax.set_ylim(ylim)
			ax.set_xscale('log')
			ax.set_yscale('log')
			ax.set_xlabel('Observed lambda [A]')
			ax.set_ylabel(ylabel)
			ax.grid(True)
			ax.set_title(plate[seli]+'-'+mjd[seli]+'-'+fiber[seli].zfill(4)+' '+
				'(zspec='+'%.3f' % zspec[seli]+', zQ='+str(zQ[seli])+', zCont='+str(zCont[seli])+')')
			ax.legend(loc=legloc,fontsize=10)
			## emission lines
			ax.plot(np.zeros(2)+3727*(1+zspec[seli]),ylim,color='0.5',alpha=0.2,zorder=1) # OII  3727
			ax.plot(np.zeros(2)+4959*(1+zspec[seli]),ylim,color='0.5',alpha=0.2,zorder=1) # OIII 4959
			ax.plot(np.zeros(2)+5007*(1+zspec[seli]),ylim,color='0.5',alpha=0.2,zorder=1) # OIII 5007
			## observed magnitudes
			tmpfs = 10
			tmpc  = 'r'
			tmpc2 = 'b'
			if (funit=='flam'):
				tmpy0  = 7e2
				tmpy1  = 5e2
				tmpy2  = 4e2
			if (funit=='fnu'):
				tmpy0 = 8e1
				tmpy1 = 6e1
				tmpy2 = 5e1
			ax.text(3.1e3,tmpy0,'Band',color=tmpc,fontsize=tmpfs)
			ax.text(3.1e3,tmpy1,'ABmag',color=tmpc,fontsize=tmpfs)
			for j in xrange(5):
				filt = ['g','r','z','w1','w2'][j]
				ax.text(lc[j],tmpy0,filt, color=tmpc,fontsize=tmpfs,ha='center')
				exec 'ax.text(lc[j],tmpy1,"%.2f"%(22.5-2.5*np.log10(f_'+filt+'[seli])), color=tmpc,fontsize=tmpfs,ha="center")'
				exec 'ax.text(lc[j],tmpy2,"%.2f"%fast_'+filt+'[seli], color=tmpc2,fontsize=tmpfs,ha="center")'
			## FAST best-fit values
			tmpfs = 10
			tmpc  = 'b'
			if (funit=='flam'):
				tmpy = 0.75
			if (funit=='fnu'):
				tmpy = 0.45
			tmpdy = 0.05
			ax.text(0.80,tmpy-0*tmpdy,'FAST output',color=tmpc,fontsize=tmpfs,transform=ax.transAxes)
			ax.text(0.80,tmpy-1*tmpdy,'lmass= ' +'%.2f'%lmass[seli],color=tmpc,fontsize=tmpfs,transform=ax.transAxes)
			ax.text(0.80,tmpy-2*tmpdy,'lage  = '+'%.2f'%lage[seli], color=tmpc,fontsize=tmpfs,transform=ax.transAxes)
			ax.text(0.80,tmpy-3*tmpdy,'ltau  = '+'%.2f'%ltau[seli], color=tmpc,fontsize=tmpfs,transform=ax.transAxes)
			ax.text(0.80,tmpy-4*tmpdy,'Av    = '+'%.2f'%Av[seli],   color=tmpc,fontsize=tmpfs,transform=ax.transAxes)
			pdf.savefig()
			plt.close()

	return True



def get_sector(RA,DEC,CHUNK):
	# chunk infos
	tiledir = get_tile_dir(CHUNK)
	hdu     = fits.open(tiledir+'geometry-'+CHUNK+'.fits')
	geomsect= hdu[1].data['sector']
	hdu     = fits.open(tiledir+'sector-'+CHUNK+'.fits')
	data    = hdu[1].data
	sector  = data['sector']
	tiles   = data['tiles']
	ntiles  = data['ntiles']
	area    = data['area']
	# objects
	nobj       = len(RA)
	mng        = mangle.Mangle(tiledir+'geometry-'+CHUNK+'.ply')
	objpolyid  = mng.get_polyids(RA,DEC)
	objsector  = np.zeros(nobj,dtype='int')-1.
	objtiles   = np.zeros((nobj,tiles.shape[1]),dtype='int')-1
	objntiles  = np.zeros(nobj,dtype='int')
	objarea    = np.zeros(nobj)
	for i in xrange(len(sector)):
		ind_i = np.where(geomsect==sector[i])[0]
		#print 'i=',i,'sector[i]=',sector[i],'ind_i=',ind_i
		for j in ind_i:
			tmp          = (objpolyid==j)
			#print '\tj=',j,'tmp[tmp].shape=',tmp[tmp].shape
			objsector[tmp] = sector[i]
			objtiles[np.arange(len(objtiles))[tmp]] = tiles[i,:]
			objntiles[tmp] = ntiles[i]
			objarea[tmp]   = area[i]
	return objsector,objtiles,objntiles,objarea



def get_sector_TSR_SSR(RA,DEC,CHUNK,TSR_SSR_FILE,OUTROOT='None'):
	# get sector
	sector,a2,a3,a4 = get_sector(RA,DEC,CHUNK)
	# reading TSR,SSR
	all_sector,a2,a3,a4,all_TSR,all_SSR = np.loadtxt(TSR_SSR_FILE,unpack=True) 
	# get TSR,SSR
	TSR = np.zeros(len(RA))-1
	SSR = np.zeros(len(RA))-1
	for sect,tsr,ssr in zip(all_sector,all_TSR,all_SSR):
		tmp = (sector==sect)
		TSR[tmp] = tsr
		SSR[tmp] = ssr
    # plot if asked
	if (OUTROOT!='None'):
		#footprint+tiles
		footprint_ra,footprint_dec,xlim,ylim = get_footprint(CHUNK)
		parfile = get_tile_dir(CHUNK)+'/tiles-'+CHUNK+'.par'
		a1,tilera,tiledec = read_tiles(parfile)
		tile_rad          = get_tile_rad()
		for quant in ['TSR','SSR']:
			fig,ax = plt.subplots(figsize=(20,10))
			# building tile ellipses
			ells = [Ellipse(xy=np.array([x,y]),width=2.*tile_rad/np.cos(y/180.*np.pi),height=2.*tile_rad,angle=0)
					for x,y in zip(tilera,tiledec)]
			# eBOSS/ELG footprint
			ax.plot(footprint_ra,footprint_dec,c='k',zorder=10,linewidth=3)
			if (quant=='TSR'):
				c   = TSR
				clab= 'Sector Tiling Success Rate [%]'
			if (quant=='SSR'):
				c   = SSR
				clab = 'Sector Spectroscopic Success Rate [%]'
			clim = [0,100]
			tmp  = (c>0)
			SC   = ax.scatter(RA[tmp],DEC[tmp],
					c=100.*c[tmp],edgecolor='none',s=2,
					vmin=clim[0],vmax=clim[1],cmap=matplotlib.cm.jet)
			for e in ells:
				e.set_color('k')
				e.set_facecolor('none')
				e.set_linewidth(1)
				e.set_alpha(0.3)
				e.set_clip_box(ax.bbox)
				ax.add_patch(e)
			# colorbar
			cbar = plt.colorbar(SC)
			cbar.set_label(clab,fontsize=20)
			cbar.ax.tick_params(labelsize=20)
			cbar.set_clim(clim)
			# plot stuff
			ax.grid(True)
			ax.set_xlabel('R.A. [deg.]',fontsize=20)
			ax.set_ylabel('DEC. [deg.]',fontsize=20)
			ax.tick_params(axis='both', which='major', labelsize=20)
			ax.set_xlim(xlim)
			ax.set_ylim(ylim)
			plt.savefig(OUTROOT+'.'+quant+'.png',bbox_inches='tight')
			plt.close()
	return TSR,SSR


def get_hpsystweight(RA,DEC,CHUNK):
	hdu  = fits.open(os.getenv('MYSAS_DIR')+'/Catalogs/ELG_systweights_hp.fits')
	data = hdu[1].data
	hpind= data['hpind']
	exec 'hpsystweight = data["hpsystweight_'+CHUNK+'"]'
	phi       = RA*np.pi/180.
	theta     = (90.-DEC)*np.pi/180.
	pixind    = hp.ang2pix(256,theta,phi,nest=False)
	systweight= np.zeros(len(RA))-1.
	for ind in np.unique(pixind):
		systweight[(pixind==ind)] = hpsystweight[(hpind==ind)]
	return systweight


def dorands(CHUNK,LATESTFITS,TSR_SSR_FILE,INFITS,OUTFITS):
	# data sector/plate_SSR
	hdu         = fits.open(LATESTFITS)
	data        = hdu[1].data
	datasector  = data['sector']
	datatile    = data['TILE']
	dataplate   = data['PLATE']
	datacartid  = data['CARTID']
	dataxfocal  = data['XFOCAL']
	datayfocal  = data['YFOCAL']
	datafiberid = data['FIBERID']
	datarsn2    = data['plate_rSN2']
	datapssr    = data['plate_SSR']
	datasng     = data['SN_MEDIAN'][:,1]
	datahasfiber= data['hasfiber']
	datazreliable= data['Z_reliable']
	# input rands
	hdu        = fits.open(INFITS)
	randra     = hdu[1].data['ra']
	randdec    = hdu[1].data['dec']
	randsector = hdu[1].data['sector']
	# TSR/SSR
	randtsr,randssr= get_sector_TSR_SSR(randra,randdec,CHUNK,TSR_SSR_FILE)
	# all tiles
	parfile  = get_tile_dir(CHUNK)+'/tiles-'+CHUNK+'.par'
	alltile, alltilera, alltiledec = read_tiles(parfile)
	# rand sectors
	randtile = np.zeros(len(randra),dtype='int')
	randsector,a2,a3,a4 = get_sector(randra,randdec,CHUNK)
	# attributing tile to rands: loop on sectors
	for sector in np.unique(datasector):
		# rands in this sector
		selrand_i  = (randsector==sector)
		nrand_i    = len(randsector[selrand_i])
		randtile_i = np.zeros(nrand_i,dtype='int')
		# data in this sector
		seldata_i   = (datasector==sector)
		datasector_i= datasector[seldata_i]
		datatile_i  = datatile[seldata_i]
		ndata_i     = float(len(datasector_i))
		# stuff
		count       = 0
		tmpn        = 0
		tmptile     = 0.
		for tile in np.unique(datatile_i):
			tmp     = (datatile_i==tile)
			ndata_j = float(len(datatile_i[tmp]))
			nrand_j = int(ndata_j/ndata_i*nrand_i)
			randtile_i[count:count+nrand_j] = tile
			count += nrand_j
			# to keep track of the tile with the most data
			if (ndata_j>tmpn):	
				tmpn   = np.copy(ndata_j)
				tmptile= np.copy(tile)		
		# if some rands still have tile=0, we attribute tmptile
		randtile_i[(randtile_i==0)] = tmptile
		# we shuffle
		np.random.shuffle(randtile_i)
		# we store in randtile
		randtile[selrand_i] = randtile_i
	print str(len(randtile[randtile==0]))+'/'+str(len(randtile))+' rands have tile=0'
	# attributing pSSR/XFOCAL/YFOCAL to rands
	randplate  = np.zeros(len(randra),dtype='int')
	randcartid = np.zeros(len(randra),dtype='int')
	randhasfiber = np.zeros(len(randra),dtype='int')
	randsng    = np.zeros(len(randra))
	randrsn2   = np.zeros(len(randra))
	randpssr   = np.zeros(len(randra))
	randxfoc   = np.zeros(len(randra))
	randyfoc   = np.zeros(len(randra))
	norm_focal = get_tile_rad() / 325. # 1.49 deg to 325 mm
	randfiberid= np.zeros(len(randra),dtype='int')
	randzreliable = np.zeros(len(randra),dtype='bool')
	for tile in np.unique(datatile):
		selrand = (randtile==tile)
		seldata = (datatile==tile)
		# plate
		randplate[selrand] = dataplate[seldata][0]
		# cartid
		randcartid[selrand]= datacartid[seldata][0]
		# sng
		randsng[selrand]   = np.random.choice(datasng[seldata],size=len(randsng[selrand]),replace=True)
		# rSN2 , pSSR
		randrsn2[selrand]  = datarsn2[seldata][0]
		randpssr[selrand]  = datapssr[seldata][0]
		# XFOCAL/YFOCAL [method correct to better than 0.4 pixel]
		## tile ra,dec info
		tmp          = (alltile.astype('int')==tile)
		tra          = alltilera[tmp][0]
		tdec         = alltiledec[tmp][0]
		tSkyCoord    = SkyCoord(ra=tra*u.deg,dec=tdec*u.deg,frame='fk5')
		aframe       = tSkyCoord.skyoffset_frame()
		## centering the randra,randdec on the tile center
		randSkyCoord = SkyCoord(ra=randra[selrand]*u.deg,dec=randdec[selrand]*u.deg,frame='fk5')
		a            = randSkyCoord.transform_to(aframe)
		## converting deg to pixel
		randxfoc[selrand] = a.lon.value / norm_focal
		randyfoc[selrand] = a.lat.value / norm_focal
		## fiberid+hasfiber+zreliable
		collist = []
		collist.append(fits.Column(name='x', format='E',array=dataxfocal[seldata]))
		collist.append(fits.Column(name='y', format='E',array=datayfocal[seldata]))
		collist.append(fits.Column(name='f', format='K',array=datafiberid[seldata]))
		collist.append(fits.Column(name='hf',format='K',array=datahasfiber[seldata]))
		collist.append(fits.Column(name='zok',format='L',array=datazreliable[seldata]))
		tmphdu = fits.BinTableHDU.from_columns(fits.ColDefs(collist))
		tmphdu.writeto('tmp.data.fits_'+pid,clobber=True)
		collist = []
		collist.append(fits.Column(name='index',format='K',array=np.arange(len(randxfoc[selrand]))))
		collist.append(fits.Column(name='x',format='E',array=randxfoc[selrand]))
		collist.append(fits.Column(name='y',format='E',array=randyfoc[selrand]))
		tmphdu = fits.BinTableHDU.from_columns(fits.ColDefs(collist))
		tmphdu.writeto('tmp.rand.fits_'+pid,clobber=True)
		tmpstr = (STILTSCMD+' tmatch2 '+
					'matcher=2d find=best2 join=all2 params=100. '+
					'in1=tmp.data.fits_'+pid+' ifmt1=fits values1="x y" '+
					'in2=tmp.rand.fits_'+pid+' ifmt2=fits values2="x y" '+
					'out=tmp.fits_'+pid+' ofmt=fits '+
					"ocmd='sort index';")
		subprocess.call(tmpstr, shell=True)
		tmphdu                 = fits.open('tmp.fits_'+pid)
		randfiberid  [selrand] = tmphdu[1].data['f']
		randhasfiber [selrand] = tmphdu[1].data['hf']
		randzreliable[selrand] = tmphdu[1].data['zok']
		subprocess.call('rm tmp.data.fits_'+pid+' tmp.rand.fits_'+pid+' tmp.fits_'+pid, shell=True)
	# creating fits columns
	newcols = [c for c in hdu[1].columns]
	newcols.append(fits.Column(name='chunk',      format='A7',array = [CHUNK for x in randra]))
	newcols.append(fits.Column(name='sector_TSR', format='E',array=randtsr))
	newcols.append(fits.Column(name='sector_SSR', format='E',array=randssr))
	newcols.append(fits.Column(name='tile',       format='I',array=randtile))
	newcols.append(fits.Column(name='plate',      format='I',array=randplate))
	newcols.append(fits.Column(name='SN_MEDIAN_G',format='E',array=randsng))
	newcols.append(fits.Column(name='plate_rSN2', format='E',array=randrsn2))
	newcols.append(fits.Column(name='plate_SSR',  format='E',array=randpssr))
	newcols.append(fits.Column(name='XFOCAL',     format='E',array=randxfoc))
	newcols.append(fits.Column(name='YFOCAL',     format='E',array=randyfoc))
	newcols.append(fits.Column(name='FIBERID',    format='I',array=randfiberid))
	newcols.append(fits.Column(name='CARTID',     format='I',array=randcartid))
	newcols.append(fits.Column(name='hasfiber',   format='I',array=randhasfiber))
	newcols.append(fits.Column(name='Z_reliable', format='L',array=randzreliable))
	newhdu = fits.BinTableHDU.from_columns(newcols)
	# cutting on obs. plates
	tmp    = (randtsr>0)
	newhdu.data = newhdu.data[tmp]
	newhdu.writeto(OUTFITS,clobber=True)
	return True



def get_sector_plate(PLATE,CHUNK,RUN2D):
	# get TILEID
	hdu     = fits.open(os.getenv('EBOSS_ROOT')+'/spectro/redux/'+RUN2D+'/platelist.fits')
	p       = hdu[1].data['PLATE']
	t       = hdu[1].data['TILEID']
	TILEID  = t[p==PLATE][0]
	# read sectors
	hdu     = fits.open(os.getenv('EBOSS_ROOT')+'/ebosstilelist/trunk/outputs/'+CHUNK+'/sector-'+CHUNK+'.fits')
	s       = hdu[1].data['SECTOR']
	t       = hdu[1].data['TILES']
	keep    = np.zeros(len(s),dtype='bool')
	for i in xrange(4):
		tmp      = (t[:,i]==TILEID)
		keep[tmp] = True
	return s[keep]



def final_merge(CHUNKLIST,RUN2D,OUTDATAFITS,OUTRANDFITS,DATACOLNAMES=None,RANDCOLNAMES=None):
	# concatenating the latest rand files [keeping RANDCOLNAMES columns]
	INFITSLIST = []
	for i in xrange(len(CHUNKLIST)):
		CHUNK = CHUNKLIST[i]
		INFITSLIST.append(os.getenv('MYELG_DIR')+'/'+CHUNK+'/cats/'+CHUNK+'.'+RUN2D+'.latest.rands.fits')
	a = append_fits(INFITSLIST,OUTRANDFITS,COLNAMES=RANDCOLNAMES)
	# concatenating the latest data files [keeping DATACOLNAMES columns]
	INFITSLIST = []
	for i in xrange(len(CHUNKLIST)):
		CHUNK = CHUNKLIST[i]
		INFITSLIST.append(os.getenv('MYELG_DIR')+'/'+CHUNK+'/cats/'+CHUNK+'.'+RUN2D+'.latest.fits')
	a = append_fits(INFITSLIST,TMPDIR+'/tmp.dataspec.fits_'+pid,COLNAMES=DATACOLNAMES)

	# reading TS
	tshdu      = fits.open(os.getenv('MYELG_DIR')+'/ELG_TS_master.fits')
	tsdata     = tshdu[1].data
	tscol      = tshdu[1].columns
	tscolnames = np.array(tscol.names)
	# reading spectro.
	sphdu      = fits.open(TMPDIR+'/tmp.dataspec.fits_'+pid)
	spdata     = sphdu[1].data
	spcol      = sphdu[1].columns
	spcolnames = np.array(spcol.names)
	# we select TS objects *not* in spectro.
	tsuniq     = np.array([d+'-'+b+'-'+str(o) for d,b,o in 
						zip(tsdata['decals_dr'],tsdata['brickname'],tsdata['decals_objid'])],
						dtype='object')
	spuniq     = np.array([d+'-'+b+'-'+str(o) for d,b,o in 
						zip(spdata['decals_dr'],spdata['brickname'],spdata['decals_objid'])],
						dtype='object')
	tsinsp     = np.in1d(tsuniq,spuniq) # objects in spectro.
	tsdata     = tsdata[~tsinsp]
	nts        = tsdata.shape
	# we duplicate the spcol structure [DATACOLNAMES columns] for the ~tsinsp subsample
	# - if the column is in TS, we take the TS values
	# - if not, we add zeros
	newcollist = []
	for c in spcol:
		if (len(tscolnames[tscolnames==c.name])==0):
			newcollist.append(fits.Column(name=c.name,format=c.format,array=np.zeros(nts)))
		else:
			newcollist.append(fits.Column(name=c.name,format=c.format,array=tsdata[c.name]))
	newcols= fits.ColDefs(newcollist)
	# we create a table with spcol columns and ts data
	newhdu = fits.BinTableHDU.from_columns(newcols)
	newhdu.writeto(TMPDIR+'/tmp.dataphot.fits_'+pid,clobber=True)
	# we append dataspectro and dataphot
	a = append_fits([TMPDIR+'/tmp.dataspec.fits_'+pid,TMPDIR+'/tmp.dataphot.fits_'+pid],OUTDATAFITS,COLNAMES=DATACOLNAMES)
	# cleaning
	subprocess.call('rm '+TMPDIR+'/tmp.dataspec.fits_'+pid+' '+TMPDIR+'/tmp.dataphot.fits_'+pid, shell=True)

