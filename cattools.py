import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u


def cut(cat, w):
	''' Trim catalog columns using a boolean array
	Copied from Julian Bautista's code and AJR doesn't get how it works
	'''

	size = cat.size
	for f in cat.__dict__.items():
		if hasattr(f[1], 'size') and f[1].size % size == 0:
			self.__dict__[f[0]] = f[1][w]
	self.size = self.RA.size

def cutrange(cat,col,min,max):
	''' 
	cut within range given catalog, column, min, and max
	seems to be super slow and didn't work
	'''

	size = cat.size
	mask = np.ones(cat.size,dtype=bool)
	for i in range(0,cat.size):
		if cat[col][i] < min or cat[col][i] > max:
			mask[i] = False
	cat = cat[mask,...]
	return cat		

def fiber_collision(cat, apply_noz=1,type='ELG'):

	'''
	Pulled out of Julian Bautista's LRG code to be adapted to generic catalogs
	'''

	print ' '
	print '========================================'
	print '==== Finding fiber-collision pairs  ===='
	print '========================================'
	print ' '

	if type == 'ELG':
		sectors = cat['sector']
	else:
		sectors = cat['SECTOR']
	unique_sectors = np.unique(sectors)
  
	cnt_legacy_removed = 0
 
	#-- making copies for testing purposes
	all_z = np.copy(cat['Z'])
	#all_imatch = N.copy(self.IMATCH)
	if type == 'ELG':
		all_imatch = np.zeros(cat.size, dtype=int)
		for i in range(0,cat.size):
			if cat[i]['hasfiber'] == 1:
				all_imatch[i] = 1
		all_dups = np.copy(cat['isdupl'])	
		all_ids = np.copy(cat['decals_objid'])	
		#np.copy(cat['hasfiber']) #ELG equivalent?
	all_weight_cp = np.ones(cat.size)
	#all_weight_noz = N.ones(self.size)

	all_pair_in_sector_good = np.zeros(cat.size, dtype=int)
	all_pair_in_sector_tot = np.zeros(cat.size, dtype=int)
	all_gal_in_sector_good = np.zeros(cat.size, dtype=int)
	all_gal_in_sector_bad = np.zeros(cat.size, dtype=int)
	all_cp_pair_over_poss = np.zeros(cat.size)
	all_cp_gal_over_poss = np.zeros(cat.size)
	all_compboss = np.zeros(cat.size) #completeness counting close pairs as observations
	all_compreal= np.zeros(cat.size) #completeness just n_withfib/n_tot
	all_cp_match = np.zeros(cat.size) #this will hold the index of the colliding galaxy
	all_index = np.zeros(cat.size) #this will hold the index of the target galaxies
	for i in range(0,cat.size):
		all_index[i] = i



	print 'Solving fiber collisions in ', unique_sectors.size, 'sectors' 
	for i, sect in enumerate(unique_sectors):
		if type == 'ELG':
			w = (sectors == sect) & (all_dups == False) #& all_imatch == True#& (all_imatch != 2) & (self.vetobits == 0)
			#w2 = (sectors == sect) & (all_dups == 1)
			w2 = (sectors == sect) & (all_dups == 0) & (cat['sector_TSR'] != 0)
		plates = np.unique(cat['PLATE'][w])

		z = all_z[w]
		index = all_ids[w]
		cp_match = all_cp_match[w]
		imatch = all_imatch[w]
		weight_cp = np.zeros(z.size)
		nspec = sum(imatch>0)
		compreal = sum(imatch)/(float(np.sum(w)))
		#if np.mean(cat['sector_TSR'][w]) == 0:
		
		all_compreal[w] = compreal
	   
		#if (imatch==0).all():
		#	continue
		if abs(1.-compreal/np.mean(cat['sector_TSR'][w2])) > .01:# == 0:
			print compreal,np.mean(cat['sector_TSR'][w2])#,np.sum(w2)
			print '  %d of %d'%(i, unique_sectors.size), \
				  '\tSector#:', sect,\
				  '\tTargets:', np.sum(w), \
				  '\tnplates:', len(plates),\
				  '\tnspec:', sum(imatch>0)#,\
			  #'\tnduplicates:', np.sum(w2))

		#-- call to spheregroup for galaxies (with new redshifts) in this sector
		#-- with linking length of 62''.
		if type == 'ELG':
			pairs = spheregroup(cat['ra'][w], cat['dec'][w], angle=62./3600)
		else:
			pairs = spheregroup(cat['RA'][w], cat['DEC'][w], angle=62./3600)
	   
   

		
		pair_in_sector_good = 0
		pair_in_sector_tot = 0
		gal_in_sector_good = 0
		gal_in_sector_bad = 0
		
		for pair in pairs:
			#for ELGs, imatch is True or False, 0 or 1
			#targets identified as close pairs get imatch = 3
			imatch1 = imatch[pair[0]]
			imatch2 = imatch[pair[1]]
			#print imatch1, imatch2

			if imatch1 != 0 or imatch2 != 0:
				pair_in_sector_tot += 1
				
			#-- if collision is resolved
			if imatch1 != 0 and imatch1 != 3 and imatch2 != 0 and imatch2 != 3:
				gal_in_sector_good += 1
				pair_in_sector_good += 1
				#print 'This collision is solved through many plates'

			#-- solving collision
			elif imatch1 == 0 and imatch2 !=0 and imatch2 != 3:
				z[pair[0]] = z[pair[1]]
				cp_match[pair[0]] = index[pair[1]]
				imatch[pair[0]] = 3
				imatch1 = 3
				weight_cp[pair[1]] += 1
				gal_in_sector_bad += 1
				#print 'Solving collision: putting right redshift on the left'

			elif imatch2 == 0 and imatch1 !=0 and imatch1 != 3:
				z[pair[1]] = z[pair[0]]
				cp_match[pair[1]] = index[pair[0]]
				imatch[pair[1]] = 3
				imatch2 = 3
				weight_cp[pair[0]] += 1
				gal_in_sector_bad += 1
				#print 'Solving collision: putting left redshift on the right'
	

		if pair_in_sector_tot > 0: 
			cp_pair_over_poss = pair_in_sector_good*1./pair_in_sector_tot
		else: 
			cp_pair_over_poss = 0
		if gal_in_sector_good+gal_in_sector_bad > 0: 
			cp_gal_over_poss = gal_in_sector_good*1. / \
							  (gal_in_sector_good + gal_in_sector_bad)
		else:
			cp_gal_over_poss = 0
		if plates.size == 1:
			cp_pair_over_poss = 0.
			cp_gal_over_poss = 0.

		#if np.mean(cat['sector_TSR'][w]) == 0:
		#print min(weight_cp),max(weight_cp)
		all_z[w] = z
		all_cp_match[w] = cp_match
		all_imatch[w] = imatch
		all_weight_cp[w] += weight_cp
		compboss = (nspec+np.sum(weight_cp))/float(np.sum(w))
		all_compboss[w] = compboss
		if abs(1.-compreal/np.mean(cat['sector_TSR'][w2])) > .01:
			print '  Close-pairs:', len(pairs)
			print '     pairs numbers :', pair_in_sector_good, \
						pair_in_sector_tot, cp_pair_over_poss 
			print '     galaxy numbers:', gal_in_sector_good,  \
						gal_in_sector_bad, cp_gal_over_poss

			print compboss,compreal
		all_pair_in_sector_good[w] = pair_in_sector_good
		all_pair_in_sector_tot[w] = pair_in_sector_tot
		all_gal_in_sector_good[w] = gal_in_sector_good
		all_gal_in_sector_bad[w] = gal_in_sector_bad
		all_cp_pair_over_poss[w] = cp_pair_over_poss
		all_cp_gal_over_poss[w] = cp_gal_over_poss

	
		#-- now look for LEGACY close pairs in this sector.
		#-- remove some according to probability cp_gal_over_poss
		if type != 'ELG': 
			ww = (sectors == sect) & (self.vetobits==0)
			imatch = all_imatch[ww]
			z = all_z[ww]
			weight_cp = all_weight_cp[ww]
			weight_noz = all_weight_noz[ww]
			ra = self.RA[ww]
			dec = self.DEC[ww]

			#-- if there is any legacy redshifts do something
			if (all_imatch[ww] == 2).any():
				pairs = Utils.spheregroup(ra, dec, angle=62./3600)
			
				for pair in pairs:
					imatch1 = imatch[pair[0]]
					imatch2 = imatch[pair[1]]

					if  (imatch1 == 2 or imatch2 == 2) and \
						 imatch1 != 3 and imatch2 != 3 and \
						 imatch1 != 0 and imatch2 != 0 and \
						 imatch1 != 8 and imatch2 != 8:
						r1 = N.random.rand()
						if r1 < cp_gal_over_poss: continue
						r2 = N.random.rand() 
						if r2 > 0.5: 
							  imatch[pair[0]] = 8
							  imatch1 = 8 # this needs to update too so the same galaxy doesn't get updated twice.
							  z[pair[0]] = z[pair[1]]
							  weight_cp[pair[1]] += weight_cp[pair[0]]
						else:
							  imatch[pair[1]] = 8
							  imatch2 = 8 # this needs to update too so the same galaxy doesn't get updated twice.
							  z[pair[1]] = z[pair[0]]
							  weight_cp[pair[0]] += weight_cp[pair[1]]
					
						cnt_legacy_removed += 1L
				#-- end of pair loop            
			#-- if any imatch==2 in this sector
			all_z[ww] = z
			all_imatch[ww] = imatch
			all_weight_cp[ww] = weight_cp


			##-- Here goes the redshift failure correction ### To be changed hopefully
			if apply_noz == 0:
				continue

			#-- find failed targets in this sector
			wz = (imatch == 7)     
			if sum(wz) == 0:
				continue
	  
			#-- find closest neighbor with imatch = 1, 4 (star) or 9 (wrong class)
			wgood = N.where( (imatch == 1) | (imatch == 4) | (imatch == 9))[0]
			if wgood.size == 0:
				print '  No good redshifts found to do redshift failure correction'
				continue
 
			print '     Failures:', sum(wz), '\t Number of targets available to apply correction', wgood.size
			id1, id2, dist = Utils.spherematch(ra[wz], dec[wz], ra[wgood], dec[wgood], \
										 angle=2.)
	   
			if dist_root!='': 
				for j in range(dist.size): 
					print>>fout, sect, dist[j]

			#-- attribute same redshift and transfer weight_cp 
			imatch[wz] = 5
			z[wz] = z[wgood[id2]]
			for i in range(id2.size):
				weight_noz[wgood[id2[i]]] += weight_cp[wz][id1[i]] 
 
			all_z[ww] = z
			all_imatch[ww] = imatch
			all_weight_noz[ww] = weight_noz
   
	##-- end loop over regions
	return all_weight_cp,all_imatch,all_cp_match,all_compreal,all_compboss
	#self.Z = all_z
	#self.IMATCH = all_imatch
	#self.WEIGHT_CP = all_weight_cp
	#self.WEIGHT_NOZ = all_weight_noz       
	#self.PAIRS_GOOD =     all_pair_in_sector_good
	#self.PAIRS_TOT =  all_pair_in_sector_tot
	#self.GALS_GOOD =     all_gal_in_sector_good
	#self.GALS_BAD =     all_gal_in_sector_bad
	#self.COMP_CP =    all_cp_pair_over_poss
	#self.COMP_GAL=    all_cp_gal_over_poss


	#if dist_root != '':
	#	fout.close()   

#class Utils:
	'''
	taken from Julian Bautista's LRG code
	'''
    #@staticmethod
def spherematch(ra1, dec1, ra2, dec2, angle=1./3600):
	''' Implementation of spherematch using SkyCoord from astropy
	Inputs
	------
	ra1, dec1, ra2, dec2: arrays 
		Coordinates to be matched
	angle: float
		Angle in degrees for defining maximum separation

	Returns
	-------
	idx1, idx2: arrays
		Index of the matching 
	distance: array
		Distance between matches
	'''

	c1 = SkyCoord(ra=ra1*u.degree, dec=dec1*u.degree)     
	c2 = SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree)
	idx, d2d, d3d = c1.match_to_catalog_sky(c2)  
	w = d2d.value <= angle
	idx[~w] = -1
	
	idx1 = N.where(w)[0]
	idx2 = idx[idx>-1] 
	distance = d2d.value[w]

	return idx1, idx2, distance

    #@staticmethod
def spheregroup(ra, dec, angle=62./3600):
	''' Gets list of index for (ra, dec) pairs at distances smaller than angle'''

	c = SkyCoord(ra=ra*u.degree, dec=dec*u.degree)
	idc1, idc2, d2d, d3d = c.search_around_sky(c, angle*u.degree)

	#-- exclude pairs of same object twice
	#ww = (d2d>0.)
	ww = (idc1!=idc2) & (d2d>0.)
	
	#-- if no collisions return empty list
	if sum(ww) == 0: 
		return []

	i1 = idc1[ww]
	i2 = idc2[ww] 
	distance = d2d[ww]

	#-- removing duplicate pairs
	pairs = [ [ii1, ii2] for ii1, ii2 in zip(i1, i2) if ii1<ii2]
	return pairs

    #@staticmethod
def maskbits(value):
	''' Get mask bit values from integer '''

	if type(value)==list or type(value)==N.ndarray:
		return [Utils.maskbits(v) for v in value]
	else:
		return [pos for pos, char in enumerate(bin(value)[::-1]) if char == '1']

    #@staticmethod
    #-- healpix map rotation 
def rotate_map(m, coord=['C', 'G'], rot=None):
	"""
	STOLEN FROM THE GREAT SPIDER COLLABORATION, and simplified to non-polarized data.
	Rotate an input map from one coordinate system to another or to place a
	particular point at centre in rotated map. 

	e.g. m = rotate_map(m, rot=[phi,90.-theta,0.])

	takes point at original theta, phi to new coord ra=dec=0

	Arguments
	---------
	m : array_like
		A single map
	coord : list of two coordinates, optional.
		Coordinates to rotate between.  Default: ['C', 'G']
	rot : scalar or sequence, optional
		Describe the rotation to apply.
		In the form (lon, lat, psi) (unit: degrees) : the point at
		longitude lon and latitude lat will be at the center of the rotated
		map. An additional rotation of angle psi around this direction is applied
	"""

	res = hp.get_nside(m)
	mi = m
	if rot is None:
		R = hp.Rotator(coord=coord, inv=False)
	else:
		R = hp.Rotator(rot=rot, inv=False)

	# rotate new coordinate system to original coordinates
	theta, phi = hp.pix2ang(res, N.arange(len(mi)))
	mtheta, mphi = R(theta, phi)
	mr = hp.get_interp_val(mi, mtheta, mphi)

	return mr

