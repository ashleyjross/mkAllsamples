'''
Copied from Julian Bautista to be used as a general script to execute eBOSS reconstruction
'''

#from ebosscat import Catalog
from recon import Recon
import argparse
from astropy.io import fits
import numpy as np
from cattools import *

dir = '/Users/ashleyross/fitsfiles/' #directory on Ashley's computer where catalogs are


'''
argument parser allows script to be run like 
> python do_recon.py -d 'datafile' -r 'randomfile' ...
'''


parser = argparse.ArgumentParser()

parser.add_argument('-reg', '--region', help='SGC or NGC',default='SGC')
parser.add_argument('-v', '--version', help='version',default='test')
parser.add_argument('-o', '--output', help='Output catalogs root name',default='rec')
parser.add_argument('-t', '--type', help='Target class',default='ELG')
parser.add_argument('--nthreads', \
    help='Number of threads', type=int, default=1)
parser.add_argument('--niter', \
    help='Number of iterations', type=int, default=3)
parser.add_argument('--nbins', \
    help='Number of bins for FFTs', type=int, default=512)
parser.add_argument('--padding', default=200., \
    help='Size in Mpc/h of the zero padding region', type=float)
parser.add_argument('--zmin', help='Redshift lower bound', type=float,default=.6)
parser.add_argument('--zmax', help='Redshift upper bound', type=float,default=1.1)
parser.add_argument('--smooth', help='Smoothing scale in Mpc/h', \
    type=float, default=15.)
#parser.add_argument('--bias', \
#    help='Estimate of the bias of the sample', type=float, required=True)
#parser.add_argument('--f', \
#    help='Estimate of the growth rate', type=float, required=True, default=0.817)
args = parser.parse_args()
print args

'''
First thing that is needed is data and randoms with ra,dec,z,weight columns
'''

if args.type == 'ELG':
	from mksimpELG import *
	mkgalELGsimp(args.region,zmin=args.zmin,zmax=args.zmax,vo=args.version)
	mkranELGsimp(args.region,vo=args.version)
	cat = fits.open(dir+args.type+args.region+args.version+'.dat.fits')[1].data
	ran = fits.open(dir+args.type+args.region+args.version+'.ran.fits')[1].data

	cat.weight = cat.WEIGHT_SYSTOT
	ran.weight = ran.WEIGHT_SYSTOT
	bias = 1.4 #get good value for this
	f = .82 #eventually calculate from Cosmo

nbins=args.nbins
nthreads=args.nthreads
padding =args.padding 
zmin=args.zmin
zmax=args.zmax
smooth=args.smooth
#bias = args.bias
#f = args.f
opt_box = 1 #optimize box dimensions

#-- selecting galaxies
#w = (cat.IMATCH==1)|(cat.IMATCH==2)|(cat.IMATCH==101)|(cat.IMATCH==102)
#w = w & ((cat.Z>=zmin)&(cat.Z<=zmax))
#cat.cut(w)
#wr = ((ran.Z>=zmin)&(ran.Z<=zmax))
#ran.cut(wr)

#cat = cutrange(cat,'Z',zmin,zmax)
#ran = cutrange(ran,'Z',zmin,zmax)

rec = Recon(cat, ran, nbins=nbins, smooth=smooth, f=f, bias=bias, \
            padding=padding, opt_box=opt_box, nthreads=nthreads)
for i in range(args.niter):
    rec.iterate(i)
rec.apply_shifts()
rec.summary()

cat.RA, cat.DEC, cat.Z = rec.get_new_radecz(rec.cat) 
ran.RA, ran.DEC, ran.Z = rec.get_new_radecz(rec.ran) 

cols = []
RAc = fits.Column(name='RA',format='D', array=cat.RA)
cols.append(RAc)
DECc = fits.Column(name='DEC',format='D', array=cat.DEC)
cols.append(DECc)
Zc = fits.Column(name='Z',format='D', array=cat.Z)
cols.append(Zc)
fkpc = fits.Column(name='WEIGHT_FKP',format='D', array=cat.WEIGHT_FKP)
cols.append(fkpc)
sysc = fits.Column(name='WEIGHT_SYSTOT',format='D', array=cat.WEIGHT_SYSTOT)
cols.append(sysc)

hdulist = fits.BinTableHDU.from_columns(cols)
header = hdulist.header
hdulist.writeto(dir+args.type+args.region+args.version+args.output+'.dat.fits', overwrite=True)

cols = []
RAc = fits.Column(name='RA',format='D', array=ran.RA)
cols.append(RAc)
DECc = fits.Column(name='DEC',format='D', array=ran.DEC)
cols.append(DECc)
Zc = fits.Column(name='Z',format='D', array=ran.Z)
cols.append(Zc)
fkpc = fits.Column(name='WEIGHT_FKP',format='D', array=ran.WEIGHT_FKP)
cols.append(fkpc)
sysc = fits.Column(name='WEIGHT_SYSTOT',format='D', array=ran.WEIGHT_SYSTOT)
cols.append(sysc)

hdulist = fits.BinTableHDU.from_columns(cols)
header = hdulist.header
hdulist.writeto(dir+args.type+args.region+args.version+args.output+'.ran.fits', overwrite=True)


#cat.export(args.output+'.dat.fits')
#ran.export(args.output+'.ran.fits')


