import numpy as np
import glob
import os
import zipfile
import urllib.request
import json
from astropy.io import fits
from astropy.utils.data import download_file
import astropy.units as u
from astropy import wcs
from astropy.wcs import WCS
from astropy.visualization import ImageNormalize, ManualInterval, LogStretch, LinearStretch, AsinhStretch
from jwst.pipeline import Detector1Pipeline
from jwst.pipeline import Spec2Pipeline
from jwst.pipeline import Spec3Pipeline
from jwst.associations import asn_from_list as afl
from jwst.associations.lib.rules_level2_base import DMSLevel2bBase
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base
from jwst.datamodels import dqflags
# individual steps
#from jwst.assign_wcs import AssignWcsStep
#from jwst.assign_wcs import nirspec
#from jwst.background import BackgroundStep
#from jwst.imprint import ImprintStep
#from jwst.msaflagopen import MSAFlagOpenStep
#from jwst.extract_2d import Extract2dStep
#from jwst.srctype import SourceTypeStep
#from jwst.wavecorr import WavecorrStep
#from jwst.flatfield import FlatFieldStep
#from jwst.pathloss import PathLossStep
#from jwst.photom import PhotomStep
#from jwst.cube_build import CubeBuildStep
#from jwst.extract_1d import Extract1dStep
# data models
from jwst import datamodels

pipestates = [0,1,1] # Corresponding to stages 1,2,3a and 3c. If you want to 
                     # use some of them, put 1, If not, write 0. For NIRSpec 
                     # simulations they start in stage 2 so pipestates[0] = 0

print('Loaded libraries. Starting stage 0!')
print(' ')

# ============================ FUNCTIONS ================================== #

def writel2asn(files,asnfile,**kwargs):
	asn = afl.asn_from_list(files,rule=DMSLevel2bBase,product_name='Level2')
	if ('bg' in kwargs):
		for bgfile in kwargs['bg']:
			asn['products'][0]['members'].append({'expname': bgfile,  \
                                         'exptype':'background'})
	_, serialized = asn.dump()
	with open(asnfile, 'w') as outfile:
		outfile.write(serialized)

def writel3asn(files,asnfile,**kwargs):
	asn = afl.asn_from_list(files,rule=DMS_Level3_Base,product_name= \
                         'Level3_nirspec_ifusim_comb_1234')
	if ('bg' in kwargs):
		for bgfile in kwargs['bg']:
			asn['products'][0]['members'].append({'expname': bgfile, \
                                         'exptype':'background'})
	_, serialized = asn.dump()
	with open(asnfile, 'w') as outfile:
		outfile.write(serialized)
        
# ========================= DIRECTORIES =================================== #

main_dir = './'
output_dir = './nirspec_files/'
stage1_dir = './nirspec_files/stage1/'
stage2_dir = './nirspec_files/stage2/'
stage3_dir = './nirspec_files/stage3/'

if not os.path.exists(output_dir):
    os.makedirs(output_dir)    
if not os.path.exists(stage1_dir):
    os.makedirs(stage1_dir)
if not os.path.exists(stage2_dir):
    os.makedirs(stage2_dir)
if not os.path.exists(stage3_dir):
    os.makedirs(stage3_dir)    

sstring = 'det*exp1.fits'    # I dont know the name of the stage0 files for
                             # NIRSpec. look for them :)   
                             
filesdet1 = sorted(glob.glob(sstring))
print('there are ' + str(len(filesdet1))+ ' files to be processed in stage 1.')
print(' ')    

# ======================== PIPELINE STAGE 1 =============================== #

if pipestates[0] == 1:
	det1 = Detector1Pipeline()  
	det1.output_dir = output_dir   #careful with this
	det1.dq_init.skip = False
	det1.saturation.skip = False
	det1.superbias.skip = False
	det1.refpix.skip = True
	det1.linearity.skip = False
	det1.persistence.skip = False
	det1.dark_current.skip = False
	det1.jump.skip = False
	det1.save_results = True
	det1(filesdet1)
print('finished stage 1!')
print(' ')

# ======================= JSON FILE 1 ===================================== #

os.chdir(output_dir)
sstring = '*rate.fits'
calfiles = np.array(sorted(glob.glob(sstring)))
writel2asn(calfiles,'l2_asn.json')
print('there are ' + str(len(calfiles))+ ' files to be processed in stage 2.')
print(' ')

lines = []
with open(r"l2_asn.json", 'r') as fp:
    lines = fp.readlines()
with open(r"l2_asn.json", 'w') as fp:
    for number, line in enumerate(lines):
        if number not in [3,4,5,7,8]:  # Removing some not needed lines
            fp.write(line)

asn_file = "l2_asn.json"
with open(asn_file) as f_obj:
    asn_data = json.load(f_obj)

# ======================== PIPELINE STAGE 2 =============================== #

if pipestates[1] == 1:
	spec2 = Spec2Pipeline()
	spec2.output_dir = main_dir
	spec2.assign_wcs.skip = False
	spec2.bkg_subtract.skip = False
	spec2.imprint_subtract.skip = False
	spec2.srctype.skip = False
	spec2.msa_flagging.skip = False
	spec2.flat_field.skip = False
	spec2.pathloss.skip = False
	spec2.photom.skip = False
	spec2.cube_build.skip = True
	spec2.extract_1d.skip = True
	spec2.save_results = True
	result = spec2(asn_file)
print('finished stage 2!')
print(' ')

# ======================= JSON FILE 2 ===================================== #

sstring = '*cal.fits'
calfiles = np.array(sorted(glob.glob(sstring)))
writel3asn(calfiles,'l3_asn.json')
print('there are ' + str(len(calfiles))+ ' files to be processed in stage 3.')
print(' ')

lines = []
with open(r"l3_asn.json", 'r') as fp:
    lines = fp.readlines()
with open(r"l3_asn.json", 'w') as fp:
    for number, line in enumerate(lines):
        if number not in [3,4,5,7,9]:  # Removing some not needed lines
            fp.write(line)

asn_file = "l3_asn.json"
with open(asn_file) as f_obj:
    asn_data = json.load(f_obj)

# ======================== PIPELINE STAGE 3 =============================== #

if pipestates[2] == 1:
    spec3 = Spec3Pipeline()
    spec3.save_results = True
    spec3.output_dir = main_dir
    spec3.outlier_detection.skip = True # skip this step for now, because the
                                        # simulations do not include outliers
    result = spec3(asn_file)
print('finished stage 3!')
print(' ')

# =========================== REUBICATE FILES ============================= #

print('Reubicating files')
print(' ')

files = sorted(glob.glob('*.fits'))

for i in range(len(files)):
    if files[i].endswith('rate.fits'):
        os.rename(files[i], 'stage1/'+files[i])
    elif files[i].endswith('cal.fits'):
        os.rename(files[i], 'stage2/'+files[i])
    elif files[i].startswith('Level3'):
        os.rename(files[i], 'stage3/'+files[i])
        
print('Data finished !!')
print(' ')
