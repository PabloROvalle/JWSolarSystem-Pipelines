import os
import glob
import shutil
import numpy as np
from astropy.io import fits
from astropy.utils.data import download_file
from jwst.associations import asn_from_list as afl
from astropy.visualization import ImageNormalize, ManualInterval, LogStretch
import json
from jwst.pipeline import calwebb_detector1
from jwst import datamodels
from jwst.associations.lib.rules_level2_base import DMSLevel2bBase
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base

import jwst
print(jwst.__version__)

pipestates = [1,1,0] 

def writel2asn(files,asnfile,**kwargs):
	asn = afl.asn_from_list(files,rule=DMSLevel2bBase,product_name='Level2')
	if ('bg' in kwargs):
		for bgfile in kwargs['bg']:
			asn['products'][0]['members'].append({'expname': bgfile, 'exptype':'background'})
	_, serialized = asn.dump()
	with open(asnfile, 'w') as outfile:
		outfile.write(serialized)

def writel3asn(files,asnfile,**kwargs):
	asn = afl.asn_from_list(files,rule=DMS_Level3_Base,product_name='Level3')
	if ('bg' in kwargs):
		for bgfile in kwargs['bg']:
			asn['products'][0]['members'].append({'expname': bgfile, 'exptype':'background'})
	_, serialized = asn.dump()
	with open(asnfile, 'w') as outfile:
		outfile.write(serialized)
        
main_dir = './'
output_dir = './nircam_files/'     # You have to paste your uncals and trapsfilled (if you have) inside this repo
stage1_dir = './nircam_files/stage1/'
stage2_dir = './nircam_files/stage2/'
stage3_dir = './nircam_files/stage3/'

if not os.path.exists(output_dir):
    os.makedirs(output_dir)    
if not os.path.exists(stage1_dir):
    os.makedirs(stage1_dir)
if not os.path.exists(stage2_dir):
    os.makedirs(stage2_dir)
if not os.path.exists(stage3_dir):
    os.makedirs(stage3_dir)

os.chdir(output_dir)

sstring = '*uncal.fits'   
                               
filesdet1 = sorted(glob.glob(sstring))
print('there are ' + str(len(filesdet1))+ ' files to be processed in stage 1.')
print(' ')

# ======================== PIPELINE STAGE 1 =============================== #

if pipestates[0] == 1:
	det1 =  jwst.pipeline.Detector1Pipeline() 
	det1.output_dir = main_dir 
    
	det1.dq_init.skip = False           #You can modify- https://jwst-pipeline.readthedocs.io/en/latest/jwst/dq_init/
	det1.saturation.skip = True        #You can modify- https://jwst-pipeline.readthedocs.io/en/latest/jwst/saturation/
	det1.superbias.skip = False         #You can modify- https://jwst-pipeline.readthedocs.io/en/latest/jwst/su[erbias/
	det1.refpix.skip = False             #You can modify- https://jwst-pipeline.readthedocs.io/en/latest/jwst/refpix/
	det1.linearity.skip = False         #You can modify- https://jwst-pipeline.readthedocs.io/en/latest/jwst/linearity/
	det1.persistence.skip = False       #You can modify- https://jwst-pipeline.readthedocs.io/en/latest/jwst/persistence/
	det1.dark_current.skip = False      #You can modify- https://jwst-pipeline.readthedocs.io/en/latest/jwst/dark_current/
	det1.jump.skip = False              #You can modify- https://jwst-pipeline.readthedocs.io/en/latest/jwst/jump/
    
	det1.save_results = True
	det1(filesdet1[0])      # Da fallo al leer la variable completa 
print('finished stage 1!')
print(' ')

sstring = '*rate.fits'
calfiles = np.array(sorted(glob.glob(sstring)))
writel2asn(calfiles,'l2_asn.json')
print('there are ' + str(len(calfiles)) + ' files to be processed in stage 2.')
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
	image2 = jwst.pipeline.Image2Pipeline()
	image2.output_dir = main_dir 
	image2.assign_wcs.skip = False        #You can modify- https://jwst-pipeline.readthedocs.io/en/latest/jwst/assign_wcs/
	image2.bkg_subtract.skip = False      #You can modify- https://jwst-pipeline.readthedocs.io/en/latest/jwst/bkg_subtract/
	image2.flat_field.skip = False        #You can modify- https://jwst-pipeline.readthedocs.io/en/latest/jwst/flat_field/
	image2.photom.skip = False            #You can modify- https://jwst-pipeline.readthedocs.io/en/latest/jwst/photom/
	image2.resample.skip = False          #You can modify- https://jwst-pipeline.readthedocs.io/en/latest/jwst/resample/
	image2.save_results = True
	result = image2(asn_file)
print('finished stage 2!')
print(' ')

sstring = '*_cal.fits'
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
    
asn_file = "l3_asn.json"
with open(asn_file) as f_obj:
    asn_data = json.load(f_obj)

# ======================== PIPELINE STAGE 3 =============================== #

if pipestates[2] == 1:
    image3 = jwst.pipeline.Image3Pipeline()
    image3.save_results = True
    image3.output_dir = output_dir
    
    # OUTLIER_DETECTION -- Arrate !!! Puede que arregle algo si lo pasaos a True (y hace skip de ese paso)
    
    image3.assign_mtwcs.skip = False       #You can modify- https://jwst-pipeline.readthedocs.io/en/latest/jwst/assign_mtwcs/
    image3.outlier_detection.skip = False   #You can modify- https://jwst-pipeline.readthedocs.io/en/latest/jwst/outlier_detection/
    image3.skymatch.skip = False           #You can modify- https://jwst-pipeline.readthedocs.io/en/latest/jwst/skymatch/
    image3.tweakreg.skip = False           #You can modify- https://jwst-pipeline.readthedocs.io/en/latest/jwst/tweakreg/
    image3.resample.skip = False           #You can modify- https://jwst-pipeline.readthedocs.io/en/latest/jwst/resample/
    image3.source_catalog.skip = False     #You can modify- https://jwst-pipeline.readthedocs.io/en/latest/jwst/source_catalog/
    
    result = image3(asn_file)
print('finished stage 3!')
print(' ')

print('Reubicating files')
print(' ')

files = sorted(glob.glob('*.fits'))

for i in range(len(files)):
    if files[i].endswith('rate.fits'):
        os.rename(files[i], 'stage1/'+files[i])
    elif files[i].endswith('rateints.fits'):
        os.rename(files[i], 'stage1/'+files[i])
    elif files[i].endswith('_cal.fits'):
        os.rename(files[i], 'stage2/'+files[i])
    elif files[i].endswith('_i2d.fits'):
        os.rename(files[i], 'stage2/'+files[i])
    elif files[i].startswith('Level3'):
        os.rename(files[i], 'stage3/'+files[i])
        
print('Data finished !!')
print(' ')
