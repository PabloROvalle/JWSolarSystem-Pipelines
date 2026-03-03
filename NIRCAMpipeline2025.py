import os
print('')
print('Running JWST NIRCam pipeline 2025.')
print('')
user_home_dir = os.path.expanduser('~')
sci_dir = os.path.join(user_home_dir, 'pipeline/nircam_files/')

# Science processing
dodet1 = True  # calwebb_detector1
doimage2 = True  # calwebb_image2
doimage3 = True  # calwebb_image3
remove_1group_uncal = True  # Remove uncal_1group folder at the end of the run.
unlimit_desat = 0   # 0=dont do it ;  1=do it only for 1 group  ; 2=do the pipeline and the 1 group desaturation 


print('############################################################')
print('')
print('PIPELINE CONFIGURATION')
print('')
if dodet1 == True:
    print('You will execute stage 1 in this run.')
elif doimage2 == True:
    print('You will execute stage 2 in this run.')
elif doimage3 == True:
    print('You will execute stage 3 in this run.')
elif dodet1 == True and doimage2 == False and doimage3 == False:
    print('You are not launching the pipeline.')
print('')

if unlimit_desat == 0:
    print('You are not using the unlimited desaturation method.')
elif unlimit_desat == 1:
    print('You are only using the unlimited desaturation method (1 single group).')
elif unlimit_desat == 2:
    print('You are using the unlimited desaturation method (1 single group) and the standard method (>2 groups).')
print('')
print('############################################################')
print('')


if (os.getenv('CRDS_PATH') is None):
    os.environ['CRDS_PATH'] = os.path.join(os.path.expanduser('~'), 'crds')
# Check whether the CRDS server URL has been set.  If not, set it.
if (os.getenv('CRDS_SERVER_URL') is None):
    os.environ['CRDS_SERVER_URL'] = 'https://jwst-crds.stsci.edu'

# Echo CRDS path in use
print(f"CRDS local filepath: {os.environ['CRDS_PATH']}")
print(f"CRDS file server: {os.environ['CRDS_SERVER_URL']}")
if os.getenv('CRDS_CONTEXT'):
    print(f"CRDS CONTEXT: {os.environ['CRDS_CONTEXT']}")
    

# Basic system utilities for interacting with files
# ----------------------General Imports------------------------------------
import glob
import time
from pathlib import Path
from astropy.io import fits
import numpy as np
from astroquery.mast import Observations
from astropy.table import Table
from astropy.coordinates import SkyCoord
import shutil
import jwst
import crds
from jwst.pipeline import Detector1Pipeline
from jwst.pipeline import Image2Pipeline
from jwst.pipeline import Image3Pipeline
from asdf import AsdfFile
from jwst import datamodels
from jwst.associations import asn_from_list  # Tools for creating association files
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base 

# Echo pipeline version and CRDS context in use
print(f"JWST Calibration Pipeline Version: {jwst.__version__}")
print(f"Using CRDS Context: {crds.get_context_name('jwst')}")


 # Define output subdirectories to keep science data products organized
uncal_dir = os.path.join(sci_dir, 'uncal')  # Uncalibrated pipeline inputs should be here
det1_dir = os.path.join(sci_dir, 'stage1')  # calwebb_detector1 pipeline outputs will go here
image2_dir = os.path.join(sci_dir, 'stage2')  # calwebb_spec2 pipeline outputs will go here
image3_dir = os.path.join(sci_dir, 'stage3')  # calwebb_spec3 pipeline outputs will go here
uncal_dir_1group = os.path.join(sci_dir, 'uncal_1group') 
det1_dir_1group = os.path.join(sci_dir, 'stage1_1group')
image2_dir_1group = os.path.join(sci_dir, 'stage2_1group') 
image3_dir_1group = os.path.join(sci_dir, 'stage3_1group')


# We need to check that the desired output directories exist, and if not
# create them
if not os.path.exists(det1_dir):
    os.makedirs(det1_dir)
if not os.path.exists(image2_dir):
    os.makedirs(image2_dir)
if not os.path.exists(image3_dir):
    os.makedirs(image3_dir)
if unlimit_desat==1 or unlimit_desat==2:
    if not os.path.exists(uncal_dir_1group):
            os.makedirs(uncal_dir_1group)
    if not os.path.exists(det1_dir_1group):
            os.makedirs(det1_dir_1group)
    if not os.path.exists(image2_dir_1group):
            os.makedirs(image2_dir_1group)
    if not os.path.exists(image3_dir_1group):
            os.makedirs(image3_dir_1group)  


for filename in glob.glob(os.path.join(uncal_dir, '*.*')):
    shutil.copy(filename, uncal_dir_1group)

print('')
if unlimit_desat==1 or unlimit_desat==2:
    sstring = os.path.join(uncal_dir_1group, 'jw*uncal.fits')
    uncal_files = np.array(sorted(glob.glob(sstring)))
    print(uncal_files)
    #uncal_files = select_ch_band_files(uncal_files, use_ch, use_band)
    print('Reducing uncalibrated files to 1 group.')
    for file in uncal_files:
        hdu = fits.open(file)
        hdr= hdu['SCI'].header
        data = hdu['SCI'].data
        first = data.shape[0]
        second = data.shape[1]
        third = data.shape[2]
        fourth = data.shape[3]
        if int(second) == 1:
            break
        else:
            print('max number of groups --> ' + str(second))
            k = 1
            print('you are creating a fits file with ' +str(k) + ' groups')
            newdata = np.zeros((first,k,third,fourth))
            for ii in range(k):
                newdata[:,ii,:,:] = data[:,ii,:,:]
            newdata = newdata.astype('uint16')
            hdu['SCI'].data = newdata
            hdu[0].header['NGROUPS'] = k
            print('IN PROCESS ...')
            hdu.writeto(file, overwrite=True)
            print('DONE!')
    print('')
    print('Finished compressing 1 group to all cubes.')
    print('')
  
# List uncal files
uncal_files = sorted(glob.glob(os.path.join(uncal_dir, '*_uncal.fits')))
sw_uncal_files = [uncfile for uncfile in uncal_files if 'long' not in uncfile]
lw_uncal_files = [uncfile for uncfile in uncal_files if 'long' in uncfile]

colnames = ('Instrument', 'Filter', 'Pupil', 'Number of Integrations', 'Number of Groups',
            'Readout pattern', 'Dither position number')
dtypes = ('S7', 'S10', 'S10', 'i4', 'i4', 'S15', 'i4')
meta_check = Table(names=(colnames), dtype=dtypes)

# Open example files and get metadata for display
if len(sw_uncal_files) > 0:
    sw_examine = datamodels.open(sw_uncal_files[0])
    sw_row = [sw_examine.meta.instrument.name, sw_examine.meta.instrument.filter,
              sw_examine.meta.instrument.pupil, sw_examine.meta.exposure.nints,
              sw_examine.meta.exposure.ngroups, sw_examine.meta.exposure.readpatt,
              sw_examine.meta.dither.position_number]
    meta_check.add_row(sw_row)

if len(lw_uncal_files) > 0:
    lw_examine = datamodels.open(lw_uncal_files[0])
    lw_row = [lw_examine.meta.instrument.name, lw_examine.meta.instrument.filter,
              lw_examine.meta.instrument.pupil, lw_examine.meta.exposure.nints,
              lw_examine.meta.exposure.ngroups, lw_examine.meta.exposure.readpatt,
              lw_examine.meta.dither.position_number]
    meta_check.add_row(lw_row)

# Print out exposure info
print('########### Meta checking ##############')
print(meta_check)
print('######################')

###### PIPELINE 2025 NIRCAM #######

##### Detector1 #######
print('')
print('#########################################')
print('STAGE 1')
print('#########################################')
print('')

if unlimit_desat==0 or unlimit_desat==2:
    det1dict = {}
    det1dict['group_scale'], det1dict['dq_init'], det1dict['saturation'] = {}, {}, {}
    det1dict['ipc'], det1dict['superbias'], det1dict['refpix'] = {}, {}, {}
    det1dict['linearity'], det1dict['persistence'], det1dict['dark_current'], = {}, {}, {}
    det1dict['charge_migration'], det1dict['jump'], det1dict['clean_flicker_noise'] = {}, {}, {}
    det1dict['ramp_fit'], det1dict['gain_scale'] = {}, {}
    det1dict['persistence']['skip'] = True
    det1dict['jump']['expand_large_events'] = True
    if dodet1:
        for uncal in uncal_files:
            print('processing stage 1 file (all groups): ', uncal)
            print('')
            rate_result = Detector1Pipeline.call(uncal, output_dir=det1_dir, steps=det1dict, save_results=True)
    else:
        print('Skipping Detector1 processing')
if unlimit_desat==1 or unlimit_desat==2:
    det1dict = {}
    det1dict['group_scale'], det1dict['dq_init'], det1dict['saturation'] = {}, {}, {}
    det1dict['ipc'], det1dict['superbias'], det1dict['refpix'] = {}, {}, {}
    det1dict['linearity'], det1dict['persistence'], det1dict['dark_current'], = {}, {}, {}
    det1dict['charge_migration'], det1dict['jump'], det1dict['clean_flicker_noise'] = {}, {}, {}
    det1dict['ramp_fit'], det1dict['gain_scale'] = {}, {}
    det1dict['persistence']['skip'] = True
    det1dict['saturation']['skip'] = True
    det1dict['jump']['expand_large_events'] = True
    if dodet1:
        for uncal in uncal_files:
            print('processing stage 1 file (1 group): ', uncal)
            print('')
            rate_result = Detector1Pipeline.call(uncal, output_dir=det1_dir_1group, steps=det1dict, save_results=True)
    else:
        print('Skipping Detector1 processing')    


   
##### Image2 #######
print('')
print('#########################################')
print('STAGE 2')
print('#########################################')
print('')

if unlimit_desat==0 or unlimit_desat==2:
    image2dict = {}
    image2dict['assign_wcs'], image2dict['flat_field'] = {}, {}
    image2dict['photom'], image2dict['resample'] = {}, {}
    sstring = os.path.join(det1_dir, 'jw*rate.fits')  # Use files from the detector1 output folder
    rate_files = sorted(glob.glob(sstring))
    rate_files = [os.path.abspath(fname) for fname in rate_files]
    print(f"Found  {len(rate_files)} science files")
    if doimage2:
        for rate in rate_files:
            print('processing stage 2 file (all groups): ', file)
            print('')
            cal_result = Image2Pipeline.call(rate, output_dir=image2_dir, steps=image2dict, save_results=True)
    else:
        print("Skipping Image2 processing.")
if unlimit_desat==1 or unlimit_desat==2:
    image2dict = {}
    image2dict['assign_wcs'], image2dict['flat_field'] = {}, {}
    image2dict['photom'], image2dict['resample'] = {}, {}
    sstring = os.path.join(det1_dir_1group, 'jw*rate.fits')  # Use files from the detector1 output folder
    rate_files = sorted(glob.glob(sstring))
    rate_files = [os.path.abspath(fname) for fname in rate_files]
    print(f"Found  {len(rate_files)} science files")
    if doimage2:
        for rate in rate_files:
            print('processing stage 2 file (1 group): ', file)
            print('')
            cal_result = Image2Pipeline.call(rate, output_dir=image2_dir_1group, steps=image2dict, save_results=True)
    else:
        print("Skipping Image2 processing.")


##### Image3 #######
print('')
print('#########################################')
print('STAGE 3')
print('#########################################')
print('')

if unlimit_desat==0 or unlimit_desat==2:
    image3dict = {}
    image3dict['assign_mtwcs'], image3dict['tweakreg'], image3dict['skymatch'] = {}, {}, {}
    image3dict['outlier_detection'], image3dict['resample'], image3dict['source_catalog'] = {}, {}, {}
    #image3dict['outlier_detection']['skip'] = True
    #image3dict['tweakreg']['abs_refcat'] = 'GAIADR3'
    sw_sstring = os.path.join(image2_dir, 'jw*nrc??_cal.fits')     # shortwave files. Detectors a1-a4, b1-b4
    lw_sstring = os.path.join(image2_dir, 'jw*nrc*long_cal.fits')  # longwave files. Detectors along, blong 
    sw_cal_files = sorted(glob.glob(sw_sstring))
    lw_cal_files = sorted(glob.glob(lw_sstring))
    sw_cal_files = [os.path.abspath(fname) for fname in sw_cal_files]
    lw_cal_files = [os.path.abspath(fname) for fname in lw_cal_files]
    print(f'Found {len(sw_cal_files)} shortwave science files to process')
    print(f'Found {len(lw_cal_files)} longwave science files to process')
    do_swimage3 = False
    if doimage3:
        if len(sw_cal_files) > 0:
            # Only create an association file if there are SW data files to process
            do_swimage3 = True
            sw_product_name = 'image3_sw'
            sw_association = asn_from_list.asn_from_list(sw_cal_files,
                                                         rule=DMS_Level3_Base,
                                                         product_name=sw_product_name)
            sw_association.data['asn_type'] = 'image3'
            program = datamodels.open(sw_cal_files[0]).meta.observation.program_number
            sw_association.data['program'] = program
            sw_asn_filename, sw_serialized = sw_association.dump(format="json")
            sw_association_im3 = os.path.join(sci_dir, sw_asn_filename)
            with open(sw_association_im3, "w") as fd:
                fd.write(sw_serialized)
    do_lwimage3 = False
    if doimage3:
        if len(lw_cal_files) > 0:
            # Only create an association file if there are SW data files to process
            do_lwimage3 = True
            lw_product_name = 'image3_lw'
            lw_association = asn_from_list.asn_from_list(lw_cal_files,
                                                         rule=DMS_Level3_Base,
                                                         product_name=lw_product_name)
            lw_association.data['asn_type'] = 'image3'
            program = datamodels.open(lw_cal_files[0]).meta.observation.program_number
            lw_association.data['program'] = program
            lw_asn_filename, lw_serialized = lw_association.dump(format="json")
            lw_association_im3 = os.path.join(sci_dir, lw_asn_filename)
            with open(lw_association_im3, "w") as fd:
                fd.write(lw_serialized)
    if doimage3 and do_lwimage3:
        lw_i2d_result = Image3Pipeline.call(lw_association_im3, output_dir=image3_dir, steps=image3dict, save_results=True)
    else:
        print('Skipping Image3 LW processing')
    if doimage3 and do_lwimage3:
        lw_i2d_file = os.path.join(image3_dir, f'{lw_product_name}_i2d.fits')
        lw_data = datamodels.open(lw_i2d_file)
        tree = {"wcs": lw_data.meta.wcs}
        wcs_file = AsdfFile(tree)
        gwcs_filename = os.path.join(image3_dir + 'lw_gwcs.asdf')
        print(f'Saving gWCS into {gwcs_filename}')
        wcs_file.write_to(gwcs_filename)
        ysize, xsize = lw_data.data.shape
    if doimage3 and do_swimage3:
        sw_i2d_result = Image3Pipeline.call(sw_association_im3, output_dir=image3_dir, steps=image3dict, save_results=True)
    else:
        print('Skipping Image3 SW processing')
        
if unlimit_desat==1 or unlimit_desat==2: 
    image3dict = {}
    image3dict['assign_mtwcs'], image3dict['tweakreg'], image3dict['skymatch'] = {}, {}, {}
    image3dict['outlier_detection'], image3dict['resample'], image3dict['source_catalog'] = {}, {}, {}
    #image3dict['outlier_detection']['skip'] = True
    #image3dict['tweakreg']['abs_refcat'] = 'GAIADR3'
    sw_sstring = os.path.join(image2_dir_1group, 'jw*nrc??_cal.fits')     # shortwave files. Detectors a1-a4, b1-b4
    lw_sstring = os.path.join(image2_dir_1group, 'jw*nrc*long_cal.fits')  # longwave files. Detectors along, blong 
    sw_cal_files = sorted(glob.glob(sw_sstring))
    lw_cal_files = sorted(glob.glob(lw_sstring))
    sw_cal_files = [os.path.abspath(fname) for fname in sw_cal_files]
    lw_cal_files = [os.path.abspath(fname) for fname in lw_cal_files]
    print(f'Found {len(sw_cal_files)} shortwave science files to process')
    print(f'Found {len(lw_cal_files)} longwave science files to process')
    do_swimage3 = False
    if doimage3:
        if len(sw_cal_files) > 0:
            # Only create an association file if there are SW data files to process
            do_swimage3 = True
            sw_product_name = 'image3_sw'
            sw_association = asn_from_list.asn_from_list(sw_cal_files,
                                                         rule=DMS_Level3_Base,
                                                         product_name=sw_product_name)
            sw_association.data['asn_type'] = 'image3'
            program = datamodels.open(sw_cal_files[0]).meta.observation.program_number
            sw_association.data['program'] = program
            sw_asn_filename, sw_serialized = sw_association.dump(format="json")
            sw_association_im3 = os.path.join(sci_dir, sw_asn_filename)
            with open(sw_association_im3, "w") as fd:
                fd.write(sw_serialized)
    do_lwimage3 = False
    if doimage3:
        if len(lw_cal_files) > 0:
            # Only create an association file if there are SW data files to process
            do_lwimage3 = True
            lw_product_name = 'image3_lw'
            lw_association = asn_from_list.asn_from_list(lw_cal_files,
                                                         rule=DMS_Level3_Base,
                                                         product_name=lw_product_name)
            lw_association.data['asn_type'] = 'image3'
            program = datamodels.open(lw_cal_files[0]).meta.observation.program_number
            lw_association.data['program'] = program
            lw_asn_filename, lw_serialized = lw_association.dump(format="json")
            lw_association_im3 = os.path.join(sci_dir, lw_asn_filename)
            with open(lw_association_im3, "w") as fd:
                fd.write(lw_serialized)
    if doimage3 and do_lwimage3:
        lw_i2d_result = Image3Pipeline.call(lw_association_im3, output_dir=image3_dir_1group, steps=image3dict, save_results=True)
    else:
        print('Skipping Image3 LW processing')
    if doimage3 and do_lwimage3:
        lw_i2d_file = os.path.join(image3_dir_1group, f'{lw_product_name}_i2d.fits')
        lw_data = datamodels.open(lw_i2d_file)
        tree = {"wcs": lw_data.meta.wcs}
        wcs_file = AsdfFile(tree)
        gwcs_filename = os.path.join(image3_dir_1group + 'lw_gwcs.asdf')
        print(f'Saving gWCS into {gwcs_filename}')
        wcs_file.write_to(gwcs_filename)
        ysize, xsize = lw_data.data.shape
    if doimage3 and do_swimage3:
        sw_i2d_result = Image3Pipeline.call(sw_association_im3, output_dir=image3_dir_1group, steps=image3dict, save_results=True)
    else:
        print('Skipping Image3 SW processing') 

if remove_1group_uncal == True:
    shutil.rmtree(uncal_dir_1group)
print('')
print('################################')
print('##### Processing Finished! #####')
print('################################')
print('')
print('')





