import os

print('')
user_home_dir = os.path.expanduser('~')
sci_dir = os.path.join(user_home_dir, 'pipeline/MRS_files/')
bg_dir = '' # If no background observation, use an empty string

# --------------------------Set Processing Steps--------------------------
# Whether or not to process only data from a given MRS band/channel (useful
# if overriding reference files)
# Note that BOTH parameters must be set in order to work
use_ch = ''  # '12' or '34'
use_band = ''  # 'SHORT', 'MEDIUM', or 'LONG'

unlimit_desat = 2   # 0=dont do it ;  1=do it only for 1 group  ; 2=do the pipeline and the 1 group desaturation 


# Science processing
dodet1 = True  # calwebb_detector1
dospec2 = False  # calwebb_spec2
dospec3 = False  # calwebb_spec3
fringing = False # Defringing step
remove_1group_uncal = True  # Remove uncal_1group folder at the end of the run.

# Background processing
dodet1bg = False  # calwebb_detector1
dospec2bg = False  # calwebb_spec2 (needed for Master Background subtraction)
master_bg = False  # Master-background subtraction in spec3 (subtract spectrum generated from the backgrounds).
pixel_bg = False  # Pixel-based background subtraction in spec2 (direct pixel subtraction).

print('############################################################')
print('')
print('PIPELINE CONFIGURATION')
print('')
if dodet1 == True:
    print('You will execute stage 1 in this run.')
elif dospec2 == True:
    print('You will execute stage 2 in this run.')
elif dospec3 == True:
    print('You will execute stage 3 in this run.')
elif dodet1 == True and dospec2 == False and dospec3 == False:
    print('You are not launching the pipeline.')
print('')
if fringing == True:
    print('You are defringing in stage 2. It will take longer.')
else:
    print('You are not defringing in stage 2.')
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

# ------------------------Set CRDS context and paths----------------------

# Set CRDS reference file context.  Leave commented-out to use the default context
# (latest reference files associated with the calibration pipeline version)
# or set a specific context here.
#%env CRDS_CONTEXT  jwst_1295.pmap

# Check whether the local CRDS cache directory has been set.
# If not, set it to the user home directory
if (os.getenv('CRDS_PATH') is None):
    os.environ['CRDS_PATH'] = os.path.join(os.path.expanduser('~'), 'crds')
# Check whether the CRDS server URL has been set.  If not, set it.
if (os.getenv('CRDS_SERVER_URL') is None):
    os.environ['CRDS_SERVER_URL'] = 'https://jwst-crds.stsci.edu'

# Echo CRDS path and context in use
print('CRDS local filepath:', os.environ['CRDS_PATH'])
print('CRDS file server:', os.environ['CRDS_SERVER_URL'])
print('')

import glob
import copy
import shutil
import time
from pathlib import Path
import numpy as np
from astropy.io import fits
from astroquery.mast import Observations
import jwst
import crds
from jwst.pipeline import Detector1Pipeline
from jwst.pipeline import Spec2Pipeline
from jwst.pipeline import Spec3Pipeline
from jwst import datamodels  # JWST datamodels
from jwst.associations import asn_from_list as afl  # Tools for creating association files
from jwst.associations.lib.rules_level2_base import DMSLevel2bBase  # Definition of a Lvl2 association file
from jwst.associations.lib.rules_level3_base import DMS_Level3_Base  # Definition of a Lvl3 association file
from jwst.stpipe import Step  # Import the wrapper class for pipeline steps

print("JWST Calibration Pipeline Version = {}".format(jwst.__version__))
print('')

# Define a convenience function to select only files of a given channel/band from an input set
def select_ch_band_files(files, use_ch, use_band):
    if ((use_ch != '') & (use_band != '')):
        keep = np.zeros(len(files))
        for ii in range(0, len(files)):
            with fits.open(files[ii]) as hdu:
                hdu.verify()
                hdr = hdu[0].header
                if ((hdr['CHANNEL'] == use_ch) & (hdr['BAND'] == use_band)):
                    keep[ii] = 1
        indx = np.where(keep == 1)
        files_culled = files[indx]
    else:
        files_culled = files
    return files_culled

def writel2asn(onescifile, bgfiles, selfcalfiles, asnfile, prodname):
    # Define the basic association of science files
    asn = afl.asn_from_list([onescifile], rule=DMSLevel2bBase, product_name=prodname)  # Wrap in array since input was single exposure
    #Channel/band configuration for this sci file
    with fits.open(onescifile) as hdu:
        hdu.verify()
        hdr = hdu[0].header
        this_channel, this_band = hdr['CHANNEL'], hdr['BAND']
    # If backgrounds were provided, find which are appropriate to this
    # channel/band and add to association
    for file in bgfiles:
        with fits.open(file) as hdu:
            hdu.verify()
            if ((hdu[0].header['CHANNEL'] == this_channel) & (hdu[0].header['BAND'] == this_band)):
                asn['products'][0]['members'].append({'expname': file, 'exptype': 'background'})
    # If provided with a list of files to use for bad pixel self-calibration, find which
    # are appropriate to this detector and add to association
    for file in selfcalfiles:
        with fits.open(file) as hdu:
            hdu.verify()
            if (hdu[0].header['CHANNEL'] == this_channel):
                asn['products'][0]['members'].append({'expname': file, 'exptype': 'selfcal'})                
    # Write the association to a json file
    _, serialized = asn.dump()
    with open(asnfile, 'w') as outfile:
        outfile.write(serialized)

def writel3asn(scifiles, bgfiles, asnfile, prodname):
    # Define the basic association of science files
    asn = afl.asn_from_list(scifiles, rule=DMS_Level3_Base, product_name=prodname)

    # Add background files to the association
    for file in bgfiles:
        asn['products'][0]['members'].append({'expname': file, 'exptype': 'background'})

    # Write the association to a json file
    _, serialized = asn.dump()
    with open(asnfile, 'w') as outfile:
        outfile.write(serialized)



# Science and bg Folders
# Define output subdirectories to keep science data products organized
uncal_dir = os.path.join(sci_dir, 'uncal')  # Uncalibrated pipeline inputs should be here
det1_dir = os.path.join(sci_dir, 'stage1')  # calwebb_detector1 pipeline outputs will go here
spec2_dir = os.path.join(sci_dir, 'stage2')  # calwebb_spec2 pipeline outputs will go here
spec3_dir = os.path.join(sci_dir, 'stage3')  # calwebb_spec3 pipeline outputs will go here
uncal_dir_1group = os.path.join(sci_dir, 'uncal_1group') 
det1_dir_1group = os.path.join(sci_dir, 'stage1_1group')
spec2_dir_1group = os.path.join(sci_dir, 'stage2_1group') 
spec3_dir_1group = os.path.join(sci_dir, 'stage3_1group')

# Output subdirectories to keep background data products organized
uncal_bgdir = os.path.join(bg_dir, 'uncal')  # Uncalibrated pipeline inputs should be here
det1_bgdir = os.path.join(bg_dir, 'stage1')  # calwebb_detector1 pipeline outputs will go here
spec2_bgdir = os.path.join(bg_dir, 'stage2')  # calwebb_spec2 pipeline outputs will go here

# We need to check that the desired output directories exist, and if not create them
if not os.path.exists(det1_dir):
    os.makedirs(det1_dir)
if not os.path.exists(spec2_dir):
    os.makedirs(spec2_dir)
if not os.path.exists(spec3_dir):
    os.makedirs(spec3_dir)
if unlimit_desat==1 or unlimit_desat==2:
    if not os.path.exists(uncal_dir_1group):
            os.makedirs(uncal_dir_1group)
    if not os.path.exists(det1_dir_1group):
            os.makedirs(det1_dir_1group)
    if not os.path.exists(spec2_dir_1group):
            os.makedirs(spec2_dir_1group)
    if not os.path.exists(spec3_dir_1group):
            os.makedirs(spec3_dir_1group)
if (bg_dir != ''):
    if not os.path.exists(det1_bgdir):
        os.makedirs(det1_bgdir)
    if not os.path.exists(spec2_bgdir):
        os.makedirs(spec2_bgdir)


# Avoiding bugs related to bg 
if (bg_dir == ''):
    dodet1bg = False
    dospec2bg = False


########## Creating 1 group uncals in case they are needed  ############

for filename in glob.glob(os.path.join(uncal_dir, '*.*')):
    shutil.copy(filename, uncal_dir_1group)


print('')
if unlimit_desat==1 or unlimit_desat==2:
    sstring = os.path.join(uncal_dir_1group, 'jw*mirifu*uncal.fits')
    uncal_files = np.array(sorted(glob.glob(sstring)))
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

###### PIPELINE 2025 #######

##### Detector1 #######
print('')
print('#########################################')
print('STAGE 1')
print('#########################################')
print('')
if unlimit_desat==0 or unlimit_desat==2:
    det1dict = {}
    det1dict['group_scale'], det1dict['dq_init'], det1dict['emicorr'], det1dict['saturation'], det1dict['ipc'] = {}, {}, {}, {}, {}
    det1dict['firstframe'], det1dict['lastframe'], det1dict['reset'], det1dict['linearity'], det1dict['rscd'] = {}, {}, {}, {}, {}
    det1dict['dark_current'], det1dict['refpix'], det1dict['charge_migration'], det1dict['jump'], det1dict['ramp_fit'] = {}, {}, {}, {}, {}
    det1dict['gain_scale'] = {}
    det1dict['firstframe']['bright_use_group1'] = True
    sstring = os.path.join(uncal_dir, 'jw*mirifu*uncal.fits')
    uncal_files = np.array(sorted(glob.glob(sstring)))
    uncal_files = select_ch_band_files(uncal_files, use_ch, use_band)
    print('Found ' + str(len(uncal_files)) + ' science input files')
    if dodet1:
        for file in uncal_files:
            print('processing stage 1 file (all groups): ', file)
            print('')
            Detector1Pipeline.call(file, steps=det1dict, save_results=True, output_dir=det1_dir)
    else:
        print('Skipping Detector1 processing for SCI data')
        print('')
        
if unlimit_desat==1 or unlimit_desat==2:
    det1dict = {}
    det1dict['group_scale'], det1dict['dq_init'], det1dict['emicorr'], det1dict['saturation'], det1dict['ipc'] = {}, {}, {}, {}, {}
    det1dict['firstframe'], det1dict['lastframe'], det1dict['reset'], det1dict['linearity'], det1dict['rscd'] = {}, {}, {}, {}, {}
    det1dict['dark_current'], det1dict['refpix'], det1dict['charge_migration'], det1dict['jump'], det1dict['ramp_fit'] = {}, {}, {}, {}, {}
    det1dict['gain_scale'] = {}
    det1dict['saturation']['skip'] = True
    det1dict['firstframe']['bright_use_group1'] = True
    sstring = os.path.join(uncal_dir_1group, 'jw*mirifu*uncal.fits')
    uncal_files = np.array(sorted(glob.glob(sstring)))
    uncal_files = select_ch_band_files(uncal_files, use_ch, use_band)
    print('Found ' + str(len(uncal_files)) + ' science input files')
    if dodet1:
        for file in uncal_files:
            print('processing stage 1 file (1 group): ', file)
            print('')
            Detector1Pipeline.call(file, steps=det1dict, save_results=True, output_dir=det1_dir_1group)
    else:
        print('Skipping Detector1 processing for SCI data')
        print('')
        
##### Spectroscopy2 #######
print('')
print('#########################################')
print('STAGE 2')
print('#########################################')
print('')
if unlimit_desat==0 or unlimit_desat==2:
    spec2dict = {}
    spec2dict['assign_wcs'], spec2dict['badpix_selfcal'], spec2dict['bkg_subtract'], spec2dict['flat_field'], spec2dict['srctype'] = {}, {}, {}, {}, {}
    spec2dict['straylight'], spec2dict['fringe'], spec2dict['photom'], spec2dict['residual_fringe'], spec2dict['pixel_replace'] = {}, {}, {}, {}, {}
    spec2dict['cube_build'], spec2dict['extract_1d'] = {}, {}
    if (pixel_bg is True):
        spec2dict['bkg_subtract']['skip'] = False
    else:
        spec2dict['bkg_subtract']['skip'] = True 
    if fringing == True:
        spec2dict['residual_fringe']['skip'] = False
    sstring = os.path.join(det1_dir, 'jw*mirifu*rate.fits')  # Use files from the detector1 output folder
    ratefiles = sorted(glob.glob(sstring))
    for ii in range(0, len(ratefiles)):
        ratefiles[ii] = os.path.abspath(ratefiles[ii])
    ratefiles = np.array(ratefiles)
    ratefiles = select_ch_band_files(ratefiles, use_ch, use_band)
    sstring = os.path.join(det1_bgdir, 'jw*mirifu*rate.fits')
    bgfiles = sorted(glob.glob(sstring))
    for ii in range(0, len(bgfiles)):
        bgfiles[ii] = os.path.abspath(bgfiles[ii])
    bgfiles = np.array(bgfiles)
    bgfiles = select_ch_band_files(bgfiles, use_ch, use_band)
    selfcalfiles = ratefiles.copy()
    selfcalfiles = np.append(selfcalfiles, bgfiles)
    print('Found ' + str(len(ratefiles)) + ' science files')
    print('Found ' + str(len(bgfiles)) + ' background files')
    print('Found ' + str(len(selfcalfiles)) + ' potential selfcal files')
    spec2dict_sci = copy.deepcopy(spec2dict)
    spec2dict_sci['cube_build']['skip'] = True
    spec2dict_sci['extract_1d']['skip'] = True
    if dospec2:
        for file in ratefiles:
            print('processing stage 2 file (all groups): ', file)
            print('')
            asnfile = os.path.join(sci_dir, 'l2asn.json')
            writel2asn(file, bgfiles, selfcalfiles, asnfile, 'Level2')
            Spec2Pipeline.call(asnfile, steps=spec2dict_sci, save_results=True, output_dir=spec2_dir)
    else:
        print('Skipping Spec2 processing for SCI data')
        print('')
        
if unlimit_desat==1 or unlimit_desat==2:
    spec2dict = {}
    spec2dict['assign_wcs'], spec2dict['badpix_selfcal'], spec2dict['bkg_subtract'], spec2dict['flat_field'], spec2dict['srctype'] = {}, {}, {}, {}, {}
    spec2dict['straylight'], spec2dict['fringe'], spec2dict['photom'], spec2dict['residual_fringe'], spec2dict['pixel_replace'] = {}, {}, {}, {}, {}
    spec2dict['cube_build'], spec2dict['extract_1d'] = {}, {}
    if (pixel_bg is True):
        spec2dict['bkg_subtract']['skip'] = False
    else:
        spec2dict['bkg_subtract']['skip'] = True 
    if fringing == True:
        spec2dict['residual_fringe']['skip'] = False
    sstring = os.path.join(det1_dir_1group, 'jw*mirifu*rate.fits')  # Use files from the detector1 output folder
    ratefiles = sorted(glob.glob(sstring))
    for ii in range(0, len(ratefiles)):
        ratefiles[ii] = os.path.abspath(ratefiles[ii])
    ratefiles = np.array(ratefiles)
    ratefiles = select_ch_band_files(ratefiles, use_ch, use_band)
    sstring = os.path.join(det1_bgdir, 'jw*mirifu*rate.fits')
    bgfiles = sorted(glob.glob(sstring))
    for ii in range(0, len(bgfiles)):
        bgfiles[ii] = os.path.abspath(bgfiles[ii])
    bgfiles = np.array(bgfiles)
    bgfiles = select_ch_band_files(bgfiles, use_ch, use_band)
    selfcalfiles = ratefiles.copy()
    selfcalfiles = np.append(selfcalfiles, bgfiles)
    print('Found ' + str(len(ratefiles)) + ' science files')
    print('Found ' + str(len(bgfiles)) + ' background files')
    print('Found ' + str(len(selfcalfiles)) + ' potential selfcal files')
    spec2dict_sci = copy.deepcopy(spec2dict)
    spec2dict_sci['cube_build']['skip'] = True
    spec2dict_sci['extract_1d']['skip'] = True
    if dospec2:
        for file in ratefiles:
            print('processing stage 2 file (1 group): ', file)
            print('')
            asnfile = os.path.join(sci_dir, 'l2asn.json')
            writel2asn(file, bgfiles, selfcalfiles, asnfile, 'Level2')
            Spec2Pipeline.call(asnfile, steps=spec2dict_sci, save_results=True, output_dir=spec2_dir_1group)
    else:
        print('Skipping Spec2 processing for SCI data')
        print('')

##### Spectroscopy3 #######
print('')
print('#########################################')
print('STAGE 3')
print('#########################################')
print('')
if unlimit_desat==0 or unlimit_desat==2:
    spec3dict = {}
    spec3dict['assign_mtwcs'], spec3dict['master_background'], spec3dict['outlier_detection'], spec3dict['mrs_imatch'], spec3dict['cube_build'] = {}, {}, {}, {}, {}
    spec3dict['pixel_replace'], spec3dict['extract_1d'], spec3dict['spectral_leak'] = {}, {}, {}
    if (master_bg is True):
        spec3dict['master_background']['skip'] = False
    else:
        spec3dict['master_background']['skip'] = True
    spec3dict['pixel_replace']['skip'] = False
    spec3dict['pixel_replace']['algorithm'] = 'mingrad'
    spec3dict['cube_build']['output_type'] = 'band'  # 'band', 'channel' (default)
    spec3dict['cube_build']['coord_system'] = 'ifualign'  # Cube rotation: 'ifualign', 'skyalign' (default)
    spec3dict['extract_1d']['ifu_autocen'] = True # Recommended
    sstring = os.path.join(spec2_dir, 'jw*mirifu*_cal.fits')
    calfiles = sorted(glob.glob(sstring))
    for ii in range(0, len(calfiles)):
        calfiles[ii] = os.path.abspath(calfiles[ii])
    calfiles = np.array(calfiles)
    calfiles = select_ch_band_files(calfiles, use_ch, use_band)
    sstring = os.path.join(spec2_bgdir, 'jw*mirifu*x1d.fits')
    bgfiles = sorted(glob.glob(sstring))
    for ii in range(0, len(bgfiles)):
        bgfiles[ii] = os.path.abspath(bgfiles[ii])
    bgfiles = np.array(bgfiles)
    bgfiles = select_ch_band_files(bgfiles, use_ch, use_band)
    print('Found ' + str(len(calfiles)) + ' science files to process')
    print('Found ' + str(len(bgfiles)) + ' background files to process')
    print('processing stage 3 files. Please wait...')
    asnfile = os.path.join(sci_dir, 'l3asn.json')
    if dospec3:
        writel3asn(calfiles, bgfiles, asnfile, 'Level3')
    if dospec3:
        print('processing stage 3 file (all groups): ', file)
        print('')
        Spec3Pipeline.call(asnfile, steps=spec3dict, save_results=True, output_dir=spec3_dir)
    else:
        print('Skipping Spec3 processing')
        print('')
        

if unlimit_desat==1 or unlimit_desat==2:
    spec3dict = {}
    spec3dict['assign_mtwcs'], spec3dict['master_background'], spec3dict['outlier_detection'], spec3dict['mrs_imatch'], spec3dict['cube_build'] = {}, {}, {}, {}, {}
    spec3dict['pixel_replace'], spec3dict['extract_1d'], spec3dict['spectral_leak'] = {}, {}, {}
    if (master_bg is True):
        spec3dict['master_background']['skip'] = False
    else:
        spec3dict['master_background']['skip'] = True
    spec3dict['pixel_replace']['skip'] = False
    spec3dict['pixel_replace']['algorithm'] = 'mingrad'
    spec3dict['cube_build']['output_type'] = 'band'  # 'band', 'channel' (default)
    spec3dict['cube_build']['coord_system'] = 'ifualign'  # Cube rotation: 'ifualign', 'skyalign' (default)
    spec3dict['extract_1d']['ifu_autocen'] = True # Recommended
    sstring = os.path.join(spec2_dir_1group, 'jw*mirifu*_cal.fits')
    calfiles = sorted(glob.glob(sstring))
    for ii in range(0, len(calfiles)):
        calfiles[ii] = os.path.abspath(calfiles[ii])
    calfiles = np.array(calfiles)
    calfiles = select_ch_band_files(calfiles, use_ch, use_band)
    sstring = os.path.join(spec2_bgdir, 'jw*mirifu*x1d.fits')
    bgfiles = sorted(glob.glob(sstring))
    for ii in range(0, len(bgfiles)):
        bgfiles[ii] = os.path.abspath(bgfiles[ii])
    bgfiles = np.array(bgfiles)
    bgfiles = select_ch_band_files(bgfiles, use_ch, use_band)
    print('Found ' + str(len(calfiles)) + ' science files to process')
    print('Found ' + str(len(bgfiles)) + ' background files to process')
    print('processing stage 3 files. Please wait...')
    asnfile = os.path.join(sci_dir, 'l3asn.json')
    if dospec3:
        writel3asn(calfiles, bgfiles, asnfile, 'Level3')
    if dospec3:
        print('processing stage 3 file (1 group): ', file)
        print('')
        Spec3Pipeline.call(asnfile, steps=spec3dict, save_results=True, output_dir=spec3_dir_1group)
    else:
        print('Skipping Spec3 processing')
        print('')
        


if remove_1group_uncal == True:
    shutil.rmtree(uncal_dir_1group)
print('')
print('################################')
print('##### Processing Finished! #####')
print('################################')
print('')
print('')
