"""
This script is written for CASA 6.2.1.7
Version0 script written by M. Benisty / S. Facchini (Feb 2st 2022)
Version 2 adapted by Stefano Facchini and Myriam Benisty (Mar 2nd 2022)
Combined version adapted to include ACA data by Jane Huang (May 19th 2022)
Version 3 includes EB alignment by Ryan Loomis and Rich Teague (Aug 2022)
"""

# Script applied on MWC 758 by Francesco Zagaria (frzagaria@mpia.de) in December 2024 using CASA 6.6.5-31

"""
Please use the local CASA 6.2.1 installation that also has numba installed
(numba is needed for the alignment):
/lustre/cv/projects/exoALMA/casa_monolithic_6.2.1/bin/casa
Do NOT use casa-alma, since that points to CASA 6.4
To run the script in parallel, do the following:
/lustre/cv/projects/exoALMA/casa_monolithic_6.2.1/bin/mpicasa -n 8 /lustre/cv/projects/exoALMA/casa_monolithic_6.2.1/bin/casa
Do not use more than 8 cores!
"""

import os
import numpy as np
import shutil
import matplotlib
#increase the number of figures that can be open before issuing a warning
matplotlib.rcParams.update({'figure.max_open_warning':80})

#path to your local copy of the gihub repository
#(make sure to do a 'git pull' to have the most up-to-date version)
github_path = '/data/beegfs/astro-storage/groups/benisty/frzagaria/mycasa_versions/'
import sys
sys.path.append(github_path)
#import alignment_default as alignment
execfile(os.path.join(github_path,'reduction_utils_py3_mpi.py'))

execfile(os.path.join(github_path,'keplerian_mask.py'))

prefix = 'MWC_758'
path_to_vis = '../'

# System properties.
incl  = 21.   # deg, from Dong et al. 2018
PA    = 62.   # deg, from Dong et al. 2018
v_sys =  5.9  # km/s, from listobs

# Whether to run tclean in parallel or not.
use_parallel = False

#Continuum imaging
LB_iteration2_cont_averaged = f'{prefix}_time_ave_continuum'

#Define line mask and noise annulus
mask_pa        = PA
mask_semimajor = 0.8
mask_semiminor = 0.8
mask_ra  = '05h30m27.535947s' #05:30:27.535947
mask_dec = '25d19m56.615500s' #25:19:56.615500

mask = f'ellipse[[{mask_ra},{mask_dec}], [{mask_semimajor:.3f}arcsec, {mask_semiminor:.3f}arcsec], {mask_pa:.1f}deg]'

noise_annulus = f"annulus[[{mask_ra}, {mask_dec}],['4.arcsec', '6.arcsec']]"

tclean_wrapper_kwargs = {
    'deconvolver':'multiscale',
    'smallscalebias':0.6,'cycleniter':300,
    'niter':10000000,'mask':mask,#'imsize':2400,
    'parallel':use_parallel,'savemodel':'modelcolumn',
    'interactive':False,'gridder':'standard',
}

robust    = [0.5,1.0,1.5]

cellsize  = np.array([0.027,0.044,0.066])/8.
threshold = np.array([6.62e-03,6.46e-03,6.90e-03])
imsize    = np.array([3600,2400,1600])

scales    = [
    [0,2,4,8,16,24,32],
    [0,2,4,8,16,24,32],
    [0,2,4,8,16,24,32],
]

for _robust,_cellsize,_threshold,_scales,_imsize in zip(robust,cellsize,threshold,scales,imsize):
    imagename = LB_iteration2_cont_averaged+f'{_robust}robust_1.0sigma_0.02gain'
    tclean_wrapper(
        vis       = path_to_vis+LB_iteration2_cont_averaged+'.ms', 
        imagename = imagename,
        threshold = f'{_threshold}mJy',
        scales    = _scales,
        robust    = _robust,
        cellsize  = f'{_cellsize}arcsec',
        gain      = 0.02,
        imsize    = _imsize,
        **tclean_wrapper_kwargs
    )
    estimate_SNR(imagename+'.image',disk_mask=mask,noise_mask=noise_annulus)
    exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)

#MWC_758_time_ave_continuum0.5robust_1.0sigma_0.02gain.image
#Beam 0.045 arcsec x 0.027 arcsec (7.99 deg)
#Flux inside disk mask: 58.61 mJy
#Peak intensity of source: 0.52 mJy/beam
#rms: 6.62e-03 mJy/beam
#Peak SNR: 79.30

#MWC_758_time_ave_continuum1.0robust_1.0sigma_0.02gain.image
#Beam 0.066 arcsec x 0.044 arcsec (8.55 deg)
#Flux inside disk mask: 58.01 mJy
#Peak intensity of source: 1.11 mJy/beam
#rms: 6.46e-03 mJy/beam
#Peak SNR: 172.33

#MWC_758_time_ave_continuum1.5robust_1.0sigma_0.02gain.image
#Beam 0.078 arcsec x 0.066 arcsec (-10.29 deg)
#Flux inside disk mask: 57.81 mJy
#Peak intensity of source: 1.80 mJy/beam
#rms: 6.91e-03 mJy/beam
#Peak SNR: 259.92

#robust: 0.5,1.0,1.5
#std:    3.0290e-01,4.3246e-01,4.3624e-01 mJy

#SO - wider spectral range
SBLB_no_ave_selfcal = f'../{prefix}_SBLB_no_ave_selfcal_time_ave.ms'
contsub_vis = f'{SBLB_no_ave_selfcal}.contsub'

vis_SO = SBLB_no_ave_selfcal[:-3]+'_SOwide.ms'
os.system(f'rm -rf {vis_SO}*')
spw_SO = '3:356~582,7:356~582,11:356~582,15:356~582,19:356~582,23:356~582,27:356~582,31:356~582,35:356~582'

split(vis=SBLB_no_ave_selfcal,outputvis=vis_SO,spw=spw_SO,datacolumn='data',keepflags=False)
listobs(vis=vis_SO,listfile=vis_SO+'.listobs.txt',overwrite=True)

split(vis=contsub_vis,outputvis=f'{vis_SO}.contsub',spw=spw_SO,datacolumn='data',keepflags=False)
listobs(vis=vis_SO+'.contsub',listfile=vis_SO+'.contsub.listobs.txt',overwrite=True)

#SO - wider spectral range final imaging
SBLB_no_ave_selfcal = f'{prefix}_SBLB_no_ave_selfcal_time_ave.ms'
vis_SO = SBLB_no_ave_selfcal[:-3]+'_SOwide.ms'

#Define line mask and noise annulus
mask_pa        = PA
mask_semimajor = 0.8
mask_semiminor = 0.8
mask_ra  = '05h30m27.535947s' #05:30:27.535947
mask_dec = '25d19m56.615500s' #25:19:56.615500

mask = f'ellipse[[{mask_ra},{mask_dec}], [{mask_semimajor:.3f}arcsec, {mask_semiminor:.3f}arcsec], {mask_pa:.1f}deg]'

noise_annulus = f"annulus[[{mask_ra}, {mask_dec}],['4.arcsec', '6.arcsec']]"

tclean_kwargs = {
    'specmode':'cube','restoringbeam':'common','nterms':1,
    'start':'-96.3km/s','width':'1.4km/s','nchan':147,
    'deconvolver':'multiscale','niter':50000,
    'outframe':'LSRK','veltype':'radio','restfreq':'219.9494420GHz',
    'usemask':'user','mask':mask,'interactive':False,
}

imagename = vis_SO[:-3]+'.contsub_image_natural_1.0sigma_uvtaper0.25arcsec_0.01gain_8.0ppb'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.image.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis       = path_to_vis+vis_SO+'.contsub',
    imagename = imagename,
    scales    = [0,2,4,8],
    cell      = '0.04825arcsec',
    gain      = 0.01,
    weighting = 'natural',
    threshold = '6.85e-01mJy',
    uvtaper   = '0.25arcsec',
    imsize    = 1200,
    **tclean_kwargs
)
estimate_SNR(imagename+'.image',disk_mask=mask,noise_mask=noise_annulus)
for _chans in [70,71,72,73,74,75]:
    rms = imstat(imagename=imagename+'.image',region=noise_annulus,chans=f'{_chans}')['rms'][0]
    print('#rms: %.2e Jy/beam' %rms)
#rms: [6.60e-04,6.32e-04,6.25e-04,6.31e-04,6.37e-04,6.32e-04] Jy/beam
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#MWC_758_SBLB_no_ave_selfcal_time_ave_SOwide.contsub_image_natural_1.0sigma_uvtaper0.25arcsec_0.01gain_8.0ppb.image
#Beam 0.428 arcsec x 0.386 arcsec (-11.74 deg)
#Flux inside disk mask: 25.60 mJy
#Peak intensity of source: 4.10 mJy/beam
#rms: 6.85e-01 mJy/beam
#Peak SNR: 5.99

#SO - final imaging
SBLB_no_ave_selfcal = f'{prefix}_SBLB_no_ave_selfcal_time_ave.ms'
vis_SO = SBLB_no_ave_selfcal[:-3]+'_SO.ms'

#Define line mask and noise annulus
mask_pa        = PA
mask_semimajor = 0.8
mask_semiminor = 0.8
mask_ra  = '05h30m27.535947s' #05:30:27.535947
mask_dec = '25d19m56.615500s' #25:19:56.615500

mask = f'ellipse[[{mask_ra},{mask_dec}], [{mask_semimajor:.3f}arcsec, {mask_semiminor:.3f}arcsec], {mask_pa:.1f}deg]'

noise_annulus = f"annulus[[{mask_ra}, {mask_dec}],['4.arcsec', '6.arcsec']]"

tclean_kwargs = {
    'specmode':'cube','restoringbeam':'common','nterms':1,
    'start':'-1.10km/s','width':'1.4km/s','nchan':10,
    'deconvolver':'multiscale','niter':50000,
    'outframe':'LSRK','veltype':'radio','restfreq':'219.9494420GHz',
    'usemask':'user','mask':mask,'interactive':False,
}

tclean_kwargs_keplmask = {
    'specmode':'cube','restoringbeam':'common','nterms':1,
    'start':'-1.10km/s','width':'1.4km/s','nchan':10,
    'deconvolver':'multiscale','niter':50000,
    'interactive':False,'usemask':'user',
    'outframe':'LSRK','veltype':'radio','restfreq':'219.9494420GHz',
}

imagename = vis_SO[:-3]+'.contsub_image_natural_1.0sigma_uvtaper0.5arcsec_0.01gain_8.0ppb_fewchans_new'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.image.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis       = path_to_vis+vis_SO+'.contsub',
    imagename = imagename,
    scales    = [0,2,4,8],
    cell      = '0.079625arcsec',
    gain      = 0.01,
    weighting = 'natural',
    threshold = '8.81e-01mJy',
    uvtaper   = '0.5arcsec',
    imsize    = 800,
    **tclean_kwargs
)
estimate_SNR(imagename+'.image',disk_mask=mask,noise_mask=noise_annulus)
for _chans in [3,4,5,6,7,8]:
    rms = imstat(imagename=imagename+'.image',region=noise_annulus,chans=f'{_chans}')['rms'][0]
    print('#rms: %.2e Jy/beam' %rms)
#rms: [8.39e-04,8.16e-04,8.32e-04,8.46e-04,8.31e-04,8.21e-04] Jy/beam
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#MWC_758_SBLB_no_ave_selfcal_time_ave_SO.contsub_image_natural_1.0sigma_uvtaper0.5arcsec_0.01gain_8.0ppb_fewchans_new.image
#Beam 0.720 arcsec x 0.637 arcsec (-6.54 deg)
#Flux inside disk mask: 23.81 mJy
#Peak intensity of source: 7.05 mJy/beam
#rms: 8.81e-01 mJy/beam
#Peak SNR: 8.00

imagename = vis_SO[:-3]+'.contsub_image_natural_1.0sigma_uvtaper0.25arcsec_0.01gain_8.0ppb_fewchans_new'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.image.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis       = path_to_vis+vis_SO+'.contsub',
    imagename = imagename,
    scales    = [0,2,4,8],
    cell      = '0.04825arcsec',
    gain      = 0.01,
    weighting = 'natural',
    threshold = '6.61e-01mJy',
    uvtaper   = '0.25arcsec',
    imsize    = 1200,
    **tclean_kwargs
)
estimate_SNR(imagename+'.image',disk_mask=mask,noise_mask=noise_annulus)
for _chans in [3,4,5,6,7,8]:
    rms = imstat(imagename=imagename+'.image',region=noise_annulus,chans=f'{_chans}')['rms'][0]
    print('#rms: %.2e Jy/beam' %rms)
#rms: [6.32e-04,6.26e-04,6.32e-04,6.37e-04,6.32e-04,6.31e-04] Jy/beam
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#MWC_758_SBLB_no_ave_selfcal_time_ave_SO.contsub_image_natural_1.0sigma_uvtaper0.25arcsec_0.01gain_8.0ppb_fewchans_new.image
#Beam 0.428 arcsec x 0.386 arcsec (-11.73 deg)
#Flux inside disk mask: 30.05 mJy
#Peak intensity of source: 4.39 mJy/beam
#rms: 6.61e-01 mJy/beam
#Peak SNR: 6.65

imagename = vis_SO[:-3]+'.contsub_image_natural_1.0sigma_uvtaper0.2arcsec_0.01gain_8.0ppb_fewchans_new'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.image.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis       = path_to_vis+vis_SO+'.contsub',
    imagename = imagename,
    scales    = [0,2,4,8],
    cell      = '0.041arcsec',
    gain      = 0.01,
    weighting = 'natural',
    threshold = '6.11e-01mJy',
    uvtaper   = '0.2arcsec',
    imsize    = 1200,
    **tclean_kwargs
)
estimate_SNR(imagename+'.image',disk_mask=mask,noise_mask=noise_annulus)
for _chans in [3,4,5,6,7,8]:
    rms = imstat(imagename=imagename+'.image',region=noise_annulus,chans=f'{_chans}')['rms'][0]
    print('#rms: %.2e Jy/beam' %rms)
#rms: [5.88e-04,5.83e-04,5.89e-04,5.90e-04,5.87e-04,5.86e-04] Jy/beam
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#MWC_758_SBLB_no_ave_selfcal_time_ave_SO.contsub_image_natural_1.0sigma_uvtaper0.2arcsec_0.01gain_8.0ppb_fewchans_new.image
#Beam 0.367 arcsec x 0.328 arcsec (-20.29 deg)
#Flux inside disk mask: 31.31 mJy
#Peak intensity of source: 3.73 mJy/beam
#rms: 6.11e-01 mJy/beam
#Peak SNR: 6.10

imagename = vis_SO[:-3]+'.contsub_image_natural_1.0sigma_uvtaper0.15arcsec_0.01gain_8.0ppb_fewchans_new'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.image.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis       = path_to_vis+vis_SO+'.contsub',
    imagename = imagename,
    scales    = [0,2,4,8],
    cell      = '0.033arcsec',
    gain      = 0.01,
    weighting = 'natural',
    threshold = '5.60e-01mJy',
    uvtaper   = '0.15arcsec',
    imsize    = 1200,
    **tclean_kwargs
)
estimate_SNR(imagename+'.image',disk_mask=mask,noise_mask=noise_annulus)
for _chans in [3,4,5,6,7,8]:
    rms = imstat(imagename=imagename+'.image',region=noise_annulus,chans=f'{_chans}')['rms'][0]
    print('#rms: %.2e Jy/beam' %rms)
#rms: [5.44e-04,5.38e-04,5.45e-04,5.41e-04,5.39e-04,5.38e-04] Jy/beam
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#MWC_758_SBLB_no_ave_selfcal_time_ave_SO.contsub_image_natural_1.0sigma_uvtaper0.15arcsec_0.01gain_8.0ppb_fewchans_new.image
#Beam 0.305 arcsec x 0.264 arcsec (-29.12 deg)
#Flux inside disk mask: 33.24 mJy
#Peak intensity of source: 3.21 mJy/beam
#rms: 5.60e-01 mJy/beam
#Peak SNR: 5.73

imagename = vis_SO[:-3]+'.contsub_image_natural_1.0sigma_uvtaper0.1arcsec_0.01gain_8.0ppb_fewchans_new'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.image.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis       = path_to_vis+vis_SO+'.contsub',
    imagename = imagename,
    scales    = [0,2,4,8],
    cell      = '0.025125arcsec',
    gain      = 0.01,
    weighting = 'natural',
    threshold = '5.08e-01mJy',
    uvtaper   = '0.1arcsec',
    imsize    = 1440,
    **tclean_kwargs
)
estimate_SNR(imagename+'.image',disk_mask=mask,noise_mask=noise_annulus)
for _chans in [3,4,5,6,7,8]:
    rms = imstat(imagename=imagename+'.image',region=noise_annulus,chans=f'{_chans}')['rms'][0]
    print('#rms: %.2e Jy/beam' %rms)
#rms: [5.00e-04,4.95e-04,5.00e-04,4.94e-04,4.92e-04,4.90e-04] Jy/beam
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#MWC_758_SBLB_no_ave_selfcal_time_ave_SO.contsub_image_natural_1.0sigma_uvtaper0.1arcsec_0.01gain_8.0ppb_fewchans_new.image
#Beam 0.247 arcsec x 0.201 arcsec (-32.14 deg)
#Flux inside disk mask: 34.47 mJy
#Peak intensity of source: 2.90 mJy/beam
#rms: 5.08e-01 mJy/beam
#Peak SNR: 5.71

imagename = vis_SO[:-3]+'.contsub_image_natural_1.0sigma_uvtaper0.05arcsec_0.01gain_8.0ppb_fewchans_new'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.image.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis       = path_to_vis+vis_SO+'.contsub',
    imagename = imagename,
    scales    = [0,2,4,8],
    cell      = '0.017625arcsec',
    gain      = 0.01,
    weighting = 'natural',
    threshold = '4.52e-01mJy',
    uvtaper   = '0.05arcsec',
    imsize    = 2000,
    **tclean_kwargs
)
estimate_SNR(imagename+'.image',disk_mask=mask,noise_mask=noise_annulus)
for _chans in [3,4,5,6,7,8]:
    rms = imstat(imagename=imagename+'.image',region=noise_annulus,chans=f'{_chans}')['rms'][0]
    print('#rms: %.2e Jy/beam' %rms)
#rms: [4.52e-04,4.48e-04,4.51e-04,4.43e-04,4.39e-04,4.36e-04] Jy/beam
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#MWC_758_SBLB_no_ave_selfcal_time_ave_SO.contsub_image_natural_1.0sigma_uvtaper0.05arcsec_0.01gain_8.0ppb_fewchans_new.image
#Beam 0.185 arcsec x 0.141 arcsec (-29.71 deg)
#Flux inside disk mask: 34.37 mJy
#Peak intensity of source: 2.58 mJy/beam
#rms: 4.52e-01 mJy/beam
#Peak SNR: 5.72

imagename = vis_SO[:-3]+'.contsub_image_natural_1.0sigma_0.01gain_8.0ppb_fewchans_new'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.image.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis       = path_to_vis+vis_SO+'.contsub',
    imagename = imagename,
    scales    = [0,2,4,8],
    cell      = '0.009arcsec',
    gain      = 0.01,
    weighting = 'natural',
    threshold = '3.82e-01mJy',
    imsize    = 4000,
    **tclean_kwargs
)
estimate_SNR(imagename+'.image',disk_mask=mask,noise_mask=noise_annulus)
for _chans in [3,4,5,6,7,8]:
    rms = imstat(imagename=imagename+'.image',region=noise_annulus,chans=f'{_chans}')['rms'][0]
    print('#rms: %.2e Jy/beam' %rms)
#rms: [3.93e-04,3.94e-04,3.92e-04,3.81e-04,3.73e-04,3.66e-04] Jy/beam
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#MWC_758_SBLB_no_ave_selfcal_time_ave_SO.contsub_image_natural_1.0sigma_0.01gain_8.0ppb_fewchans_new.image
#Beam 0.091 arcsec x 0.072 arcsec (-18.61 deg)
#Flux inside disk mask: 37.31 mJy
#Peak intensity of source: 2.01 mJy/beam
#rms: 3.82e-01 mJy/beam
#Peak SNR: 5.25

imagename = vis_SO[:-3]+'.contsub_image_natural_1.0sigma_uvtaper0.5arcsec_0.01gain_8.0ppb_fewchans_new_keplmask'
make_mask(image=imagename[:-9]+'.image', inc=19, PA=240., mstar=1.4, dist=156., vlsr=5.9e3, nbeams=1.0, r_min=0.0, r_max=1.0)
keplmask = imagename[:-9]+'.mask.image'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.image.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis       = path_to_vis+vis_SO+'.contsub',
    imagename = imagename,
    scales    = [0,2],
    cell      = '0.079625arcsec',
    gain      = 0.01,
    weighting = 'natural',
    threshold = '8.83e-01mJy',
    uvtaper   = '0.5arcsec',
    imsize    = 800,
    mask      = keplmask,
    **tclean_kwargs_keplmask
)
estimate_SNR(imagename+'.image',disk_mask=mask,noise_mask=noise_annulus)
for _chans in [3,4,5,6,7,8]:
    rms = imstat(imagename=imagename+'.image',region=noise_annulus,chans=f'{_chans}')['rms'][0]
    print('#rms: %.2e Jy/beam' %rms)
#rms: [8.40e-04,8.15e-04,8.28e-04,8.47e-04,8.30e-04,8.23e-04] Jy/beam
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#MWC_758_SBLB_no_ave_selfcal_time_ave_SO.contsub_image_natural_1.0sigma_uvtaper0.5arcsec_0.01gain_8.0ppb_fewchans_new_keplmask.image
#Beam 0.720 arcsec x 0.637 arcsec (-6.54 deg)
#Flux inside disk mask: 24.06 mJy
#Peak intensity of source: 7.00 mJy/beam
#rms: 8.83e-01 mJy/beam
#Peak SNR: 7.93

imagename = vis_SO[:-3]+'.contsub_image_natural_1.0sigma_uvtaper0.25arcsec_0.01gain_8.0ppb_fewchans_new_keplmask'
make_mask(image=imagename[:-9]+'.image', inc=19, PA=240., mstar=1.4, dist=156., vlsr=5.9e3, nbeams=1.0, r_min=0.0, r_max=1.0)
keplmask = imagename[:-9]+'.mask.image'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.image.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis       = path_to_vis+vis_SO+'.contsub',
    imagename = imagename,
    scales    = [0,2],
    cell      = '0.04825arcsec',
    gain      = 0.01,
    weighting = 'natural',
    threshold = '6.61e-01mJy',
    uvtaper   = '0.25arcsec',
    imsize    = 1200,
    mask      = keplmask,
    **tclean_kwargs_keplmask
)
estimate_SNR(imagename+'.image',disk_mask=mask,noise_mask=noise_annulus)
for _chans in [3,4,5,6,7,8]:
    rms = imstat(imagename=imagename+'.image',region=noise_annulus,chans=f'{_chans}')['rms'][0]
    print('#rms: %.2e Jy/beam' %rms)
#rms: [6.32e-04,6.25e-04,6.28e-04,6.34e-04,6.31e-04,6.32e-04] Jy/beam
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#MWC_758_SBLB_no_ave_selfcal_time_ave_SO.contsub_image_natural_1.0sigma_uvtaper0.25arcsec_0.01gain_8.0ppb_fewchans_new_keplmask.image
#Beam 0.428 arcsec x 0.386 arcsec (-11.73 deg)
#Flux inside disk mask: 30.43 mJy
#Peak intensity of source: 4.28 mJy/beam
#rms: 6.61e-01 mJy/beam
#Peak SNR: 6.48

imagename = vis_SO[:-3]+'.contsub_image_natural_1.0sigma_uvtaper0.2arcsec_0.01gain_8.0ppb_fewchans_new_keplmask'
make_mask(image=imagename[:-9]+'.image', inc=19, PA=240., mstar=1.4, dist=156., vlsr=5.9e3, nbeams=1.0, r_min=0.0, r_max=1.0)
keplmask = imagename[:-9]+'.mask.image'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.image.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis       = path_to_vis+vis_SO+'.contsub',
    imagename = imagename,
    scales    = [0,2],
    cell      = '0.041arcsec',
    gain      = 0.01,
    weighting = 'natural',
    threshold = '6.12e-01mJy',
    uvtaper   = '0.2arcsec',
    imsize    = 1200,
    mask      = keplmask,
    **tclean_kwargs_keplmask
)
estimate_SNR(imagename+'.image',disk_mask=mask,noise_mask=noise_annulus)
for _chans in [3,4,5,6,7,8]:
    rms = imstat(imagename=imagename+'.image',region=noise_annulus,chans=f'{_chans}')['rms'][0]
    print('#rms: %.2e Jy/beam' %rms)
#rms: [5.89e-04,5.83e-04,5.85e-04,5.87e-04,5.86e-04,5.88e-04] Jy/beam
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#MWC_758_SBLB_no_ave_selfcal_time_ave_SO.contsub_image_natural_1.0sigma_uvtaper0.2arcsec_0.01gain_8.0ppb_fewchans_new_keplmask.image
#Beam 0.367 arcsec x 0.328 arcsec (-20.29 deg)
#Flux inside disk mask: 31.71 mJy
#Peak intensity of source: 3.64 mJy/beam
#rms: 6.12e-01 mJy/beam
#Peak SNR: 5.96

imagename = vis_SO[:-3]+'.contsub_image_natural_1.0sigma_uvtaper0.15arcsec_0.01gain_8.0ppb_fewchans_new_keplmask'
make_mask(image=imagename[:-9]+'.image', inc=19, PA=240., mstar=1.4, dist=156., vlsr=5.9e3, nbeams=1.0, r_min=0.0, r_max=1.0)
keplmask = imagename[:-9]+'.mask.image'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.image.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis       = path_to_vis+vis_SO+'.contsub',
    imagename = imagename,
    scales    = [0,2],
    cell      = '0.033arcsec',
    gain      = 0.01,
    weighting = 'natural',
    threshold = '5.60e-01mJy',
    uvtaper   = '0.15arcsec',
    imsize    = 1200,
    mask      = keplmask,
    **tclean_kwargs_keplmask
)
estimate_SNR(imagename+'.image',disk_mask=mask,noise_mask=noise_annulus)
for _chans in [3,4,5,6,7,8]:
    rms = imstat(imagename=imagename+'.image',region=noise_annulus,chans=f'{_chans}')['rms'][0]
    print('#rms: %.2e Jy/beam' %rms)
#rms: [5.44e-04,5.38e-04,5.42e-04,5.39e-04,5.38e-04,5.39e-04] Jy/beam
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#MWC_758_SBLB_no_ave_selfcal_time_ave_SO.contsub_image_natural_1.0sigma_uvtaper0.15arcsec_0.01gain_8.0ppb_fewchans_new_keplmask.image
#Beam 0.305 arcsec x 0.264 arcsec (-29.12 deg)
#Flux inside disk mask: 30.99 mJy
#Peak intensity of source: 3.13 mJy/beam
#rms: 5.60e-01 mJy/beam
#Peak SNR: 5.59

imagename = vis_SO[:-3]+'.contsub_image_natural_1.0sigma_uvtaper0.1arcsec_0.01gain_8.0ppb_fewchans_new_keplmask'
make_mask(image=imagename[:-9]+'.image', inc=19, PA=240., mstar=1.4, dist=156., vlsr=5.9e3, nbeams=1.0, r_min=0.0, r_max=1.0)
keplmask = imagename[:-9]+'.mask.image'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.image.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis       = path_to_vis+vis_SO+'.contsub',
    imagename = imagename,
    scales    = [0,2],
    cell      = '0.025125arcsec',
    gain      = 0.01,
    weighting = 'natural',
    threshold = '5.09e-01mJy',
    uvtaper   = '0.1arcsec',
    imsize    = 1440,
    mask      = keplmask,
    **tclean_kwargs_keplmask
)
estimate_SNR(imagename+'.image',disk_mask=mask,noise_mask=noise_annulus)
for _chans in [3,4,5,6,7,8]:
    rms = imstat(imagename=imagename+'.image',region=noise_annulus,chans=f'{_chans}')['rms'][0]
    print('#rms: %.2e Jy/beam' %rms)
#rms: [5.01e-04,4.95e-04,4.98e-04,4.92e-04,4.92e-04,4.91e-04] Jy/beam
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#MWC_758_SBLB_no_ave_selfcal_time_ave_SO.contsub_image_natural_1.0sigma_uvtaper0.1arcsec_0.01gain_8.0ppb_fewchans_new_keplmask.image
#Beam 0.247 arcsec x 0.201 arcsec (-32.14 deg)
#Flux inside disk mask: 28.39 mJy
#Peak intensity of source: 2.91 mJy/beam
#rms: 5.09e-01 mJy/beam
#Peak SNR: 5.72

imagename = vis_SO[:-3]+'.contsub_image_natural_1.0sigma_uvtaper0.05arcsec_0.01gain_8.0ppb_fewchans_new_keplmask'
make_mask(image=imagename[:-9]+'.image', inc=19, PA=240., mstar=1.4, dist=156., vlsr=5.9e3, nbeams=1.0, r_min=0.0, r_max=1.0)
keplmask = imagename[:-9]+'.mask.image'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.image.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis       = path_to_vis+vis_SO+'.contsub',
    imagename = imagename,
    scales    = [0,2,4],
    cell      = '0.017625arcsec',
    gain      = 0.01,
    weighting = 'natural',
    threshold = '4.52e-01mJy',
    uvtaper   = '0.05arcsec',
    imsize    = 2000,
    mask      = keplmask,
    **tclean_kwargs_keplmask
)
estimate_SNR(imagename+'.image',disk_mask=mask,noise_mask=noise_annulus)
for _chans in [3,4,5,6,7,8]:
    rms = imstat(imagename=imagename+'.image',region=noise_annulus,chans=f'{_chans}')['rms'][0]
    print('#rms: %.2e Jy/beam' %rms)
#rms: [4.53e-04,4.48e-04,4.49e-04,4.42e-04,4.39e-04,4.37e-04] Jy/beam
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#MWC_758_SBLB_no_ave_selfcal_time_ave_SO.contsub_image_natural_1.0sigma_uvtaper0.05arcsec_0.01gain_8.0ppb_fewchans_new_keplmask.image
#Beam 0.185 arcsec x 0.141 arcsec (-29.71 deg)
#Flux inside disk mask: 17.22 mJy
#Peak intensity of source: 2.71 mJy/beam
#rms: 4.52e-01 mJy/beam
#Peak SNR: 5.99

imagename = vis_SO[:-3]+'.contsub_image_natural_1.0sigma_0.01gain_8.0ppb_fewchans_new_keplmask'
make_mask(image=imagename[:-9]+'.image', inc=19, PA=240., mstar=1.4, dist=156., vlsr=5.9e3, nbeams=1.25, r_min=0.0, r_max=1.0)
keplmask = imagename[:-9]+'.mask.image'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.image.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis       = path_to_vis+vis_SO+'.contsub',
    imagename = imagename,
    scales    = [0,2,4],
    cell      = '0.009arcsec',
    gain      = 0.01,
    weighting = 'natural',
    threshold = '3.83e-01mJy',
    imsize    = 4000,
    mask      = keplmask,
    **tclean_kwargs_keplmask
)
estimate_SNR(imagename+'.image',disk_mask=mask,noise_mask=noise_annulus)
for _chans in [3,4,5,6,7,8]:
    rms = imstat(imagename=imagename+'.image',region=noise_annulus,chans=f'{_chans}')['rms'][0]
    print('#rms: %.2e Jy/beam' %rms)
#rms: [3.93e-04,3.94e-04,3.92e-04,3.81e-04,3.73e-04,3.67e-04] Jy/beam
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#MWC_758_SBLB_no_ave_selfcal_time_ave_SO.contsub_image_natural_1.0sigma_0.01gain_8.0ppb_fewchans_new_keplmask.image
#Beam 0.091 arcsec x 0.072 arcsec (-18.61 deg)
#Flux inside disk mask: -14.30 mJy
#Peak intensity of source: 2.09 mJy/beam
#rms: 3.83e-01 mJy/beam
#Peak SNR: 5.45

imagenames = [
    'natural_1.0sigma_uvtaper0.5arcsec_0.01gain_8.0ppb_fewchans_new',
    'natural_1.0sigma_uvtaper0.25arcsec_0.01gain_8.0ppb_fewchans_new',
    'natural_1.0sigma_uvtaper0.2arcsec_0.01gain_8.0ppb_fewchans_new',
    'natural_1.0sigma_uvtaper0.15arcsec_0.01gain_8.0ppb_fewchans_new',
    'natural_1.0sigma_uvtaper0.1arcsec_0.01gain_8.0ppb_fewchans_new',
    'natural_1.0sigma_uvtaper0.05arcsec_0.01gain_8.0ppb_fewchans_new',
    'natural_1.0sigma_0.01gain_8.0ppb_fewchans_new',
    'natural_1.0sigma_uvtaper0.5arcsec_0.01gain_8.0ppb_fewchans_new_keplmask',
    'natural_1.0sigma_uvtaper0.25arcsec_0.01gain_8.0ppb_fewchans_new_keplmask',
    'natural_1.0sigma_uvtaper0.2arcsec_0.01gain_8.0ppb_fewchans_new_keplmask',
    'natural_1.0sigma_uvtaper0.15arcsec_0.01gain_8.0ppb_fewchans_new_keplmask',
    'natural_1.0sigma_uvtaper0.1arcsec_0.01gain_8.0ppb_fewchans_new_keplmask',
    'natural_1.0sigma_uvtaper0.05arcsec_0.01gain_8.0ppb_fewchans_new_keplmask',
    'natural_1.0sigma_0.01gain_8.0ppb_fewchans_new_keplmask',
]

for _fname in imagenames:
    imagename = f'./integrated_intensity_8.0ppb_fewchans_new/MWC_758_SBLB_no_ave_selfcal_time_ave_SO.contsub_image_{_fname}.image_M0_noclip.fits'
    rms = imstat(imagename=imagename,region=noise_annulus)['rms'][0] #Jy/beam m/s
    print('#noclip rms_ii: %.2e Jy/beam m/s' %rms)

    imagename = f'./integrated_intensity_8.0ppb_fewchans_new/MWC_758_SBLB_no_ave_selfcal_time_ave_SO.contsub_image_{_fname}.image_M0_bestchans.fits'
    rms = imstat(imagename=imagename,region=noise_annulus)['rms'][0] #Jy/beam m/s
    print('#bestchans rms_ii: %.2e Jy/beam m/s' %rms)

#rms_noclip:       5.28e+00,3.95e+00,3.64e+00,3.31e+00,2.98e+00,2.62e+00,2.15e+00,5.30e+00,3.95e+00,3.63e+00,3.31e+00,2.98e+00,2.62e+00,2.15e+00 Jy/beam m/s
#rms_bestchans_37: 3.83e+00,2.83e+00,2.61e+00,2.38e+00,2.15e+00,1.90e+00,1.57e+00,3.84e+00,2.83e+00,2.60e+00,2.38e+00,2.15e+00,1.90e+00,1.57e+00 Jy/beam m/s
#rms_bsetchans_28: 4.64e+00,3.48e+00,3.21e+00,2.92e+00,2.64e+00,2.32e+00,1.90e+00,4.65e+00,3.48e+00,3.20e+00,2.38e+00,2.63e+00,2.31e+00,1.90e+00 Jy/beam m/s

for _fname in imagenames:
    imagename = f'./peak_intensity_8.0ppb_fewchans_new/MWC_758_SBLB_no_ave_selfcal_time_ave_SO.contsub_image_{_fname}.image_M8_noclip.fits'
    rms = imstat(imagename=imagename,region=noise_annulus)['rms'][0] #Jy/beam
    print('#noclip rms_pi: %.2e Jy/beam' %rms)

    imagename = f'./peak_intensity_8.0ppb_fewchans_new/MWC_758_SBLB_no_ave_selfcal_time_ave_SO.contsub_image_{_fname}.image_M8_bestchans.fits'
    rms = imstat(imagename=imagename,region=noise_annulus)['rms'][0] #Jy/beam
    print('#bestchans rms_pi: %.2e Jy/beam' %rms)

#rms_noclip:       1.40e-03,1.04e-03,9.58e-04,8.77e-04,7.97e-04,7.10e-04,6.05e-04,1.40e-03,1.04e-03,9.59e-04,8.78e-04,7.98e-04,7.11e-04,6.06e-04 Jy/beam
#rms_bestchans_37: 1.06e-03,7.93e-04,7.37e-04,6.80e-04,6.25e-04,5.66e-04,4.96e-04,1.06e-03,7.91e-04,7.36e-04,6.80e-04,6.25e-04,5.66e-04,4.96e-04 Jy/beam
#rms_bestchans_28: 1.17e-03,8.82e-04,8.21e-04,7.58e-04,6.95e-04,6.27e-04,5.47e-04,1.17e-03,8.81e-04,8.21e-04,7.58e-04,6.96e-04,6.28e-04,5.47e-04 Jy/beam