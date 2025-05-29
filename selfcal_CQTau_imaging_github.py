"""
This script is written for CASA 6.2.1.7
Version0 script written by M. Benisty / S. Facchini (Feb 2st 2022)
Version 2 adapted by Stefano Facchini and Myriam Benisty (Mar 2nd 2022)
Combined version adapted to include ACA data by Jane Huang (May 19th 2022)
Version 3 includes EB alignment by Ryan Loomis and Rich Teague (Aug 2022)
"""

# Script applied on CQ Tau by Francesco Zagaria (frzagaria@mpia.de) in December 2024 using CASA 6.6.5-31

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
# import alignment_default as alignment
execfile(os.path.join(github_path,'reduction_utils_py3_mpi.py'))

execfile(os.path.join(github_path,'keplerian_mask.py'))

prefix = 'CQ_Tau'
path_to_vis = '../'

# System properties.
incl  = 35.  # deg, from Ubeira-Gabellini et al. 2019
PA    = 55.  # deg, from Ubeira-Gabellini et al. 2019
v_sys =  6.2 # km/s, from listobs

# Whether to run tclean in parallel or not.
use_parallel = False

#Continuum imaging
LB_iteration2_cont_averaged = f'{prefix}_time_ave_continuum'

#Define mask
mask_pa        = PA
mask_semimajor = 0.7#0.85
mask_semiminor = 0.6#0.72
mask_ra  = '05h35m58.472722s'
mask_dec = '24d44m53.613450s'

mask = f'ellipse[[{mask_ra},{mask_dec}], [{mask_semimajor:.3f}arcsec, {mask_semiminor:.3f}arcsec], {mask_pa:.1f}deg]'

noise_annulus = f"annulus[[{mask_ra}, {mask_dec}],['4.arcsec', '6.arcsec']]"

tclean_wrapper_kwargs = {
    'deconvolver':'multiscale',
    'smallscalebias':0.6,'cycleniter':300,
    'niter':1000000,'mask':mask,'imsize':2400,
    'parallel':use_parallel,'savemodel':'modelcolumn',
    'interactive':False,'gridder':'standard',
}

robust    = [-1.0,-0.5,0.0,0.5,1.0,1.5]

cellsize  = np.array([0.031,0.033,0.040,0.058,0.088,0.099])/8.
threshold = np.array([3.15e-02,2.06e-02,1.43e-02,1.19e-02,1.36e-02,1.51e-02])

scales    = [
    [0,2,4,8,16,24,32], #added larger scales 'cause imagss look rather patchy and the models are not as continuous as for 0.0 and 0.5
    [0,2,4,8,16,24,32],
    [0,2,4,8,16,24],
    [0,2,4,8,16,24],
    [0,2,4,8,16,20], #zero-pixel scale leads to negative single-pixel components (maybe go away with more conservative gain)
    [0,2,4,8,16,20], #zero-pixel scale leads to negative single-pixel components (maybe go away with more conservative gain)
]

for _robust,_cellsize,_threshold,_scales in zip(robust,cellsize,threshold,scales):
    imagename = LB_iteration2_cont_averaged+f'{_robust}robust_1.0sigma_0.02gain'
    tclean_wrapper(
        vis       = path_to_vis+LB_iteration2_cont_averaged+'.ms', 
        imagename = imagename,
        threshold = f'{_threshold}mJy',
        scales    = _scales,
        robust    = _robust,
        cellsize  = f'{_cellsize}arcsec',
        gain      = 0.02,
        **tclean_wrapper_kwargs
    )
    estimate_SNR(imagename+'.image',disk_mask=mask,noise_mask=noise_annulus)
    exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)

#CQ_Tau_time_ave_continuum-1.0robust_1.0sigma_0.02gain.image
#Beam 0.045 arcsec x 0.031 arcsec (-12.58 deg)
#Flux inside disk mask: 145.89 mJy
#Peak intensity of source: 0.92 mJy/beam
#rms: 3.15e-02 mJy/beam
#Peak SNR: 29.32

#CQ_Tau_time_ave_continuum-0.5robust_1.0sigma_0.02gain.image
#Beam 0.047 arcsec x 0.033 arcsec (-10.94 deg)
#Flux inside disk mask: 146.68 mJy
#Peak intensity of source: 0.95 mJy/beam
#rms: 2.06e-02 mJy/beam
#Peak SNR: 46.07

#CQ_Tau_time_ave_continuum0.0robust_1.0sigma_0.02gain.image
#Beam 0.055 arcsec x 0.040 arcsec (-4.06 deg)
#Flux inside disk mask: 146.80 mJy
#Peak intensity of source: 1.29 mJy/beam
#rms: 1.43e-02 mJy/beam
#Peak SNR: 90.04

#CQ_Tau_time_ave_continuum0.5robust_1.0sigma_0.02gain.image
#Beam 0.079 arcsec x 0.058 arcsec (11.88 deg)
#Flux inside disk mask: 147.11 mJy
#Peak intensity of source: 2.56 mJy/beam
#rms: 1.19e-02 mJy/beam
#Peak SNR: 215.10

#CQ_Tau_time_ave_continuum1.0robust_1.0sigma_0.02gain.image
#Beam 0.106 arcsec x 0.088 arcsec (8.22 deg)
#Flux inside disk mask: 146.90 mJy
#Peak intensity of source: 4.86 mJy/beam
#rms: 1.36e-02 mJy/beam
#Peak SNR: 358.95

#CQ_Tau_time_ave_continuum1.5robust_1.0sigma_0.02gain.image
#Beam 0.119 arcsec x 0.099 arcsec (1.10 deg)
#Flux inside disk mask: 146.84 mJy
#Peak intensity of source: 5.95 mJy/beam
#rms: 1.51e-02 mJy/beam
#Peak SNR: 392.57

#robust: -1.0,-0.5,0.0,0.5,1.0,1.5
#std:    4.7706e-01,2.7093e-01,3.7384e-01,2.1975e-01,2.1435e-01,1.9137e-01 mJy

#SO - wider spectral range
SBLB_no_ave_selfcal = f'../{prefix}_SBLB_no_ave_selfcal_time_ave.ms'
contsub_vis = f'{SBLB_no_ave_selfcal}.contsub'

vis_SO = SBLB_no_ave_selfcal[:-3]+'_SOwide.ms'
os.system(f'rm -rf {vis_SO}*')
spw_SO = '5:763~959,6:0~253,14:748~959,15:0~238,22:770~959,23:0~261,30:770~959,31:0~261'

split(vis=SBLB_no_ave_selfcal,outputvis=vis_SO,spw=spw_SO,datacolumn='data',keepflags=False)
listobs(vis=vis_SO,listfile=vis_SO+'.listobs.txt',overwrite=True)

split(vis=contsub_vis,outputvis=f'{vis_SO}.contsub',spw=spw_SO,datacolumn='data',keepflags=False)
listobs(vis=vis_SO+'.contsub',listfile=vis_SO+'.contsub.listobs.txt',overwrite=True)

#SO - wider spectral range final imaging
SBLB_no_ave_selfcal = f'{prefix}_SBLB_no_ave_selfcal_time_ave.ms'
vis_SO = SBLB_no_ave_selfcal[:-3]+'_SOwide.ms'

#Define mask
mask_pa        = PA
mask_semimajor = 0.7
mask_semiminor = 0.6
mask_ra  = '05h35m58.472722s'#05:35:58.472722
mask_dec = '24d44m53.613450s'#24:44:53.613450

mask = f'ellipse[[{mask_ra},{mask_dec}], [{mask_semimajor:.3f}arcsec, {mask_semiminor:.3f}arcsec], {mask_pa:.1f}deg]'

noise_annulus = f"annulus[[{mask_ra}, {mask_dec}],['4.arcsec', '6.arcsec']]"

tclean_kwargs = {
    'specmode':'cube','restoringbeam':'common','nterms':1,
    'start':'-118.4km/s','width':'0.7km/s','nchan':358,
    'deconvolver':'multiscale','niter':50000,
    'outframe':'LSRK','veltype':'radio','restfreq':'219.9494420GHz',
    'usemask':'user','mask':mask,'interactive':False,
}

imagename = vis_SO[:-3]+'.contsub_image_natural_1.0sigma_uvtaper0.05arcsec_0.01gain_8.0ppb_fewchans'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.image.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis       = path_to_vis+vis_SO+'.contsub',
    imagename = imagename,
    scales    = [0,2,4,8],
    cell      = '0.0215arcsec',
    gain      = 0.01,
    weighting = 'natural',
    threshold = '8.85e-01mJy',
    uvtaper   = '0.05arcsec',
    imsize    = 1600,
    **tclean_kwargs
)
estimate_SNR(imagename+'.image',disk_mask=mask,noise_mask=noise_annulus)
for _chans in [170,171,172,173,174,182,183,184,185,186]:
    rms = imstat(imagename=imagename+'.image',region=noise_annulus,chans=f'{_chans}')['rms'][0]
    print('#rms: %.2e Jy/beam' %rms)
#rms: [8.67e-04,8.66e-04,8.83e-04,8.76e-04,8.65e-04,9.18e-04,9.03e-04,8.90e-04,8.93e-04,8.77e-04] Jy/beam
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#CQ_Tau_SBLB_no_ave_selfcal_time_ave_SOwide.contsub_image_natural_1.0sigma_uvtaper0.05arcsec_0.01gain_8.0ppb_fewchans.image
#Beam 0.211 arcsec x 0.172 arcsec (-10.79 deg)
#Flux inside disk mask: 142.26 mJy
#Peak intensity of source: 8.88 mJy/beam
#rms: 8.85e-01 mJy/beam
#Peak SNR: 10.03

#SO - final imaging
SBLB_no_ave_selfcal = f'{prefix}_SBLB_no_ave_selfcal_time_ave.ms'
vis_SO = SBLB_no_ave_selfcal[:-3]+'_SO.ms'

#Define mask
mask_pa        = PA
mask_semimajor = 0.7
mask_semiminor = 0.6
mask_ra  = '05h35m58.472722s'#05:35:58.472722
mask_dec = '24d44m53.613450s'#24:44:53.613450

mask = f'ellipse[[{mask_ra},{mask_dec}], [{mask_semimajor:.3f}arcsec, {mask_semiminor:.3f}arcsec], {mask_pa:.1f}deg]'

noise_annulus = f"annulus[[{mask_ra}, {mask_dec}],['4.arcsec', '6.arcsec']]"

tclean_kwargs = {
    'specmode':'cube','restoringbeam':'common','nterms':1,
    'start':'-2.20km/s','width':'0.7km/s','nchan':24,
    'deconvolver':'multiscale','niter':50000,
    'outframe':'LSRK','veltype':'radio','restfreq':'219.9494420GHz',
    'usemask':'user','mask':mask,'interactive':False,
}

tclean_kwargs_keplmask = {
    'specmode':'cube','restoringbeam':'common','nterms':1,
    'start':'-2.20km/s','width':'0.7km/s','nchan':24,
    'deconvolver':'multiscale','niter':50000,
    'interactive':False,'usemask':'user',
    'outframe':'LSRK','veltype':'radio','restfreq':'219.9494420GHz',
}

imagename = vis_SO[:-3]+'.contsub_image_natural_1.0sigma_uvtaper0.1arcsec_0.01gain_8.0ppb_fewchans_new'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.image.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis       = path_to_vis+vis_SO+'.contsub',
    imagename = imagename,
    scales    = [0,2,4,8],
    cell      = '0.021875arcsec',
    gain      = 0.01,
    weighting = 'natural',
    threshold = '8.28e-01mJy',
    uvtaper   = '0.1arcsec',
    imsize    = 1200,
    **tclean_kwargs
)
estimate_SNR(imagename+'.image',disk_mask=mask,noise_mask=noise_annulus)
for _chans in [4,5,6,7,8,16,17,18,19,20]:
    rms = imstat(imagename=imagename+'.image',region=noise_annulus,chans=f'{_chans}')['rms'][0]
    print('#rms: %.2e Jy/beam' %rms)
#rms: [8.11e-04,8.12e-04,8.21e-04,8.14e-04,8.10e-04,8.51e-04,8.41e-04,8.34e-04,8.39e-04,8.26e-04] Jy/beam
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#CQ_Tau_SBLB_no_ave_selfcal_time_ave_SO.contsub_image_natural_1.0sigma_uvtaper0.1arcsec_0.01gain_8.0ppb_fewchans_new.image
#Beam 0.203 arcsec x 0.175 arcsec (-15.08 deg)
#Flux inside disk mask: 89.91 mJy
#Peak intensity of source: 8.88 mJy/beam
#rms: 8.28e-01 mJy/beam
#Peak SNR: 10.72

imagename = vis_SO[:-3]+'.contsub_image_natural_1.0sigma_uvtaper0.05arcsec_0.01gain_8.0ppb_fewchans_new'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.image.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis       = path_to_vis+vis_SO+'.contsub',
    imagename = imagename,
    scales    = [0,2,4,8],
    cell      = '0.016625arcsec',
    gain      = 0.01,
    weighting = 'natural',
    threshold = '7.44e-01mJy',
    uvtaper   = '0.05arcsec',
    imsize    = 1600,
    **tclean_kwargs
)
estimate_SNR(imagename+'.image',disk_mask=mask,noise_mask=noise_annulus)
for _chans in [4,5,6,7,8,16,17,18,19,20]:
    rms = imstat(imagename=imagename+'.image',region=noise_annulus,chans=f'{_chans}')['rms'][0]
    print('#rms: %.2e Jy/beam' %rms)
#rms: [7.32e-04,7.34e-04,7.38e-04,7.31e-04,7.29e-04,7.55e-04,7.53e-04,7.52e-04,7.59e-04,7.48e-04] Jy/beam
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#CQ_Tau_SBLB_no_ave_selfcal_time_ave_SO.contsub_image_natural_1.0sigma_uvtaper0.05arcsec_0.01gain_8.0ppb_fewchans_new.image
#Beam 0.159 arcsec x 0.133 arcsec (-6.76 deg)
#Flux inside disk mask: 93.74 mJy
#Peak intensity of source: 7.67 mJy/beam
#rms: 7.44e-01 mJy/beam
#Peak SNR: 10.31

imagename = vis_SO[:-3]+'.contsub_image_natural_1.0sigma_0.01gain_8.0ppb_fewchans_new'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.image.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis       = path_to_vis+vis_SO+'.contsub',
    imagename = imagename,
    scales    = [0,2,4,8],
    cell      = '0.013125arcsec',
    gain      = 0.01,
    weighting = 'natural',
    threshold = '7.03e-01mJy',
    imsize    = 2000,
    **tclean_kwargs
)
estimate_SNR(imagename+'.image',disk_mask=mask,noise_mask=noise_annulus)
for _chans in [4,5,6,7,8,16,17,18,19,20]:
    rms = imstat(imagename=imagename+'.image',region=noise_annulus,chans=f'{_chans}')['rms'][0]
    print('#rms: %.2e Jy/beam' %rms)
#rms: [6.94e-04,6.95e-04,6.96e-04,6.89e-04,6.89e-04,7.09e-04,7.11e-04,7.14e-04,7.22e-04,7.14e-04] Jy/beam
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#CQ_Tau_SBLB_no_ave_selfcal_time_ave_SO.contsub_image_natural_1.0sigma_0.01gain_8.0ppb_fewchans_new.image
#Beam 0.126 arcsec x 0.105 arcsec (-2.15 deg)
#Flux inside disk mask: 105.22 mJy
#Peak intensity of source: 6.38 mJy/beam
#rms: 7.03e-01 mJy/beam
#Peak SNR: 9.07

imagename = vis_SO[:-3]+'.contsub_image_1.0robust_1.0sigma_0.01gain_8.0ppb_fewchans_new'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.image.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis       = path_to_vis+vis_SO+'.contsub',
    imagename = imagename,
    scales    = [0,2,4,8],
    cell      = '0.011625arcsec',
    gain      = 0.01,
    weighting = 'briggs',
    robust    = 1.0,
    threshold = '7.14e-01mJy',
    imsize    = 2400,
    **tclean_kwargs
)
estimate_SNR(imagename+'.image',disk_mask=mask,noise_mask=noise_annulus)
for _chans in [4,5,6,7,8,16,17,18,19,20]:
    rms = imstat(imagename=imagename+'.image',region=noise_annulus,chans=f'{_chans}')['rms'][0]
    print('#rms: %.2e Jy/beam' %rms)
#rms: [7.07e-04,7.08e-04,7.05e-04,6.97e-04,6.99e-04,7.17e-04,7.19e-04,7.27e-04,7.37e-04,7.29e-04] Jy/beam
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#CQ_Tau_SBLB_no_ave_selfcal_time_ave_SO.contsub_image_1.0robust_1.0sigma_0.01gain_8.0ppb_fewchans_new.image
#Beam 0.111 arcsec x 0.093 arcsec (0.08 deg)
#Flux inside disk mask: 115.84 mJy
#Peak intensity of source: 5.70 mJy/beam
#rms: 7.14e-01 mJy/beam
#Peak SNR: 7.98

imagename = vis_SO[:-3]+'.contsub_image_0.8robust_1.0sigma_0.01gain_8.0ppb_fewchans_new'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.image.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis       = path_to_vis+vis_SO+'.contsub',
    imagename = imagename,
    scales    = [0,2,4,8],
    cell      = '0.0105arcsec',
    gain      = 0.01,
    weighting = 'briggs',
    robust    = 0.8,
    threshold = '7.32e-01mJy',
    imsize    = 2700,
    **tclean_kwargs
)
estimate_SNR(imagename+'.image',disk_mask=mask,noise_mask=noise_annulus)
for _chans in [4,5,6,7,8,16,17,18,19,20]:
    rms = imstat(imagename=imagename+'.image',region=noise_annulus,chans=f'{_chans}')['rms'][0]
    print('#rms: %.2e Jy/beam' %rms)
#rms: [7.26e-04,7.27e-04,7.24e-04,7.16e-04,7.17e-04,7.34e-04,7.38e-04,7.45e-04,7.56e-04,7.49e-04] Jy/beam
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#CQ_Tau_SBLB_no_ave_selfcal_time_ave_SO.contsub_image_0.8robust_1.0sigma_0.01gain_8.0ppb_fewchans_new.image
#Beam 0.102 arcsec x 0.084 arcsec (0.17 deg)
#Flux inside disk mask: 123.66 mJy
#Peak intensity of source: 5.37 mJy/beam
#rms: 7.32e-01 mJy/beam
#Peak SNR: 7.34

imagename = vis_SO[:-3]+'.contsub_image_0.7robust_1.0sigma_0.01gain_8.0ppb_fewchans_new'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.image.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis       = path_to_vis+vis_SO+'.contsub',
    imagename = imagename,
    scales    = [0,2,4,8],
    cell      = '0.00975arcsec',
    gain      = 0.01,
    weighting = 'briggs',
    robust    = 0.7,
    threshold = '7.47e-01mJy',
    imsize    = 2880,
    **tclean_kwargs
)
estimate_SNR(imagename+'.image',disk_mask=mask,noise_mask=noise_annulus)
for _chans in [4,5,6,7,8,16,17,18,19,20]:
    rms = imstat(imagename=imagename+'.image',region=noise_annulus,chans=f'{_chans}')['rms'][0]
    print('#rms: %.2e Jy/beam' %rms)
#rms: [7.40e-04,7.42e-04,7.36e-04,7.27e-04,7.30e-04,7.47e-04,7.51e-04,7.62e-04,7.76e-04,7.68e-04] Jy/beam
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#CQ_Tau_SBLB_no_ave_selfcal_time_ave_SO.contsub_image_0.7robust_1.0sigma_0.01gain_8.0ppb_fewchans_new.image
#Beam 0.097 arcsec x 0.078 arcsec (0.18 deg)
#Flux inside disk mask: 125.71 mJy
#Peak intensity of source: 5.00 mJy/beam
#rms: 7.47e-01 mJy/beam
#Peak SNR: 6.69

imagename = vis_SO[:-3]+'.contsub_image_0.5robust_1.0sigma_0.01gain_8.0ppb_fewchans_new'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.image.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis       = path_to_vis+vis_SO+'.contsub',
    imagename = imagename,
    scales    = [0,2,4,8],
    cell      = '0.0085arcsec',
    gain      = 0.01,
    weighting = 'briggs',
    robust    = 0.5,
    threshold = '7.93e-01mJy',
    imsize    = 3600,
    **tclean_kwargs
)
estimate_SNR(imagename+'.image',disk_mask=mask,noise_mask=noise_annulus)
for _chans in [4,5,6,7,8,16,17,18,19,20]:
    rms = imstat(imagename=imagename+'.image',region=noise_annulus,chans=f'{_chans}')['rms'][0]
    print('#rms: %.2e Jy/beam' %rms)
#rms: [7.86e-04,7.84e-04,7.80e-04,7.74e-04,7.77e-04,7.93e-04,8.00e-04,8.13e-04,8.24e-04,8.14e-04] Jy/beam
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#CQ_Tau_SBLB_no_ave_selfcal_time_ave_SO.contsub_image_0.5robust_1.0sigma_0.01gain_8.0ppb_fewchans_new.image
#Beam 0.088 arcsec x 0.067 arcsec (0.15 deg)
#Flux inside disk mask: 114.53 mJy
#Peak intensity of source: 4.44 mJy/beam
#rms: 7.93e-01 mJy/beam
#Peak SNR: 5.60

imagename = vis_SO[:-3]+'.contsub_image_natural_1.0sigma_uvtaper0.1arcsec_0.01gain_8.0ppb_fewchans_new_keplmask'
make_mask(image=imagename[:-9]+'.image', inc=-36.2, PA=235., mstar=1.4, dist=149., vlsr=6.2e3, nbeams=1.25, r_min=0.0, r_max=1.0)
keplmask = imagename[:-9]+'.mask.image'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.image.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis       = path_to_vis+vis_SO+'.contsub',
    imagename = imagename,
    scales    = [0,2,4],
    cell      = '0.021875arcsec',
    gain      = 0.01,
    weighting = 'natural',
    threshold = '8.29e-01mJy',
    uvtaper   = '0.1arcsec',
    mask      = keplmask,
    imsize    = 1200,
    **tclean_kwargs_keplmask
)
estimate_SNR(imagename+'.image',disk_mask=mask,noise_mask=noise_annulus)
for _chans in [4,5,6,7,8,16,17,18,19,20]:
    rms = imstat(imagename=imagename+'.image',region=noise_annulus,chans=f'{_chans}')['rms'][0]
    print('#rms: %.2e Jy/beam' %rms)
#rms: [8.14e-04,8.15e-04,8.24e-04,8.16e-04,8.11e-04,8.51e-04,8.42e-04,8.36e-04,8.42e-04,8.28e-04] Jy/beam
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#CQ_Tau_SBLB_no_ave_selfcal_time_ave_SO.contsub_image_natural_1.0sigma_uvtaper0.1arcsec_0.01gain_8.0ppb_fewchans_new_keplmask.image
#Beam 0.203 arcsec x 0.175 arcsec (-15.08 deg)
#Flux inside disk mask: 92.90 mJy
#Peak intensity of source: 8.94 mJy/beam
#rms: 8.29e-01 mJy/beam
#Peak SNR: 10.78

imagename = vis_SO[:-3]+'.contsub_image_natural_1.0sigma_uvtaper0.05arcsec_0.01gain_8.0ppb_fewchans_new_keplmask'
make_mask(image=imagename[:-9]+'.image', inc=-36.2, PA=235., mstar=1.4, dist=149., vlsr=6.2e3, nbeams=1.25, r_min=0.0, r_max=1.0)
keplmask = imagename[:-9]+'.mask.image'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.image.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis       = path_to_vis+vis_SO+'.contsub',
    imagename = imagename,
    scales    = [0,2,4],
    cell      = '0.016625arcsec',
    gain      = 0.01,
    weighting = 'natural',
    threshold = '7.45e-01mJy',
    uvtaper   = '0.05arcsec',
    mask      = keplmask,
    imsize    = 1600,
    **tclean_kwargs_keplmask
)
estimate_SNR(imagename+'.image',disk_mask=mask,noise_mask=noise_annulus)
for _chans in [4,5,6,7,8,16,17,18,19,20]:
    rms = imstat(imagename=imagename+'.image',region=noise_annulus,chans=f'{_chans}')['rms'][0]
    print('#rms: %.2e Jy/beam' %rms)
#rms: [7.34e-04,7.36e-04,7.40e-04,7.32e-04,7.30e-04,7.56e-04,7.53e-04,7.53e-04,7.61e-04,7.50e-04] Jy/beam
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#CQ_Tau_SBLB_no_ave_selfcal_time_ave_SO.contsub_image_natural_1.0sigma_uvtaper0.05arcsec_0.01gain_8.0ppb_fewchans_new_keplmask.image
#Beam 0.159 arcsec x 0.133 arcsec (-6.76 deg)
#Flux inside disk mask: 98.51 mJy
#Peak intensity of source: 7.69 mJy/beam
#rms: 7.45e-01 mJy/beam
#Peak SNR: 10.33

imagename = vis_SO[:-3]+'.contsub_image_natural_1.0sigma_0.01gain_8.0ppb_fewchans_new_keplmask'
make_mask(image=imagename[:-9]+'.image', inc=-36.2, PA=235., mstar=1.4, dist=149., vlsr=6.2e3, nbeams=1.25, r_min=0.0, r_max=1.0)
keplmask = imagename[:-9]+'.mask.image'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.image.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis       = path_to_vis+vis_SO+'.contsub',
    imagename = imagename,
    scales    = [0,2,4],
    cell      = '0.013125arcsec',
    gain      = 0.01,
    weighting = 'natural',
    threshold = '7.04e-01mJy',
    mask      = keplmask,
    imsize    = 2000,
    **tclean_kwargs_keplmask
)
estimate_SNR(imagename+'.image',disk_mask=mask,noise_mask=noise_annulus)
for _chans in [4,5,6,7,8,16,17,18,19,20]:
    rms = imstat(imagename=imagename+'.image',region=noise_annulus,chans=f'{_chans}')['rms'][0]
    print('#rms: %.2e Jy/beam' %rms)
#rms: [6.95e-04,6.96e-04,6.97e-04,6.90e-04,6.89e-04,7.09e-04,7.11e-04,7.15e-04,7.24e-04,7.15e-04] Jy/beam
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#CQ_Tau_SBLB_no_ave_selfcal_time_ave_SO.contsub_image_natural_1.0sigma_0.01gain_8.0ppb_fewchans_new_keplmask.image
#Beam 0.126 arcsec x 0.105 arcsec (-2.15 deg)
#Flux inside disk mask: 101.62 mJy
#Peak intensity of source: 6.44 mJy/beam
#rms: 7.04e-01 mJy/beam
#Peak SNR: 9.15

imagename = vis_SO[:-3]+'.contsub_image_1.0robust_1.0sigma_0.01gain_8.0ppb_fewchans_new_keplmask'
make_mask(image=imagename[:-9]+'.image', inc=-36.2, PA=235., mstar=1.4, dist=149., vlsr=6.2e3, nbeams=1.25, r_min=0.0, r_max=1.0)
keplmask = imagename[:-9]+'.mask.image'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.image.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis       = path_to_vis+vis_SO+'.contsub',
    imagename = imagename,
    scales    = [0,2,4],
    cell      = '0.011625arcsec',
    gain      = 0.01,
    weighting = 'briggs',
    robust    = 1.0,
    threshold = '7.14e-01mJy',
    mask      = keplmask,
    imsize    = 2400,
    **tclean_kwargs_keplmask
)
estimate_SNR(imagename+'.image',disk_mask=mask,noise_mask=noise_annulus)
for _chans in [4,5,6,7,8,16,17,18,19,20]:
    rms = imstat(imagename=imagename+'.image',region=noise_annulus,chans=f'{_chans}')['rms'][0]
    print('#rms: %.2e Jy/beam' %rms)
#rms: [7.08e-04,7.09e-04,7.06e-04,6.98e-04,6.99e-04,7.17e-04,7.20e-04,7.28e-04,7.38e-04,7.31e-04] Jy/beam
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#CQ_Tau_SBLB_no_ave_selfcal_time_ave_SO.contsub_image_1.0robust_1.0sigma_0.01gain_8.0ppb_fewchans_new_keplmask.image
#Beam 0.111 arcsec x 0.093 arcsec (0.08 deg)
#Flux inside disk mask: 116.75 mJy
#Peak intensity of source: 5.68 mJy/beam
#rms: 7.14e-01 mJy/beam
#Peak SNR: 7.95

imagename = vis_SO[:-3]+'.contsub_image_0.8robust_1.0sigma_0.01gain_8.0ppb_fewchans_new_keplmask'
make_mask(image=imagename[:-9]+'.image', inc=-36.2, PA=235., mstar=1.4, dist=149., vlsr=6.2e3, nbeams=1.25, r_min=0.0, r_max=1.0)
keplmask = imagename[:-9]+'.mask.image'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.image.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis       = path_to_vis+vis_SO+'.contsub',
    imagename = imagename,
    scales    = [0,2,4],
    cell      = '0.0105arcsec',
    gain      = 0.01,
    weighting = 'briggs',
    robust    = 0.8,
    threshold = '7.32e-01mJy',
    mask      = keplmask,
    imsize    = 2880,
    **tclean_kwargs_keplmask
)
estimate_SNR(imagename+'.image',disk_mask=mask,noise_mask=noise_annulus)
for _chans in [4,5,6,7,8,16,17,18,19,20]:
    rms = imstat(imagename=imagename+'.image',region=noise_annulus,chans=f'{_chans}')['rms'][0]
    print('#rms: %.2e Jy/beam' %rms)
#rms: [7.27e-04,7.26e-04,7.23e-04,7.16e-04,7.17e-04,7.34e-04,7.37e-04,7.48e-04,7.59e-04,7.49e-04] Jy/beam
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#CQ_Tau_SBLB_no_ave_selfcal_time_ave_SO.contsub_image_0.8robust_1.0sigma_0.01gain_8.0ppb_fewchans_new_keplmask.image
#Beam 0.102 arcsec x 0.084 arcsec (0.17 deg)
#Flux inside disk mask: 115.43 mJy
#Peak intensity of source: 5.01 mJy/beam
#rms: 7.32e-01 mJy/beam
#Peak SNR: 6.84

imagename = vis_SO[:-3]+'.contsub_image_0.7robust_1.0sigma_0.01gain_8.0ppb_fewchans_new_keplmask'
make_mask(image=imagename[:-9]+'.image', inc=-36.2, PA=235., mstar=1.4, dist=149., vlsr=6.2e3, nbeams=2.0, r_min=0.0, r_max=1.0)
keplmask = imagename[:-9]+'.mask.image'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.image.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis       = path_to_vis+vis_SO+'.contsub',
    imagename = imagename,
    scales    = [0,2,4],
    cell      = '0.00975arcsec',
    gain      = 0.01,
    weighting = 'briggs',
    robust    = 0.7,
    threshold = '7.48e-01mJy',
    mask      = keplmask,
    imsize    = 2880,
    **tclean_kwargs_keplmask
)
estimate_SNR(imagename+'.image',disk_mask=mask,noise_mask=noise_annulus)
for _chans in [4,5,6,7,8,16,17,18,19,20]:
    rms = imstat(imagename=imagename+'.image',region=noise_annulus,chans=f'{_chans}')['rms'][0]
    print('#rms: %.2e Jy/beam' %rms)
#rms: [7.42e-04,7.43e-04,7.37e-04,7.28e-04,7.31e-04,7.47e-04,7.52e-04,7.63e-04,7.77e-04,7.69e-04] Jy/beam
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#CQ_Tau_SBLB_no_ave_selfcal_time_ave_SO.contsub_image_0.7robust_1.0sigma_0.01gain_8.0ppb_fewchans_new_keplmask.image
#Beam 0.097 arcsec x 0.078 arcsec (0.18 deg)
#Flux inside disk mask: 130.40 mJy
#Peak intensity of source: 4.80 mJy/beam
#rms: 7.48e-01 mJy/beam
#Peak SNR: 6.42

imagename = vis_SO[:-3]+'.contsub_image_0.5robust_1.0sigma_0.01gain_8.0ppb_fewchans_new_keplmask'
make_mask(image=imagename[:-9]+'.image', inc=-36.2, PA=235., mstar=1.4, dist=149., vlsr=6.2e3, nbeams=2.0, r_min=0.0, r_max=1.0)
keplmask = imagename[:-9]+'.mask.image'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.image.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis       = path_to_vis+vis_SO+'.contsub',
    imagename = imagename,
    scales    = [0,2],
    cell      = '0.008375arcsec',
    gain      = 0.01,
    weighting = 'briggs',
    robust    = 0.5,
    threshold = '7.94e-01mJy',
    mask      = keplmask,
    imsize    = 3600,
    **tclean_kwargs_keplmask
)
estimate_SNR(imagename+'.image',disk_mask=mask,noise_mask=noise_annulus)
for _chans in [4,5,6,7,8,16,17,18,19,20]:
    rms = imstat(imagename=imagename+'.image',region=noise_annulus,chans=f'{_chans}')['rms'][0]
    print('#rms: %.2e Jy/beam' %rms)
#rms: [7.91e-04,7.86e-04,7.78e-04,7.73e-04,7.78e-04,7.92e-04,7.99e-04,8.11e-04,8.25e-04,8.17e-04] Jy/beam
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#CQ_Tau_SBLB_no_ave_selfcal_time_ave_SO.contsub_image_0.5robust_1.0sigma_0.01gain_8.0ppb_fewchans_new_keplmask.image
#Beam 0.088 arcsec x 0.068 arcsec (0.15 deg)
#Flux inside disk mask: 109.03 mJy
#Peak intensity of source: 4.97 mJy/beam
#rms: 7.94e-01 mJy/beam
#Peak SNR: 6.27

imagenames = [
    'natural_1.0sigma_uvtaper0.1arcsec_0.01gain_8.0ppb_fewchans_new',
    'natural_1.0sigma_uvtaper0.05arcsec_0.01gain_8.0ppb_fewchans_new',
    'natural_1.0sigma_0.01gain_8.0ppb_fewchans_new',
    '1.0robust_1.0sigma_0.01gain_8.0ppb_fewchans_new',
    '0.8robust_1.0sigma_0.01gain_8.0ppb_fewchans_new',
    '0.7robust_1.0sigma_0.01gain_8.0ppb_fewchans_new',
    '0.5robust_1.0sigma_0.01gain_8.0ppb_fewchans_new',
    'natural_1.0sigma_uvtaper0.1arcsec_0.01gain_8.0ppb_fewchans_new_keplmask',
    'natural_1.0sigma_uvtaper0.05arcsec_0.01gain_8.0ppb_fewchans_new_keplmask',
    'natural_1.0sigma_0.01gain_8.0ppb_fewchans_new_keplmask',
    '1.0robust_1.0sigma_0.01gain_8.0ppb_fewchans_new_keplmask',
    '0.8robust_1.0sigma_0.01gain_8.0ppb_fewchans_new_keplmask',
    '0.7robust_1.0sigma_0.01gain_8.0ppb_fewchans_new_keplmask',
    '0.5robust_1.0sigma_0.01gain_8.0ppb_fewchans_new_keplmask',
]

for _fname in imagenames:
    imagename = f'./integrated_intensity_8.0ppb_fewchans_new/CQ_Tau_SBLB_no_ave_selfcal_time_ave_SO.contsub_image_{_fname}.image_M0.fits'
    rms = imstat(imagename=imagename,region=noise_annulus)['rms'][0] #Jy/beam m/s
    print('#rms_ii: %.2e Jy/beam m/s' %rms)

#rms_bestchans: 4.67e+00,4.22e+00,3.99e+00,4.07e+00,4.17e+00,4.26e+00,4.52e+00,4.68e+00,4.22e+00,4.00e+00,4.07e+00,4.18e+00,4.26e+00,4.54e+00 Jy/beam m/s

for _fname in imagenames:
    imagename = f'./peak_intensity_8.0ppb_fewchans_new/CQ_Tau_SBLB_no_ave_selfcal_time_ave_SO.contsub_image_{_fname}.image_M8.fits'
    rms = imstat(imagename=imagename,region=noise_annulus)['rms'][0] #Jy/beam
    print('#rms_pi: %.2e Jy/beam' %rms)

#rms_bestchans: 1.52e-03,1.37e-03,1.29e-03,1.31e-03,1.34e-03,1.37e-03,1.45e-03,1.53e-03,1.37e-03,1.29e-03,1.31e-03,1.34e-03,1.37e-03,1.45e-03 Jy/beam