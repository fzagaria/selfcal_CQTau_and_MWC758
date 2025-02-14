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
matplotlib.rcParams.update({'figure.max_open_warning':120})

#path to your local copy of the gihub repository
#(make sure to do a 'git pull' to have the most up-to-date version)
github_path = '/data/beegfs/astro-storage/groups/benisty/frzagaria/mycasa_versions/'
import sys
sys.path.append(github_path)
import alignment_default as alignment
execfile(os.path.join(github_path,'reduction_utils_py3_mpi.py'))

prefix = 'MWC_758'

# System properties.
incl  = 21.   # deg, from Dong et al. 2018
PA    = 62.   # deg, from Dong et al. 2018
v_sys =  5.9  # km/s, from listobs

# Whether to run tclean in parallel or not.
use_parallel = False

data_folderpath = '/data/beegfs/astro-storage/groups/benisty/frzagaria/SO_detections/MWC758/MWC758_Band6_data/'

# data_folderpath = '/lustre/cv/projects/exoALMA/ALMA_PL_calibrated_data/PDS_66'
# # TM_names = {'LB':'TM1','SB':'TM2'}
# PL_calibrated_vis = {key:os.path.join(data_folderpath,f'{TM_name}/calibrated_final.ms')
#                     for key,TM_name in TM_names.items()}
#
# for baseline_key,TM_name in TM_names.items():
#     listobs(
#         vis=PL_calibrated_vis[baseline_key],
#         listfile=f'{prefix}_{TM_name}_calibrated_final.ms.txt',
#         overwrite=True,
#     )

#2017.1.00940.S: ang. res. = 0.128", spec. res. = 1.468 km/s, cont. sens. = 0.0199 mJy/beam (SB0 to 3)
path_2017_1_00940_S_SBs = './2017.1.00940.S/science_goal.uid___A001_X1273_X416/group.uid___A001_X1273_X417/member.uid___A001_X1273_X418/calibrated/'
data_2017_1_00940_S_SBs = ['uid___A002_Xc7b4ea_X18e9.ms.split.cal','uid___A002_Xc805c1_X18d7.ms.split.cal','uid___A002_Xc8592e_Xbe8f.ms.split.cal','uid___A002_Xc86fe5_X24ab.ms.split.cal']

#2017.1.00940.S: ang. res. = 0.022", spec. res. = 1.468 km/s, cont. sens. = 0.0166 mJy/beam (LB0 to 4)
path_2017_1_00940_S_LBs = './2017.1.00940.S/science_goal.uid___A001_X1273_X416/group.uid___A001_X1273_X417/member.uid___A001_X1273_X414/calibrated/'
data_2017_1_00940_S_LBs = ['uid___A002_Xc59134_X24e9.ms.split.cal','uid___A002_Xc59134_X2987.ms.split.cal','uid___A002_Xc59134_X642c.ms.split.cal','uid___A002_Xc5b7d7_X2014.ms.split.cal','uid___A002_Xc5b7d7_X63c7.ms.split.cal']

PL_calibrated_path = [path_2017_1_00940_S_SBs,path_2017_1_00940_S_LBs]

PL_calibrated_data = [data_2017_1_00940_S_SBs,data_2017_1_00940_S_LBs]

PL_calibrated_name = [['SB_EB0','SB_EB1','SB_EB2','SB_EB3'],['LB_EB0','LB_EB1','LB_EB2','LB_EB3','LB_EB4']]

field_name    = ['MWC758', 'MWC758']

spw_name      = ['19,23,25,27', '19,23,25,27']

#Create initial .ms files to work on.
#Splitting also for spws is crucial, otherwise the pipeline-calibrated .ms retain the full spw structure (with 10s of IDs) that avg_cont doesn't like!
for _path, _data, _name, _field, _spw in zip(PL_calibrated_path, PL_calibrated_data, PL_calibrated_name, field_name, spw_name):
    for i in range(len(_data)):
        listobs(
            vis       = os.path.join(os.path.join(data_folderpath,_path),_data[i]),
            listfile  = f'{data_folderpath}'+f'{_path}'+f'{_data[i]}'+'.listobs.txt',
            overwrite = True,
        )

        os.system(f'rm -rf {prefix}_{_name[i]}.ms')
        split(
            vis        = os.path.join(os.path.join(data_folderpath,_path),_data[i]),
            outputvis  = f'{prefix}_{_name[i]}.ms',
            datacolumn = 'DATA',
            field      = _field,
            spw        = _spw, #in some scans of the SBs the target was observed for atm calibration, good to keep?
                #Choose the spwIDs corresponidng to the 'source' entry in listobs already at this stage.
                #It is not compulsory, but otherwise the avg_cont function retains the full spw structure of your .ms table
                #i.e., say you have up to spw 25, with only 17,19,23,25 for your target 'source', 
                #when you only give avg_cont a maxchanwidth, it will get a width_array for spw 0 to 25
                #sth that the split function within avg_cont doesn't like (only the target 'source' are required).
                #Of course, if you skip this step, you can provide directly the contspws and width_array to avg_cont.
            intent     = 'OBSERVE_TARGET#ON_SOURCE',
            keepflags  = False, 
        )

        listobs(
            vis       = f'{prefix}_{_name[i]}.ms',
            listfile  = f'{prefix}_{_name[i]}.ms.listobs.txt',
            overwrite = True,
        )

# spw_sets = {
#     'SB':['25,27,29,31','93,95,97,99'],
#     'LB':['25,27,29,31','92,94,96,98','112,114,116,118','132,134,136,138']
# }
# fields = {'SB':'PDS_66','LB':'PDS_66'}
#
# for baseline_key,vis in PL_calibrated_vis.items():
#     split_all_obs(msfile=vis,nametemplate=f'{prefix}_{baseline_key}_EB')
#
# number_of_EBs = {baseline_key:len(spw_set) for baseline_key,spw_set in spw_sets.items()}
#
# for baseline_key,n_EB in number_of_EBs.items():
#     for i in range(n_EB):
#         vis = f'{prefix}_{baseline_key}_EB{i}.ms'
#         listobs(
#             vis=vis,
#             listfile=f'{vis}.txt',
#             overwrite=True,
#         )

#Check that spws and field numbering is what you expect them to be now from the listobs files
#In particular, make sure that there are 4 spws and that there is a single field with field ID = 0

# rest_freq_SiS_1211    = 217.8176630e9 #J=12-11
# rest_freq_SiS_1312    = 235.9613634e9 #J=13-12

rest_freq_H2CO_918_919 = 216.5686510e9 #Spw_0 (216.52 to 218.50)
rest_freq_DCN_F21      = 217.2384000e9 #J=3-2
rest_freq_DCN_F22      = 217.2386307e9 #J=3-2
rest_freq_H2CO_303_202 = 218.2221920e9
rest_freq_H2CO_322_221 = 218.4756320e9
rest_freq_H2CO_321_220 = 218.7600660e9 #Spw_1 (218.53 to 220.40)
rest_freq_C18O         = 219.5603541e9 #J=2-1
rest_freq_SO           = 219.9494420e9 #3Sigma 6(5)-5(4)
rest_freq_13CO         = 220.3986842e9 #J=2-1 ignoring fine structure
rest_freq_12CO         = 230.5380000e9 #J=2-1 Spw2 (230.53 to 232.40)
# continuum                              #Spw3 (233.21 to 235.19)

#SBs and LBs
#0: 233.35, 235.10 continuum
#1: 216.60, 218.40 rest_freq_H2CO_918_919,rest_freq_DCN_F21,rest_freq_DCN_F22,rest_freq_H2CO_303_202,rest_freq_H2CO_322_221
#2: 230.50, 232.40 rest_freq_12CO
#3: 218.50, 220.40 rest_freq_H2CO_322_221,rest_freq_H2CO_321_220,rest_freq_C18O,rest_freq_SO,rest_freq_13CO

number_of_EBs  = {'SB':4,'LB':5}

data_params_LB = {
    f'LB{i}': {
        'vis':          f'{prefix}_LB_EB{i}.ms',
        'name':         f'LB_EB{i}',
        'field':        'MWC758',
        'line_spws':    np.array([1,1,1,1,1, 2, 3,3,3,3,3]), #list of spws containing lines
        'line_freqs':   np.array([
            rest_freq_H2CO_918_919,rest_freq_DCN_F21,rest_freq_DCN_F22,rest_freq_H2CO_303_202,rest_freq_H2CO_322_221,
            rest_freq_12CO,
            rest_freq_H2CO_322_221,rest_freq_H2CO_321_220,rest_freq_C18O,rest_freq_SO,rest_freq_13CO
        ]), #frequencies (Hz) corresponding to line_spws
        'spwcont_forplot': '0',
        'cont_spws':       '0,1,2,3',
        'width_array':     [2,2,32,32],
    } for i in range(number_of_EBs['LB']) #Phase RMS 13.935, 15.003, 21.429, 32.285, 34.154
}

data_params_SB = {
    f'SB{i}': {
        'vis':          f'{prefix}_SB_EB{i}.ms',
        'name':         f'SB_EB{i}',
        'field':        'MWC758',
        'line_spws':    np.array([1,1,1,1,1, 2, 3,3,3,3,3]), #list of spws containing lines
        'line_freqs':   np.array([
            rest_freq_H2CO_918_919,rest_freq_DCN_F21,rest_freq_DCN_F22,rest_freq_H2CO_303_202,rest_freq_H2CO_322_221,
            rest_freq_12CO,
            rest_freq_H2CO_322_221,rest_freq_H2CO_321_220,rest_freq_C18O,rest_freq_SO,rest_freq_13CO
        ]), #frequencies (Hz) corresponding to line_spws
        'spwcont_forplot': '0',
        'cont_spws':       '0,1,2,3',
        'width_array':     [8,8,192,192],
    } for i in range(number_of_EBs['SB']) #Phase RMS 38.807, 20.869, 24.399, 16.723
}

data_params = data_params_LB.copy()
data_params.update(data_params_SB)

figures_foldername = 'figures'
os.mkdir(figures_foldername)

def get_figures_folderpath(foldername):
    return os.path.join(figures_foldername,foldername)

def make_figures_folder(folderpath):
    if os.path.isdir(folderpath):
        print(f'going to deleted folder {folderpath} and its content')
        valid_answer = 'yes'
        answer = input(f'to confirm, type \'{valid_answer}\': ')
        if answer == valid_answer:
            shutil.rmtree(folderpath)
        else:
            print('aborting')
            return
    os.mkdir(folderpath)
    return folderpath

preselfcal_amp_figures_folder = get_figures_folderpath('1_preselfcal_amp_figures')
make_figures_folder(preselfcal_amp_figures_folder)

#adjust these plot ranges according to your data, check from weblogs or QA2 report
plotranges = {'SB':[0,3150,0,0.1], 'LB':[0,16200,0,0.1]} #xmin,xmax,ymin,ymax 

for params in data_params.values():
    
    plot_filename   = prefix+'_'+params['name']+'_amp-v-freq_preselfcal.png'
    plotms(
        vis         = params['vis'],
        xaxis       = 'freq',
        yaxis       = 'amp',
        field       = params['field'],
        ydatacolumn = 'data',
        avgtime     = '1e8',
        avgscan     = True,
        avgbaseline = False,
        iteraxis    = 'spw',
        coloraxis   = 'corr',
        showgui     = False,
        exprange    = 'all',
        plotfile    = os.path.join(preselfcal_amp_figures_folder,plot_filename),
        freqframe   = 'LSRK',
        # yselfscale  = True,
        showatm     = True,
        overwrite   = True,
    )

    plot_filename   = prefix+'_'+params['name']+'_amp-v-chan_preselfcal.png'
    plotms(
        vis         = params['vis'],
        xaxis       = 'channel',
        yaxis       = 'amp',
        field       = params['field'],
        ydatacolumn = 'data',
        avgtime     = '1e8',
        avgscan     = True,
        avgbaseline = True,
        iteraxis    = 'spw',
            #There is a bug in plotms that plots the data wrongly if there is unequal flagging
            #between the polarisations
            #the bug occurs if we plot amp and we average over baselines
            #thus we color by polarization such that
            #we easily identify this issue
        coloraxis   = 'corr',
        showgui     = False,
        exprange    = 'all',
        plotfile    = os.path.join(preselfcal_amp_figures_folder,plot_filename),
        # yselfscale  = True,
        # showatm     = True,
        overwrite   = True,
    )
    
    baseline_key, _ = params['name'].split('_')
    plot_filename   = prefix+'_'+params['name']+'_amp-v-uvdist_cont_spw_preselfcal.png'
    plotms(
        vis         = params['vis'],
        xaxis       = 'UVdist',
        yaxis       = 'amp',
        spw         = params['spwcont_forplot'],
        field       = params['field'],
        ydatacolumn = 'data',
        avgtime     = '1e8',
        avgscan     = True,
        avgchannel  = '3840',
        coloraxis   = 'corr',
        showgui     = False,
        plotrange   = plotranges[baseline_key],
        plotfile    = os.path.join(preselfcal_amp_figures_folder,plot_filename),
        # yselfscale  = True,
        overwrite   = True,
    )

#For PDS 66 everything looks fine, but I'll also try the other plots with single polarizations
#
#For CQ Tau, the plotms polarisation bug seems indeed there in LB EB1 spw3, so let's plot
#polarizations separately to verify
#for corr in ('XX','YY'):
#    plot_filename = prefix+'_'+data_params['LB1']['name']+f'_chan-v-amp_preselfcal_pol{corr}_spw3.png'
#    plotms(
#        vis         = data_params['LB1']['vis'],
#        xaxis       = 'channel',
#        yaxis       = 'amp',
#        field       = data_params['LB1']['field'],
#        ydatacolumn = 'data',
#        avgtime     = '1e8',
#        avgscan     = True,
#        avgbaseline = True,
#        spw         = '3',
#        correlation = corr,
#        showgui     = False,
#        exprange    = 'all',
#        plotfile    = os.path.join(preselfcal_amp_figures_folder,plot_filename)
#    )
#when plotting polarizations separately, everything looks fine also here

#Additional step: flag the line-contaminated channels, but do not average spectrally.
#Then, plot amplitude vs frequency to check if flagging worked fine.
for params in data_params.values():
    os.system(f'rm -rf '+prefix+'_flagtest_'+params['name']+'_initcont.ms')

    baseline_key, _ = params['name'].split('_')
    _, idx_key      = params['name'].split('EB')
    contspws        = params['cont_spws']

    flagchannels_string = get_flagchannels(ms_dict=params,output_prefix=prefix,velocity_range=np.array([-15.,15.]) + v_sys)
    # Flagchannels input string for LB_EB0: '1:123~125, 1:80~82, 1:80~82, 1:17~19, 1:1~3, 2:0~17, 3:1919~1919, 3:1676~1698, 3:856~879, 3:458~480, 3:0~20'
    # Flagchannels input string for LB_EB1: '1:123~125, 1:80~82, 1:80~82, 1:17~19, 1:1~3, 2:0~17, 3:1919~1919, 3:1676~1698, 3:856~879, 3:458~480, 3:0~20'
    # Flagchannels input string for LB_EB2: '1:123~125, 1:80~82, 1:80~82, 1:17~19, 1:1~3, 2:0~17, 3:1919~1919, 3:1676~1698, 3:856~879, 3:458~480, 3:0~20'
    # Flagchannels input string for LB_EB3: '1:123~125, 1:80~82, 1:80~82, 1:17~19, 1:1~3, 2:0~17, 3:1919~1919, 3:1676~1698, 3:856~879, 3:458~480, 3:0~20'
    # Flagchannels input string for LB_EB4: '1:123~125, 1:80~82, 1:80~82, 1:17~19, 1:1~3, 2:0~17, 3:1919~1919, 3:1676~1698, 3:856~879, 3:458~480, 3:0~20'
    # Flagchannels input string for SB_EB0: '1:123~125, 1:80~82, 1:80~82, 1:17~19, 1:1~3, 2:0~16, 3:1919~1919, 3:1675~1698, 3:856~878, 3:457~480, 3:0~20'
    # Flagchannels input string for SB_EB1: '1:123~125, 1:80~82, 1:80~82, 1:17~19, 1:1~3, 2:0~16, 3:1919~1919, 3:1675~1697, 3:856~878, 3:457~480, 3:0~20'
    # Flagchannels input string for SB_EB2: '1:123~125, 1:80~82, 1:80~82, 1:17~19, 1:1~3, 2:0~16, 3:1919~1919, 3:1675~1697, 3:855~878, 3:457~480, 3:0~20'
    # Flagchannels input string for SB_EB3: '1:123~125, 1:80~82, 1:80~82, 1:17~19, 1:1~3, 2:0~16, 3:1919~1919, 3:1675~1697, 3:855~878, 3:457~480, 3:0~20'

    avg_cont(ms_dict=params,output_prefix=prefix+'_flagtest',flagchannels=flagchannels_string,contspws=contspws,width_array=np.ones(len(contspws.split(",")),dtype=int))

    listobs(
        vis       = prefix+'_flagtest_'+params['name']+'_initcont.ms',
        listfile  = prefix+'_flagtest_'+params['name']+'_initcont.ms.listobs.txt',
        overwrite = True,
    )

    plot_filename   = prefix+'_flagtest_'+params['name']+'_amp-v-freq_preselfcal.png'
    plotms(
        vis         = prefix+'_flagtest_'+params['name']+'_initcont.ms',
        xaxis       = 'freq',
        yaxis       = 'amp',
        field       = params['field'],
        ydatacolumn = 'data',
        avgtime     = '1e8',
        avgscan     = True,
        avgbaseline = False,
        iteraxis    = 'spw',
        coloraxis   = 'corr',
        showgui     = False,
        exprange    = 'all',
        plotfile    = os.path.join(preselfcal_amp_figures_folder,plot_filename),
        freqframe   = 'LSRK',
        # yselfscale  = True,
        showatm     = True,
        overwrite   = True,
    )

#Check the maximum channel width you can average over to avoid bandwidth smearing. Depends on the maximum baseline and channel frequency.
#Chose a conservative reduction of the peak response to a point source at the edge of the primary beam of <1%
def delta_freq_smearing(max_baseline,ref_freq,reduction_fac=0.01,antenna_diam=12.):
    """
    Returns maximum channel size to avoid bandwidth smearing by (reduction_fac)*100% at the edge of the primary beam in MHz.
    Parameters:
    max_baseline:  maximum baseline data (m) use au.getBaselineLenghts(msfile)
    ref_freq:      reference channel frequency (GHz)
    reduction_fac: response reduction factor, default is 1%
    antenna_diam:  antenna diameter (m), default is 12 m for ALMA
    Returns:
    maximum channel width consistent with <(reduction_fac)*100% bandwidth smearing at the edge of the primary beam in MHz.
    """
    R = 1. - reduction_fac
    beta_max = np.sqrt(1./R**2 - 1.)
    
    return beta_max * 2 * np.sqrt(np.log(2)) * antenna_diam * ref_freq * 1e3 / max_baseline

max_chan_width = []
for params in data_params.values():
    baseline_key, _ = params['name'].split('_')
    _, idx_key      = params['name'].split('EB')
    contspws        = np.array(params['cont_spws'].split(","),dtype=int)
    
    chanfreqs = []
    for spw in contspws:
        tb.open(params['vis']+'/SPECTRAL_WINDOW')
        chanfreqs.append(np.amin(tb.getcol('CHAN_FREQ', startrow=spw, nrow=1)[:,0]))
        tb.close()

    ref_freq     = np.amin(chanfreqs)
    max_baseline = au.getBaselineLengths(params['vis'])[-1][1]
    
    max_chan_width.append(delta_freq_smearing(max_baseline=max_baseline,ref_freq=ref_freq/1e9))

print(max_chan_width) #38.064, 38.064, 38.064, 38.064, 38.064, 196.085, 244.922, 244.922, 244.922
print(f'These datasets require a maximum channel width of {np.amin(max_chan_width)} MHz to avoid bandwidth smearing.')

#LBs
#0,  128, 15625.000kHz   2 31.25MHz 
#1,  128, 15625.000kHz   2 31.25MHz
#2, 1920,   976.562kHz  32 31.25MHz
#3, 1920,   976.562kHz  32 31.25MHz

#SBs
#0,  128, 15625.000kHz   8 125.00MHz   8 125.00MHz
#1,  128, 15625.000kHz   8 125.00MHz   8 125.00MHz
#2, 1920,   976.562kHz 192 187.50MHz 240 234.37MHz
#3, 1920,   976.562kHz 192 187.50MHz 240 234.37MHz

#Do the averaging (using the width_array corresponding to the minimum between the maximum channel width above and 250 MHz, listed above)
for params in data_params.values():
    os.system(f'rm -rf '+prefix+'_'+params['name']+'_initcont.ms')

    flagchannels_string = get_flagchannels(ms_dict=params,output_prefix=prefix,velocity_range=np.array([-15.,15.]) + v_sys)
    #Double-check that the channels idenfied are at the center of the spws, due to a potential issue with the data_desc_id key in the ms table of some programs

    avg_cont(ms_dict=params,output_prefix=prefix,flagchannels=flagchannels_string,contspws=params['cont_spws'],width_array=params['width_array'])
        #In principle you can just provide a maxchanwidth as in exoALMA, However, some of the channels can be dropped if 
        #the number of channels is not an integer multiple of the width_array computed by avg_cont for the given maxchanwidth. e.g.,

        #Artificially dropping some of the channels because width_array averaging is by 6 channels which is not a factor of 128. 
        #Don't like it, provide width_array manually set on the max_chan_width just printed (sad but different for every EB...)
        #WARN    MSTransformManager::dropNonUniformWidthChannels       Not enough input channels to populate output channel n# 21 from SPW 0.
        #WARN    MSTransformManager::dropNonUniformWidthChannels+      The channel will be dropped in order to have an uniform grid.
        #SpwID  Name                                      #Chans   Frame   Ch0(MHz)  ChanWid(kHz)  TotBW(kHz) CtrFreq(MHz) BBC Num  Corrs  
        #0      X339763914#ALMA_RB_07#BB_1#SW-01#FULL_RES    128   TOPO  356222.334     15625.000   2000000.0 357214.5212        1  XX  YY
        #TotBW(kHz) CtrFreq(MHz) BBC Num  Corrs  
        #0      X339763914#ALMA_RB_07#BB_1#SW-01#FULL_RES     21   TOPO  356261.396     93750.000   1968750.0 357198.8962        1  XX  YY
  
        #WARN    MSTransformManager::checkCorrelatorPreaveraging       The data has already been preaveraged by the correlator but further smoothing or averaging has been requested.
        #LB_EB3,4,5,6,7,8,14,15,16, SB_EB4,5,6

for baseline_key,n_EB in number_of_EBs.items():
    for i in range(n_EB):
        vis = f'{prefix}_{baseline_key}_EB{i}_initcont.ms'
        listobs(
            vis       = vis,
            listfile  = f'{vis}.listobs.txt',
            overwrite = True,
        )

# for params in data_params.values():
#     flagchannels_string = get_flagchannels(ms_dict=params,output_prefix=prefix,
#                                            velocity_range=np.array([-15.,15.])+v_sys)
#     avg_cont(ms_dict=params,output_prefix=prefix,flagchannels=flagchannels_string,
#              contspws=params['cont_spws'],width_array=params['width_array'])
#
# Flagchannels input string for LB_EB0: '0:836~3004, 2:794~3043, 3:787~3055'
# Flagchannels input string for LB_EB1: '0:837~3005, 2:794~3043, 3:787~3055'
# Flagchannels input string for LB_EB2: '0:837~3005, 2:794~3043, 3:787~3055'
# Flagchannels input string for LB_EB3: '0:835~3003, 2:796~3045, 3:784~3051'
# Flagchannels input string for SB_EB0: '0:837~3005, 2:794~3043, 3:788~3056'
# Flagchannels input string for SB_EB1: '0:837~3005, 2:794~3043, 3:788~3056'

preselfcal_initcont_amp_folder = get_figures_folderpath('2_preselfcal_initcont_amp_figures')
make_figures_folder(preselfcal_initcont_amp_folder)

uv_ranges = {'LB':'405~430m','SB':'125~150m'}

for params in data_params.values():

    plot_filename   = prefix+'_'+params['name']+'_amp-v-freq_initcont_preselfcal.png'
    plotms(
        vis         = prefix+'_'+params['name']+'_initcont.ms',
        xaxis       = 'freq',
        yaxis       = 'amp',
        field       = params['field'],
        ydatacolumn = 'data',
        avgtime     = '1e8',
        avgscan     = True,
        avgbaseline = False,
        iteraxis    = 'spw',
        coloraxis   = 'corr',
        showgui     = False,
        exprange    = 'all',
        plotfile    = os.path.join(preselfcal_initcont_amp_folder,plot_filename),
        freqframe   = 'LSRK',
        # yselfscale  = True,
        showatm     = True,
        overwrite   = True,
    )

    baseline_key, _ = params['name'].split('_')
    plot_filename   = prefix+'_'+params['name']+'_amp-v-uvdist_initcont_preselfcal.png'
    plotms(
        vis         = prefix+'_'+params['name']+'_initcont.ms',
        xaxis       = 'UVdist',
        yaxis       = 'amp',
        # spw         = params['spwcont_forplot'],
        field       = params['field'],
        ydatacolumn = 'data',
        avgtime     = '1e8',
        avgscan     = True,
        avgchannel  = '3840',
        iteraxis    = 'spw',
        coloraxis   = 'corr',
        # default: coloraxis   = 'spw',
        showgui     = False,
        exprange    = 'all',
        plotrange   = plotranges[baseline_key],
        plotfile    = os.path.join(preselfcal_initcont_amp_folder,plot_filename),
        # yselfscale  = True,
        overwrite   = True,
    )

    baseline_key, _ = params['name'].split('_')
    plot_filename   = prefix+'_'+params['name']+'_amp-v-time_initcont_preselfcal.png'
    plotms(
        vis         = prefix+'_'+params['name']+'_initcont.ms',
        xaxis       = 'time',
        yaxis       = 'amp',
        # default: avgspw      = True,
            #For some obscure reason plotms doesn't like to load the cache when averaging over all the spw 
            #(both in vector and scalar mode and when reload is enforced) it just plots an empty canvas and breaks. Chose to iterate over spw instead.

            #SEVERE  plotms::::casa   Exception Reported: Error in plotms: <_InactiveRpcError of RPC that terminated with:
            #SEVERE  plotms::::casa+  status = StatusCode.INVALID_ARGUMENT
            #SEVERE  plotms::::casa+  details = "plot parameters failed"
            #SEVERE  plotms::::casa+  debug_error_string = "UNKNOWN:Error received from peer  {grpc_message:"plot parameters failed", grpc_status:3, created_time:"2024-10-30T06:19:55.017628422+01:00"}"
            #SEVERE  plotms::::casa+ >
        uvrange     = uv_ranges[baseline_key],
        avgchannel  = '3840',
        avgbaseline = True,
        iteraxis    = 'spw',
        coloraxis   = 'corr',
        # correlation = corr,
        showgui     = False,
        exprange    = 'all',
        plotfile    = os.path.join(preselfcal_initcont_amp_folder,plot_filename),
        # yselfscale  = True,
        overwrite   = True,
    )
#for SB, we see "waterfall features" (times where the amp suddently sharply decreases), we will try to fix this with self-cal

#Define simple masks and clean scales for imaging
mask_pa        = PA  #position angle of mask in degrees
mask_semimajor = 0.8 #semimajor axis of mask in arcsec
mask_semiminor = 0.8 #mask_semimajor*np.cos(incl/180.*np.pi) #semiminor axis of mask in arcsec
mask_ra        = '05h30m27.535947s' #taken from listobs
mask_dec       = '25d19m56.615500s' #taken from listobs

mask = f'ellipse[[{mask_ra},{mask_dec}], [{mask_semimajor:.3f}arcsec, {mask_semiminor:.3f}arcsec], {mask_pa:.1f}deg]'
noise_annulus = f"annulus[[{mask_ra}, {mask_dec}],['4.arcsec', '6.arcsec']]"

# Cellsize: ~beam/8
cellsize =  {'LB':'0.002arcsec','SB':'0.014arcsec'}

# Image size: ~primary beam 1.22*lam/A = 32'' with A=12m (19 arcsec)
imsize = {'LB':20000,'SB':3600}
scales = {'LB':[0,2,4,8,16,24],'SB':[0,2,4,8,12]} #32 seems too big of a scale for the LBs, same for 16 for the SBs

# Threshold: ~6 sigma for phase cal [1.35e-02,1.52e-02,1.96e-02,2.53e-02,3.27e-02],[3.57e-02,4.02e-02,4.46e-02,3.54e-02]
thresholds = {
    'LB':['0.0804mJy','0.0912mJy','0.1176mJy','0.1524mJy','0.1962mJy'],
    'SB':['0.2232mJy','0.2394mJy','0.2712mJy','0.2112mJy']
} #clean down to 6 sigma for phase cal

image_png_plot_sizes = [3,10] #sizes in arcsec of the zoomed and overview plots of the pngs

preselfcal_images_png_folder = get_figures_folderpath('3_preselfcal_images')
make_figures_folder(preselfcal_images_png_folder)

for baseline_key,params in zip(('LB','SB'),(data_params_LB,data_params_SB)):
    for EB_key,p in params.items():
        os.system(f'rm -rf '+prefix+'_'+p['name']+'_initcont_image*')
        _, idx_key    = p['name'].split('EB')

        imagename = prefix+'_'+p['name']+'_initcont_image'
        tclean_wrapper(
            vis            = prefix+'_'+p['name']+'_initcont.ms',
            imagename      = imagename,
            # deconvolver    = 'hogbom',
            deconvolver    = 'multiscale',
            scales         = scales[baseline_key],
            smallscalebias = 0.6,                  #Default from Cornwell et al. (2008) and in CASA 5.5 (biases to smaller scales)
            gain           = 0.3,                  #Default in DSHARP and exoALMA
            cycleniter     = 300,                  #Default in DSHARP and exoALMA
            niter          = 1000000,
            mask           = mask,
            threshold      = thresholds[baseline_key][int(idx_key)],
            cellsize       = cellsize[baseline_key],
            imsize         = imsize[baseline_key],
            parallel       = use_parallel,
            savemodel      = 'modelcolumn',
        )
        estimate_SNR(f'{imagename}.image',disk_mask=mask,noise_mask=noise_annulus)
        rms = imstat(imagename=f'{imagename}.image',region=noise_annulus)['rms'][0]
        p['rms'] = rms
        generate_image_png(
            image=f'{imagename}.image',plot_sizes=image_png_plot_sizes,
            color_scale_limits=[-3*rms,10*rms],
            save_folder=preselfcal_images_png_folder
        )

#MWC_758_LB_EB0_initcont_image.image
#Beam 0.030 arcsec x 0.019 arcsec (-12.89 deg)
#Flux inside disk mask: 5.12 mJy
#Peak intensity of source: 0.22 mJy/beam
#rms: 1.33e-02 mJy/beam
#Peak SNR: 16.50

#MWC_758_LB_EB1_initcont_image.image
#Beam 0.036 arcsec x 0.019 arcsec (-36.40 deg)
#Flux inside disk mask: 4.12 mJy
#Peak intensity of source: 0.19 mJy/beam
#rms: 1.50e-02 mJy/beam
#Peak SNR: 12.98

#MWC_758_LB_EB2_initcont_image.image
#Beam 0.030 arcsec x 0.019 arcsec (14.36 deg)
#Flux inside disk mask: 1.45 mJy
#Peak intensity of source: 0.22 mJy/beam
#rms: 1.93e-02 mJy/beam
#Peak SNR: 11.13

#MWC_758_LB_EB3_initcont_image.image
#Beam 0.039 arcsec x 0.022 arcsec (16.03 deg)
#Flux inside disk mask: 13.87 mJy
#Peak intensity of source: 0.38 mJy/beam
#rms: 2.52e-02 mJy/beam
#Peak SNR: 14.96

#MWC_758_LB_EB4_initcont_image.image
#Beam 0.035 arcsec x 0.023 arcsec (-7.27 deg)
#Flux inside disk mask: 13.87 mJy
#Peak intensity of source: 0.36 mJy/beam
#rms: 3.24e-02 mJy/beam
#Peak SNR: 11.01

#MWC_758_SB_EB0_initcont_image.image
#Beam 0.155 arcsec x 0.108 arcsec (-11.28 deg)
#Flux inside disk mask: 47.27 mJy
#Peak intensity of source: 3.14 mJy/beam
#rms: 3.72e-02 mJy/beam
#Peak SNR: 84.19

#MWC_758_SB_EB1_initcont_image.image
#Beam 0.234 arcsec x 0.153 arcsec (-41.10 deg)
#Flux inside disk mask: 58.29 mJy
#Peak intensity of source: 5.65 mJy/beam
#rms: 3.97e-02 mJy/beam
#Peak SNR: 142.22

#MWC_758_SB_EB2_initcont_image.image
#Beam 0.201 arcsec x 0.157 arcsec (-31.14 deg)
#Flux inside disk mask: 53.48 mJy
#Peak intensity of source: 5.23 mJy/beam
#rms: 4.43e-02 mJy/beam
#Peak SNR: 118.14

#MWC_758_SB_EB3_initcont_image.image
#Beam 0.201 arcsec x 0.157 arcsec (-28.54 deg)
#Flux inside disk mask: 57.32 mJy
#Peak intensity of source: 5.61 mJy/beam
#rms: 3.47e-02 mJy/beam
#Peak SNR: 161.66

#Selfcal individual EBs

#Self-calibration parameters
single_EB_contspws = '0~3'
single_EB_spw_mapping = [0,0,0,0]

individual_EB_selfcal_shift_folder = get_figures_folderpath('4_individual_EB_selfcal_and_shift_figures')
make_figures_folder(individual_EB_selfcal_shift_folder)

#One round of phase-only self-cal
for params in data_params.values():
    
    vis = prefix+'_'+params['name']+'_initcont.ms'
    single_EB_p1 = prefix+'_'+params['name']+'_initcont.p1'
    
    os.system(f'rm -rf {single_EB_p1}')
    gaincal(vis=vis,caltable=single_EB_p1,gaintype='T',spw=single_EB_contspws,combine='scan,spw',calmode='p',solint='inf',minsnr=4,minblperant=3)

    #Print calibration png file
    plotfilename = prefix+'_'+params['name']+'_initcont_gain_p1_phase_vs_time.png'
    # plotms(single_EB_p1,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,plotfile=os.path.join(individual_EB_selfcal_shift_folder,plotfilename))
    plotms(
        single_EB_p1,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,
        plotfile=os.path.join(individual_EB_selfcal_shift_folder,plotfilename),
        customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',
        #iteraxis='spw',exprange='all'
    )
    if params['name'] == 'LB_EB0':
        flagdata(vis=prefix+'_'+params['name']+'_initcont.p1',mode='manual',antenna='DA54,DV24')

        plotfilename = prefix+'_'+params['name']+'_initcont_gain_p1_phase_vs_time_flagged.png'
        plotms(
            single_EB_p1,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,
            plotfile=os.path.join(individual_EB_selfcal_shift_folder,plotfilename),
            customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',
            #iteraxis='spw',exprange='all'
        )

    #Apply the solutions
    os.system(f'rm -rf '+prefix+'_'+params['name']+'_initcont_selfcal.ms')
    applycal(vis=vis,spw=single_EB_contspws,spwmap=single_EB_spw_mapping,gaintable=[single_EB_p1],interp='linearPD',applymode='calonly',calwt=True)
    split(vis=vis,outputvis=prefix+'_'+params['name']+'_initcont_selfcal.ms',datacolumn='corrected')

#LB_EB0:  3 of 46 solutions flagged due to SNR < 4 in spw=0 at 2017/10/10/09:36:15.1
#LB_EB1:  8 of 44 solutions flagged due to SNR < 4 in spw=0 at 2017/10/10/11:02:38.7
#LB_EB2: 11 of 43 solutions flagged due to SNR < 4 in spw=0 at 2017/10/11/08:04:17.3
#LB_EB3: 11 of 46 solutions flagged due to SNR < 4 in spw=0 at 2017/10/15/07:58:16.4
#LB_EB4: 13 of 44 solutions flagged due to SNR < 4 in spw=0 at 2017/10/16/09:14:17.5

#Image the self-calibrated EBs with same cellsize and imsize for ratio
for baseline_key,params in zip(('LB','SB'),(data_params_LB,data_params_SB)):
    for EB_key,p in params.items():
        os.system(f'rm -rf '+prefix+'_'+p['name']+'_initcont_selfcal_image*')
        _, idx_key    = p['name'].split('EB')

        imagename = prefix+'_'+p['name']+'_initcont_selfcal_image'
        tclean_wrapper(
            vis            = prefix+'_'+p['name']+'_initcont_selfcal.ms',
            imagename      = imagename,
            # deconvolver    = 'hogbom',
            deconvolver    = 'multiscale',
            scales         = scales[baseline_key],
            smallscalebias = 0.6,                  #Default from Cornwell et al. (2008) and in CASA 5.5 (biases to smaller scales)
            gain           = 0.3,                  #Default in DSHARP and exoALMA
            cycleniter     = 300,                  #Default in DSHARP and exoALMA
            niter          = 1000000,
            mask           = mask,
            threshold      = thresholds[baseline_key][int(idx_key)],
            cellsize       = cellsize[baseline_key],
            imsize         = imsize[baseline_key],
            parallel       = use_parallel,
            savemodel      = 'modelcolumn',
        )
        estimate_SNR(f'{imagename}.image',disk_mask=mask,noise_mask=noise_annulus)
        rms = imstat(imagename=f'{imagename}.image',region=noise_annulus)['rms'][0]
        p['rms'] = rms
        generate_image_png(
            image=f'{imagename}.image',plot_sizes=image_png_plot_sizes,
            color_scale_limits=[-3*rms,10*rms],
            save_folder=individual_EB_selfcal_shift_folder
        )

#MWC_758_LB_EB0_initcont_selfcal_image.image
#Beam 0.030 arcsec x 0.019 arcsec (-12.89 deg)
#Flux inside disk mask: 4.82 mJy
#Peak intensity of source: 0.24 mJy/beam
#rms: 1.33e-02 mJy/beam
#Peak SNR: 18.23

#MWC_758_LB_EB1_initcont_selfcal_image.image
#Beam 0.036 arcsec x 0.019 arcsec (-36.40 deg)
#Flux inside disk mask: 4.18 mJy
#Peak intensity of source: 0.21 mJy/beam
#rms: 1.50e-02 mJy/beam
#Peak SNR: 13.80

#MWC_758_LB_EB2_initcont_selfcal_image.image
#Beam 0.030 arcsec x 0.019 arcsec (14.36 deg)
#Flux inside disk mask: 1.99 mJy
#Peak intensity of source: 0.25 mJy/beam
#rms: 1.93e-02 mJy/beam
#Peak SNR: 13.05

#MWC_758_LB_EB3_initcont_selfcal_image.image
#Beam 0.039 arcsec x 0.022 arcsec (16.03 deg)
#Flux inside disk mask: 14.08 mJy
#Peak intensity of source: 0.37 mJy/beam
#rms: 2.51e-02 mJy/beam
#Peak SNR: 14.59

#MWC_758_LB_EB4_initcont_selfcal_image.image
#Beam 0.035 arcsec x 0.023 arcsec (-7.27 deg)
#Flux inside disk mask: 13.57 mJy
#Peak intensity of source: 0.39 mJy/beam
#rms: 3.24e-02 mJy/beam
#Peak SNR: 12.15

#MWC_758_SB_EB0_initcont_selfcal_image.image
#Beam 0.155 arcsec x 0.108 arcsec (-11.28 deg)
#Flux inside disk mask: 48.08 mJy
#Peak intensity of source: 3.21 mJy/beam
#rms: 3.34e-02 mJy/beam
#Peak SNR: 96.11

#MWC_758_SB_EB1_initcont_selfcal_image.image
#Beam 0.234 arcsec x 0.153 arcsec (-41.10 deg)
#Flux inside disk mask: 58.66 mJy
#Peak intensity of source: 5.75 mJy/beam
#rms: 3.49e-02 mJy/beam
#Peak SNR: 164.88

#MWC_758_SB_EB2_initcont_selfcal_image.image
#Beam 0.201 arcsec x 0.157 arcsec (-31.14 deg)
#Flux inside disk mask: 54.80 mJy
#Peak intensity of source: 5.42 mJy/beam
#rms: 3.34e-02 mJy/beam
#Peak SNR: 162.22

#MWC_758_SB_EB3_initcont_selfcal_image.image
#Beam 0.201 arcsec x 0.157 arcsec (-28.55 deg)
#Flux inside disk mask: 57.91 mJy
#Peak intensity of source: 5.74 mJy/beam
#rms: 2.67e-02 mJy/beam
#Peak SNR: 215.26

#For the images with same pixel scales, compute the image ratios
for baseline_key,params in zip(('LB','SB'),(data_params_LB,data_params_SB)):
    for p in params.values():
        imagename = prefix+'_'+p['name']+'_initcont_selfcal_image'
        #ratio image, to be compared to the ratio image after alignment
        ref_image = f'{prefix}_{baseline_key}_EB0_initcont_selfcal_image.image'
        ref_rms = params[f'{baseline_key}0']['rms']
        ratio_image = imagename+'.ratio'
        os.system(f'rm -rf {ratio_image}')
        immath(imagename=[ref_image,imagename+'.image'],mode='evalexpr',outfile=ratio_image,expr=f'iif(IM0 > 3*{ref_rms}, IM1/IM0, 0)')
        generate_image_png(
            ratio_image,plot_sizes=[2*mask_semimajor,2*mask_semimajor],
            color_scale_limits=[0.5,1.5],image_units='ratio',
            save_folder=individual_EB_selfcal_shift_folder
        )

#Align data (go from *initcont_selfcal.ms to *initcont_shift.ms)
#Select the LB EB to act as the reference (usually the best SNR one)
reference_for_LB_alignment = f'{prefix}_LB_EB0_initcont_selfcal.ms'
assert 'LB' in reference_for_LB_alignment, 'you need to choose an LB EB for alignment of LB'

alignment_offsets = {}

#All the other EBs will be aligned to the reference EB...
#We also include the reference EB itself to make sure the coordinate changes are copied over
offset_LB_EBs = ['{}_{}_initcont_selfcal.ms'.format(prefix, params['name']) for params in data_params_LB.values()]

#Select the continuum spw with the large bandwidth
continuum_spw_id = 0
#Find the relative offsets and update the phase centers for all offset_EBs.
#NOTE: We might want to play around with the values of npix and cell_size.
npix      = 1024
cell_size = 0.008

for params in data_params_LB.values():
    offset_ms    = prefix+'_'+params['name']+'_initcont_selfcal.ms'
    plotfilename = 'uv_overlap_'+params['name']+'.png'
    if offset_ms == reference_for_LB_alignment:
        #for some reason the fitter fails when computing the offset of an EB to itself, so we skip the ref EB
        continue
    
    offset = alignment.find_offset(
        reference_ms=reference_for_LB_alignment,
        offset_ms=offset_ms,
        npix=npix,cell_size=cell_size,spwid=continuum_spw_id,plot_uv_grid=False,
        uv_grid_plot_filename=os.path.join(individual_EB_selfcal_shift_folder,plotfilename)
    )
    print(f'#Offset for {offset_ms}: ',offset)
    #Offset for MWC_758_LB_EB1_initcont_selfcal.ms:  [ 0.00272863  0.01096198]
    #Offset for MWC_758_LB_EB2_initcont_selfcal.ms:  [-0.00732089  0.00034065]
    #Offset for MWC_758_LB_EB3_initcont_selfcal.ms:  [-0.0128694  -0.00634683]
    #Offset for MWC_758_LB_EB4_initcont_selfcal.ms:  [-0.0021153  -0.01235212]
    
    for _npix in [256,512,1024,2048]:
        plotfilename = 'uv_overlap_'+params['name']+f'_{_npix}pixels.png'
        offset = alignment.find_offset(
            reference_ms=reference_for_LB_alignment,
            offset_ms=offset_ms,
            npix=_npix,cell_size=cell_size,spwid=continuum_spw_id,plot_uv_grid=True,
            uv_grid_plot_filename=os.path.join(individual_EB_selfcal_shift_folder,plotfilename)
        )
        print(f'#Offset for {offset_ms}: ',offset)
        #A bit concerning:
        #Offset for MWC_758_LB_EB1_initcont_selfcal.ms:  [-0.00524202 -0.00139445]
        #Offset for MWC_758_LB_EB1_initcont_selfcal.ms:  [-0.01065405 -0.00067   ]
        #Offset for MWC_758_LB_EB1_initcont_selfcal.ms:  [ 0.00272863  0.01096198]
        #Offset for MWC_758_LB_EB1_initcont_selfcal.ms:  [-0.01164673  0.02087499]

        #Offset for MWC_758_LB_EB2_initcont_selfcal.ms:  [-0.00662535  0.00404133]
        #Offset for MWC_758_LB_EB2_initcont_selfcal.ms:  [-0.02170253  0.00325939]
        #Offset for MWC_758_LB_EB2_initcont_selfcal.ms:  [-0.00732089  0.00034065]
        #Offset for MWC_758_LB_EB2_initcont_selfcal.ms:  [-0.01210021 -0.00571096]

        #Offset for MWC_758_LB_EB3_initcont_selfcal.ms:  [-0.00269578 -0.00693527]
        #Offset for MWC_758_LB_EB3_initcont_selfcal.ms:  [-0.00551198  0.00030559]
        #Offset for MWC_758_LB_EB3_initcont_selfcal.ms:  [-0.0128694  -0.00634683]
        #Offset for MWC_758_LB_EB3_initcont_selfcal.ms:  [-0.01751555 -0.03452715]

        #Offset for MWC_758_LB_EB4_initcont_selfcal.ms:  [-0.00544589 -0.00729134]
        #Offset for MWC_758_LB_EB4_initcont_selfcal.ms:  [ 0.00225751 -0.01001635]
        #Offset for MWC_758_LB_EB4_initcont_selfcal.ms:  [-0.0021153  -0.01235212]
        #Offset for MWC_758_LB_EB4_initcont_selfcal.ms:  [-0.00198791 -0.00768387]

alignment.align_measurement_sets(
    reference_ms=reference_for_LB_alignment,align_ms=offset_LB_EBs,
    npix=npix,cell_size=cell_size,spwid=continuum_spw_id
)
#New coordinates for MWC_758_LB_EB0_initcont_selfcal.ms no shift, reference MS.
#New coordinates for MWC_758_LB_EB1_initcont_selfcal.ms requires a shift of [ 0.0027286, 0.01096200]
#New coordinates for MWC_758_LB_EB2_initcont_selfcal.ms requires a shift of [-0.0073209, 0.00034065]
#New coordinates for MWC_758_LB_EB3_initcont_selfcal.ms requires a shift of [-0.0128690,-0.00634680]
#New coordinates for MWC_758_LB_EB4_initcont_selfcal.ms requires a shift of [-0.0021153,-0.01235200]

#Insert offsets from the alignment output
alignment_offsets['LB_EB0'] = [ 0.0,       0.0       ]
alignment_offsets['LB_EB1'] = [ 0.0027286, 0.01096200]
alignment_offsets['LB_EB2'] = [-0.0073209, 0.00034065]
alignment_offsets['LB_EB3'] = [-0.0128690,-0.00634680]
alignment_offsets['LB_EB4'] = [-0.0021153,-0.01235200]

shifted_LB_EBs = [EB.replace('.ms','_shift.ms') for EB in offset_LB_EBs]

#To check if alignment worked, calculate shift again and verify that shifts are small (i.e. a fraction of the cell size):
for shifted_ms in shifted_LB_EBs:
    if shifted_ms == reference_for_LB_alignment.replace('.ms','_shift.ms'):
        #For some reason the fitter fails when computing the offset of an EB to itself, so we skip the ref EB
        continue
    offset = alignment.find_offset(
        reference_ms=reference_for_LB_alignment,
        offset_ms=shifted_ms,npix=npix,plot_uv_grid=False,
        cell_size=cell_size,spwid=continuum_spw_id
    )
    print(f'#Offset for {shifted_ms}: ',offset)
    #Offset for MWC_758_LB_EB1_initcont_selfcal_shift.ms:  [ 6.28430839e-05 -4.49028662e-05]
    #Offset for MWC_758_LB_EB2_initcont_selfcal_shift.ms:  [-2.84107951e-05  2.59390592e-05]
    #Offset for MWC_758_LB_EB3_initcont_selfcal_shift.ms:  [ 4.36137441e-05 -2.18412985e-05]
    #Offset for MWC_758_LB_EB4_initcont_selfcal_shift.ms:  [ 9.17612073e-06 -4.94944473e-05]

#Merge shifted LB EBs for aligning SB EBs
LB_concat_shifted = f'{prefix}_LB_concat_shifted.ms'
os.system(f'rm -rf {LB_concat_shifted}')
concat(
    vis=shifted_LB_EBs,concatvis=LB_concat_shifted,
    dirtol='0.1arcsec',freqtol='2.0GHz',
    copypointing=False
)
listobs(vis=LB_concat_shifted,listfile=f'{LB_concat_shifted}.listobs.txt',overwrite=True)

#Align SB EBs to concat shifted LB EBs
reference_for_SB_alignment = LB_concat_shifted

offset_SB_EBs = ['{}_{}_initcont_selfcal.ms'.format(prefix, params['name']) for params in data_params_SB.values()]

for params in data_params_SB.values():
    offset_ms    = prefix+'_'+params['name']+'_initcont_selfcal.ms'

    plotfilename = 'uv_overlap_'+params['name']+'.png'
    if offset_ms == reference_for_SB_alignment:
        #For some reason the fitter fails when computing the offset of an EB to itself, so we skip the ref EB
        continue
    offset = alignment.find_offset(
        reference_ms=reference_for_SB_alignment,
        offset_ms=offset_ms,
        npix=npix,cell_size=cell_size,spwid=continuum_spw_id,plot_uv_grid=False,
        uv_grid_plot_filename=os.path.join(individual_EB_selfcal_shift_folder,plotfilename)
    )
    print(f'#Offset for {offset_ms}: ',offset)
    #Offset for MWC_758_SB_EB0_initcont_selfcal.ms:  [ 0.00028434 -0.01428973]
    #Offset for MWC_758_SB_EB1_initcont_selfcal.ms:  [ 0.01114162 -0.0224137 ]
    #Offset for MWC_758_SB_EB2_initcont_selfcal.ms:  [-0.02411884 -0.00031381]
    #Offset for MWC_758_SB_EB3_initcont_selfcal.ms:  [-0.00655696 -0.01639861]

    for _npix in [256,512,1024,2048]:
        plotfilename = 'uv_overlap_'+params['name']+f'_{_npix}pixels.png'
        offset = alignment.find_offset(
            reference_ms=reference_for_SB_alignment,
            offset_ms=offset_ms,
            npix=_npix,cell_size=cell_size,spwid=continuum_spw_id,plot_uv_grid=True,
            uv_grid_plot_filename=os.path.join(individual_EB_selfcal_shift_folder,plotfilename)
        )
        print(f'#Offset for {offset_ms}: ',offset)
        
        #Offset for MWC_758_SB_EB0_initcont_selfcal.ms:  [ 0.00068351 -0.01348155]
        #Offset for MWC_758_SB_EB0_initcont_selfcal.ms:  [-0.00171016 -0.01055136]
        #Offset for MWC_758_SB_EB0_initcont_selfcal.ms:  [ 0.00028434 -0.01428973]
        #Offset for MWC_758_SB_EB0_initcont_selfcal.ms:  [-0.00287896 -0.01167131]

        #Offset for MWC_758_SB_EB1_initcont_selfcal.ms:  [ 0.01060743 -0.0282687 ]
        #Offset for MWC_758_SB_EB1_initcont_selfcal.ms:  [ 0.00864346 -0.02217845]
        #Offset for MWC_758_SB_EB1_initcont_selfcal.ms:  [ 0.01114162 -0.0224137 ]
        #Offset for MWC_758_SB_EB1_initcont_selfcal.ms:  [ 0.00349056 -0.02467989]

        #Offset for MWC_758_SB_EB2_initcont_selfcal.ms:  [-0.01976164 -0.00235033]
        #Offset for MWC_758_SB_EB2_initcont_selfcal.ms:  [-0.02467678  0.00180465]
        #Offset for MWC_758_SB_EB2_initcont_selfcal.ms:  [-0.02411884 -0.00031381]
        #Offset for MWC_758_SB_EB2_initcont_selfcal.ms:  [-0.02269373  0.00197626]

        #Offset for MWC_758_SB_EB3_initcont_selfcal.ms:  [-0.00352746 -0.01399474]
        #Offset for MWC_758_SB_EB3_initcont_selfcal.ms:  [-0.00672193 -0.01163956]
        #Offset for MWC_758_SB_EB3_initcont_selfcal.ms:  [-0.00655696 -0.01639861]
        #Offset for MWC_758_SB_EB3_initcont_selfcal.ms:  [-0.006268   -0.01442109]
   
alignment.align_measurement_sets(
    reference_ms=reference_for_SB_alignment,align_ms=offset_SB_EBs,
    npix=npix,cell_size=cell_size,spwid=continuum_spw_id
)
#New coordinates for MWC_758_SB_EB0_initcont_selfcal.ms requires a shift of [ 0.00028434,-0.01429   ]
#New coordinates for MWC_758_SB_EB1_initcont_selfcal.ms requires a shift of [ 0.011142,  -0.022414  ]
#New coordinates for MWC_758_SB_EB2_initcont_selfcal.ms requires a shift of [-0.024119,  -0.00031381]
#New coordinates for MWC_758_SB_EB3_initcont_selfcal.ms requires a shift of [-0.006557,  -0.016399  ]

alignment_offsets['SB_EB0'] = [ 0.00028434,-0.01429   ]
alignment_offsets['SB_EB1'] = [ 0.011142,  -0.022414  ]
alignment_offsets['SB_EB2'] = [-0.024119,  -0.00031381]
alignment_offsets['SB_EB3'] = [-0.006557,  -0.016399  ]

shifted_SB_EBs = [EB.replace('.ms','_shift.ms') for EB in offset_SB_EBs]

#Check by calculating offset again
for params in data_params_SB.values():
    shifted_ms = prefix+'_'+params['name']+'_initcont_selfcal_shift.ms'
    offset = alignment.find_offset(
        reference_ms=reference_for_SB_alignment,
        offset_ms=shifted_ms,npix=npix,plot_uv_grid=False,
        cell_size=cell_size,spwid=continuum_spw_id
    )
    print(f'#Offset for {shifted_ms}: ',offset)
    #Offset for MWC_758_SB_EB0_initcont_selfcal_shift.ms:  [-3.79661973e-05 -6.92752483e-05]
    #Offset for MWC_758_SB_EB1_initcont_selfcal_shift.ms:  [ 1.43485533e-05 -7.31554965e-06]
    #Offset for MWC_758_SB_EB2_initcont_selfcal_shift.ms:  [-3.45129141e-05 -2.58725160e-05]
    #Offset for MWC_758_SB_EB3_initcont_selfcal_shift.ms:  [-7.40714058e-05 -8.27564886e-05]

#Remove the '_selfcal' part of the names to match the naming convention below.
for shifted_EB in shifted_LB_EBs+shifted_SB_EBs:
    os.system('mv {} {}'.format(shifted_EB, shifted_EB.replace('_selfcal', '')))

#Check that the images are indeed aligned after the shift
for baseline_key,params in zip(('LB','SB'),(data_params_LB,data_params_SB)):
    for EB_key,p in params.items():
        os.system(f'rm -rf '+prefix+'_'+p['name']+'_initcont_shift_image*')
        _, idx_key    = p['name'].split('EB')

        imagename = prefix+'_'+p['name']+'_initcont_shift_image'
        tclean_wrapper(
            vis            = prefix+'_'+p['name']+'_initcont_shift.ms',
            imagename      = imagename,
            # deconvolver    = 'hogbom',
            deconvolver    = 'multiscale',
            scales         = scales[baseline_key],
            smallscalebias = 0.6,                  #Default from Cornwell et al. (2008) and in CASA 5.5 (biases to smaller scales)
            gain           = 0.3,                  #Default in DSHARP and exoALMA
            cycleniter     = 300,                  #Default in DSHARP and exoALMA
            niter          = 1000000,
            mask           = mask,
            threshold      = thresholds[baseline_key][int(idx_key)],
            cellsize       = cellsize[baseline_key],
            imsize         = imsize[baseline_key],
            parallel       = use_parallel,
            savemodel      = 'modelcolumn',
        )
        estimate_SNR(f'{imagename}.image',disk_mask=mask,noise_mask=noise_annulus)
        rms = imstat(imagename=f'{imagename}.image',region=noise_annulus)['rms'][0]
        p['rms'] = rms
        generate_image_png(
            image=f'{imagename}.image',plot_sizes=image_png_plot_sizes,
            color_scale_limits=[-3*rms,10*rms],
            save_folder=individual_EB_selfcal_shift_folder
        )

#MWC_758_LB_EB0_initcont_shift_1024pixs_image.image
#Beam 0.030 arcsec x 0.019 arcsec (-12.89 deg)
#Flux inside disk mask: 5.05 mJy
#Peak intensity of source: 0.24 mJy/beam
#rms: 1.33e-02 mJy/beam
#Peak SNR: 18.39

#MWC_758_LB_EB1_initcont_shift_1024pixs_image.image
#Beam 0.036 arcsec x 0.019 arcsec (-36.40 deg)
#Flux inside disk mask: 3.88 mJy
#Peak intensity of source: 0.21 mJy/beam
#rms: 1.50e-02 mJy/beam
#Peak SNR: 13.84

#MWC_758_LB_EB2_initcont_shift_1024pixs_image.image
#Beam 0.030 arcsec x 0.019 arcsec (14.36 deg)
#Flux inside disk mask: 1.81 mJy
#Peak intensity of source: 0.25 mJy/beam
#rms: 1.93e-02 mJy/beam
#Peak SNR: 13.11

#MWC_758_LB_EB3_initcont_shift_1024pixs_image.image
#Beam 0.039 arcsec x 0.022 arcsec (16.03 deg)
#Flux inside disk mask: 13.70 mJy
#Peak intensity of source: 0.36 mJy/beam
#rms: 2.51e-02 mJy/beam
#Peak SNR: 14.51

#MWC_758_LB_EB4_initcont_shift_1024pixs_image.image
#Beam 0.035 arcsec x 0.023 arcsec (-7.27 deg)
#Flux inside disk mask: 13.48 mJy
#Peak intensity of source: 0.38 mJy/beam
#rms: 3.24e-02 mJy/beam
#Peak SNR: 11.57

#MWC_758_SB_EB0_initcont_shift_1024pixs_image.image
#Beam 0.155 arcsec x 0.108 arcsec (-11.28 deg)
#Flux inside disk mask: 48.03 mJy
#Peak intensity of source: 3.21 mJy/beam
#rms: 3.35e-02 mJy/beam
#Peak SNR: 95.95

#MWC_758_SB_EB1_initcont_shift_1024pixs_image.image
#Beam 0.234 arcsec x 0.153 arcsec (-41.10 deg)
#Flux inside disk mask: 58.62 mJy
#Peak intensity of source: 5.76 mJy/beam
#rms: 3.49e-02 mJy/beam
#Peak SNR: 164.84

#MWC_758_SB_EB2_initcont_shift_1024pixs_image.image
#Beam 0.201 arcsec x 0.157 arcsec (-31.14 deg)
#Flux inside disk mask: 54.97 mJy
#Peak intensity of source: 5.41 mJy/beam
#rms: 3.35e-02 mJy/beam
#Peak SNR: 161.32

#MWC_758_SB_EB3_initcont_shift_1024pixs_image.image
#Beam 0.201 arcsec x 0.157 arcsec (-28.55 deg)
#Flux inside disk mask: 57.92 mJy
#Peak intensity of source: 5.72 mJy/beam
#rms: 2.67e-02 mJy/beam
#Peak SNR: 214.05

for baseline_key,params in zip(('LB','SB'),(data_params_LB,data_params_SB)):
    for p in params.values():
        imagename = prefix+'_'+p['name']+'_initcont_shift_image'
        #ratio image, to be compared to the ratio image after alignment
        ref_image = f'{prefix}_{baseline_key}_EB0_initcont_shift_image.image'
        ref_rms = params[f'{baseline_key}0']['rms']
        ratio_image = imagename+'.ratio'
        os.system(f'rm -rf {ratio_image}')
        immath(imagename=[ref_image,imagename+'.image'],mode='evalexpr',outfile=ratio_image,expr=f'iif(IM0 > 3*{ref_rms}, IM1/IM0, 0)')
        generate_image_png(
            ratio_image,plot_sizes=[2*mask_semimajor,2*mask_semimajor],
            color_scale_limits=[0.5,1.5],image_units='ratio',
            save_folder=individual_EB_selfcal_shift_folder
        )
#By checking the different images in casaviewer, the alignment of some EBs improved, for others LB_EB2 and SB_EB2 it didn't

for i,params in enumerate(data_params_LB.values()):
    print('#'+params['name'])
    mask = f"ellipse[[{mask_ra},{mask_dec}], [{mask_semimajor:.3f}arcsec, {mask_semiminor:.3f}arcsec], {mask_pa:.1f}deg]"
    
    fit_gaussian(prefix+'_'+params['name']+'_initcont_selfcal_image.image',region=mask)

    #LB_EB0, 05h30m27.519563s +25d19m57.07381s
    #Pixel coordinates of peak: x = 10111.062 y = 10229.153
    #LB_EB1, 05h30m27.518385s +25d19m57.06236s
    #Pixel coordinates of peak: x = 10119.050 y = 10223.432
    #LB_EB2, 05h30m27.518407s +25d19m57.07056s
    #Pixel coordinates of peak: x = 10118.897 y = 10227.562
    #LB_EB3, 05h30m27.518816s +25d19m57.06485s
    #Pixel coordinates of peak: x = 10116.091 y = 10224.847
    #LB_EB4, 05h30m27.519260s +25d19m57.05768s
    #Pixel coordinates of peak: x = 10113.079 y = 10221.300
    
    fit_gaussian(prefix+'_'+params['name']+'_initcont_shift_image.image',region=mask)

    #LB_EB0, 05h30m27.520036s +25d19m57.06502s
    #Pixel coordinates of peak: x = 10120.681 y = 10234.042
    #LB_EB1, 05h30m27.518671s +25d19m57.04264s
    #Pixel coordinates of peak: x = 10129.940 y = 10222.849
    #LB_EB2, 05h30m27.519467s +25d19m57.06156s
    #Pixel coordinates of peak: x = 10124.542 y = 10232.309
    #LB_EB3, 05h30m27.520220s +25d19m57.06252s
    #Pixel coordinates of peak: x = 10119.438 y = 10232.790
    #LB_EB4, 05h30m27.519858s +25d19m57.06106s
    #Pixel coordinates of peak: x = 10121.890 y = 10232.059
 
    #Differences are huge (but also Gaussians on such a low SNR asymmetric disc... don't know how reliable...)!

for i,params in enumerate(data_params_SB.values()):
    print('#'+params['name'])
    mask = f"ellipse[[{mask_ra},{mask_dec}], [{mask_semimajor:.3f}arcsec, {mask_semiminor:.3f}arcsec], {mask_pa:.1f}deg]"

    fit_gaussian(prefix+'_'+params['name']+'_initcont_selfcal_image.image',region=mask)

    #SB_EB0, 05h30m27.523115s +25d19m56.85150s
    #Pixel coordinates of peak: x = 1812.284 y = 1817.167
    #SB_EB1, 05h30m27.524016s +25d19m56.79284s
    #Pixel coordinates of peak: x = 1811.384 y = 1813.020
    #SB_EB2, 05h30m27.518972s +25d19m56.94008s
    #Pixel coordinates of peak: x = 1816.234 y = 1823.590
    #SB_EB3, 05h30m27.522873s +25d19m56.80906s
    #Pixel coordinates of peak: x = 1812.454 y = 1814.237

    fit_gaussian(prefix+'_'+params['name']+'_initcont_shift_image.image',region=mask)

    #SB_EB0, 05h30m27.523645s +25d19m56.87042s
    #Pixel coordinates of peak: x = 1813.745 y = 1819.534
    #SB_EB1, 05h30m27.523252s +25d19m56.83625s
    #Pixel coordinates of peak: x = 1814.126 y = 1817.094
    #SB_EB2, 05h30m27.522646s +25d19m56.90281s
    #Pixel coordinates of peak: x = 1814.713 y = 1821.848
    #SB_EB3, 05h30m27.524017s +25d19m56.83294s
    #Pixel coordinates of peak: x = 1813.385 y = 1816.857

    #Differences are sub-pixels in the x-direction but significant in the y-direction!

#Chose not to apply the shift because it artificially increases the offset between EBs (regardless of npix, and number of spws included)
#likely due to the very low SNR of the LB EBs

#Now that everything is aligned, we inspect the flux calibration
for params in data_params.values():
    msfile = prefix+'_'+params['name']+'_initcont_selfcal.ms' #msfile = prefix+'_'+params['name']+'_initcont_shift.ms'
    #Export MS contents into numpy save files
    export_MS(msfile)

#Plot deprojected visibility profiles for all data together
list_npz_files = []
for baseline_key,n_EB in number_of_EBs.items():
    list_npz_files += [f'{prefix}_{baseline_key}_EB{i}_initcont_selfcal.vis.npz' for i in range(n_EB)] #list_npz_files += [f'{prefix}_{baseline_key}_EB{i}_initcont_shift.vis.npz' for i in range(n_EB)]

deprojected_vis_profiles_folder = get_figures_folderpath('5_deprojected_vis_profiles')
make_figures_folder(deprojected_vis_profiles_folder)

plot_deprojected(
    filelist=list_npz_files,
    fluxscale=[1.]*(number_of_EBs['LB']+number_of_EBs['SB']),
    PA=PA,incl=incl,show_err=True,
    plot_label=os.path.join(deprojected_vis_profiles_folder,f'{prefix}_flux_scale_EB_preselfcal.png')
)

flux_comparison_folder = get_figures_folderpath('6_flux_comparisons')
make_figures_folder(flux_comparison_folder)

##We choose an EB to compare the flux scaling; if possible the best SB EB.
##Since PDS 66 shows decoherence in the SB EBs (waterfalls), I chose as reference 
##the LB EB with the highest SNR (as for the alignment), that is, LB EB1.

#There are traces of waterfall features, but they're rather mild. Go for SB_EB3 that shows
#less signs of decoherence, a higher SNR and lower phase RMS (16.723 deg) than the others.
#The phase RMS of LB_EB0 and LB_EB1 is slightly lower (13.935, 15.003) but their SNR
#is too small to be good enough flux scale references

#Decoherence in SB_EB1 is clear from plot of the flux ratios.
#Flux offset of SB_EB0 is clear from the deprojected visibilities.

flux_ref_EB = 'SB_EB3' 

for params in data_params.values():
    plot_label = os.path.join(flux_comparison_folder,'flux_comparison_'+params['name']+f'_to_{flux_ref_EB}.png')
    estimate_flux_scale(
        reference=f'{prefix}_{flux_ref_EB}_initcont_selfcal.vis.npz',#reference=f'{prefix}_{flux_ref_EB}_initcont_shift.vis.npz',
        comparison=prefix+'_'+params['name']+'_initcont_selfcal.vis.npz',#comparison=prefix+'_'+params['name']+'_initcont_shift.vis.npz',
        incl=incl,PA=PA,plot_label=plot_label
    )

#The ratio of the fluxes of MWC_758_LB_EB0_initcont_selfcal.vis.npz to
#MWC_758_SB_EB3_initcont_selfcal.vis.npz is 0.95742
#The scaling factor for gencal is 0.978 for your comparison measurement
#The error on the weighted mean ratio is 4.902e-03, although it's likely that
#the weights in the measurement sets are too off by some constant factor

#The ratio of the fluxes of MWC_758_LB_EB1_initcont_selfcal.vis.npz to
#MWC_758_SB_EB3_initcont_selfcal.vis.npz is 0.93723
#The scaling factor for gencal is 0.968 for your comparison measurement
#The error on the weighted mean ratio is 4.833e-03, although it's likely that
#the weights in the measurement sets are too off by some constant factor

#The ratio of the fluxes of MWC_758_LB_EB2_initcont_selfcal.vis.npz to
#MWC_758_SB_EB3_initcont_selfcal.vis.npz is 0.90976
#The scaling factor for gencal is 0.954 for your comparison measurement
#The error on the weighted mean ratio is 6.064e-03, although it's likely that
#the weights in the measurement sets are too off by some constant factor

#The ratio of the fluxes of MWC_758_LB_EB3_initcont_selfcal.vis.npz to
#MWC_758_SB_EB3_initcont_selfcal.vis.npz is 1.01642
#The scaling factor for gencal is 1.008 for your comparison measurement
#The error on the weighted mean ratio is 4.064e-03, although it's likely that
#the weights in the measurement sets are too off by some constant factor

#The ratio of the fluxes of MWC_758_LB_EB4_initcont_selfcal.vis.npz to
#MWC_758_SB_EB3_initcont_selfcal.vis.npz is 1.11297
#The scaling factor for gencal is 1.055 for your comparison measurement
#The error on the weighted mean ratio is 5.086e-03, although it's likely that
#the weights in the measurement sets are too off by some constant factor

#The ratio of the fluxes of MWC_758_SB_EB0_initcont_selfcal.vis.npz to
#MWC_758_SB_EB3_initcont_selfcal.vis.npz is 0.88668
#The scaling factor for gencal is 0.942 for your comparison measurement
#The error on the weighted mean ratio is 1.119e-03, although it's likely that
#the weights in the measurement sets are too off by some constant factor

#The ratio of the fluxes of MWC_758_SB_EB1_initcont_selfcal.vis.npz to
#MWC_758_SB_EB3_initcont_selfcal.vis.npz is 0.98104
#The scaling factor for gencal is 0.990 for your comparison measurement
#The error on the weighted mean ratio is 7.696e-04, although it's likely that
#the weights in the measurement sets are too off by some constant factor

#The ratio of the fluxes of MWC_758_SB_EB2_initcont_selfcal.vis.npz to
#MWC_758_SB_EB3_initcont_selfcal.vis.npz is 0.93757
#The scaling factor for gencal is 0.968 for your comparison measurement
#The error on the weighted mean ratio is 8.283e-04, although it's likely that
#the weights in the measurement sets are too off by some constant factor

#The ratio of the fluxes of MWC_758_SB_EB3_initcont_selfcal.vis.npz to
#MWC_758_SB_EB3_initcont_selfcal.vis.npz is 1.00000
#The scaling factor for gencal is 1.000 for your comparison measurement
#The error on the weighted mean ratio is 6.559e-04, although it's likely that
#the weights in the measurement sets are too off by some constant factor

#The flux offsets are all >4%, for LB_EB1 and LB_EB2 some strange trends in amp_vs_time could be fixed by phase-selfcal
#Similarly for SB_EB1 and SB_EB3. It's tough to assess how robust these offsets are due to the uv-depence of the ratio
#
#Selfcal SBs first and the SBs + LBs in phase-only mode and reassess the flux offsets. Only then decide if scaling is needed

#We correct flux differences >4%
#
#The general procedure is as follows:
#
#If you need to re-scale fluxes, use command:
#rescale_flux(vis=prefix+'_LB_EB0_initcont_shift.ms', gencalparameter=[1.044])
#
#1) check flux scaling, and if flux differences are >4%, re-scale the fluxes UNLESS
#    THERE IS CLEAR DE-COHERENCE (e.g. seen as a systematic decrease of scaling with UV distance)
#2) do not scale those EBs with phase decoherence. Proceed with self-cal following the steps below 
#   to the point that those EBs are self calibrated
#3) check flux scaling again
#4) if no flux scaling has been applied in 1), and you now still see a flux offset, then
#    apply the flux scaling determined in step 3) to the non-self caled EBs
#    and repeat the self-cal. The additional code two do the second iteration of steps 1-3 is
#    already in this template. It is essentially a copy/paste of the first iteration,
#    but it saves all output with different filenames. If you use that additional code,
#    make sure you use the right files when proceeding with the script after self-cal.
#
#If all SB EBs suffer from phase decoherence and cannot be used as the reference
#for flux scaling:
#1) concat ACA data and self-cal them in phase
#2) concat them to SB data and self-cal them in phase
#3) check flux offsets
#4) if there are flux offsets, go back to 1), but before re-starting the process apply
#    the corrections to the non-self caled ACA and SB data, and re-do steps 1) - 3).
#    For this you should add additional code (essentially copy/paste the code from the first
#    iteration) that repeats steps 1-3, but saves all output with different filenames. Remember
#    to use the right filenames when you continue the script after step 6
#5) Check flux offsets of LB EBs to a correct SB EB
#6) concat LB data to ACA+SB and continue with script
#
# Notes:
# Check intervals with exoALMA paper: 360,120,60,30,18s
# Check which flux scale to apply (before or after phase-selfcal?): after selfcal as in the PDS 66 calibration script
# In Myriam's script on J1604 they also ran amplitude selfcal to see how it improved the flux scaling.
# I will go for 360,120,60,20s,    SNR3/ANT4, scan/spw for the SB and
#               360,120,60,30,18s, SNR2/ANT4, scan/spw for the LBs
#
# For MWC 758, I see strange amp_vs_time in the SB EBs and descending trends in the flux ratios, indicating decoherence. 
# Therefore, selfcal (phase-only) should fix it. I'll proceed by self-calibrating the SBs without rescaling
# Then I'll check the flux offsets at the end of the process...

#Begin of SB self-cal - iteration 1
#For phase self-cal, clean down to ~6sigma
SB_selfcal_folder = get_figures_folderpath('7_selfcal_SB_figures')
make_figures_folder(SB_selfcal_folder)

SB_cont_p0 = prefix+'_SB_contp0'
os.system('rm -rf %s.ms*' %SB_cont_p0)
concat(
    vis=[f'{prefix}_SB_EB{i}_initcont_selfcal.ms' for i in range(number_of_EBs['SB'])],#vis=[f'{prefix}_SB_EB{i}_initcont_shift.ms' for i in range(number_of_EBs['SB'])],
    concatvis=SB_cont_p0+'.ms',dirtol='0.1arcsec',copypointing=False
)
listobs(vis=SB_cont_p0+'.ms',listfile=SB_cont_p0+'.ms.listobs.txt',overwrite=True)

#Define new SB mask using new center read from listobs
mask_pa        = PA
mask_semimajor = 0.8
mask_semiminor = 0.8
mask_ra  = '05h30m27.535799s' 
mask_dec = '25d19m56.611160s'

SB_mask= f"ellipse[[{mask_ra},{mask_dec}], [{mask_semimajor:.3f}arcsec, {mask_semiminor:.3f}arcsec], {mask_pa:.1f}deg]"

noise_annulus_SB = f"annulus[[{mask_ra}, {mask_dec}],['4.arcsec', '6.arcsec']]"

SB_tclean_wrapper_kwargs = {
    'deconvolver':'multiscale','scales':[0,2,4,8,12],
    'smallscalebias':0.6,'gain':0.3,'cycleniter':300,
    'niter':1000000,'mask':SB_mask,
    'cellsize':'0.014arcsec','imsize':3600,
    'parallel':use_parallel,'savemodel':'modelcolumn',
    'robust':0.5,'interactive':False,
    'gridder':'standard',
}

tclean_wrapper(
    vis       = SB_cont_p0+'.ms',
    imagename = SB_cont_p0,
    threshold = '0.0960mJy', 
    **SB_tclean_wrapper_kwargs
)
estimate_SNR(SB_cont_p0+'.image',disk_mask=SB_mask,noise_mask=noise_annulus_SB)
rms_SB = imstat(imagename=SB_cont_p0+'.image',region=noise_annulus_SB)['rms'][0]
generate_image_png(
    SB_cont_p0+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_SB,10*rms_SB],
    save_folder=SB_selfcal_folder
)
#MWC_758_SB_contp0.image
#Beam 0.193 arcsec x 0.146 arcsec (-28.31 deg)
#Flux inside disk mask: 56.86 mJy
#Peak intensity of source: 5.03 mJy/beam
#rms: 1.58e-02 mJy/beam
#Peak SNR: 317.43

#Look for references antennas from weblog, and pick the first that are listed, overlapping with all EBs
#SB_EB0: DV05, DV08, DA64, DA60
#SB_EB1: DA49, DV23, DA64, DV06 
#SB_EB2: DV08, DA60, DV15, DV01
#SB_EB3: DA47, DV08, DV15, DA60

#Get station numbers
for ref_ant in ('DV08','DA64','DV15','DA60'):
    get_station_numbers(SB_cont_p0+'.ms',ref_ant)
    #Observation ID 0: DV08@A042
    #Observation ID 1: DV08@A042
    #Observation ID 2: DV08@A042
    #Observation ID 3: DV08@A042
    #Observation ID 0: DA64@A015
    #Observation ID 1: DA64@A015
    #Observation ID 2: DA64@A015
    #Observation ID 3: DA64@A015
    #Observation ID 0: DV15@A130
    #Observation ID 1: DV15@A047
    #Observation ID 2: DV15@A047
    #Observation ID 3: DV15@A047
    #Observation ID 0: DA60@A043
    #Observation ID 1: DA60@A043
    #Observation ID 2: DA60@A043
    #Observation ID 3: DA60@A043

SB_refant      = 'DV08@A042, DA64@A015, DV15@A047, DA60@A043'

SB_contspws    = '0~15'
SB_spw_mapping = [
    0,0,0,0,
    4,4,4,4,
    8,8,8,8,
    12,12,12,12
]

#First round of phase-only self-cal
#NOTE: you need .p1 instead of _p1 in the caltable name if you want flagdata to work (i.e., flagging problematic antennas non interactively)...
SB_p1 = prefix+'_SB.p1'
os.system('rm -rf '+SB_p1)
gaincal(
    vis=SB_cont_p0+'.ms',caltable=SB_p1,
    gaintype='G',combine='scan,spw',calmode='p',solint='inf',
    spw=SB_contspws,refant=SB_refant,
    minsnr=3.,minblperant=4
)

#Inspect gain tables interactively and decide whether to manually flag something
#plotms(SB_p1,xaxis='time',yaxis='GainPhase',iteraxis='spw',showgui=True)
#Print calibration png file
#plotms(
#    SB_p1,xaxis='time',yaxis='GainPhase',iteraxis='spw',exprange='all',
#    overwrite=True,showgui=False,
#    plotfile=os.path.join(SB_selfcal_folder,f'{prefix}_SB_p1_gain_phase_vs_time.png')
#)

#Inspect gain tables and decide whether to flag something
plotms(
    SB_p1,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(SB_selfcal_folder,f'{prefix}_SB_p1_gain_phase_vs_time.png'),
)

#Apply the calibration gains and split-off the corrected ms
applycal(
    vis=SB_cont_p0+'.ms',spw=SB_contspws,spwmap=SB_spw_mapping,
    gaintable=[SB_p1],interp='linearPD',calwt=True,applymode='calonly'
)
SB_cont_p1 = SB_cont_p0.replace('p0','p1')
os.system('rm -rf %s.ms*' %SB_cont_p1)
split(vis=SB_cont_p0+'.ms',outputvis=SB_cont_p1+'.ms',datacolumn='corrected')

#CLEAN, again with a ~6sigma threshold
tclean_wrapper(
    vis       = SB_cont_p1+'.ms',
    imagename = SB_cont_p1,
    threshold = '0.0918mJy',
    **SB_tclean_wrapper_kwargs
)
estimate_SNR(SB_cont_p1+'.image',disk_mask=SB_mask,noise_mask=noise_annulus_SB)
generate_image_png(
    SB_cont_p1+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_SB,10*rms_SB],
    save_folder=SB_selfcal_folder
)
#MWC_758_SB_contp1.image
#Beam 0.193 arcsec x 0.146 arcsec (-28.31 deg)
#Flux inside disk mask: 56.81 mJy
#Peak intensity of source: 5.06 mJy/beam
#rms: 1.53e-02 mJy/beam
#Peak SNR: 330.85

#Re-run the first step of phase-only selfcal with infinite integration interval to improve the model and get a better EB alignment
#One step is enough, there is not much improvement later on.
#
#Step .p1_bis
SB_p1_bis = SB_p1.replace('p1','p1_bis')
os.system('rm -rf '+SB_p1_bis)
gaincal(
    vis=SB_cont_p1+'.ms',caltable=SB_p1_bis,
    gaintype='T',combine='scan,spw',calmode='p',solint='inf',
    spw=SB_contspws,refant=SB_refant,
    minsnr=3.,minblperant=4
)

#Inspect gain tables interactively and decide whether to manually flag something
#plotms(SB_p1_bis,xaxis='time',yaxis='GainPhase',iteraxis='spw',showgui=True)
#Print calibration png file
#plotms(
#    SB_p1_bis,xaxis='time',yaxis='GainPhase',iteraxis='spw',exprange='all',
#    overwrite=True,showgui=False,
#    plotfile=os.path.join(SB_selfcal_folder,f'{prefix}_SB_p1_bis_gain_phase_vs_time.png')
#)

#Inspect gain tables and decide whether to flag something
plotms(
    SB_p1_bis,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(SB_selfcal_folder,f'{prefix}_SB_p1_bis_gain_phase_vs_time.png'),
)

#Apply the calibration gains and split-off the corrected ms
applycal(
    vis=SB_cont_p1+'.ms',spw=SB_contspws,spwmap=SB_spw_mapping,
    gaintable=[SB_p1_bis],interp='linearPD',calwt=True,applymode='calonly'
)
SB_cont_p1_bis = SB_cont_p1.replace('p1','p1_bis')
os.system('rm -rf %s.ms*' %SB_cont_p1_bis)
split(vis=SB_cont_p1+'.ms',outputvis=SB_cont_p1_bis+'.ms',datacolumn='corrected')

#CLEAN, again with a ~6sigma threshold
tclean_wrapper(
    vis       = SB_cont_p1_bis+'.ms',
    imagename = SB_cont_p1_bis,
    threshold = '0.0918mJy',
    **SB_tclean_wrapper_kwargs
)
estimate_SNR(SB_cont_p1_bis+'.image',disk_mask=SB_mask,noise_mask=noise_annulus_SB)
generate_image_png(
    SB_cont_p1_bis+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_SB,10*rms_SB],
    save_folder=SB_selfcal_folder
)
#MWC_758_SB_contp1_bis.image
#Beam 0.193 arcsec x 0.146 arcsec (-28.31 deg)
#Flux inside disk mask: 56.82 mJy
#Peak intensity of source: 5.06 mJy/beam
#rms: 1.53e-02 mJy/beam
#Peak SNR: 330.75

#Second round of phase-only self-cal
SB_p2 = SB_p1_bis.replace('p1_bis','p2')
os.system('rm -rf '+SB_p2)
gaincal(
    vis=SB_cont_p1_bis+'.ms',caltable=SB_p2,
    gaintype='T',combine='scan,spw',calmode='p',solint='360s',
    spw=SB_contspws,refant=SB_refant,
    minsnr=3.,minblperant=4
)
#4 of 43 solutions flagged due to SNR < 3 in spw=0 at 2017/12/09/04:53:08.3
#1 of 43 solutions flagged due to SNR < 3 in spw=0 at 2017/12/09/05:02:34.0
#4 of 43 solutions flagged due to SNR < 3 in spw=0 at 2017/12/09/05:14:50.5
#4 of 43 solutions flagged due to SNR < 3 in spw=0 at 2017/12/09/05:23:12.0
#1 of 45 solutions flagged due to SNR < 3 in spw=4 at 2017/12/17/06:11:47.4
#1 of 46 solutions flagged due to SNR < 3 in spw=8 at 2017/12/27/04:10:33.4
#1 of 46 solutions flagged due to SNR < 3 in spw=8 at 2017/12/27/04:31:22.8
#1 of 46 solutions flagged due to SNR < 3 in spw=8 at 2017/12/27/04:35:39.6
#1 of 46 solutions flagged due to SNR < 3 in spw=8 at 2017/12/27/04:39:01.8

#Inspect gain tables interactively and decide whether to manually flag something
#plotms(SB_p2,xaxis='time',yaxis='GainPhase',iteraxis='spw',showgui=True)
#Print calibration png file
#plotms(
#    SB_p2,xaxis='time',yaxis='GainPhase',iteraxis='spw',exprange='all',
#    overwrite=True,showgui=False,
#    plotfile=os.path.join(SB_selfcal_folder,f'{prefix}_SB_p2_gain_phase_vs_time.png')
#)

#Inspect gain tables and decide whether to flag something
plotms(
    SB_p2,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(SB_selfcal_folder,f'{prefix}_SB_p2_gain_phase_vs_time.png'),
)

#Flag problematic antennas
flagdata(vis=SB_p2,mode='manual',spw='4', timerange='2017/12/17/06:00:40~2017/12/17/06:01:00',antenna='DV13')
flagdata(vis=SB_p2,mode='manual',spw='8', timerange='2017/12/27/04:10:00~2017/12/27/04:11:00',antenna='DA42,DA44,DV13')
flagdata(vis=SB_p2,mode='manual',spw='8', timerange='2017/12/27/04:35:00~2017/12/27/04:36:00',antenna='DA44')
flagdata(vis=SB_p2,mode='manual',spw='8', timerange='2017/12/27/04:38:00~2017/12/27/04:40:00',antenna='DA44')

#Inspect gain tables to check if flagging worked
plotms(
    SB_p2,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(SB_selfcal_folder,f'{prefix}_SB_p2_gain_phase_vs_time_flagged.png'),
)

#Apply the calibration gains and split-off the corrected ms
applycal(
    vis=SB_cont_p1_bis+'.ms',spw=SB_contspws,spwmap=SB_spw_mapping,
    gaintable=[SB_p2],interp='linearPD',calwt=True,applymode='calonly'
)
SB_cont_p2 = SB_cont_p1_bis.replace('p1_bis','p2')
os.system('rm -rf %s.ms*' %SB_cont_p2)
split(vis=SB_cont_p1_bis+'.ms',outputvis=SB_cont_p2+'.ms',datacolumn='corrected')

#CLEAN, again with a ~6sigma threshold
tclean_wrapper(
    vis       = SB_cont_p2+'.ms',
    imagename = SB_cont_p2,
    threshold = '0.0792mJy',
    **SB_tclean_wrapper_kwargs
)
estimate_SNR(SB_cont_p2+'.image',disk_mask=SB_mask,noise_mask=noise_annulus_SB)
generate_image_png(
    SB_cont_p2+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_SB,10*rms_SB],
    save_folder=SB_selfcal_folder
)
#MWC_758_SB_contp2.image
#Beam 0.193 arcsec x 0.146 arcsec (-28.31 deg)
#Flux inside disk mask: 57.35 mJy
#Peak intensity of source: 5.34 mJy/beam
#rms: 1.31e-02 mJy/beam
#Peak SNR: 406.88

#Third round of phase-only self-cal
SB_p3 = SB_p2.replace('p2','p3')
os.system('rm -rf '+SB_p3)
gaincal(
    vis=SB_cont_p2+'.ms',caltable=SB_p3,
    gaintype='T',combine='scan,spw',calmode='p',solint='120s', #diff from exoALMA, kept combine='scan' because some scans for the SB_EBs~60sec
    spw=SB_contspws,refant=SB_refant,
    minsnr=3.,minblperant=4
)
#1 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:40:28.5
#1 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:42:29.5
#1 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:44:31.2
#1 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:49:49.4
#2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:51:51.1
#4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:53:08.3
#1 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:58:55.3
#3 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:00:56.3
#1 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:02:58.0
#6 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:04:14.3
#1 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:09:31.4
#3 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:11:32.4
#2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:13:34.1
#5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:14:50.5
#1 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:17:53.0
#1 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:21:55.7
#3 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:23:12.0
#2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:25:41.7
#2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:30:05.7
#1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:49:53.8
#1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:51:57.5
#3 of 43 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:54:57.1
#1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/03:56:58.6
#2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:10:33.4
#8 of 45 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:22:11.4
#1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:31:22.8
#1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:35:43.5
#1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:39:01.8
#1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:40:44.5
#1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:57:48.9
#3 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:09:28.0

#Inspect gain tables interactively and decide whether to manually flag something
#plotms(SB_p3,xaxis='time',yaxis='GainPhase',iteraxis='spw',showgui=True)
#Print calibration png file
#plotms(
#    SB_p3,xaxis='time',yaxis='GainPhase',iteraxis='spw',exprange='all',
#    overwrite=True,showgui=False,
#    plotfile=os.path.join(SB_selfcal_folder,f'{prefix}_SB_p3_gain_phase_vs_time.png')
#)

#Inspect gain tables and decide whether to flag something
plotms(
    SB_p3,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(SB_selfcal_folder,f'{prefix}_SB_p3_gain_phase_vs_time.png'),
)

#Flag problematic antennas
flagdata(vis=SB_p3,mode='manual',spw='0', timerange='2017/12/09/04:58:30~2017/12/09/04:59:00',antenna='DV15')
flagdata(vis=SB_p3,mode='manual',spw='0', timerange='2017/12/09/05:04:00~2017/12/09/05:04:30',antenna='DA45,DA49,PM02,DV10')
flagdata(vis=SB_p3,mode='manual',spw='0', timerange='2017/12/09/05:11:00~2017/12/09/05:12:00',antenna='DA46')
flagdata(vis=SB_p3,mode='manual',spw='0', timerange='2017/12/09/05:13:00~2017/12/09/05:14:00',antenna='DA55')
flagdata(vis=SB_p3,mode='manual',spw='0', timerange='2017/12/09/05:17:30~2017/12/09/05:18:00',antenna='DV10')
flagdata(vis=SB_p3,mode='manual',spw='0', timerange='2017/12/09/05:19:30~2017/12/09/05:20:00',antenna='DA46')
flagdata(vis=SB_p3,mode='manual',spw='0', timerange='2017/12/09/05:23:00~2017/12/09/05:23:20',antenna='DA55')
flagdata(vis=SB_p3,mode='manual',spw='0', timerange='2017/12/09/05:25:30~2017/12/09/05:26:00',antenna='DA51')

flagdata(vis=SB_p3,mode='manual',spw='4', timerange='2017/12/17/05:53:30~2017/12/17/05:54:00',antenna='DA43,DA44')
flagdata(vis=SB_p3,mode='manual',spw='4', timerange='2017/12/17/05:54:30~2017/12/17/05:55:00',antenna='DA41,DA42,DA43,DA44')
flagdata(vis=SB_p3,mode='manual',spw='4', timerange='2017/12/17/05:58:30~2017/12/17/05:59:00',antenna='DA44,DA45,DA50,DA55,DV07,DV13')
flagdata(vis=SB_p3,mode='manual',spw='4', timerange='2017/12/17/06:00:30~2017/12/17/06:01:00',antenna='DV13')
flagdata(vis=SB_p3,mode='manual',spw='4', timerange='2017/12/17/06:02:30~2017/12/17/06:03:00',antenna='DV07,DV13')
flagdata(vis=SB_p3,mode='manual',spw='4', timerange='2017/12/17/06:04:00~2017/12/17/06:04:20',antenna='DV13')
flagdata(vis=SB_p3,mode='manual',spw='4', timerange='2017/12/17/06:08:20~2017/12/17/06:08:40',antenna='DA55,DV07,DV11')
flagdata(vis=SB_p3,mode='manual',spw='4', timerange='2017/12/17/06:11:40~2017/12/17/06:12:00',antenna='DV11')

flagdata(vis=SB_p3,mode='manual',spw='8', timerange='2017/12/27/04:09:00~2017/12/27/04:09:30',antenna='DA42,DV13')
flagdata(vis=SB_p3,mode='manual',spw='8', timerange='2017/12/27/04:10:30~2017/12/27/04:11:00',antenna='DA42,DA44,DV13')
flagdata(vis=SB_p3,mode='manual',spw='8', timerange='2017/12/27/04:33:30~2017/12/27/04:34:00',antenna='DA44')
flagdata(vis=SB_p3,mode='manual',spw='8', timerange='2017/12/27/04:35:00~2017/12/27/04:36:00',antenna='DA44')
flagdata(vis=SB_p3,mode='manual',spw='8', timerange='2017/12/27/04:37:00~2017/12/27/04:38:00',antenna='DA44')
flagdata(vis=SB_p3,mode='manual',spw='8', timerange='2017/12/27/04:39:00~2017/12/27/04:40:00',antenna='DA44')

flagdata(vis=SB_p3,mode='manual',spw='12',timerange='2017/12/28/03:57:30~2017/12/28/03:58:00',antenna='DA55')
flagdata(vis=SB_p3,mode='manual',spw='12',timerange='2017/12/28/04:04:00~2017/12/28/04:04:30',antenna='DA44,DA55,DA46')
flagdata(vis=SB_p3,mode='manual',spw='12',timerange='2017/12/28/04:09:00~2017/12/28/04:09:30',antenna='DA42,DV13')
flagdata(vis=SB_p3,mode='manual',spw='12',timerange='2017/12/28/04:15:00~2017/12/28/04:15:30',antenna='DA44,DV11,DV13')
flagdata(vis=SB_p3,mode='manual',spw='12',timerange='2017/12/28/04:20:30~2017/12/28/04:21:00',antenna='DA44,DV11,DV13')
flagdata(vis=SB_p3,mode='manual',spw='12',timerange='2017/12/28/04:30:00~2017/12/28/04:30:30',antenna='DV13')

#Inspect gain tables to check if flagging worked
plotms(
    SB_p3,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(SB_selfcal_folder,f'{prefix}_SB_p3_gain_phase_vs_time_flagged.png'),
)

#Apply the calibration gains and split-off the corrected ms
applycal(
    vis=SB_cont_p2+'.ms',spw=SB_contspws,spwmap=SB_spw_mapping,
    gaintable=[SB_p3],interp='linearPD',calwt=True,applymode='calonly'
)
SB_cont_p3 = SB_cont_p2.replace('p2','p3')
os.system('rm -rf %s.ms*' %SB_cont_p3)
split(vis=SB_cont_p2+'.ms',outputvis=SB_cont_p3+'.ms',datacolumn='corrected')

#CLEAN, again with a ~6sigma threshold
tclean_wrapper(
    vis       = SB_cont_p3+'.ms',
    imagename = SB_cont_p3,
    threshold = '0.0768mJy',
    **SB_tclean_wrapper_kwargs
)
estimate_SNR(SB_cont_p3+'.image',disk_mask=SB_mask,noise_mask=noise_annulus_SB)
generate_image_png(
    SB_cont_p3+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_SB,10*rms_SB],
    save_folder=SB_selfcal_folder
)
#MWC_758_SB_contp3.image
#Beam 0.193 arcsec x 0.146 arcsec (-28.31 deg)
#Flux inside disk mask: 57.63 mJy
#Peak intensity of source: 5.52 mJy/beam
#rms: 1.28e-02 mJy/beam
#Peak SNR: 432.31

#Fourth round of phase-only self-cal
SB_p4 = SB_p3.replace('p3','p4')
os.system('rm -rf '+SB_p4)
gaincal(
    vis=SB_cont_p3+'.ms',caltable=SB_p4,
    gaintype='T',combine='spw',calmode='p',solint='60s',
    spw=SB_contspws,refant=SB_refant,
    minsnr=3.,minblperant=4
)
#2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:39:58.6
#2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:40:59.0
#2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:41:59.5
#2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:43:00.0
#1 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:44:00.5
#1 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:45:02.5
#2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:48:19.3
#1 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:49:19.8
#1 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:50:20.3
#3 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:51:20.7
#2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:52:22.8
#4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:53:08.3
#2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:55:02.5
#4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:55:48.0
#4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:58:25.0
#2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:59:25.5
#4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:00:25.9
#3 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:01:26.4
#2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:02:26.9
#5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:03:28.9
#5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:04:14.3
#3 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:09:01.2
#2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:10:01.6
#5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:11:02.1
#1 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:12:02.6
#3 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:13:03.1
#2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:14:05.1
#5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:14:50.5
#2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:17:22.7
#3 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:18:23.2
#2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:19:23.7
#2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:21:24.6
#2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:22:26.6
#2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:23:12.0
#3 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:25:41.7
#3 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:28:50.3
#5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:30:36.2
#1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:35:40.8
#2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:50:24.3
#1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:51:25.3
#1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:52:27.4
#1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:54:45.0
#1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:11:47.4
#2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:02:54.0
#1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:07:45.5
#1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:12:47.6
#1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:14:37.8
#2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:17:39.2
#1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:21:59.3
#1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:31:22.8
#1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:34:13.0
#1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:35:13.5
#1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:36:14.0
#1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:37:14.3
#1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:38:16.3
#1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:39:01.8
#1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:40:29.5
#1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:41:14.9
#1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:57:48.9
#1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:09:16.6

#Inspect gain tables interactively and decide whether to manually flag something
#plotms(SB_p4,xaxis='time',yaxis='GainPhase',iteraxis='spw',showgui=True)
#Print calibration png file
#plotms(
#    SB_p4,xaxis='time',yaxis='GainPhase',iteraxis='spw',exprange='all',
#    overwrite=True,showgui=False,
#    plotfile=os.path.join(SB_selfcal_folder,f'{prefix}_SB_p4_gain_phase_vs_time.png')
#)

#Inspect gain tables and decide whether to flag something
plotms(
    SB_p4,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(SB_selfcal_folder,f'{prefix}_SB_p4_gain_phase_vs_time.png'),
)

#Flag problematic antennas
flagdata(vis=SB_p4,mode='manual',spw='0', timerange='2017/12/09/04:52:20~2017/12/09/04:52:25',antenna='DV10')
flagdata(vis=SB_p4,mode='manual',spw='0', timerange='2017/12/09/04:54:30~2017/12/09/04:55:30',antenna='DV15')
flagdata(vis=SB_p4,mode='manual',spw='0', timerange='2017/12/09/04:55:30~2017/12/09/04:56:00',antenna='DV11')
flagdata(vis=SB_p4,mode='manual',spw='0', timerange='2017/12/09/05:04:00~2017/12/09/05:04:20',antenna='DA45,DA50,PM02,DA49,DV10')
flagdata(vis=SB_p4,mode='manual',spw='0', timerange='2017/12/09/05:14:00~2017/12/09/05:14:20',antenna='DA55')
flagdata(vis=SB_p4,mode='manual',spw='0', timerange='2017/12/09/05:17:20~2017/12/09/05:17:40',antenna='DV10')
flagdata(vis=SB_p4,mode='manual',spw='0', timerange='2017/12/09/05:20:20~2017/12/09/05:20:40',antenna='DA46')
flagdata(vis=SB_p4,mode='manual',spw='0', timerange='2017/12/09/05:21:20~2017/12/09/05:21:40',antenna='DA46')
flagdata(vis=SB_p4,mode='manual',spw='0', timerange='2017/12/09/05:22:00~2017/12/09/05:23:00',antenna='DA46')
flagdata(vis=SB_p4,mode='manual',spw='0', timerange='2017/12/09/05:23:10~2017/12/09/05:23:20',antenna='DA55,DV10')
flagdata(vis=SB_p4,mode='manual',spw='0', timerange='2017/12/09/05:25:00~2017/12/09/05:26:00',antenna='DA51,DA46')

flagdata(vis=SB_p4,mode='manual',spw='4', timerange='2017/12/17/05:49:20~2017/12/17/05:49:30',antenna='DA55')
flagdata(vis=SB_p4,mode='manual',spw='4', timerange='2017/12/17/05:53:10~2017/12/17/05:53:30',antenna='DA55')
flagdata(vis=SB_p4,mode='manual',spw='4', timerange='2017/12/17/05:53:30~2017/12/17/05:54:00',antenna='DA43,DA44')
flagdata(vis=SB_p4,mode='manual',spw='4', timerange='2017/12/17/05:54:30~2017/12/17/05:55:00',antenna='DA41,DA42,DA43,DA44,DA55,DV13')
flagdata(vis=SB_p4,mode='manual',spw='4', timerange='2017/12/17/05:58:10~2017/12/17/05:58:20',antenna='DA44,DA45,DA55,DA50,DA55,DV07,DV13')
flagdata(vis=SB_p4,mode='manual',spw='4', timerange='2017/12/17/05:59:10~2017/12/17/05:59:20',antenna='DA44,DA45,DA55,DA50,DA55,DV07,DV13')
flagdata(vis=SB_p4,mode='manual',spw='4', timerange='2017/12/17/06:00:10~2017/12/17/06:00:20',antenna='DV13')
flagdata(vis=SB_p4,mode='manual',spw='4', timerange='2017/12/17/06:01:20~2017/12/17/06:01:30',antenna='DV07,DV13')
flagdata(vis=SB_p4,mode='manual',spw='4', timerange='2017/12/17/06:02:20~2017/12/17/06:02:30',antenna='DV07,DV13')
flagdata(vis=SB_p4,mode='manual',spw='4', timerange='2017/12/17/06:03:20~2017/12/17/06:03:30',antenna='DV07,DV13')
flagdata(vis=SB_p4,mode='manual',spw='4', timerange='2017/12/17/06:04:00~2017/12/17/06:04:10',antenna='DV13')
flagdata(vis=SB_p4,mode='manual',spw='4', timerange='2017/12/17/06:08:40~2017/12/17/06:09:00',antenna='DA55,DV07')
flagdata(vis=SB_p4,mode='manual',spw='4', timerange='2017/12/17/06:09:40~2017/12/17/06:10:00',antenna='DA55,DV07')

flagdata(vis=SB_p4,mode='manual',spw='8', timerange='2017/12/27/04:08:30~2017/12/27/04:09:00',antenna='DA42,DV13')
flagdata(vis=SB_p4,mode='manual',spw='8', timerange='2017/12/27/04:09:30~2017/12/27/04:10:00',antenna='DA42,DA44,DV11,DV13')
flagdata(vis=SB_p4,mode='manual',spw='8', timerange='2017/12/27/04:10:30~2017/12/27/04:11:00',antenna='DA42,DA44,DV13')
flagdata(vis=SB_p4,mode='manual',spw='8', timerange='2017/12/27/04:25:30~2017/12/27/04:25:40',antenna='DV11')
flagdata(vis=SB_p4,mode='manual',spw='8', timerange='2017/12/27/04:33:00~2017/12/27/04:33:20',antenna='DA44')
flagdata(vis=SB_p4,mode='manual',spw='8', timerange='2017/12/27/04:34:00~2017/12/27/04:34:20',antenna='DA44')
flagdata(vis=SB_p4,mode='manual',spw='8', timerange='2017/12/27/04:35:00~2017/12/27/04:35:20',antenna='DA44')
flagdata(vis=SB_p4,mode='manual',spw='8', timerange='2017/12/27/04:36:00~2017/12/27/04:36:20',antenna='DA44')
flagdata(vis=SB_p4,mode='manual',spw='8', timerange='2017/12/27/04:37:00~2017/12/27/04:37:20',antenna='DA44')
flagdata(vis=SB_p4,mode='manual',spw='8', timerange='2017/12/27/04:38:00~2017/12/27/04:38:20',antenna='DA44')
flagdata(vis=SB_p4,mode='manual',spw='8', timerange='2017/12/27/04:39:00~2017/12/27/04:39:20',antenna='DA44')

flagdata(vis=SB_p4,mode='manual',spw='12',timerange='2017/12/28/03:46:20~2017/12/28/03:46:30',antenna='DA46')
flagdata(vis=SB_p4,mode='manual',spw='12',timerange='2017/12/28/03:56:00~2017/12/28/03:56:20',antenna='DV11')
flagdata(vis=SB_p4,mode='manual',spw='12',timerange='2017/12/28/03:57:00~2017/12/28/03:57:10',antenna='DA55,DV11')
flagdata(vis=SB_p4,mode='manual',spw='12',timerange='2017/12/28/03:59:00~2017/12/28/03:59:30',antenna='DA55')
flagdata(vis=SB_p4,mode='manual',spw='12',timerange='2017/12/28/04:00:00~2017/12/28/04:00:20',antenna='DA44,DA55,DV07')
flagdata(vis=SB_p4,mode='manual',spw='12',timerange='2017/12/28/04:02:50~2017/12/28/04:03:00',antenna='DA44,DV07,DA46')
flagdata(vis=SB_p4,mode='manual',spw='12',timerange='2017/12/28/04:03:40~2017/12/28/04:04:00',antenna='DA44,DA55,DA46')
flagdata(vis=SB_p4,mode='manual',spw='12',timerange='2017/12/28/04:04:40~2017/12/28/04:05:00',antenna='DA44')
flagdata(vis=SB_p4,mode='manual',spw='12',timerange='2017/12/28/04:09:10~2017/12/28/04:09:20',antenna='DA42,DV13')
flagdata(vis=SB_p4,mode='manual',spw='12',timerange='2017/12/28/04:21:20~2017/12/28/04:21:40',antenna='DA44,DV11,DV13')
flagdata(vis=SB_p4,mode='manual',spw='12',timerange='2017/12/28/04:22:20~2017/12/28/04:22:40',antenna='DV11,DV13')
flagdata(vis=SB_p4,mode='manual',spw='12',timerange='2017/12/28/04:30:10~2017/12/28/04:30:30',antenna='DV13')

#Inspect gain tables to check if flagging worked
plotms(
    SB_p4,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(SB_selfcal_folder,f'{prefix}_SB_p4_gain_phase_vs_time_flagged.png'),
)

#Apply the calibration gains and split-off the corrected ms
applycal(
    vis=SB_cont_p3+'.ms',spw=SB_contspws,spwmap=SB_spw_mapping,
    gaintable=[SB_p4],interp='linearPD',calwt=True,applymode='calonly'
)
SB_cont_p4 = SB_cont_p3.replace('p3','p4')
os.system('rm -rf %s.ms*' %SB_cont_p4)
split(vis=SB_cont_p3+'.ms',outputvis=SB_cont_p4+'.ms',datacolumn='corrected')

#CLEAN, again with a ~6sigma threshold
tclean_wrapper(
    vis       = SB_cont_p4+'.ms',
    imagename = SB_cont_p4,
    threshold = '0.0768mJy',
    **SB_tclean_wrapper_kwargs
)
estimate_SNR(SB_cont_p4+'.image',disk_mask=SB_mask,noise_mask=noise_annulus_SB)
generate_image_png(
    SB_cont_p4+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_SB,10*rms_SB],
    save_folder=SB_selfcal_folder
)
#MWC_758_SB_contp4.image
#Beam 0.193 arcsec x 0.146 arcsec (-28.31 deg)
#Flux inside disk mask: 57.80 mJy
#Peak intensity of source: 5.60 mJy/beam
#rms: 1.28e-02 mJy/beam
#Peak SNR: 437.52

#Fifth round of phase-only selfcal
SB_p5 = SB_p4.replace('p4','p5')
os.system('rm -rf '+SB_p5)
gaincal(
    vis=SB_cont_p4+'.ms',caltable=SB_p5,
    gaintype='T',combine='spw',calmode='p',solint='20s',
    spw=SB_contspws,refant=SB_refant,
    minsnr=3.,minblperant=4
)
# 4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:39:40.5
# 4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:40:01.7
# 5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:40:19.8
# 3 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:40:41.0
# 4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:41:02.2
# 5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:41:20.3
# 6 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:41:41.5
# 5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:42:02.6
# 3 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:42:20.8
# 5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:42:41.9
# 4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:43:03.1
# 4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:43:21.3
# 4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:43:39.4
# 2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:44:00.6
# 2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:44:21.7
# 3 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:44:41.4
# 4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:44:59.6
# 2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:45:20.7
# 2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:47:00.8
# 4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:47:22.0
# 5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:47:40.2
# 5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:48:01.3
# 4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:48:22.5
# 6 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:48:40.6
# 3 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:49:01.8
# 7 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:49:23.0
# 5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:49:41.1
# 2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:50:02.3
# 5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:50:23.5
# 6 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:50:41.6
# 2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:50:59.7
# 4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:51:20.9
# 9 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:51:42.1
# 5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:52:01.8
# 5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:52:19.9
# 5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:52:41.1
# 7 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:53:02.2
# 8 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:53:17.4
# 4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:54:44.5
# 6 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:55:05.7
# 7 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:55:23.8
# 6 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:55:45.0
#13 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:56:00.1
# 5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:58:06.8
# 7 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:58:28.0
# 8 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:58:46.1
# 4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:59:07.3
# 6 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:59:28.5
# 7 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:59:46.6
# 6 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:00:07.8
# 7 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:00:28.9
# 6 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:00:47.1
# 4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:01:08.3
# 7 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:01:29.4
# 4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:01:47.6
# 9 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:02:05.7
# 2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:02:26.9
# 8 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:02:48.0
# 9 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:03:07.7
# 8 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:03:25.9
# 8 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:03:47.0
# 9 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:04:08.2
#10 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:04:23.3
# 5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:08:43.0
# 7 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:09:04.2
# 6 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:09:22.3
# 4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:09:43.5
# 6 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:10:04.7
# 4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:10:22.8
# 5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:10:44.0
# 6 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:11:05.1
# 6 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:11:23.3
# 5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:11:44.5
# 9 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:12:05.6
# 6 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:12:23.8
# 9 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:12:41.9
# 4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:13:03.1
# 6 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:13:24.2
# 6 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:13:43.9
# 4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:14:02.1
# 6 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:14:23.2
# 8 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:14:44.4
#10 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:14:59.5
# 3 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:17:04.6
# 9 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:17:25.7
# 3 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:17:43.9
# 5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:18:05.0
# 7 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:18:26.2
# 7 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:18:44.4
# 5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:19:05.5
# 9 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:19:26.7
# 6 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:19:44.8
# 6 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:20:06.0
# 5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:20:27.2
# 8 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:20:45.3
# 7 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:21:03.5
# 5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:21:24.6
# 3 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:21:45.8
# 6 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:22:05.5
# 4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:22:23.6
# 4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:22:44.8
# 3 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:23:06.0
#10 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:23:21.1
# 3 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:25:23.5
# 6 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:25:44.7
# 7 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:26:02.8
# 6 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:27:31.8
# 5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:27:53.0
# 9 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:28:11.1
# 5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:28:32.3
# 4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:28:53.5
# 7 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:29:11.6
# 4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:29:32.8
# 6 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:29:53.9
# 5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:30:12.1
# 3 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:30:33.2
#12 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:30:48.4
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:29:33.5
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:29:54.6
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:30:12.8
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:30:33.9
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:30:55.1
# 4 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:31:13.2
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:31:34.4
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:31:55.6
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:32:13.7
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:32:34.9
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:32:56.1
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:33:32.3
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:33:53.5
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:34:14.7
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:34:34.4
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:34:52.6
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:35:13.7
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:35:34.7
# 3 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:35:50.1
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:37:12.1
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:37:33.3
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:37:51.4
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:38:12.6
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:38:33.7
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:38:51.9
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:39:13.1
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:39:34.2
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:39:52.2
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:40:34.6
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:40:52.8
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:41:11.0
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:41:53.3
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:42:13.0
# 3 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:42:31.2
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:42:52.3
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:43:13.5
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:43:28.5
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:44:30.9
# 6 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:45:46.5
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:47:26.4
# 3 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:47:44.6
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:48:05.8
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:48:26.9
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:48:45.1
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:49:27.3
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:49:45.6
# 4 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:50:06.6
# 4 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:50:27.9
# 4 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:50:45.9
# 4 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:51:04.3
# 3 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:51:25.2
# 4 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:51:46.3
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:52:06.3
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:52:24.3
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:52:45.5
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:53:06.7
# 5 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:53:21.7
# 3 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:54:02.8
# 3 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:54:21.0
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:54:42.2
# 3 of 43 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:54:57.1
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:58:00.6
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:58:21.8
# 3 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:58:39.9
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:59:01.1
# 4 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:59:22.3
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:59:40.4
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:00:01.6
# 3 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:00:22.6
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:00:40.7
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:01:23.2
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:01:41.0
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:01:59.6
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:02:20.6
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:02:41.9
# 3 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:03:01.4
# 3 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:03:19.5
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:04:02.2
# 3 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:04:17.3
# 3 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:06:01.1
# 3 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:06:19.2
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:06:40.4
# 4 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:07:01.5
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:07:19.7
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:07:40.5
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:08:02.5
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:08:20.6
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:08:41.4
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:09:02.5
# 3 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:09:20.6
# 3 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:09:59.6
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:10:21.1
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:10:40.8
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:10:58.9
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:11:20.4
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:11:41.1
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:11:56.4
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:12:58.9
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:13:20.1
# 3 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:13:38.2
# 6 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:14:14.5
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:15:34.4
# 3 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:15:55.6
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:16:13.7
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/03:56:46.6
# 4 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/03:57:07.7
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/03:57:25.9
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/03:57:47.0
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/03:58:08.2
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/03:58:26.3
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/03:58:47.5
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/03:59:08.7
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:00:27.3
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:00:45.5
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:01:06.6
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:01:27.8
# 4 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:01:47.5
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:02:05.6
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:02:26.8
# 4 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:02:47.9
# 4 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:03:03.1
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:04:26.0
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:04:47.1
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:05:05.3
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:05:47.6
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:06:05.8
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:06:26.9
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:06:48.1
# 5 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:07:06.2
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:07:27.4
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:07:48.6
# 4 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:08:06.7
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:08:24.9
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:09:07.2
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:09:26.9
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:09:45.0
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:10:06.2
# 4 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:10:27.4
# 4 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:10:42.5
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:11:44.1
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:12:05.2
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:12:23.4
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:12:44.5
# 4 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:12:59.7
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:14:19.6
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:14:40.8
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:14:58.9
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:15:20.1
# 4 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:15:41.3
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:15:59.4
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:16:20.6
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:16:41.8
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:16:59.9
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:17:21.1
# 4 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:17:42.3
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:18:00.4
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:18:18.5
# 4 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:18:39.7
# 4 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:19:00.9
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:19:20.5
# 4 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:19:38.7
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:19:59.9
# 4 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:20:36.1
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:20:56.0
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:21:17.0
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:21:35.1
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:21:56.3
# 8 of 45 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:22:11.4
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:25:15.4
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:25:36.5
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:25:54.7
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:26:15.8
# 4 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:26:37.0
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:26:55.1
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:27:16.3
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:27:37.5
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:28:16.8
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:28:38.0
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:28:56.1
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:29:35.4
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:29:56.6
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:30:16.5
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:30:34.4
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:30:55.6
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:31:16.5
# 4 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:31:31.8
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:32:54.5
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:33:15.6
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:33:33.8
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:33:54.9
# 4 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:34:16.1
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:34:34.3
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:34:55.4
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:35:16.6
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:35:34.7
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:35:55.9
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:36:17.1
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:36:35.2
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:36:53.4
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:37:14.6
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:37:35.4
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:37:55.5
# 4 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:38:13.4
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:38:34.7
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:38:55.9
# 4 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:39:11.1
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:40:11.4
# 4 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:40:32.5
# 5 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:40:50.7
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:41:11.8
# 6 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:41:27.0
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:43:06.9
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:43:25.1
# 2 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:44:03.0
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:45:03.6
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:45:24.7
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:45:42.9
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:47:25.7
# 2 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:48:02.0
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:48:23.1
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:48:44.3
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:49:22.1
# 2 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:50:04.5
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:52:02.7
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:52:20.8
# 2 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:54:03.7
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:54:21.8
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:55:22.3
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:55:40.4
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:56:01.6
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:57:00.6
# 2 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:57:42.9
# 2 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:57:58.0
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:59:00.0
# 2 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:59:21.1
# 3 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:00:15.6
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:01:36.3
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:03:58.4
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:05:35.1
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:06:55.3
# 2 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:07:52.8
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:09:13.5
# 3 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:09:28.0
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:13:11.1
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:14:54.0
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:16:12.6
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:17:13.1
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:20:31.2
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:20:49.3
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:21:31.7
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:22:11.0
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:24:51.5
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:25:11.3
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:25:28.9
# 2 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:26:26.5
# 3 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:28:43.0
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:30:22.5
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:30:40.6

#Inspect gain tables interactively and decide whether to manually flag something
#plotms(SB_p5,xaxis='time',yaxis='GainPhase',iteraxis='spw',showgui=True)
#Print calibration png file
#plotms(
#    SB_p5,xaxis='time',yaxis='GainPhase',iteraxis='spw',exprange='all',
#    overwrite=True,showgui=False,
#    plotfile=os.path.join(SB_selfcal_folder,f'{prefix}_SB_p5_gain_phase_vs_time.png')
#)

#Inspect gain tables and decide whether to flag something
plotms(
    SB_p5,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(SB_selfcal_folder,f'{prefix}_SB_p5_gain_phase_vs_time.png'),
)

#Flag problematic antennas
flagdata(vis=SB_p5,mode='manual',spw='0', timerange='2017/12/09/04:40:35~2017/12/09/04:40:45',antenna='DV10')
flagdata(vis=SB_p5,mode='manual',spw='0', timerange='2017/12/09/04:41:00~2017/12/09/04:41:05',antenna='DV11')
flagdata(vis=SB_p5,mode='manual',spw='0', timerange='2017/12/09/04:47:39~2017/12/09/04:47:41',antenna='DA46')
flagdata(vis=SB_p5,mode='manual',spw='0', timerange='2017/12/09/04:54:40~2017/12/09/04:54:50',antenna='DV12,DV11,DV15')
flagdata(vis=SB_p5,mode='manual',spw='0', timerange='2017/12/09/05:02:05~2017/12/09/05:02:10',antenna='DA45,DA50,DA51')
flagdata(vis=SB_p5,mode='manual',spw='0', timerange='2017/12/09/05:02:45~2017/12/09/05:02:50',antenna='DA45')
flagdata(vis=SB_p5,mode='manual',spw='0', timerange='2017/12/09/05:03:20~2017/12/09/05:03:30',antenna='DA49')
flagdata(vis=SB_p5,mode='manual',spw='0', timerange='2017/12/09/05:04:20~2017/12/09/05:04:25',antenna='DA45,DA49,DA46,DA47')
flagdata(vis=SB_p5,mode='manual',spw='0', timerange='2017/12/09/05:08:40~2017/12/09/05:08:45',antenna='DA45,DA49,DA50,DV07,PM02')
flagdata(vis=SB_p5,mode='manual',spw='0', timerange='2017/12/09/05:09:00~2017/12/09/05:09:05',antenna='PM02')
flagdata(vis=SB_p5,mode='manual',spw='0', timerange='2017/12/09/05:17:40~2017/12/09/05:17:50',antenna='DV10')
flagdata(vis=SB_p5,mode='manual',spw='0', timerange='2017/12/09/05:20:20~2017/12/09/05:20:40',antenna='DA46')
flagdata(vis=SB_p5,mode='manual',spw='0', timerange='2017/12/09/05:20:40~2017/12/09/05:20:50',antenna='DA50')
flagdata(vis=SB_p5,mode='manual',spw='0', timerange='2017/12/09/05:22:20~2017/12/09/05:22:25',antenna='DA46')
flagdata(vis=SB_p5,mode='manual',spw='0', timerange='2017/12/09/05:22:40~2017/12/09/05:22:50',antenna='DA55')
flagdata(vis=SB_p5,mode='manual',spw='0', timerange='2017/12/09/05:25:20~2017/12/09/05:25:40',antenna='DA51')
flagdata(vis=SB_p5,mode='manual',spw='0', timerange='2017/12/09/05:25:40~2017/12/09/05:25:50',antenna='DA51')
flagdata(vis=SB_p5,mode='manual',spw='0', timerange='2017/12/09/05:26:00~2017/12/09/05:26:10',antenna='DA51')
flagdata(vis=SB_p5,mode='manual',spw='0', timerange='2017/12/09/05:30:30~2017/12/09/05:30:35',antenna='DA46')

flagdata(vis=SB_p5,mode='manual',spw='4', timerange='2017/12/17/05:49:05~2017/12/17/05:49:10',antenna='DA55')
flagdata(vis=SB_p5,mode='manual',spw='4', timerange='2017/12/17/05:53:40~2017/12/17/05:53:45',antenna='DA43,DA44')
flagdata(vis=SB_p5,mode='manual',spw='4', timerange='2017/12/17/05:54:00~2017/12/17/05:54:05',antenna='DA43,DA44')
flagdata(vis=SB_p5,mode='manual',spw='4', timerange='2017/12/17/05:54:20~2017/12/17/05:54:25',antenna='DA43')
flagdata(vis=SB_p5,mode='manual',spw='4', timerange='2017/12/17/05:54:40~2017/12/17/05:54:45',antenna='DA41,DA42,DA43,DA44,DA55,DV13')
flagdata(vis=SB_p5,mode='manual',spw='4', timerange='2017/12/17/05:54:55~2017/12/17/05:55:00',antenna='DA41,DA42,DA43,DA44')
flagdata(vis=SB_p5,mode='manual',spw='4', timerange='2017/12/17/05:58:00~2017/12/17/05:58:05',antenna='DA44,DV07,DV13')
flagdata(vis=SB_p5,mode='manual',spw='4', timerange='2017/12/17/05:58:20~2017/12/17/05:58:25',antenna='DA44,DA45,DA50,DA55,DV07,DV13')
flagdata(vis=SB_p5,mode='manual',spw='4', timerange='2017/12/17/05:58:20~2017/12/17/05:58:25',antenna='DA44,DA45,DA50,DA55,DV07,DV13')
flagdata(vis=SB_p5,mode='manual',spw='4', timerange='2017/12/17/05:58:35~2017/12/17/05:58:40',antenna='DA45,DA50,DV07,DV13')
flagdata(vis=SB_p5,mode='manual',spw='4', timerange='2017/12/17/05:59:00~2017/12/17/05:59:05',antenna='DA44,DA45,DA50,DV07,DV13')
flagdata(vis=SB_p5,mode='manual',spw='4', timerange='2017/12/17/05:59:20~2017/12/17/05:59:30',antenna='DA45,DA50,DV07')
flagdata(vis=SB_p5,mode='manual',spw='4', timerange='2017/12/17/05:59:40~2017/12/17/05:59:45',antenna='DV07,DV13')
flagdata(vis=SB_p5,mode='manual',spw='4', timerange='2017/12/17/06:00:00~2017/12/17/06:00:10',antenna='DV13')
flagdata(vis=SB_p5,mode='manual',spw='4', timerange='2017/12/17/06:00:20~2017/12/17/06:00:25',antenna='DV13')
flagdata(vis=SB_p5,mode='manual',spw='4', timerange='2017/12/17/06:00:30~2017/12/17/06:00:45',antenna='DV13')
flagdata(vis=SB_p5,mode='manual',spw='4', timerange='2017/12/17/06:01:00~2017/12/17/06:01:05',antenna='DV13')
flagdata(vis=SB_p5,mode='manual',spw='4', timerange='2017/12/17/06:01:20~2017/12/17/06:01:25',antenna='DV07,DV13')
flagdata(vis=SB_p5,mode='manual',spw='4', timerange='2017/12/17/06:01:40~2017/12/17/06:01:45',antenna='DV07,DV13')
flagdata(vis=SB_p5,mode='manual',spw='4', timerange='2017/12/17/06:01:55~2017/12/17/06:02:00',antenna='DV07,DV13')
flagdata(vis=SB_p5,mode='manual',spw='4', timerange='2017/12/17/06:02:10~2017/12/17/06:02:30',antenna='DV07,DV13')
flagdata(vis=SB_p5,mode='manual',spw='4', timerange='2017/12/17/06:02:40~2017/12/17/06:02:50',antenna='DV07,DV13')
flagdata(vis=SB_p5,mode='manual',spw='4', timerange='2017/12/17/06:03:00~2017/12/17/06:03:10',antenna='DV07,DV13')
flagdata(vis=SB_p5,mode='manual',spw='4', timerange='2017/12/17/06:03:10~2017/12/17/06:03:20',antenna='DV13')
flagdata(vis=SB_p5,mode='manual',spw='4', timerange='2017/12/17/06:03:40~2017/12/17/06:03:50',antenna='DV07,DV13')
flagdata(vis=SB_p5,mode='manual',spw='4', timerange='2017/12/17/06:04:00~2017/12/17/06:04:10',antenna='DV13')
flagdata(vis=SB_p5,mode='manual',spw='4', timerange='2017/12/17/06:08:40~2017/12/17/06:08:50',antenna='DA55,DV07')
flagdata(vis=SB_p5,mode='manual',spw='4', timerange='2017/12/17/06:09:00~2017/12/17/06:09:10',antenna='DA55,DV07')
flagdata(vis=SB_p5,mode='manual',spw='4', timerange='2017/12/17/06:09:20~2017/12/17/06:09:30',antenna='DV07')
flagdata(vis=SB_p5,mode='manual',spw='4', timerange='2017/12/17/06:09:30~2017/12/17/06:09:40',antenna='DA45,DV07')
flagdata(vis=SB_p5,mode='manual',spw='4', timerange='2017/12/17/06:04:00~2017/12/17/06:04:10',antenna='DV13')
flagdata(vis=SB_p5,mode='manual',spw='4', timerange='2017/12/17/06:09:50~2017/12/17/06:10:00',antenna='DV07')
flagdata(vis=SB_p5,mode='manual',spw='4', timerange='2017/12/17/06:10:20~2017/12/17/06:10:25',antenna='DV07')
flagdata(vis=SB_p5,mode='manual',spw='4', timerange='2017/12/17/06:11:50~2017/12/17/06:12:00',antenna='DV11')

flagdata(vis=SB_p5,mode='manual',spw='8', timerange='2017/12/27/04:07:40~2017/12/27/04:07:50',antenna='DA42')
flagdata(vis=SB_p5,mode='manual',spw='8', timerange='2017/12/27/04:08:00~2017/12/27/04:08:10',antenna='DA42')
flagdata(vis=SB_p5,mode='manual',spw='8', timerange='2017/12/27/04:08:20~2017/12/27/04:08:30',antenna='DA42,DV13')
flagdata(vis=SB_p5,mode='manual',spw='8', timerange='2017/12/27/04:08:40~2017/12/27/04:08:50',antenna='DA42,DV13')
flagdata(vis=SB_p5,mode='manual',spw='8', timerange='2017/12/27/04:09:00~2017/12/27/04:09:10',antenna='DA42,DV13')
flagdata(vis=SB_p5,mode='manual',spw='8', timerange='2017/12/27/04:09:20~2017/12/27/04:09:30',antenna='DA42,DV13')
flagdata(vis=SB_p5,mode='manual',spw='8', timerange='2017/12/27/04:09:40~2017/12/27/04:09:50',antenna='DA42,DA44,DV13')
flagdata(vis=SB_p5,mode='manual',spw='8', timerange='2017/12/27/04:10:00~2017/12/27/04:10:10',antenna='DA42,DA44,DV13')
flagdata(vis=SB_p5,mode='manual',spw='8', timerange='2017/12/27/04:10:20~2017/12/27/04:10:30',antenna='DA42,DA44,DA55,DV13')
flagdata(vis=SB_p5,mode='manual',spw='8', timerange='2017/12/27/04:10:40~2017/12/27/04:10:50',antenna='DA42,DA44,DV13')
flagdata(vis=SB_p5,mode='manual',spw='8', timerange='2017/12/27/04:11:40~2017/12/27/04:11:50',antenna='DA42,DA44,DV13')
flagdata(vis=SB_p5,mode='manual',spw='8', timerange='2017/12/27/04:25:30~2017/12/27/04:25:40',antenna='DV11')
flagdata(vis=SB_p5,mode='manual',spw='8', timerange='2017/12/27/04:25:50~2017/12/27/04:26:00',antenna='DA55')
flagdata(vis=SB_p5,mode='manual',spw='8', timerange='2017/12/27/04:29:50~2017/12/27/04:30:00',antenna='DV07,DV11')
flagdata(vis=SB_p5,mode='manual',spw='8', timerange='2017/12/27/04:33:10~2017/12/27/04:33:20',antenna='DA44')
flagdata(vis=SB_p5,mode='manual',spw='8', timerange='2017/12/27/04:33:30~2017/12/27/04:33:40',antenna='DA44,DA45')
flagdata(vis=SB_p5,mode='manual',spw='8', timerange='2017/12/27/04:34:10~2017/12/27/04:34:20',antenna='DA44')
flagdata(vis=SB_p5,mode='manual',spw='8', timerange='2017/12/27/04:34:50~2017/12/27/04:35:00',antenna='DA44')
flagdata(vis=SB_p5,mode='manual',spw='8', timerange='2017/12/27/04:35:30~2017/12/27/04:35:40',antenna='DA44')
flagdata(vis=SB_p5,mode='manual',spw='8', timerange='2017/12/27/04:35:50~2017/12/27/04:36:00',antenna='DA44')
flagdata(vis=SB_p5,mode='manual',spw='8', timerange='2017/12/27/04:36:10~2017/12/27/04:36:20',antenna='DA44')
flagdata(vis=SB_p5,mode='manual',spw='8', timerange='2017/12/27/04:36:30~2017/12/27/04:36:40',antenna='DA44')
flagdata(vis=SB_p5,mode='manual',spw='8', timerange='2017/12/27/04:37:10~2017/12/27/04:37:20',antenna='DA44')
flagdata(vis=SB_p5,mode='manual',spw='8', timerange='2017/12/27/04:37:30~2017/12/27/04:37:40',antenna='DA44')
flagdata(vis=SB_p5,mode='manual',spw='8', timerange='2017/12/27/04:37:30~2017/12/27/04:37:40',antenna='DA44')
flagdata(vis=SB_p5,mode='manual',spw='8', timerange='2017/12/27/04:37:50~2017/12/27/04:38:00',antenna='DA44')
flagdata(vis=SB_p5,mode='manual',spw='8', timerange='2017/12/27/04:38:30~2017/12/27/04:38:40',antenna='DA44,DV13')
flagdata(vis=SB_p5,mode='manual',spw='8', timerange='2017/12/27/04:38:50~2017/12/27/04:39:00',antenna='DA44')
flagdata(vis=SB_p5,mode='manual',spw='8', timerange='2017/12/27/04:39:10~2017/12/27/04:39:20',antenna='DA44')
flagdata(vis=SB_p5,mode='manual',spw='8', timerange='2017/12/27/04:41:10~2017/12/27/04:41:20',antenna='DV11')

flagdata(vis=SB_p5,mode='manual',spw='12',timerange='2017/12/28/03:45:40~2017/12/28/03:45:45',antenna='DA55')
flagdata(vis=SB_p5,mode='manual',spw='12',timerange='2017/12/28/03:46:40~2017/12/28/03:46:50',antenna='DA55,DA46')
flagdata(vis=SB_p5,mode='manual',spw='12',timerange='2017/12/28/03:49:00~2017/12/28/03:49:10',antenna='DA44')
flagdata(vis=SB_p5,mode='manual',spw='12',timerange='2017/12/28/03:54:40~2017/12/28/03:54:45',antenna='DA55')
flagdata(vis=SB_p5,mode='manual',spw='12',timerange='2017/12/28/03:55:20~2017/12/28/03:55:30',antenna='DA55')
flagdata(vis=SB_p5,mode='manual',spw='12',timerange='2017/12/28/03:56:20~2017/12/28/03:56:30',antenna='DA55,DV11')
flagdata(vis=SB_p5,mode='manual',spw='12',timerange='2017/12/28/03:56:40~2017/12/28/03:56:50',antenna='DA55,DV11')
flagdata(vis=SB_p5,mode='manual',spw='12',timerange='2017/12/28/03:57:20~2017/12/28/03:57:25',antenna='DA55')
flagdata(vis=SB_p5,mode='manual',spw='12',timerange='2017/12/28/03:58:50~2017/12/28/03:59:00',antenna='DA45,DA50,DA55')
flagdata(vis=SB_p5,mode='manual',spw='12',timerange='2017/12/28/03:59:10~2017/12/28/03:59:30',antenna='DA50')
flagdata(vis=SB_p5,mode='manual',spw='12',timerange='2017/12/28/03:59:50~2017/12/28/04:00:10',antenna='DA44,DA55')
flagdata(vis=SB_p5,mode='manual',spw='12',timerange='2017/12/28/04:00:10~2017/12/28/04:00:20',antenna='DA42,DA44,DV07')
flagdata(vis=SB_p5,mode='manual',spw='12',timerange='2017/12/28/04:01:30~2017/12/28/04:01:40',antenna='DA44,DA56,DV07,DA46')
flagdata(vis=SB_p5,mode='manual',spw='12',timerange='2017/12/28/04:01:50~2017/12/28/04:02:00',antenna='DV07')
flagdata(vis=SB_p5,mode='manual',spw='12',timerange='2017/12/28/04:02:10~2017/12/28/04:02:20',antenna='DA55,DA48,DA51')
flagdata(vis=SB_p5,mode='manual',spw='12',timerange='2017/12/28/04:03:10~2017/12/28/04:03:20',antenna='DA44')
flagdata(vis=SB_p5,mode='manual',spw='12',timerange='2017/12/28/04:03:20~2017/12/28/04:03:40',antenna='DA44,DA55,DA46')
flagdata(vis=SB_p5,mode='manual',spw='12',timerange='2017/12/28/04:03:50~2017/12/28/04:04:00',antenna='DA44,DA46')
flagdata(vis=SB_p5,mode='manual',spw='12',timerange='2017/12/28/04:04:10~2017/12/28/04:04:20',antenna='DA44')
flagdata(vis=SB_p5,mode='manual',spw='12',timerange='2017/12/28/04:04:30~2017/12/28/04:04:40',antenna='DA44')
flagdata(vis=SB_p5,mode='manual',spw='12',timerange='2017/12/28/04:04:50~2017/12/28/04:05:00',antenna='DA44')
flagdata(vis=SB_p5,mode='manual',spw='12',timerange='2017/12/28/04:05:10~2017/12/28/04:05:20',antenna='DA44')
flagdata(vis=SB_p5,mode='manual',spw='12',timerange='2017/12/28/04:09:10~2017/12/28/04:09:20',antenna='DA42,DV11,DV13')
flagdata(vis=SB_p5,mode='manual',spw='12',timerange='2017/12/28/04:09:20~2017/12/28/04:09:30',antenna='DA42,DV13')
flagdata(vis=SB_p5,mode='manual',spw='12',timerange='2017/12/28/04:12:20~2017/12/28/04:12:40',antenna='DA42,DA55,DV13')
flagdata(vis=SB_p5,mode='manual',spw='12',timerange='2017/12/28/04:12:50~2017/12/28/04:13:00',antenna='DA55')
flagdata(vis=SB_p5,mode='manual',spw='12',timerange='2017/12/28/04:13:10~2017/12/28/04:13:20',antenna='DV11')
flagdata(vis=SB_p5,mode='manual',spw='12',timerange='2017/12/28/04:13:30~2017/12/28/04:13:35',antenna='DV11')
flagdata(vis=SB_p5,mode='manual',spw='12',timerange='2017/12/28/04:13:50~2017/12/28/04:14:00',antenna='DV11')
flagdata(vis=SB_p5,mode='manual',spw='12',timerange='2017/12/28/04:14:10~2017/12/28/04:14:20',antenna='DV11')
flagdata(vis=SB_p5,mode='manual',spw='12',timerange='2017/12/28/04:17:10~2017/12/28/04:17:15',antenna='DA44,DV13')
flagdata(vis=SB_p5,mode='manual',spw='12',timerange='2017/12/28/04:17:30~2017/12/28/04:17:40',antenna='DA44,DA55,DV02')
flagdata(vis=SB_p5,mode='manual',spw='12',timerange='2017/12/28/04:18:10~2017/12/28/04:18:20',antenna='DA44,DA55')
flagdata(vis=SB_p5,mode='manual',spw='12',timerange='2017/12/28/04:20:40~2017/12/28/04:20:50',antenna='DA44,DV13,DA51')
flagdata(vis=SB_p5,mode='manual',spw='12',timerange='2017/12/28/04:21:10~2017/12/28/04:21:20',antenna='DA44,DV13')
flagdata(vis=SB_p5,mode='manual',spw='12',timerange='2017/12/28/04:21:30~2017/12/28/04:21:40',antenna='DA44,DV13')
flagdata(vis=SB_p5,mode='manual',spw='12',timerange='2017/12/28/04:21:45~2017/12/28/04:21:55',antenna='DV11,DV13')
flagdata(vis=SB_p5,mode='manual',spw='12',timerange='2017/12/28/04:22:05~2017/12/28/04:22:15',antenna='DV11,DV13')
flagdata(vis=SB_p5,mode='manual',spw='12',timerange='2017/12/28/04:22:30~2017/12/28/04:22:40',antenna='DA55,DV13')
flagdata(vis=SB_p5,mode='manual',spw='12',timerange='2017/12/28/04:22:45~2017/12/28/04:22:55',antenna='DV11')
flagdata(vis=SB_p5,mode='manual',spw='12',timerange='2017/12/28/04:28:00~2017/12/28/04:28:10',antenna='DV11')
flagdata(vis=SB_p5,mode='manual',spw='12',timerange='2017/12/28/04:28:40~2017/12/28/04:28:50',antenna='DA44')
flagdata(vis=SB_p5,mode='manual',spw='12',timerange='2017/12/28/04:30:00~2017/12/28/04:30:10',antenna='DV13')
flagdata(vis=SB_p5,mode='manual',spw='12',timerange='2017/12/28/04:30:20~2017/12/28/04:30:30',antenna='DV13')
flagdata(vis=SB_p5,mode='manual',spw='12',timerange='2017/12/28/04:30:35~2017/12/28/04:30:45',antenna='DV13')

#Inspect gain tables to check if flagging worked
plotms(
    SB_p5,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(SB_selfcal_folder,f'{prefix}_SB_p5_gain_phase_vs_time_flagged.png'),
)

#Apply the calibration gains and split-off the corrected ms
applycal(
    vis=SB_cont_p4+'.ms',spw=SB_contspws,spwmap=SB_spw_mapping,
    gaintable=[SB_p5],interp='linearPD',calwt=True,applymode='calonly'
)
SB_cont_p5 = SB_cont_p4.replace('p4','p5')
os.system('rm -rf %s.ms*' %SB_cont_p5)
split(vis=SB_cont_p4+'.ms',outputvis=SB_cont_p5+'.ms',datacolumn='corrected')

#CLEAN, again with a ~6sigma threshold
tclean_wrapper(
    vis       = SB_cont_p5+'.ms',
    imagename = SB_cont_p5,
    threshold = '0.0768mJy',
    **SB_tclean_wrapper_kwargs
)
estimate_SNR(SB_cont_p5+'.image',disk_mask=SB_mask,noise_mask=noise_annulus_SB)
generate_image_png(
    SB_cont_p5+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_SB,10*rms_SB],
    save_folder=SB_selfcal_folder
)
#MWC_758_SB_contp5.image
#Beam 0.193 arcsec x 0.146 arcsec (-28.31 deg)
#Flux inside disk mask: 58.12 mJy
#Peak intensity of source: 5.67 mJy/beam
#rms: 1.28e-02 mJy/beam
#Peak SNR: 442.55

#For each step of the self-cal, check how it improved things
non_self_caled_SB_vis = SB_cont_p0
self_caled_SB_visibilities = {
    'p1':SB_cont_p1,
    'p2':SB_cont_p2,
    'p3':SB_cont_p3,
    'p4':SB_cont_p4,
    'p5':SB_cont_p5
}

#SB_EBs = ('EB0','EB1','EB2','EB3')
#SB_EB_spws = ('0,1,2,3','4,5,6,7','8,9,10,11','12,13,14,15') #fill out by referring to listobs output
#
#for self_cal_step,self_caled_vis in self_caled_SB_visibilities.items():
#    for EB_key,spw in zip(SB_EBs,SB_EB_spws):
#
#        print('selfcal_step: {}'.format(self_cal_step)+', EB: {}'.format(EB_key))
#
#        nametemplate = f'{prefix}_SB_{EB_key}_{self_cal_step}_compare_amp_vs_time'
#        visibilities = [self_caled_vis+'.ms',non_self_caled_SB_vis+'.ms']
#
#        plot_amp_vs_time_comparison(
#            nametemplate=nametemplate,visibilities=visibilities,spw=spw,
#            uvrange=uv_ranges['SB'],output_folder=SB_selfcal_folder
#        )
#
# Essentially need to extract it spw by spw (avgspw doesn't work on casa 6.6.5.1...)
# But since the waterfall features were very mild in the pre-selfcal files chose to avoid to save time and space

all_SB_visibilities = self_caled_SB_visibilities.copy()
all_SB_visibilities['p0'] = SB_cont_p0

SB_flux_ref_EB = 3 #this is SB_EB3

for self_cal_step,vis_name in all_SB_visibilities.items():
    #Split out SB EBs
    vis_ms = vis_name+'.ms'
    nametemplate = vis_ms.replace('.ms','_EB')
    split_all_obs(msfile=vis_ms,nametemplate=nametemplate) 
    #Saving observation 0 of MWC_758_SB_contp0.ms to MWC_758_SB_contp0_EB0.ms
    #Saving observation 1 of MWC_758_SB_contp0.ms to MWC_758_SB_contp0_EB1.ms
    #...
    #Saving observation 0 of MWC_758_SB_contp5.ms to MWC_758_SB_contp5_EB0.ms
    #Saving observation 1 of MWC_758_SB_contp5.ms to MWC_758_SB_contp5_EB1.ms

    #If differences in the time intervals corresponding to the different scans... split by scan!
    #for i in range(number_of_EBs['SB']):
    #    nametemplate = f'{prefix}_SB_cont{self_cal_step}_EB{i}'
    #    os.system(f'rm -rf {nametemplate}*')
    #    split(
    #        vis        = f'{vis}.ms', 
    #        outputvis  = f'{vis}_EB{i}.ms', 
    #        datacolumn = 'data',
    #        keepflags  = False,
    #        scan       = scans[i]
    #    )
    #    listobs(vis=f'{vis}_EB{i}.ms',listfile=f'{vis}_EB{i}.ms.listobs.txt',overwrite=True)

    exported_ms = []
    for i in range(number_of_EBs['SB']):
        EB_vis = f'{nametemplate}{i}.ms'
        listobs(vis=EB_vis,listfile=EB_vis+'.listobs.txt',overwrite=True)
        #Export MS contents into numpy save files
        export_MS(EB_vis)
        exported_ms.append(EB_vis.replace('.ms','.vis.npz'))
        #Measurement set exported to MWC_758_SB_contp0_EB0.vis.npz
        #Measurement set exported to MWC_758_SB_contp0_EB1.vis.npz
        #...
        #Measurement set exported to MWC_758_SB_contp5_EB0.vis.npz
        #Measurement set exported to MWC_758_SB_contp5_EB1.vis.npz

    for i,exp_ms in enumerate(exported_ms):
        # png_filename = f'flux_comparison_SB_EB{i}_{self_cal_step}_to_{flux_ref_EB}.png'
        png_filename = f'flux_comparison_SB_EB{i}_{self_cal_step}_to_EB{SB_flux_ref_EB}.png'
        plot_label = os.path.join(SB_selfcal_folder,png_filename)
        estimate_flux_scale(
            #reference=f'{prefix}_{flux_ref_EB}_initcont.vis.npz',#reference=f'{prefix}_{flux_ref_EB}_initcont_shift.vis.npz',
            reference=f'{nametemplate}{SB_flux_ref_EB}.vis.npz',
            comparison=exp_ms,incl=incl,PA=PA,plot_label=plot_label
        )

    fluxscale = [1.,]*number_of_EBs['SB']
    plot_label = os.path.join(SB_selfcal_folder,f'deprojected_vis_profiles_SB_{self_cal_step}.png')
    plot_deprojected(filelist=exported_ms,fluxscale=fluxscale,PA=PA,incl=incl,show_err=True,plot_label=plot_label)

#ratio          = [0.88988,0.89305,0.89593,0.90156,0.90895,0.92577] #MWC_758_SB_contp0...p5_EB0.vis.npz vs MWC_758_SB_EB3_initcont.vis.npz
#scaling_factor = [0.943  ,0.945  ,0.947  ,0.950  ,0.953  ,0.962  ]
#ratio          = [0.98340,0.97806,0.98762,0.99364,0.99671,0.99951] #MWC_758_SB_contp0...p5_EB1.vis.npz vs MWC_758_SB_EB3_initcont.vis.npz
#scaling_factor = [0.992  ,0.989  ,0.994  ,0.997  ,0.998  ,1.000  ]
#ratio          = [0.93677,0.96069,0.95945,0.96593,0.96917,0.97148] #MWC_758_SB_contp0...p5_EB2.vis.npz vs MWC_758_SB_EB3_initcont.vis.npz
#scaling_factor = [0.968  ,0.980  ,0.980  ,0.983  ,0.984  ,0.986  ]

# SB_EB0 substantial improvement but flux offsets >4%
# SB_EB1,2 show flux offsets <4% no flux-scaling and re-run SB selfcal

#Begin of SB+LB self-cal - iteration 1
#For phase self-cal, clean down to ~6sigma
LB_selfcal_folder = get_figures_folderpath('8_selfcal_SBLB_figures')
make_figures_folder(LB_selfcal_folder)

LB_cont_p0 = prefix+'_SBLB_contp0'
os.system('rm -rf %s.ms*' %LB_cont_p0)
concat(
    vis=[SB_cont_p5+'.ms']+[f'{prefix}_LB_EB{i}_initcont_selfcal.ms' for i in range(number_of_EBs['LB'])],#vis=[SB_cont_p5+'.ms']+[f'{prefix}_LB_EB{i}_initcont_shift.ms' for i in range(number_of_EBs['LB'])],
    concatvis=LB_cont_p0+'.ms',dirtol='0.1arcsec',copypointing=False
)
listobs(vis=LB_cont_p0+'.ms',listfile=LB_cont_p0+'.ms.listobs.txt',overwrite=True)
#2024-12-23 09:24:44     SEVERE  getcell::TIME   Exception Reported: TableProxy::getCell: no such row
#2024-12-23 09:24:47     WARN    concat::::casa  Some but not all of the input MSs are lacking a populated POINTING table:
#    0: MWC_758_SB_contp5.ms
#    The joint dataset will not have a valid POINTING table.
#2024-12-23 09:24:47     WARN    concat::::casa  The setup of the input MSs is not fully consistent. The concatenation may fail
#    and/or the affected columns may contain partially only default data.
#2024-12-23 09:24:47     WARN    concat::::casa  {'MWC_758_LB_EB0_initcont.ms': {'Main': {'present_a': True, 'present_b': True, 'missingcol_a': ['CORRECTED_DATA'], 'missingcol_b': []}}, 'MWC_758_LB_EB1_initcont.ms': {'Main': {'present_a': True, 'present_b': True, 'missingcol_a': ['CORRECTED_DATA'], 'missingcol_b': []}}, 'MWC_758_LB_EB2_initcont.ms': {'Main': {'present_a': True, 'present_b': True, 'missingcol_a': ['CORRECTED_DATA'], 'missingcol_b': []}}, 'MWC_758_LB_EB3_initcont.ms': {'Main': {'present_a': True, 'present_b': True, 'missingcol_a': ['CORRECTED_DATA'], 'missingcol_b': []}}, 'MWC_758_LB_EB4_initcont.ms': {'Main': {'present_a': True, 'present_b': True, 'missingcol_a': ['CORRECTED_DATA'], 'missingcol_b': []}}}

#Define new LB mask using new center read from listobs
mask_pa        = PA
mask_semimajor = 0.8
mask_semiminor = 0.8
mask_ra  = '05h30m27.535947s' 
mask_dec = '25d19m56.615500s'

LB_mask = f'ellipse[[{mask_ra},{mask_dec}], [{mask_semimajor:.3f}arcsec, {mask_semiminor:.3f}arcsec], {mask_pa:.1f}deg]'

noise_annulus_LB = f"annulus[[{mask_ra}, {mask_dec}],['4.arcsec', '6.arcsec']]"

LB_tclean_wrapper_kwargs = {
    'deconvolver':'multiscale','scales':[0,2,4,8,16,24,32],
    'smallscalebias':0.6,'gain':0.3,'cycleniter':300,
    'niter':1000000,'mask':LB_mask,
    'cellsize':'0.005arcsec','imsize':8000,
    'parallel':use_parallel,'savemodel':'modelcolumn',
    'robust':1.0,'interactive':False, #diff from exoALMA, chose robust=1.0 to increase a bit the peak SNR and aid self-cal
    'gridder':'standard',
}

tclean_wrapper(
    vis       = LB_cont_p0+'.ms', 
    imagename = LB_cont_p0,
    threshold = '0.0408mJy',
    **LB_tclean_wrapper_kwargs
)
estimate_SNR(LB_cont_p0+'.image',disk_mask=LB_mask,noise_mask=noise_annulus_LB)
rms_LB = imstat(imagename=LB_cont_p0+'.image',region=noise_annulus_LB)['rms'][0]
generate_image_png(
    LB_cont_p0+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_LB,10*rms_LB],
    save_folder=LB_selfcal_folder
)
#MWC_758_SBLB_contp0.image
#Beam 0.065 arcsec x 0.045 arcsec (7.20 deg)
#Flux inside disk mask: 67.09 mJy
#Peak intensity of source: 1.17 mJy/beam
#rms: 6.81e-03 mJy/beam
#Peak SNR: 171.13

#Look for references antennas from the weblog, and pick the first ones that are listed, overlapping with all EBs
#Since we are combined SB EBs to LB EBs, ref antennas for both SB and LB are needed
#SBs:
#SB_EB0: DV05, DV08, DA64, DA60
#SB_EB1: DA49, DV23, DA64, DV06 
#SB_EB2: DV08, DA60, DV15, DV01
#SB_EB3: DA47, DV08, DV15, DA60
#LBs:
#LB_EB0: PM02, DA61, DA50, DA62
#LB_EB1: PM02, DA50, DA61, DA62
#LB_EB2: DA50, DA61, PM01, DV02
#LB_EB3: DV09, DA50, DV06, DV25
#LB_EB4: DV09, DA50, DV07, DV06
#
#Get station numbers
#For SB we used: SB_refant = 'DV08@A042, DA64@A015, DV15@A047, DA60@A043'
#Bear in mind that here ID0=LB0,ID1=LB1,ID2=LB2,ID3=LB3,ID4=LB4,ID5=SB0,ID6=SB1,ID7=SB2,ID8=SB3
for ref_ant in ('PM02','DA50','DV09','DA61', 'DV08','DA64','DV15','DA60'):
    get_station_numbers(LB_cont_p0+'.ms',ref_ant)
    #Observation ID 0: PM02@A111
    #Observation ID 1: PM02@A111
    #Observation ID 2: PM02@A111
    #Observation ID 3: PM02@A111
    #Observation ID 4: PM02@A111
    #Observation ID 5: PM02@A111
    #Observation ID 6: PM02@T701
    #Observation ID 0: DA50@A108
    #Observation ID 1: DA50@A108
    #Observation ID 2: DA50@A108
    #Observation ID 3: DA50@A108
    #Observation ID 4: DA50@A108
    #Observation ID 5: DA50@A108
    #Observation ID 6: DA50@A108
    #Observation ID 7: DA50@A108
    #Observation ID 8: DA50@A108
    #Observation ID 3: DV09@A007
    #Observation ID 4: DV09@A007
    #Observation ID 5: DV09@A074
    #Observation ID 6: DV09@A074
    #Observation ID 7: DV09@A074
    #Observation ID 8: DV09@A074
    #Observation ID 0: DA61@A015
    #Observation ID 1: DA61@A015
    #Observation ID 2: DA61@A015
    #Observation ID 3: DA61@A015
    #Observation ID 4: DA61@A015
    #Observation ID 5: DA61@A089
    #Observation ID 6: DA61@A089
    #Observation ID 7: DA61@A089
    #Observation ID 8: DA61@A089
    #Observation ID 0: DV08@P411
    #Observation ID 1: DV08@P411
    #Observation ID 2: DV08@P411
    #Observation ID 4: DV08@P411
    #Observation ID 5: DV08@A042
    #Observation ID 6: DV08@A042
    #Observation ID 7: DV08@A042
    #Observation ID 8: DV08@A042
    #Observation ID 0: DA64@P410
    #Observation ID 1: DA64@P410
    #Observation ID 2: DA64@P410
    #Observation ID 3: DA64@P410
    #Observation ID 4: DA64@P410
    #Observation ID 5: DA64@A015
    #Observation ID 6: DA64@A015
    #Observation ID 7: DA64@A015
    #Observation ID 8: DA64@A015
    #Observation ID 0: DV15@A130
    #Observation ID 1: DV15@A130
    #Observation ID 2: DV15@A130
    #Observation ID 3: DV15@A130
    #Observation ID 4: DV15@A130
    #Observation ID 5: DV15@A130
    #Observation ID 6: DV15@A047
    #Observation ID 7: DV15@A047
    #Observation ID 8: DV15@A047
    #Observation ID 0: DA60@S301
    #Observation ID 1: DA60@S301
    #Observation ID 2: DA60@S301
    #Observation ID 3: DA60@S301
    #Observation ID 4: DA60@S301
    #Observation ID 5: DA60@A043
    #Observation ID 6: DA60@A043
    #Observation ID 7: DA60@A043
    #Observation ID 8: DA60@A043

LB_refant = 'PM02@A111, DA50@A108, DV09@A007, DA61@A015, DV08@A042, DA64@A015, DV15@A047, DA60@A043'

LB_contspws    = '0~35'
LB_spw_mapping = [
    0,0,0,0,
    4,4,4,4,
    8,8,8,8,
    12,12,12,12,
    16,16,16,16,
    20,20,20,20,
    24,24,24,24,
    28,28,28,28,
    32,32,32,32,
]

#First round of phase-only self-cal
#NOTE: you need .p1 instead of _p1 in the caltable name if you want flagdata to work (i.e., flagging problematic antennas non interactively)...
LB_p1 = prefix+'_SBLB.p1'
os.system('rm -rf '+LB_p1)
#If you get many flagged solutions, change gaintype to 'T'
gaincal(
    vis=LB_cont_p0+'.ms',caltable=LB_p1,
    gaintype='G',combine='scan,spw',calmode='p',solint='inf',
    spw=LB_contspws,refant=LB_refant,
    minsnr=3.,minblperant=4 #diff from exoALMA, chose minsnr=3 otherwise too many too odd gains (same for the following steps)
)
#13 of 92 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:37:00.7
#15 of 88 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:03:27.8
#19 of 86 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:04:22.3
#23 of 92 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:58:17.8
#27 of 88 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:14:41.2

#Inspect gain tables interactively and decide whether to manually flag something
#plotms(LB_p1,xaxis='time',yaxis='GainPhase',iteraxis='spw',showgui=True)
#Print calibration png file
#plotms(
#    LB_p1,xaxis='time', yaxis='GainPhase',iteraxis='spw',exprange='all',
#    overwrite=True,showgui=False,
#    plotfile=os.path.join(LB_selfcal_folder,f'{prefix}_LB_p1_gain_phase_vs_time.png')
#)

#Inspect gain tables and decide whether to flag something
plotms(
    LB_p1,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(LB_selfcal_folder,f'{prefix}_LB_p1_gain_phase_vs_time.png'),
)

#Flag problematic antennas
flagdata(vis=LB_p1,mode='manual',spw='0', antenna='DA45,DA54,DV24')
flagdata(vis=LB_p1,mode='manual',spw='4', antenna='DA45,DA47,DA54,DV24')
flagdata(vis=LB_p1,mode='manual',spw='8', antenna='DA45,DA47')
flagdata(vis=LB_p1,mode='manual',spw='12',antenna='DA63')
flagdata(vis=LB_p1,mode='manual',spw='20',antenna='DV15,DV10')
flagdata(vis=LB_p1,mode='manual',spw='32',antenna='DA44,DV11')

#Inspect gain tables to check if flagging worked
plotms(
    LB_p1,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(LB_selfcal_folder,f'{prefix}_LB_p1_gain_phase_vs_time_flagged.png'),
)

#Apply the calibration gains and split-off the corrected ms
applycal(
    vis=LB_cont_p0+'.ms',spw=LB_contspws,spwmap=LB_spw_mapping,
    gaintable=[LB_p1],interp='linearPD',calwt=True,applymode='calonly'
)
LB_cont_p1 = LB_cont_p0.replace('p0','p1')
os.system('rm -rf %s.ms*' %LB_cont_p1)
split(vis=LB_cont_p0+'.ms',outputvis=LB_cont_p1+'.ms',datacolumn='corrected')

#CLEAN, again with a ~6sigma threshold
tclean_wrapper(
    vis       = LB_cont_p1+'.ms',
    imagename = LB_cont_p1,
    threshold = '0.0399mJy',
    **LB_tclean_wrapper_kwargs
)
estimate_SNR(LB_cont_p1+'.image',disk_mask=LB_mask,noise_mask=noise_annulus_LB)
generate_image_png(
    LB_cont_p1+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_LB,10*rms_LB],
    save_folder=LB_selfcal_folder
)
#MWC_758_SBLB_contp1.image
#Beam 0.065 arcsec x 0.045 arcsec (7.20 deg)
#Flux inside disk mask: 66.16 mJy
#Peak intensity of source: 1.16 mJy/beam
#rms: 6.67e-03 mJy/beam
#Peak SNR: 173.27

#Re-run the first step of phase-only selfcal with infinite integration interval to improve the model and get a better EB alignment
#One step is enough, there is not much improvement later on.
#
#Step .p1_bis
LB_p1_bis = LB_p1.replace('p1','p1_bis')
os.system('rm -rf '+LB_p1_bis)
gaincal(
    vis=LB_cont_p1+'.ms',caltable=LB_p1_bis,
    gaintype='T',combine='scan,spw',calmode='p',solint='inf',
    spw=LB_contspws,refant=LB_refant,
    minsnr=3.,minblperant=4
)
# 1 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:36:59.4
# 5 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:03:28.1
# 9 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:04:22.2
# 9 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:58:17.6
#11 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:14:43.4

#Inspect gain tables interactively and decide whether to manually flag something
#plotms(LB_p1_bis,xaxis='time',yaxis='GainPhase',iteraxis='spw',showgui=True)
#Print calibration png file
#plotms(
#    LB_p1_bis,xaxis='time',yaxis='GainPhase',iteraxis='spw',exprange='all',
#    overwrite=True,showgui=False,
#    plotfile=os.path.join(LB_selfcal_folder,f'{prefix}_LB_p1_bis_gain_phase_vs_time.png')
#)

#Inspect gain tables and decide whether to flag something
plotms(
    LB_p1_bis,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(LB_selfcal_folder,f'{prefix}_LB_p1_bis_gain_phase_vs_time.png'),
)

#Flag problematic antennas
flagdata(vis=LB_p1_bis,mode='manual',spw='4', antenna='DA54,DV24')
flagdata(vis=LB_p1_bis,mode='manual',spw='8', antenna='DA45')
flagdata(vis=LB_p1_bis,mode='manual',spw='12',antenna='DA59,DA63')
flagdata(vis=LB_p1_bis,mode='manual',spw='16',antenna='DA45')
flagdata(vis=LB_p1_bis,mode='manual',spw='20',antenna='DV10,DV15')
flagdata(vis=LB_p1_bis,mode='manual',spw='32',antenna='DA44,DV11')

#Inspect gain tables to check if flagging worked
plotms(
    LB_p1_bis,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(LB_selfcal_folder,f'{prefix}_LB_p1_bis_gain_phase_vs_time_flagged.png'),
)

#Apply the calibration gains and split-off the corrected ms
applycal(
    vis=LB_cont_p1+'.ms',spw=LB_contspws,spwmap=LB_spw_mapping,
    gaintable=[LB_p1_bis],interp='linearPD',calwt=True,applymode='calonly'
)

LB_cont_p1_bis = LB_cont_p1.replace('p1','p1_bis')
os.system('rm -rf %s.ms*' %LB_cont_p1_bis)
split(vis=LB_cont_p1+'.ms',outputvis=LB_cont_p1_bis+'.ms',datacolumn='corrected')

#CLEAN, again with a ~6sigma threshold
tclean_wrapper(
    vis       = LB_cont_p1_bis+'.ms',
    imagename = LB_cont_p1_bis,
    threshold = '0.0399mJy',
    **LB_tclean_wrapper_kwargs
)
estimate_SNR(LB_cont_p1_bis+'.image',disk_mask=LB_mask,noise_mask=noise_annulus_LB)
generate_image_png(
    LB_cont_p1_bis+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_LB,10*rms_LB],
    save_folder=LB_selfcal_folder
)
#MWC_758_SBLB_contp1_bis.image
#Beam 0.065 arcsec x 0.045 arcsec (7.20 deg)
#Flux inside disk mask: 66.28 mJy
#Peak intensity of source: 1.16 mJy/beam
#rms: 6.64e-03 mJy/beam
#Peak SNR: 174.72

#Second round of phase-only self-cal
LB_p2 = LB_p1_bis.replace('p1_bis','p2')
os.system('rm -rf '+LB_p2)
#Pietro: solint='360s' resulted in "mismatched frequencies" error, so I slighlty changed it
gaincal(
    vis=LB_cont_p1_bis+'.ms',caltable=LB_p2,
    gaintype='T',combine='scan,spw',calmode='p',solint='360s',
    spw=LB_contspws,refant=LB_refant,
    minsnr=3.,minblperant=4
)

#Inspect gain tables interactively and decide whether to manually flag something
#plotms(LB_p2,xaxis='time',yaxis='GainPhase',iteraxis='spw',showgui=True)
#Print calibration png file
#plotms(
#    LB_p2,xaxis='time', yaxis='GainPhase',iteraxis='spw',exprange='all',
#    overwrite=True,showgui=False,
#    plotfile=os.path.join(LB_selfcal_folder,f'{prefix}_LB_p2_gain_phase_vs_time.png')
#)

#Inspect gain tables and decide whether to flag something
plotms(
    LB_p2,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(LB_selfcal_folder,f'{prefix}_LB_p2_gain_phase_vs_time.png'),
)

#Flag problematic antennas
flagdata(vis=LB_p2,mode='manual',spw='0', timerange='2017/10/10/09:27:20~2017/10/10/09:27:30',antenna='DV20')

flagdata(vis=LB_p2,mode='manual',spw='4', timerange='2017/10/10/10:38:20~2017/10/10/10:38:30',antenna='DA54,DA60')
flagdata(vis=LB_p2,mode='manual',spw='4', timerange='2017/10/10/11:16:10~2017/10/10/11:16:20',antenna='DA59')
flagdata(vis=LB_p2,mode='manual',spw='4', timerange='2017/10/10/11:39:20~2017/10/10/11:39:30',antenna='DV24')

flagdata(vis=LB_p2,mode='manual',spw='8', timerange='2017/10/11/07:36:50~2017/10/11/07:37:00',antenna='DA44,DA64,DV04')
flagdata(vis=LB_p2,mode='manual',spw='8', timerange='2017/10/11/07:43:30~2017/10/11/07:43:40',antenna='DA63,DA64,DV04')
flagdata(vis=LB_p2,mode='manual',spw='8', timerange='2017/10/11/08:02:30~2017/10/11/08:02:40',antenna='DA45')
flagdata(vis=LB_p2,mode='manual',spw='8', timerange='2017/10/11/08:27:50~2017/10/11/08:28:00',antenna='DA52')
flagdata(vis=LB_p2,mode='manual',spw='8', timerange='2017/10/11/08:31:20~2017/10/11/08:31:30',antenna='DV08')

flagdata(vis=LB_p2,mode='manual',spw='12',timerange='2017/10/15/07:34:20~2017/10/15/07:34:30',antenna='DA58,DV12')
flagdata(vis=LB_p2,mode='manual',spw='12',timerange='2017/10/15/07:41:00~2017/10/15/07:41:10',antenna='DV22')
flagdata(vis=LB_p2,mode='manual',spw='12',timerange='2017/10/15/08:12:10~2017/10/15/08:12:20',antenna='DV05,DV22')

flagdata(vis=LB_p2,mode='manual',spw='16',timerange='2017/10/16/09:06:30~2017/10/16/09:06:40',antenna='DA56')
flagdata(vis=LB_p2,mode='manual',spw='16',timerange='2017/10/16/09:09:50~2017/10/16/09:10:00',antenna='PM04')
flagdata(vis=LB_p2,mode='manual',spw='16',timerange='2017/10/16/09:16:30~2017/10/16/09:16:40',antenna='DA59,DV23')
flagdata(vis=LB_p2,mode='manual',spw='16',timerange='2017/10/16/09:31:10~2017/10/16/09:31:20',antenna='DA51,DA55,DV23,PM04')

flagdata(vis=LB_p2,mode='manual',spw='20',timerange='2017/12/09/04:49:30~2017/12/09/04:50:00',antenna='DV15')
flagdata(vis=LB_p2,mode='manual',spw='20',timerange='2017/12/09/04:57:30~2017/12/09/04:58:00',antenna='DV11')
flagdata(vis=LB_p2,mode='manual',spw='20',timerange='2017/12/09/05:11:30~2017/12/09/05:12:00',antenna='DA46,DA55,DV10,DV15')
flagdata(vis=LB_p2,mode='manual',spw='20',timerange='2017/12/09/05:19:40~2017/12/09/05:20:00',antenna='DV10')
flagdata(vis=LB_p2,mode='manual',spw='20',timerange='2017/12/09/05:23:10~2017/12/09/05:23:20',antenna='DA46')

flagdata(vis=LB_p2,mode='manual',spw='24',timerange='2017/12/17/06:00:30~2017/12/17/06:01:00',antenna='DV07,DV13')
flagdata(vis=LB_p2,mode='manual',spw='24',timerange='2017/12/17/06:04:00~2017/12/17/06:04:20',antenna='DV13')

flagdata(vis=LB_p2,mode='manual',spw='28',timerange='2017/12/27/04:10:30~2017/12/27/04:10:40',antenna='DA42,DA44,DV13')
flagdata(vis=LB_p2,mode='manual',spw='28',timerange='2017/12/27/04:35:30~2017/12/27/04:35:40',antenna='DA44')
flagdata(vis=LB_p2,mode='manual',spw='28',timerange='2017/12/27/04:39:00~2017/12/27/04:39:10',antenna='DA44')

flagdata(vis=LB_p2,mode='manual',spw='32',timerange='2017/12/28/03:44:50~2017/12/28/03:45:00',antenna='DA44')
flagdata(vis=LB_p2,mode='manual',spw='32',timerange='2017/12/28/04:02:00~2017/12/28/04:02:20',antenna='DA44')
flagdata(vis=LB_p2,mode='manual',spw='32',timerange='2017/12/28/04:18:30~2017/12/28/04:19:00',antenna='DA55')
flagdata(vis=LB_p2,mode='manual',spw='32',timerange='2017/12/28/04:22:50~2017/12/28/04:23:00',antenna='DV11,DV13')
flagdata(vis=LB_p2,mode='manual',spw='32',timerange='2017/12/28/04:26:10~2017/12/28/04:26:20',antenna='DV11')
flagdata(vis=LB_p2,mode='manual',spw='32',timerange='2017/12/28/04:28:40~2017/12/28/04:29:00',antenna='DV13')

#Inspect gain tables to check if flagging worked
plotms(
    LB_p2,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(LB_selfcal_folder,f'{prefix}_LB_p2_gain_phase_vs_time_flagged.png'),
)

#Apply the calibration gains and split-off the corrected ms
applycal(
    vis=LB_cont_p1_bis+'.ms',spw=LB_contspws,spwmap=LB_spw_mapping,
    gaintable=[LB_p2],interp='linearPD',calwt=True,applymode='calonly'
)
LB_cont_p2 = LB_cont_p1_bis.replace('p1_bis','p2')
os.system('rm -rf %s.ms*' %LB_cont_p2)
split(vis=LB_cont_p1_bis+'.ms',outputvis=LB_cont_p2+'.ms',datacolumn='corrected')

#CLEAN, again with a ~6sigma threshold
tclean_wrapper(
    vis       = LB_cont_p2+'.ms',
    imagename = LB_cont_p2,
    threshold = '0.0396mJy',
    **LB_tclean_wrapper_kwargs
)
estimate_SNR(LB_cont_p2+'.image',disk_mask=LB_mask,noise_mask=noise_annulus_LB)
generate_image_png(
    LB_cont_p2+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_LB,10*rms_LB],
    save_folder=LB_selfcal_folder
)
#MWC_758_SBLB_contp2.image
#Beam 0.065 arcsec x 0.045 arcsec (7.20 deg)
#Flux inside disk mask: 66.47 mJy
#Peak intensity of source: 1.17 mJy/beam
#rms: 6.60e-03 mJy/beam
#Peak SNR: 177.94

#Third round of phase-only self-cal
LB_p3 = LB_p2.replace('p2','p3')
os.system('rm -rf '+LB_p3)
gaincal(
    vis=LB_cont_p2+'.ms',caltable=LB_p3,
    gaintype='T',combine='scan,spw',calmode='p',solint='120s', #diff from exoALMA, kept combine='scan' because scans for LB_EB0~45s
    spw=LB_contspws,refant=LB_refant,
    minsnr=3.,minblperant=4
)
#21 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:01:11.1
#23 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:03:33.2
#27 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:04:42.9
#19 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:06:15.8
#26 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:07:25.3
#19 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:08:50.5
#19 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:11:51.7
#27 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:13:02.5
#23 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:14:28.5
#24 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:15:46.0
#21 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:17:20.4
#28 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:18:30.2
#20 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:19:31.5
#22 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:22:13.6
#26 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:23:22.0
#19 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:24:54.4
#23 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:26:04.0
#19 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:27:36.5
#22 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:28:46.1
#23 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:29:44.7
#20 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:32:57.6
#20 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:35:19.0
#28 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:36:28.4
#20 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:38:00.9
#24 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:39:10.2
#21 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:40:35.4
#20 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:43:36.8
#23 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:44:46.4
#20 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:46:08.3
#20 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:47:30.5
#21 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:49:05.9
#26 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:50:14.0
#19 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:51:14.7
#21 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:53:54.3
#25 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:55:04.2
#20 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:56:36.4
#25 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:57:45.5
#21 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:59:17.6
#27 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:00:26.7
#21 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:01:25.1
#21 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:04:37.0
#21 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:06:58.0
#23 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:08:07.1
#19 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:09:38.8
#22 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:10:48.1
#18 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:12:12.7
#20 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:14:29.3
#22 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:27:57.6
#27 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:28:15.7
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:29:43.9
#24 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:30:52.2
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:32:23.4
#26 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:33:32.7
#16 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:35:04.5
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:36:41.8
#18 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:39:16.9
#26 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:40:26.4
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:41:49.1
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:43:18.2
#16 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:45:00.9
#22 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:46:09.8
#17 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:47:10.5
#18 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:49:49.0
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:50:58.3
#18 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:52:29.6
#24 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:53:38.8
#17 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:55:10.4
#25 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:56:19.3
#23 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:57:23.4
#16 of 43 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:00:29.6
#18 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:02:50.1
#26 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:03:58.9
#16 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:05:30.4
#21 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:06:39.4
#17 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:08:04.2
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:11:05.2
#22 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:12:14.1
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:13:32.3
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:14:59.2
#18 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:16:34.8
#21 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:17:43.7
#15 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:18:44.3
#17 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:21:25.1
#23 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:22:34.8
#14 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:24:06.2
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:25:15.3
#17 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:26:46.8
#23 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:27:55.8
#18 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:28:54.0
#18 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:32:08.3
#18 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:34:28.7
#23 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:35:37.6
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:37:09.7
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:38:18.1
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:39:42.4
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:41:49.3
#29 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:29:39.7
#23 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:29:57.8
#19 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:31:38.9
#30 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:32:49.3
#21 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:34:22.6
#26 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:35:32.7
#19 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:37:05.8
#22 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:38:32.7
#22 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:41:32.1
#27 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:42:13.2
#22 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:43:14.6
#19 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:45:51.6
#20 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:48:41.5
#30 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:49:24.0
#21 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:52:14.8
#29 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:52:56.4
#22 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:54:55.8
#22 of 42 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:55:39.6
#20 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:57:42.1
#25 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:58:22.8
#22 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:59:22.3
#18 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:02:39.4
#22 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:05:30.8
#23 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:06:11.8
#24 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:08:13.5
#25 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:08:55.0
#18 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:10:09.5
#19 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:13:03.7
#26 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:14:13.8
#19 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:15:35.3
#23 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:16:58.1
#19 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:18:34.2
#29 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:19:44.2
#20 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:20:45.5
#19 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:23:28.4
#28 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:24:38.4
#16 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:26:11.5
#29 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:27:21.5
#19 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:28:51.8
#24 of 41 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:30:04.5
#18 of 38 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:31:03.5
#21 of 36 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:34:18.3
#19 of 36 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:36:42.1
#27 of 36 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:37:49.6
#21 of 36 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:39:21.1
#30 of 35 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:40:32.4
#20 of 35 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:41:46.8
#21 of 35 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:44:22.7
#22 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:27:25.4
#31 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:27:43.6
#21 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:29:40.1
#26 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:30:21.7
#20 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:32:22.7
#22 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:33:04.9
#20 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:35:06.5
#18 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:36:05.4
#19 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:39:05.3
#23 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:39:46.2
#21 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:40:47.8
#18 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:43:26.2
#18 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:46:15.0
#29 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:46:58.8
#18 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:49:47.9
#23 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:50:29.9
#16 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:52:30.6
#27 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:53:13.1
#21 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:55:14.2
#34 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:55:56.3
#25 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:56:55.7
#17 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:00:09.8
#22 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:02:30.1
#36 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:03:42.9
#22 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:05:16.8
#28 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:06:26.0
#20 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:07:39.6
#20 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:10:33.7
#26 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:11:43.8
#24 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:12:45.1
#21 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:15:21.3
#18 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:17:42.5
#28 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:18:53.5
#21 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:21:18.4
#32 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:22:26.6
#19 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:23:59.5
#30 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:25:09.5
#20 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:26:42.5
#25 of 44 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:27:52.3
#24 of 44 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:28:51.3
#18 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:32:09.0
#18 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:34:29.9
#25 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:35:39.8
#20 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:37:12.5
#35 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:38:22.4
#21 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:39:36.8
#21 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:41:57.5
#21 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:43:03.0
#16 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:44:46.4
#24 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:45:56.4
#15 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:47:29.3
#29 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:48:38.6
#16 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:50:11.4
#22 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:51:39.5
#18 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:54:13.8
#25 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:55:23.6
#17 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:56:44.4
#29 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:58:09.0
#15 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:59:44.1
#26 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:00:54.1
#22 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:01:55.0
#16 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:04:42.2
#26 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:05:51.8
#20 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:07:23.8
#26 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:08:33.7
#19 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:10:06.0
#24 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:11:15.5
#24 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:12:22.3
#17 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:15:37.2
#14 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:17:58.6
#30 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:19:07.7
#18 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:20:39.9
#25 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:21:49.4
#19 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:23:03.2
#18 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:26:01.4
#33 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:27:11.2
#21 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:28:31.1
#26 of 43 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:29:56.2
#17 of 42 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:31:31.8
#29 of 41 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:32:40.5
#19 of 40 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:33:41.2
#20 of 40 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:36:31.1
#26 of 39 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:37:37.6
#16 of 40 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:39:09.8
#24 of 39 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:40:18.9
#15 of 37 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:41:51.0
#27 of 37 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:43:00.2
#21 of 40 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:43:59.2
#18 of 40 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:47:20.9
#17 of 39 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:49:39.8
#30 of 38 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:50:50.8
#16 of 39 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:52:23.7
#27 of 37 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:53:31.7
#18 of 38 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:54:45.5
#16 of 37 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:57:36.2
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:40:07.6
# 1 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:43:45.3
# 1 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:45:08.6
# 1 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:49:49.4
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:51:51.1
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:53:08.3
# 1 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:55:17.4
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:58:55.3
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:00:56.3
# 1 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:02:58.0
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:04:14.3
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:09:31.5
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:11:32.4
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:13:34.1
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:14:50.5
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:19:53.9
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:21:55.7
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:23:12.0
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:25:41.7
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:28:19.5
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:30:05.7
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:49:53.8
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:51:57.5
# 3 of 43 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:54:57.1
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:11:47.4
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/03:56:58.6
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:10:33.4
# 8 of 45 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:22:11.4
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:31:22.8
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:39:01.8
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:40:44.5
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/03:57:48.9
# 3 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/04:09:28.0

#Inspect gain tables interactively and decide whether to manually flag something
#plotms(LB_p3,xaxis='time',yaxis='GainPhase',iteraxis='spw',showgui=True)
#Print calibration png file
#plotms(
#    LB_p3,xaxis='time', yaxis='GainPhase',iteraxis='spw',exprange='all',
#    overwrite=True,showgui=False,
#    plotfile=os.path.join(LB_selfcal_folder,f'{prefix}_LB_p3_gain_phase_vs_time.png')
#)

#Inspect gain tables and decide whether to flag something
plotms(
    LB_p3,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(LB_selfcal_folder,f'{prefix}_LB_p3_gain_phase_vs_time.png'),
)

#Flag problematic antennas
flagdata(vis=LB_p3,mode='manual',spw='0', timerange='2017/10/10/09:06:10~2017/10/10/09:06:20',antenna='DA60,DV03')
flagdata(vis=LB_p3,mode='manual',spw='0', timerange='2017/10/10/09:07:20~2017/10/10/09:07:30',antenna='DV05')
flagdata(vis=LB_p3,mode='manual',spw='0', timerange='2017/10/10/09:13:00~2017/10/10/09:13:10',antenna='DV16')
flagdata(vis=LB_p3,mode='manual',spw='0', timerange='2017/10/10/09:18:25~2017/10/10/09:18:35',antenna='DA41,DA42,DV20')
flagdata(vis=LB_p3,mode='manual',spw='0', timerange='2017/10/10/09:23:20~2017/10/10/09:23:30',antenna='DA47')
flagdata(vis=LB_p3,mode='manual',spw='0', timerange='2017/10/10/09:26:00~2017/10/10/09:26:10',antenna='DA60,PM03')
flagdata(vis=LB_p3,mode='manual',spw='0', timerange='2017/10/10/09:27:30~2017/10/10/09:27:40',antenna='DA60,DV15')
flagdata(vis=LB_p3,mode='manual',spw='0', timerange='2017/10/10/09:28:40~2017/10/10/09:28:50',antenna='DA45')
flagdata(vis=LB_p3,mode='manual',spw='0', timerange='2017/10/10/09:32:50~2017/10/10/09:33:00',antenna='DA47,DA52')
flagdata(vis=LB_p3,mode='manual',spw='0', timerange='2017/10/10/09:36:20~2017/10/10/09:36:30',antenna='PM03')
flagdata(vis=LB_p3,mode='manual',spw='0', timerange='2017/10/10/09:47:30~2017/10/10/09:47:40',antenna='DA65,DV15,DV22')
flagdata(vis=LB_p3,mode='manual',spw='0', timerange='2017/10/10/09:50:10~2017/10/10/09:50:20',antenna='DV03,DA42,DV13')
flagdata(vis=LB_p3,mode='manual',spw='0', timerange='2017/10/10/09:51:10~2017/10/10/09:51:20',antenna='DA65,DV03,DV22')
flagdata(vis=LB_p3,mode='manual',spw='0', timerange='2017/10/10/10:01:20~2017/10/10/10:01:30',antenna='DA54,DV22')
flagdata(vis=LB_p3,mode='manual',spw='0', timerange='2017/10/10/10:10:40~2017/10/10/10:10:50',antenna='DA58,DA64')

flagdata(vis=LB_p3,mode='manual',spw='4', timerange='2017/10/10/10:28:10~2017/10/10/10:28:20',antenna='DA45,DA60')
flagdata(vis=LB_p3,mode='manual',spw='4', timerange='2017/10/10/10:30:50~2017/10/10/10:31:00',antenna='DA60')
flagdata(vis=LB_p3,mode='manual',spw='4', timerange='2017/10/10/10:33:30~2017/10/10/10:33:40',antenna='DA52')
flagdata(vis=LB_p3,mode='manual',spw='4', timerange='2017/10/10/10:35:00~2017/10/10/10:35:10',antenna='DA54,DV24')
flagdata(vis=LB_p3,mode='manual',spw='4', timerange='2017/10/10/10:36:40~2017/10/10/10:36:50',antenna='DA54,DA60')
flagdata(vis=LB_p3,mode='manual',spw='4', timerange='2017/10/10/10:40:20~2017/10/10/10:40:30',antenna='DV03')
flagdata(vis=LB_p3,mode='manual',spw='4', timerange='2017/10/10/10:46:00~2017/10/10/10:46:10',antenna='DA48,DV08')
flagdata(vis=LB_p3,mode='manual',spw='4', timerange='2017/10/10/10:52:20~2017/10/10/10:52:30',antenna='DA60,DV15')
flagdata(vis=LB_p3,mode='manual',spw='4', timerange='2017/10/10/10:53:30~2017/10/10/10:53:40',antenna='DA44')
flagdata(vis=LB_p3,mode='manual',spw='4', timerange='2017/10/10/11:06:30~2017/10/10/11:06:40',antenna='DV15')
flagdata(vis=LB_p3,mode='manual',spw='4', timerange='2017/10/10/11:12:10~2017/10/10/11:12:20',antenna='DV01')
flagdata(vis=LB_p3,mode='manual',spw='4', timerange='2017/10/10/11:24:00~2017/10/10/11:24:10',antenna='DA65,DV05')
flagdata(vis=LB_p3,mode='manual',spw='4', timerange='2017/10/10/11:25:10~2017/10/10/11:25:20',antenna='DV12')
flagdata(vis=LB_p3,mode='manual',spw='4', timerange='2017/10/10/11:28:50~2017/10/10/11:29:00',antenna='DA43,DA64,PM04')
flagdata(vis=LB_p3,mode='manual',spw='4', timerange='2017/10/10/11:32:00~2017/10/10/11:32:10',antenna='DA65,DV13')
flagdata(vis=LB_p3,mode='manual',spw='4', timerange='2017/10/10/11:35:30~2017/10/10/11:35:40',antenna='DA45,DV05')
flagdata(vis=LB_p3,mode='manual',spw='4', timerange='2017/10/10/11:38:10~2017/10/10/11:38:20',antenna='DV22')

flagdata(vis=LB_p3,mode='manual',spw='8', timerange='2017/10/11/07:31:30~2017/10/11/07:31:40',antenna='DV20')
flagdata(vis=LB_p3,mode='manual',spw='8', timerange='2017/10/11/07:34:20~2017/10/11/07:34:30',antenna='DV13')
flagdata(vis=LB_p3,mode='manual',spw='8', timerange='2017/10/11/07:35:30~2017/10/11/07:35:40',antenna='DA44')
flagdata(vis=LB_p3,mode='manual',spw='8', timerange='2017/10/11/07:37:00~2017/10/11/07:37:10',antenna='DA44,DA56,DV15,PM04')
flagdata(vis=LB_p3,mode='manual',spw='8', timerange='2017/10/11/07:38:30~2017/10/11/07:38:40',antenna='DA65')
flagdata(vis=LB_p3,mode='manual',spw='8', timerange='2017/10/11/07:43:10~2017/10/11/07:43:20',antenna='DV03')
flagdata(vis=LB_p3,mode='manual',spw='8', timerange='2017/10/11/07:58:20~2017/10/11/07:58:30',antenna='DV16')
flagdata(vis=LB_p3,mode='manual',spw='8', timerange='2017/10/11/08:19:40~2017/10/11/08:19:50',antenna='DA58,DA59,DV17')
flagdata(vis=LB_p3,mode='manual',spw='8', timerange='2017/10/11/08:23:20~2017/10/11/08:23:30',antenna='DA45')
flagdata(vis=LB_p3,mode='manual',spw='8', timerange='2017/10/11/08:24:30~2017/10/11/08:24:40',antenna='DV22')
flagdata(vis=LB_p3,mode='manual',spw='8', timerange='2017/10/11/08:27:20~2017/10/11/08:27:30',antenna='DV16')
flagdata(vis=LB_p3,mode='manual',spw='8', timerange='2017/10/11/08:31:00~2017/10/11/08:31:10',antenna='DV03')
flagdata(vis=LB_p3,mode='manual',spw='8', timerange='2017/10/11/08:36:40~2017/10/11/08:36:50',antenna='DA47')
flagdata(vis=LB_p3,mode='manual',spw='8', timerange='2017/10/11/08:37:40~2017/10/11/08:37:50',antenna='DA58,DV08')

flagdata(vis=LB_p3,mode='manual',spw='12',timerange='2017/10/15/07:27:20~2017/10/15/07:27:30',antenna='PM04')
flagdata(vis=LB_p3,mode='manual',spw='12',timerange='2017/10/15/07:39:00~2017/10/15/07:39:10',antenna='DV03')
flagdata(vis=LB_p3,mode='manual',spw='12',timerange='2017/10/15/07:39:40~2017/10/15/07:39:50',antenna='DA54,DA57')
flagdata(vis=LB_p3,mode='manual',spw='12',timerange='2017/10/15/07:46:45~2017/10/15/07:47:05',antenna='DA57')
flagdata(vis=LB_p3,mode='manual',spw='12',timerange='2017/10/15/07:49:40~2017/10/15/07:49:50',antenna='DA52')
flagdata(vis=LB_p3,mode='manual',spw='12',timerange='2017/10/15/08:03:40~2017/10/15/08:03:50',antenna='DV06')
flagdata(vis=LB_p3,mode='manual',spw='12',timerange='2017/10/15/08:06:20~2017/10/15/08:06:40',antenna='DA63')
flagdata(vis=LB_p3,mode='manual',spw='12',timerange='2017/10/15/08:10:30~2017/10/15/08:10:40',antenna='DV05')
flagdata(vis=LB_p3,mode='manual',spw='12',timerange='2017/10/15/08:18:50~2017/10/15/08:19:00',antenna='DA43,DV25,DV14')
flagdata(vis=LB_p3,mode='manual',spw='12',timerange='2017/10/15/08:21:10~2017/10/15/08:21:20',antenna='DV15')
flagdata(vis=LB_p3,mode='manual',spw='12',timerange='2017/10/15/08:27:50~2017/10/15/08:28:00',antenna='DA52')

flagdata(vis=LB_p3,mode='manual',spw='16',timerange='2017/10/16/08:45:50~2017/10/16/08:46:00',antenna='DA60')
flagdata(vis=LB_p3,mode='manual',spw='16',timerange='2017/10/16/08:48:30~2017/10/16/08:48:40',antenna='DA42')
flagdata(vis=LB_p3,mode='manual',spw='16',timerange='2017/10/16/08:50:10~2017/10/16/08:50:20',antenna='DA45,DA55,DA62')
flagdata(vis=LB_p3,mode='manual',spw='16',timerange='2017/10/16/08:55:20~2017/10/16/08:55:30',antenna='DA64,DV04')
flagdata(vis=LB_p3,mode='manual',spw='16',timerange='2017/10/16/09:00:50~2017/10/16/09:01:00',antenna='DA44,DV08')
flagdata(vis=LB_p3,mode='manual',spw='16',timerange='2017/10/16/09:05:50~2017/10/16/09:06:00',antenna='DA64,DV15')
flagdata(vis=LB_p3,mode='manual',spw='16',timerange='2017/10/16/09:17:50~2017/10/16/09:18:00',antenna='DA45,DA63')
flagdata(vis=LB_p3,mode='manual',spw='16',timerange='2017/10/16/09:21:40~2017/10/16/09:21:50',antenna='DA60')
flagdata(vis=LB_p3,mode='manual',spw='16',timerange='2017/10/16/09:23:00~2017/10/16/09:23:10',antenna='DV10,DV20')
flagdata(vis=LB_p3,mode='manual',spw='16',timerange='2017/10/16/09:27:10~2017/10/16/09:27:20',antenna='DV15')
flagdata(vis=LB_p3,mode='manual',spw='16',timerange='2017/10/16/09:31:30~2017/10/16/09:31:40',antenna='DA51,DA55,DV07')
flagdata(vis=LB_p3,mode='manual',spw='16',timerange='2017/10/16/09:33:40~2017/10/16/09:33:50',antenna='DA42,DV01')
flagdata(vis=LB_p3,mode='manual',spw='16',timerange='2017/10/16/09:35:00~2017/10/16/10:05:00',antenna='')

flagdata(vis=LB_p3,mode='manual',spw='20',timerange='2017/12/09/04:40:00~2017/12/09/04:40:10',antenna='DV15')
flagdata(vis=LB_p3,mode='manual',spw='20',timerange='2017/12/09/04:51:40~2017/12/09/04:52:00',antenna='DV10,DV11,DV15')
flagdata(vis=LB_p3,mode='manual',spw='20',timerange='2017/12/09/04:55:00~2017/12/09/04:55:20',antenna='DV10,DV15')
flagdata(vis=LB_p3,mode='manual',spw='20',timerange='2017/12/09/04:58:30~2017/12/09/04:59:00',antenna='DV11,DV15')
flagdata(vis=LB_p3,mode='manual',spw='20',timerange='2017/12/09/05:00:30~2017/12/09/05:01:00',antenna='DV10,DV11,DV15')
flagdata(vis=LB_p3,mode='manual',spw='20',timerange='2017/12/09/05:04:00~2017/12/09/05:04:30',antenna='DA49,DA47')
flagdata(vis=LB_p3,mode='manual',spw='20',timerange='2017/12/09/05:11:30~2017/12/09/05:11:40',antenna='DA46,DV15')
flagdata(vis=LB_p3,mode='manual',spw='20',timerange='2017/12/09/05:13:30~2017/12/09/05:13:40',antenna='DA55,DV10')
flagdata(vis=LB_p3,mode='manual',spw='20',timerange='2017/12/09/05:17:40~2017/12/09/05:18:00',antenna='DV10')
flagdata(vis=LB_p3,mode='manual',spw='20',timerange='2017/12/09/05:23:00~2017/12/09/05:23:20',antenna='DA46,DV10')
flagdata(vis=LB_p3,mode='manual',spw='20',timerange='2017/12/09/05:25:40~2017/12/09/05:26:00',antenna='DA46,DA51,DA55,DV15')

flagdata(vis=LB_p3,mode='manual',spw='24',timerange='2017/12/17/05:53:30~2017/12/17/05:54:00',antenna='DA43,DA44')
flagdata(vis=LB_p3,mode='manual',spw='24',timerange='2017/12/17/05:54:30~2017/12/17/05:55:00',antenna='DA41,DA42,DA43,DA44,DV11,DV13')
flagdata(vis=LB_p3,mode='manual',spw='24',timerange='2017/12/17/05:58:40~2017/12/17/05:59:00',antenna='DA55,DV13,DV07')
flagdata(vis=LB_p3,mode='manual',spw='24',timerange='2017/12/17/06:00:40~2017/12/17/06:01:00',antenna='DV13')
flagdata(vis=LB_p3,mode='manual',spw='24',timerange='2017/12/17/06:02:40~2017/12/17/06:03:00',antenna='DV13,DV07')
flagdata(vis=LB_p3,mode='manual',spw='24',timerange='2017/12/17/06:04:00~2017/12/17/06:04:20',antenna='DV13')

flagdata(vis=LB_p3,mode='manual',spw='28',timerange='2017/12/27/04:09:00~2017/12/27/04:09:30',antenna='DA42,DV13')
flagdata(vis=LB_p3,mode='manual',spw='28',timerange='2017/12/27/04:10:30~2017/12/27/04:11:00',antenna='DA42,DA44,DV13')
flagdata(vis=LB_p3,mode='manual',spw='28',timerange='2017/12/27/04:33:30~2017/12/27/04:34:00',antenna='DA44,DV11')
flagdata(vis=LB_p3,mode='manual',spw='28',timerange='2017/12/27/04:35:30~2017/12/27/04:36:00',antenna='DA44,DV11')
flagdata(vis=LB_p3,mode='manual',spw='28',timerange='2017/12/27/04:37:30~2017/12/27/04:38:00',antenna='DA44,DV11')
flagdata(vis=LB_p3,mode='manual',spw='28',timerange='2017/12/27/04:38:50~2017/12/27/04:39:10',antenna='DA44,DV11')
flagdata(vis=LB_p3,mode='manual',spw='28',timerange='2017/12/27/04:40:40~2017/12/27/04:40:50',antenna='DV11')

flagdata(vis=LB_p3,mode='manual',spw='32',timerange='2017/12/28/03:56:30~2017/12/28/03:56:40',antenna='DA55')
flagdata(vis=LB_p3,mode='manual',spw='32',timerange='2017/12/28/03:57:40~2017/12/28/03:57:50',antenna='DA55')
flagdata(vis=LB_p3,mode='manual',spw='32',timerange='2017/12/28/03:59:30~2017/12/28/03:59:40',antenna='DA55')
flagdata(vis=LB_p3,mode='manual',spw='32',timerange='2017/12/28/04:04:20~2017/12/28/04:04:30',antenna='DA44,DA46,DA55')
flagdata(vis=LB_p3,mode='manual',spw='32',timerange='2017/12/28/04:09:20~2017/12/28/04:09:30',antenna='DA42,DV11,DV13')
flagdata(vis=LB_p3,mode='manual',spw='32',timerange='2017/12/28/04:18:30~2017/12/28/04:18:40',antenna='DA55')
flagdata(vis=LB_p3,mode='manual',spw='32',timerange='2017/12/28/04:20:50~2017/12/28/04:21:00',antenna='DA44,DV11,DV13')
flagdata(vis=LB_p3,mode='manual',spw='32',timerange='2017/12/28/04:30:10~2017/12/28/04:30:20',antenna='DV13')

#Inspect gain tables to check if flagging worked
plotms(
    LB_p3,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(LB_selfcal_folder,f'{prefix}_LB_p3_gain_phase_vs_time_flagged.png'),
)

#Apply the calibration gains and split-off the corrected ms
applycal(
    vis=LB_cont_p2+'.ms',spw=LB_contspws,spwmap=LB_spw_mapping,
    gaintable=[LB_p3],interp='linearPD',calwt=True,applymode='calonly'
)
LB_cont_p3 = LB_cont_p2.replace('p2','p3')
os.system('rm -rf %s.ms*' %LB_cont_p3)
split(vis=LB_cont_p2+'.ms',outputvis=LB_cont_p3+'.ms',datacolumn='corrected')

#CLEAN, again with a ~6sigma threshold
tclean_wrapper(
    vis       = LB_cont_p3+'.ms',
    imagename = LB_cont_p3,
    threshold = '0.0395mJy',
    **LB_tclean_wrapper_kwargs
)
estimate_SNR(LB_cont_p3+'.image',disk_mask=LB_mask,noise_mask=noise_annulus_LB)
generate_image_png(
    LB_cont_p3+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_LB,10*rms_LB],
    save_folder=LB_selfcal_folder
)
#MWC_758_SBLB_contp3.image
#Beam 0.065 arcsec x 0.045 arcsec (7.20 deg)
#Flux inside disk mask: 66.51 mJy
#Peak intensity of source: 1.18 mJy/beam
#rms: 6.56e-03 mJy/beam
#Peak SNR: 179.16

#Fourth round of phase-only self-cal
LB_p4 = LB_p3.replace('p3','p4')
os.system('rm -rf '+LB_p4)
gaincal(
    vis=LB_cont_p3+'.ms',caltable=LB_p4,
    gaintype='T',combine='scan,spw',calmode='p',solint='60s', #diff from exoALMA, kept combine='scan' because scans for LB_EB0~30s
    spw=LB_contspws,refant=LB_refant,
    minsnr=3.,minblperant=4
)
# 22 of 45 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/09:00:28.1
# 21 of 45 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/09:01:40.1
# 22 of 46 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/09:03:01.7
# 23 of 46 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/09:04:22.7
# 21 of 46 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/09:05:44.0
# 21 of 46 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/09:07:05.1
# 18 of 46 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/09:08:26.5
# 23 of 46 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/09:09:26.3
# 19 of 46 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/09:11:21.2
# 21 of 45 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/09:12:42.4
# 22 of 45 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/09:13:54.2
# 24 of 46 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/09:15:36.9
# 22 of 46 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/09:16:48.9
# 21 of 46 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/09:18:10.2
# 22 of 45 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/09:19:31.5
# 21 of 45 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/09:21:40.8
# 22 of 46 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/09:23:02.0
# 22 of 46 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/09:24:22.9
# 21 of 46 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/09:25:43.9
# 23 of 46 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/09:27:04.9
# 22 of 46 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/09:28:25.9
# 23 of 46 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/09:29:44.6
# 22 of 46 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/09:32:14.5
# 21 of 46 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/09:33:26.9
# 21 of 46 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/09:34:47.3
# 20 of 46 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/09:36:08.3
# 19 of 46 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/09:37:29.3
# 18 of 46 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/09:38:50.4
# 22 of 46 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/09:40:11.0
# 23 of 46 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/09:41:11.4
# 20 of 46 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/09:43:05.4
# 20 of 46 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/09:44:26.4
# 22 of 46 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/09:45:38.4
# 21 of 45 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/09:47:21.8
# 21 of 45 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/09:48:33.1
# 22 of 46 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/09:49:54.0
# 22 of 46 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/09:51:14.8
# 22 of 46 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/09:53:23.0
# 22 of 46 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/09:54:44.6
# 22 of 46 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/09:56:05.2
# 20 of 46 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/09:57:25.8
# 22 of 46 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/09:58:46.4
# 20 of 46 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/10:00:07.1
# 21 of 46 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/10:01:25.1
# 21 of 46 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/10:03:54.3
# 22 of 46 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/10:05:06.1
# 19 of 46 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/10:06:26.7
# 20 of 46 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/10:07:47.2
# 20 of 45 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/10:09:07.6
# 19 of 45 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/10:10:28.1
# 18 of 45 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/10:11:48.6
# 21 of 45 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/10:12:48.7
# 19 of 46 solutions flagged due to SNR < 3 in spw=0 at 2017/10/10/10:14:29.3
# 20 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/10:27:57.6
# 29 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/10:28:15.7
# 21 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/10:29:12.4
# 20 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/10:30:32.3
# 18 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/10:31:52.4
# 19 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/10:33:13.3
# 21 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/10:34:33.2
# 20 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/10:35:53.4
# 19 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/10:36:52.9
# 20 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/10:38:45.8
# 17 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/10:40:06.9
# 20 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/10:41:17.9
# 18 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/10:43:09.6
# 19 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/10:44:29.8
# 18 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/10:45:50.0
# 17 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/10:47:10.6
# 18 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/10:49:17.9
# 17 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/10:50:38.4
# 19 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/10:51:58.6
# 19 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/10:53:19.2
# 17 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/10:54:39.2
# 17 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/10:55:59.5
# 20 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/10:57:23.4
# 16 of 43 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/10:59:47.0
# 16 of 43 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/11:00:58.6
# 16 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/11:02:18.8
# 20 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/11:03:38.8
# 17 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/11:04:59.2
# 19 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/11:06:19.4
# 19 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/11:07:39.7
# 20 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/11:08:40.4
# 19 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/11:10:33.7
# 25 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/11:11:54.5
# 20 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/11:13:05.3
# 18 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/11:14:52.4
# 18 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/11:16:03.6
# 16 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/11:17:23.9
# 18 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/11:18:44.3
# 16 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/11:20:54.2
# 16 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/11:22:15.1
# 16 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/11:23:35.1
# 16 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/11:24:55.4
# 20 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/11:26:15.7
# 18 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/11:27:35.9
# 18 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/11:28:54.0
# 21 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/11:31:25.8
# 17 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/11:32:37.7
# 16 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/11:33:57.6
# 17 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/11:35:17.8
# 20 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/11:36:38.4
# 20 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/11:37:58.4
# 15 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/11:39:18.5
# 18 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/11:40:17.8
# 22 of 44 solutions flagged due to SNR < 3 in spw=4 at 2017/10/10/11:41:49.3
# 33 of 43 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/07:29:39.7
# 21 of 43 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/07:29:57.8
# 21 of 43 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/07:31:08.4
# 18 of 43 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/07:32:30.1
# 21 of 43 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/07:33:52.1
# 22 of 43 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/07:35:13.7
# 21 of 43 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/07:36:35.2
# 19 of 43 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/07:37:57.0
# 24 of 43 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/07:38:47.4
# 20 of 43 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/07:40:31.4
# 21 of 43 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/07:41:53.1
# 24 of 43 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/07:43:14.6
# 21 of 43 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/07:45:08.2
# 18 of 43 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/07:46:20.7
# 21 of 43 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/07:47:42.4
# 24 of 43 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/07:49:04.0
# 21 of 43 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/07:51:14.6
# 21 of 43 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/07:52:36.2
# 16 of 43 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/07:53:57.9
# 19 of 42 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/07:55:19.5
# 19 of 42 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/07:56:41.1
# 20 of 43 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/07:58:02.7
# 21 of 43 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/07:59:22.3
# 21 of 43 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/08:01:55.9
# 22 of 43 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/08:03:08.7
# 18 of 43 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/08:04:30.0
# 22 of 43 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/08:05:51.7
# 20 of 43 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/08:07:13.3
# 21 of 43 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/08:08:34.8
# 19 of 43 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/08:09:56.7
# 26 of 43 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/08:10:47.5
# 22 of 43 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/08:12:32.1
# 20 of 43 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/08:13:53.7
# 17 of 43 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/08:15:06.1
# 18 of 43 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/08:16:50.1
# 18 of 43 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/08:18:02.6
# 20 of 43 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/08:19:24.0
# 18 of 43 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/08:20:45.5
# 21 of 43 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/08:22:56.8
# 19 of 43 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/08:24:18.3
# 20 of 43 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/08:25:39.9
# 24 of 43 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/08:27:01.3
# 21 of 43 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/08:28:22.8
# 19 of 41 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/08:29:44.3
# 19 of 38 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/08:31:03.5
# 21 of 34 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/08:33:34.3
# 22 of 35 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/08:34:46.7
# 21 of 35 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/08:36:08.1
# 20 of 36 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/08:37:29.4
# 21 of 36 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/08:38:50.8
# 23 of 35 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/08:40:12.2
# 22 of 35 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/08:41:33.8
# 27 of 35 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/08:42:25.2
# 20 of 35 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/08:44:07.6
# 23 of 34 solutions flagged due to SNR < 3 in spw=8 at 2017/10/11/08:45:10.9
# 24 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:27:27.5
# 20 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:28:40.0
# 20 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:30:01.6
# 16 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:31:23.2
# 20 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:32:44.8
# 18 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:34:06.4
# 19 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:35:28.0
# 22 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:36:18.8
# 18 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:38:04.6
# 19 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:39:26.0
# 23 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:40:47.8
# 24 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:42:42.9
# 19 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:43:55.4
# 20 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:45:17.2
# 18 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:46:38.7
# 20 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:48:48.1
# 19 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:50:09.8
# 18 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:51:31.3
# 17 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:52:52.9
# 17 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:54:14.5
# 22 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:55:36.1
# 24 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:56:55.8
# 23 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:59:27.3
# 19 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:00:39.7
# 22 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:02:01.1
# 21 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:03:22.5
# 20 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:04:44.1
# 23 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:06:05.8
# 21 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:07:27.2
# 31 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:08:18.0
# 22 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:10:02.2
# 20 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:11:23.7
# 24 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:12:45.1
# 27 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:14:38.2
# 23 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:15:50.5
# 18 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:17:11.9
# 24 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:18:33.4
# 21 of 44 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:20:45.0
# 24 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:22:06.5
# 23 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:23:27.9
# 22 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:24:49.3
# 22 of 44 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:26:10.7
# 23 of 44 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:27:32.1
# 21 of 44 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:28:51.3
# 25 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:31:24.8
# 22 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:32:37.0
# 21 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:33:58.3
# 21 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:35:19.6
# 25 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:36:40.9
# 26 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:38:02.3
# 24 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:39:23.9
# 31 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:40:14.9
# 18 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:41:57.5
# 21 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:43:03.0
# 19 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:44:15.1
# 17 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:45:36.2
# 18 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:46:57.7
# 19 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:48:18.4
# 18 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:49:39.6
# 17 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:51:00.7
# 22 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:51:53.9
# 20 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:53:42.4
# 20 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:55:03.7
# 18 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:56:15.4
# 19 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:58:01.2
# 18 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:59:12.9
# 21 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:00:34.0
# 22 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:01:55.0
# 16 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:04:10.7
# 16 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:05:31.7
# 21 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:06:52.5
# 19 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:08:13.6
# 20 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:09:34.4
# 19 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:10:55.3
# 22 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:12:22.3
# 22 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:14:54.1
# 19 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:16:06.0
# 17 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:17:27.0
# 18 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:18:47.6
# 19 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:20:08.4
# 19 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:21:29.2
# 19 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:22:50.1
# 30 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:23:42.3
# 15 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:25:30.3
# 23 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:26:51.4
# 22 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:28:03.1
# 21 of 43 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:29:48.2
# 17 of 41 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:30:59.7
# 19 of 41 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:32:20.4
# 20 of 40 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:33:41.2
# 18 of 38 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:35:56.8
# 16 of 39 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:37:17.4
# 17 of 39 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:38:38.2
# 18 of 39 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:39:58.8
# 17 of 37 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:41:19.5
# 19 of 37 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:42:40.1
# 26 of 40 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:43:59.2
# 23 of 40 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:46:38.4
# 18 of 40 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:47:49.9
# 20 of 39 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:49:10.3
# 17 of 38 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:50:30.7
# 18 of 37 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:51:51.6
# 16 of 37 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:53:11.7
# 19 of 38 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:54:32.1
# 26 of 38 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:55:25.6
# 17 of 37 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:57:12.1
# 21 of 34 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:58:23.8
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:39:46.5
# 1 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:40:34.9
# 1 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:41:35.3
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:42:35.8
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:43:36.3
# 1 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:44:37.7
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:45:20.7
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:47:18.8
# 1 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:48:19.3
# 1 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:49:19.8
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:50:20.3
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:51:20.7
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:52:22.8
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:53:08.3
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:55:02.5
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:55:48.0
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:58:25.0
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:00:25.9
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:01:26.4
# 1 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:02:26.9
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:03:28.9
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:04:14.3
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:09:01.2
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:10:01.6
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:11:02.1
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:12:02.6
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:13:03.1
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:14:05.1
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:14:50.5
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:17:22.7
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:18:23.2
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:19:23.7
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:21:24.6
# 1 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:22:26.6
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:23:12.0
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:25:41.7
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:27:49.8
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:28:50.3
# 1 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:29:50.8
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:30:36.3
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:34:52.4
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:50:24.3
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:51:25.3
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:52:27.4
# 3 of 43 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:54:57.1
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:11:47.4
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/03:56:58.6
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:02:47.9
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:07:45.5
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:10:33.4
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:12:47.6
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:17:39.2
# 8 of 45 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:22:11.4
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:31:22.8
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:35:13.5
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:36:14.0
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:37:14.3
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:38:16.3
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:39:01.8
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:40:29.5
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:41:14.9
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/03:57:48.9
# 3 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/04:09:28.0

#Inspect gain tables interactively and decide whether to manually flag something
#plotms(LB_p4,xaxis='time',yaxis='GainPhase',iteraxis='spw',showgui=True)
#Print calibration png file
#plotms(
#    LB_p4,xaxis='time', yaxis='GainPhase',iteraxis='spw',exprange='all',
#    overwrite=True,showgui=False,
#    plotfile=os.path.join(LB_selfcal_folder,f'{prefix}_LB_p4_gain_phase_vs_time.png')
#)

#Inspect gain tables and decide whether to flag something
plotms(
    LB_p4,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(LB_selfcal_folder,f'{prefix}_LB_p4_gain_phase_vs_time.png'),
)

#Flag problematic antennas
flagdata(vis=LB_p4,mode='manual',spw='0', timerange='2017/10/10/09:07:00~2017/10/10/09:07:10',antenna='DA48')
flagdata(vis=LB_p4,mode='manual',spw='0', timerange='2017/10/10/09:08:20~2017/10/10/09:08:30',antenna='DA58')
flagdata(vis=LB_p4,mode='manual',spw='0', timerange='2017/10/10/09:11:20~2017/10/10/09:11:30',antenna='DV16,DV17')
flagdata(vis=LB_p4,mode='manual',spw='0', timerange='2017/10/10/09:27:00~2017/10/10/09:27:10',antenna='DV15')
flagdata(vis=LB_p4,mode='manual',spw='0', timerange='2017/10/10/09:28:20~2017/10/10/09:28:30',antenna='DA60,DV15')
flagdata(vis=LB_p4,mode='manual',spw='0', timerange='2017/10/10/09:29:40~2017/10/10/09:29:50',antenna='DV16')
flagdata(vis=LB_p4,mode='manual',spw='0', timerange='2017/10/10/09:37:20~2017/10/10/09:37:30',antenna='DV15')
flagdata(vis=LB_p4,mode='manual',spw='0', timerange='2017/10/10/09:49:50~2017/10/10/09:50:00',antenna='DA42,DV13')
flagdata(vis=LB_p4,mode='manual',spw='0', timerange='2017/10/10/09:51:10~2017/10/10/09:51:20',antenna='DV22')
flagdata(vis=LB_p4,mode='manual',spw='0', timerange='2017/10/10/10:00:00~2017/10/10/10:00:10',antenna='DV20')
flagdata(vis=LB_p4,mode='manual',spw='0', timerange='2017/10/10/10:01:20~2017/10/10/10:01:30',antenna='DV22')
flagdata(vis=LB_p4,mode='manual',spw='0', timerange='2017/10/10/10:06:20~2017/10/10/10:06:30',antenna='DA63')
flagdata(vis=LB_p4,mode='manual',spw='0', timerange='2017/10/10/10:14:20~2017/10/10/10:14:30',antenna='DV05')

flagdata(vis=LB_p4,mode='manual',spw='4', timerange='2017/10/10/10:27:50~2017/10/10/10:28:00',antenna='DV22')
flagdata(vis=LB_p4,mode='manual',spw='4', timerange='2017/10/10/10:28:10~2017/10/10/10:28:20',antenna='DA60,PM04')
flagdata(vis=LB_p4,mode='manual',spw='4', timerange='2017/10/10/10:35:50~2017/10/10/10:36:00',antenna='DA60')
flagdata(vis=LB_p4,mode='manual',spw='4', timerange='2017/10/10/10:36:50~2017/10/10/10:37:00',antenna='DA60,DA65')
flagdata(vis=LB_p4,mode='manual',spw='4', timerange='2017/10/10/10:40:00~2017/10/10/10:40:10',antenna='DA53')
flagdata(vis=LB_p4,mode='manual',spw='4', timerange='2017/10/10/10:41:10~2017/10/10/10:41:20',antenna='DA64')
flagdata(vis=LB_p4,mode='manual',spw='4', timerange='2017/10/10/10:50:30~2017/10/10/10:50:40',antenna='DA54,DV12')
flagdata(vis=LB_p4,mode='manual',spw='4', timerange='2017/10/10/10:57:20~2017/10/10/10:57:30',antenna='DV20')
flagdata(vis=LB_p4,mode='manual',spw='4', timerange='2017/10/10/11:11:50~2017/10/10/11:12:00',antenna='DA41,DV01,DV13')
flagdata(vis=LB_p4,mode='manual',spw='4', timerange='2017/10/10/11:13:00~2017/10/10/11:13:10',antenna='DA53,DV01')
flagdata(vis=LB_p4,mode='manual',spw='4', timerange='2017/10/10/11:17:20~2017/10/10/11:17:30',antenna='DA48,DA59')
flagdata(vis=LB_p4,mode='manual',spw='4', timerange='2017/10/10/11:18:40~2017/10/10/11:18:50',antenna='DA58')
flagdata(vis=LB_p4,mode='manual',spw='4', timerange='2017/10/10/11:22:10~2017/10/10/11:22:20',antenna='DA52,DA60,DV20')
flagdata(vis=LB_p4,mode='manual',spw='4', timerange='2017/10/10/11:23:30~2017/10/10/11:23:40',antenna='DV12')
flagdata(vis=LB_p4,mode='manual',spw='4', timerange='2017/10/10/11:35:10~2017/10/10/11:35:20',antenna='DA47,DV05')
flagdata(vis=LB_p4,mode='manual',spw='4', timerange='2017/10/10/11:36:30~2017/10/10/11:36:40',antenna='DV22')
flagdata(vis=LB_p4,mode='manual',spw='4', timerange='2017/10/10/11:39:10~2017/10/10/11:39:20',antenna='DA48,DA54')
flagdata(vis=LB_p4,mode='manual',spw='4', timerange='2017/10/10/11:41:40~2017/10/10/11:41:50',antenna='DA44')

flagdata(vis=LB_p4,mode='manual',spw='8', timerange='2017/10/11/07:29:30~2017/10/11/07:29:40',antenna='DA65')
flagdata(vis=LB_p4,mode='manual',spw='8', timerange='2017/10/11/07:29:50~2017/10/11/07:30:00',antenna='DV08')
flagdata(vis=LB_p4,mode='manual',spw='8', timerange='2017/10/11/07:33:50~2017/10/11/07:34:00',antenna='DA53,DV13')
flagdata(vis=LB_p4,mode='manual',spw='8', timerange='2017/10/11/07:36:30~2017/10/11/07:36:40',antenna='DA44,DV15')
flagdata(vis=LB_p4,mode='manual',spw='8', timerange='2017/10/11/07:37:50~2017/10/11/07:38:00',antenna='DA44,DA58,DA65,DV15,PM04')
flagdata(vis=LB_p4,mode='manual',spw='8', timerange='2017/10/11/07:38:40~2017/10/11/07:38:50',antenna='DV08,DV20')
flagdata(vis=LB_p4,mode='manual',spw='8', timerange='2017/10/11/07:53:50~2017/10/11/07:54:00',antenna='DA45')
flagdata(vis=LB_p4,mode='manual',spw='8', timerange='2017/10/11/07:58:00~2017/10/11/07:58:10',antenna='DA56,PM03')
flagdata(vis=LB_p4,mode='manual',spw='8', timerange='2017/10/11/08:07:10~2017/10/11/08:07:20',antenna='DA42,DA44,DV15,PM04')
flagdata(vis=LB_p4,mode='manual',spw='8', timerange='2017/10/11/08:08:30~2017/10/11/08:08:40',antenna='DV22')
flagdata(vis=LB_p4,mode='manual',spw='8', timerange='2017/10/11/08:18:00~2017/10/11/08:18:10',antenna='DA54')
flagdata(vis=LB_p4,mode='manual',spw='8', timerange='2017/10/11/08:19:20~2017/10/11/08:19:30',antenna='DA58')
flagdata(vis=LB_p4,mode='manual',spw='8', timerange='2017/10/11/08:28:20~2017/10/11/08:28:30',antenna='DA43,DA46,DA63')
flagdata(vis=LB_p4,mode='manual',spw='8', timerange='2017/10/11/08:29:40~2017/10/11/08:29:50',antenna='DA44,DA53,DV13,DV15')
flagdata(vis=LB_p4,mode='manual',spw='8', timerange='2017/10/11/08:31:00~2017/10/11/08:31:10',antenna='DA56')
flagdata(vis=LB_p4,mode='manual',spw='8', timerange='2017/10/11/08:33:30~2017/10/11/08:33:40',antenna='DV23')
flagdata(vis=LB_p4,mode='manual',spw='8', timerange='2017/10/11/08:38:50~2017/10/11/08:39:00',antenna='DV20')
flagdata(vis=LB_p4,mode='manual',spw='8', timerange='2017/10/11/08:42:20~2017/10/11/08:42:30',antenna='DA47,PM01')
flagdata(vis=LB_p4,mode='manual',spw='8', timerange='2017/10/11/08:45:10~2017/10/11/08:45:20',antenna='DA45')

flagdata(vis=LB_p4,mode='manual',spw='12',timerange='2017/10/15/07:28:35~2017/10/15/07:28:45',antenna='DA43,PM03')
flagdata(vis=LB_p4,mode='manual',spw='12',timerange='2017/10/15/07:30:00~2017/10/15/07:30:10',antenna='DA52')
flagdata(vis=LB_p4,mode='manual',spw='12',timerange='2017/10/15/07:31:20~2017/10/15/07:31:30',antenna='DV22')
flagdata(vis=LB_p4,mode='manual',spw='12',timerange='2017/10/15/07:32:40~2017/10/15/07:32:50',antenna='DV15')
flagdata(vis=LB_p4,mode='manual',spw='12',timerange='2017/10/15/07:39:20~2017/10/15/07:39:30',antenna='DV10')
flagdata(vis=LB_p4,mode='manual',spw='12',timerange='2017/10/15/07:43:50~2017/10/15/07:44:00',antenna='DV15')
flagdata(vis=LB_p4,mode='manual',spw='12',timerange='2017/10/15/08:00:30~2017/10/15/08:00:40',antenna='DA63')
flagdata(vis=LB_p4,mode='manual',spw='12',timerange='2017/10/15/08:02:00~2017/10/15/08:02:10',antenna='DA53,DA64')
flagdata(vis=LB_p4,mode='manual',spw='12',timerange='2017/10/15/08:03:20~2017/10/15/08:03:30',antenna='DA44')
flagdata(vis=LB_p4,mode='manual',spw='12',timerange='2017/10/15/08:04:40~2017/10/15/08:04:50',antenna='DA53,DV13')
flagdata(vis=LB_p4,mode='manual',spw='12',timerange='2017/10/15/08:06:00~2017/10/15/08:06:10',antenna='DA42,DA46,DA63')
flagdata(vis=LB_p4,mode='manual',spw='12',timerange='2017/10/15/08:08:10~2017/10/15/08:08:20',antenna='DV04')
flagdata(vis=LB_p4,mode='manual',spw='12',timerange='2017/10/15/08:10:00~2017/10/15/08:10:10',antenna='DA65,DV15')
flagdata(vis=LB_p4,mode='manual',spw='12',timerange='2017/10/15/08:11:20~2017/10/15/08:11:30',antenna='DV05')
flagdata(vis=LB_p4,mode='manual',spw='12',timerange='2017/10/15/08:17:10~2017/10/15/08:17:20',antenna='DA64')
flagdata(vis=LB_p4,mode='manual',spw='12',timerange='2017/10/15/08:22:00~2017/10/15/08:22:10',antenna='DV12')
flagdata(vis=LB_p4,mode='manual',spw='12',timerange='2017/10/15/08:28:50~2017/10/15/08:29:00',antenna='DA52')
flagdata(vis=LB_p4,mode='manual',spw='12',timerange='2017/10/15/08:31:20~2017/10/15/08:31:30',antenna='DV24')

flagdata(vis=LB_p4,mode='manual',spw='16',timerange='2017/10/16/08:48:10~2017/10/16/08:48:20',antenna='DA42')
flagdata(vis=LB_p4,mode='manual',spw='16',timerange='2017/10/16/08:50:55~2017/10/16/08:51:05',antenna='DA45')
flagdata(vis=LB_p4,mode='manual',spw='16',timerange='2017/10/16/08:51:50~2017/10/16/08:52:00',antenna='DA45')
flagdata(vis=LB_p4,mode='manual',spw='16',timerange='2017/10/16/08:53:40~2017/10/16/08:53:50',antenna='DV24')
flagdata(vis=LB_p4,mode='manual',spw='16',timerange='2017/10/16/09:18:40~2017/10/16/09:18:50',antenna='DV08')
flagdata(vis=LB_p4,mode='manual',spw='16',timerange='2017/10/16/09:21:20~2017/10/16/09:21:30',antenna='DA53,DA64')
flagdata(vis=LB_p4,mode='manual',spw='16',timerange='2017/10/16/09:23:40~2017/10/16/09:23:50',antenna='DA53')
flagdata(vis=LB_p4,mode='manual',spw='16',timerange='2017/10/16/09:25:25~2017/10/16/09:25:35',antenna='DA51,DA59,DV20')
flagdata(vis=LB_p4,mode='manual',spw='16',timerange='2017/10/16/09:28:00~2017/10/16/09:28:10',antenna='DA64')
flagdata(vis=LB_p4,mode='manual',spw='16',timerange='2017/10/16/09:29:40~2017/10/16/09:29:50',antenna='DA52')
flagdata(vis=LB_p4,mode='manual',spw='16',timerange='2017/10/16/09:32:20~2017/10/16/09:32:30',antenna='DA45,DV10')
flagdata(vis=LB_p4,mode='manual',spw='16',timerange='2017/10/16/09:33:40~2017/10/16/09:33:50',antenna='DA42,DA45,DV01')
flagdata(vis=LB_p4,mode='manual',spw='16',timerange='2017/10/16/09:35:00~2017/10/16/10:05:00',antenna='')

flagdata(vis=LB_p4,mode='manual',spw='20',timerange='2017/12/09/04:49:15~2017/12/09/04:49:25',antenna='DV15')
flagdata(vis=LB_p4,mode='manual',spw='20',timerange='2017/12/09/04:50:15~2017/12/09/04:50:25',antenna='DV15')
flagdata(vis=LB_p4,mode='manual',spw='20',timerange='2017/12/09/04:51:15~2017/12/09/04:51:25',antenna='DV15')
flagdata(vis=LB_p4,mode='manual',spw='20',timerange='2017/12/09/04:52:15~2017/12/09/04:52:25',antenna='DV11')
flagdata(vis=LB_p4,mode='manual',spw='20',timerange='2017/12/09/04:55:00~2017/12/09/04:55:10',antenna='DV15')
flagdata(vis=LB_p4,mode='manual',spw='20',timerange='2017/12/09/04:55:45~2017/12/09/04:55:50',antenna='DV11')
flagdata(vis=LB_p4,mode='manual',spw='20',timerange='2017/12/09/04:59:20~2017/12/09/04:59:30',antenna='DV11,DV15')
flagdata(vis=LB_p4,mode='manual',spw='20',timerange='2017/12/09/05:00:20~2017/12/09/05:00:30',antenna='DV10')
flagdata(vis=LB_p4,mode='manual',spw='20',timerange='2017/12/09/05:01:25~2017/12/09/05:01:30',antenna='DV10')
flagdata(vis=LB_p4,mode='manual',spw='20',timerange='2017/12/09/05:02:20~2017/12/09/05:02:30',antenna='DV10')
flagdata(vis=LB_p4,mode='manual',spw='20',timerange='2017/12/09/05:04:10~2017/12/09/05:04:20',antenna='DA49,DA47')
flagdata(vis=LB_p4,mode='manual',spw='20',timerange='2017/12/09/05:10:00~2017/12/09/05:10:05',antenna='DV10')
flagdata(vis=LB_p4,mode='manual',spw='20',timerange='2017/12/09/05:11:00~2017/12/09/05:11:05',antenna='DV10')
flagdata(vis=LB_p4,mode='manual',spw='20',timerange='2017/12/09/05:12:00~2017/12/09/05:12:10',antenna='DA46')
flagdata(vis=LB_p4,mode='manual',spw='20',timerange='2017/12/09/05:14:00~2017/12/09/05:14:10',antenna='DA46,DA55,DV11')
flagdata(vis=LB_p4,mode='manual',spw='20',timerange='2017/12/09/05:17:20~2017/12/09/05:17:25',antenna='DV10')
flagdata(vis=LB_p4,mode='manual',spw='20',timerange='2017/12/09/05:20:20~2017/12/09/05:20:25',antenna='DA46,DV15')
flagdata(vis=LB_p4,mode='manual',spw='20',timerange='2017/12/09/05:21:20~2017/12/09/05:21:25',antenna='DA46')
flagdata(vis=LB_p4,mode='manual',spw='20',timerange='2017/12/09/05:22:25~2017/12/09/05:22:30',antenna='DA46,DA55')
flagdata(vis=LB_p4,mode='manual',spw='20',timerange='2017/12/09/05:23:10~2017/12/09/05:23:25',antenna='DA46')
flagdata(vis=LB_p4,mode='manual',spw='20',timerange='2017/12/09/05:25:40~2017/12/09/05:25:45',antenna='DA46,DA51,DA55')
flagdata(vis=LB_p4,mode='manual',spw='20',timerange='2017/12/09/05:27:45~2017/12/09/05:27:50',antenna='DA46,DA51,DV10')
flagdata(vis=LB_p4,mode='manual',spw='20',timerange='2017/12/09/05:30:35~2017/12/09/05:40:00',antenna='DA46')

flagdata(vis=LB_p4,mode='manual',spw='24',timerange='2017/12/17/05:49:20~2017/12/17/05:49:30',antenna='DA55')
flagdata(vis=LB_p4,mode='manual',spw='24',timerange='2017/12/17/05:54:20~2017/12/17/05:54:30',antenna='DA41,DA42,DA43,DA44,DA45@W205,DA55,DV13')
flagdata(vis=LB_p4,mode='manual',spw='24',timerange='2017/12/17/05:54:50~2017/12/17/05:55:00',antenna='DA41,DA42,DA43,DA44')
flagdata(vis=LB_p4,mode='manual',spw='24',timerange='2017/12/17/05:58:15~2017/12/17/05:58:20',antenna='DV13,DA44,DA55,DV07')
flagdata(vis=LB_p4,mode='manual',spw='24',timerange='2017/12/17/05:59:15~2017/12/17/05:59:20',antenna='DA55,DV13,DV07')
flagdata(vis=LB_p4,mode='manual',spw='24',timerange='2017/12/17/06:00:15~2017/12/17/06:00:20',antenna='DV13')
flagdata(vis=LB_p4,mode='manual',spw='24',timerange='2017/12/17/06:01:20~2017/12/17/06:01:25',antenna='DV13,DV07')
flagdata(vis=LB_p4,mode='manual',spw='24',timerange='2017/12/17/06:02:20~2017/12/17/06:02:25',antenna='DV13,DV07')
flagdata(vis=LB_p4,mode='manual',spw='24',timerange='2017/12/17/06:03:20~2017/12/17/06:03:25',antenna='DV13,DV07')
flagdata(vis=LB_p4,mode='manual',spw='24',timerange='2017/12/17/06:04:05~2017/12/17/06:04:10',antenna='DV13')

flagdata(vis=LB_p4,mode='manual',spw='28',timerange='2017/12/27/04:08:45~2017/12/27/04:08:50',antenna='DA42,DV13')
flagdata(vis=LB_p4,mode='manual',spw='28',timerange='2017/12/27/04:09:45~2017/12/27/04:09:50',antenna='DA42,DA44,DV11,DV13')
flagdata(vis=LB_p4,mode='manual',spw='28',timerange='2017/12/27/04:10:30~2017/12/27/04:10:35',antenna='DA42,DA44,DV13')
flagdata(vis=LB_p4,mode='manual',spw='28',timerange='2017/12/27/04:12:00~2017/12/27/04:12:10',antenna='DA44')
flagdata(vis=LB_p4,mode='manual',spw='28',timerange='2017/12/27/04:21:30~2017/12/27/04:21:40',antenna='DA55')
flagdata(vis=LB_p4,mode='manual',spw='28',timerange='2017/12/27/04:25:25~2017/12/27/04:25:35',antenna='DV11')
flagdata(vis=LB_p4,mode='manual',spw='28',timerange='2017/12/27/04:26:30~2017/12/27/04:26:40',antenna='DV11')
flagdata(vis=LB_p4,mode='manual',spw='28',timerange='2017/12/27/04:33:10~2017/12/27/04:33:20',antenna='DA44')
flagdata(vis=LB_p4,mode='manual',spw='28',timerange='2017/12/27/04:34:10~2017/12/27/04:34:20',antenna='DA44,DV11')
flagdata(vis=LB_p4,mode='manual',spw='28',timerange='2017/12/27/04:35:00~2017/12/27/04:35:20',antenna='DA45,DA44')
flagdata(vis=LB_p4,mode='manual',spw='28',timerange='2017/12/27/04:36:00~2017/12/27/04:36:20',antenna='DA44')
flagdata(vis=LB_p4,mode='manual',spw='28',timerange='2017/12/27/04:37:00~2017/12/27/04:37:20',antenna='DA44')
flagdata(vis=LB_p4,mode='manual',spw='28',timerange='2017/12/27/04:38:00~2017/12/27/04:38:20',antenna='DA44')
flagdata(vis=LB_p4,mode='manual',spw='28',timerange='2017/12/27/04:39:00~2017/12/27/04:39:20',antenna='DA44')

flagdata(vis=LB_p4,mode='manual',spw='32',timerange='2017/12/28/03:45:00~2017/12/28/03:45:10',antenna='DA55')
flagdata(vis=LB_p4,mode='manual',spw='32',timerange='2017/12/28/03:56:00~2017/12/28/03:56:10',antenna='DV11')
flagdata(vis=LB_p4,mode='manual',spw='32',timerange='2017/12/28/03:57:00~2017/12/28/03:57:10',antenna='DA55')
flagdata(vis=LB_p4,mode='manual',spw='32',timerange='2017/12/28/03:59:10~2017/12/28/03:59:20',antenna='DA55')
flagdata(vis=LB_p4,mode='manual',spw='32',timerange='2017/12/28/04:00:00~2017/12/28/04:00:10',antenna='DA55,DA44')
flagdata(vis=LB_p4,mode='manual',spw='32',timerange='2017/12/28/04:02:50~2017/12/28/04:03:00',antenna='DA44')
flagdata(vis=LB_p4,mode='manual',spw='32',timerange='2017/12/28/04:03:50~2017/12/28/04:04:00',antenna='DA44,DA55,DA46')
flagdata(vis=LB_p4,mode='manual',spw='32',timerange='2017/12/28/04:04:50~2017/12/28/04:05:00',antenna='DA44')
flagdata(vis=LB_p4,mode='manual',spw='32',timerange='2017/12/28/04:09:20~2017/12/28/04:09:30',antenna='DA42,DV13')
flagdata(vis=LB_p4,mode='manual',spw='32',timerange='2017/12/28/04:12:45~2017/12/28/04:12:55',antenna='DA42,DA55,DV11,DV13')
flagdata(vis=LB_p4,mode='manual',spw='32',timerange='2017/12/28/04:17:50~2017/12/28/04:18:00',antenna='DA55')
flagdata(vis=LB_p4,mode='manual',spw='32',timerange='2017/12/28/04:18:35~2017/12/28/04:18:45',antenna='DA55')
flagdata(vis=LB_p4,mode='manual',spw='32',timerange='2017/12/28/04:20:25~2017/12/28/04:20:35',antenna='DV11')
flagdata(vis=LB_p4,mode='manual',spw='32',timerange='2017/12/28/04:21:25~2017/12/28/04:21:35',antenna='DV11,DV13,DA44')
flagdata(vis=LB_p4,mode='manual',spw='32',timerange='2017/12/28/04:22:25~2017/12/28/04:22:35',antenna='DV11,DV13')
flagdata(vis=LB_p4,mode='manual',spw='32',timerange='2017/12/28/04:30:10~2017/12/28/04:30:30',antenna='DV13')

#Inspect gain tables to check if flagging worked
plotms(
    LB_p4,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(LB_selfcal_folder,f'{prefix}_LB_p4_gain_phase_vs_time_flagged.png'),
)

#Apply the calibration gains and split-off the corrected ms
applycal(
    vis=LB_cont_p3+'.ms',spw=LB_contspws,spwmap=LB_spw_mapping,
    gaintable=[LB_p4],interp='linearPD',calwt=True,applymode='calonly'
)
LB_cont_p4 = LB_cont_p3.replace('p3','p4')
os.system('rm -rf %s.ms*' %LB_cont_p4)
split(vis=LB_cont_p3+'.ms',outputvis=LB_cont_p4+'.ms',datacolumn='corrected')

#CLEAN, again with a ~6sigma threshold
tclean_wrapper(
    vis       = LB_cont_p4+'.ms',
    imagename = LB_cont_p4,
    threshold = '0.0393mJy',
    **LB_tclean_wrapper_kwargs
)
estimate_SNR(LB_cont_p4+'.image',disk_mask=LB_mask,noise_mask=noise_annulus_LB)
generate_image_png(
    LB_cont_p4+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_LB,10*rms_LB],
    save_folder=LB_selfcal_folder
)
#MWC_758_SBLB_contp4.image
#Beam 0.065 arcsec x 0.045 arcsec (7.20 deg)
#Flux inside disk mask: 66.22 mJy
#Peak intensity of source: 1.22 mJy/beam
#rms: 6.52e-03 mJy/beam
#Peak SNR: 186.57
"""
#Fifth round of phase-only self-cal
LB_p5 = LB_p4.replace('p4','p5')
os.system('rm -rf '+LB_p5)
gaincal(
    vis=LB_cont_p4+'.ms',caltable=LB_p5,
    gaintype='T',combine='spw',calmode='p',solint='30s',
    spw=LB_contspws,refant=LB_refant,
    minsnr=3.,minblperant=4
)
#22 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:00:25.1
#34 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:00:43.2
#23 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:01:28.2
#22 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:01:55.2
#23 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:02:49.6
#25 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:03:16.6
#23 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:04:10.6
#24 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:04:37.8
#22 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:05:32.0
#22 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:05:59.0
#19 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:06:53.0
#23 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:07:20.3
#22 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:08:14.3
#22 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:08:41.5
#23 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:09:23.2
#38 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:09:41.4
#20 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:11:09.1
#26 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:11:36.3
#20 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:12:30.3
#21 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:12:57.5
#21 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:13:51.4
#28 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:14:09.3
#22 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:15:33.9
#30 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:15:52.0
#21 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:16:36.9
#24 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:17:04.1
#20 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:17:58.0
#25 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:18:25.2
#22 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:19:19.4
#23 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:19:46.3
#24 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:21:28.7
#25 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:21:55.9
#22 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:22:49.9
#24 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:23:16.9
#24 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:24:10.8
#26 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:24:38.0
#23 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:25:31.8
#20 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:25:59.0
#24 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:26:52.8
#23 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:27:20.0
#21 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:28:13.8
#21 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:28:41.0
#21 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:29:28.8
#25 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:30:00.4
#24 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:32:11.5
#35 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:32:29.6
#24 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:33:14.8
#23 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:33:41.5
#22 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:34:35.2
#22 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:35:02.4
#22 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:35:56.2
#23 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:36:23.3
#20 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:37:17.2
#19 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:37:44.3
#21 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:38:38.4
#21 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:39:05.2
#19 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:39:58.9
#22 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:40:26.1
#22 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:41:08.4
#27 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:41:26.5
#23 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:42:53.3
#23 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:43:20.5
#21 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:44:14.3
#23 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:44:41.4
#21 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:45:35.3
#32 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:45:53.1
#23 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:47:18.8
#30 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:47:36.6
#21 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:48:21.0
#21 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:48:48.2
#23 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:49:41.9
#25 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:50:09.0
#23 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:51:02.6
#22 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:51:29.7
#20 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:53:11.3
#22 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:53:38.1
#21 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:54:32.5
#21 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:54:59.1
#21 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:55:53.1
#23 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:56:19.8
#22 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:57:13.7
#23 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:57:40.4
#24 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:58:34.3
#22 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:59:01.1
#20 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:59:55.0
#22 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:00:21.7
#22 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:01:09.6
#22 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:01:39.6
#20 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:03:51.3
#32 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:04:09.5
#23 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:04:54.1
#19 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:05:20.8
#21 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:06:14.6
#22 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:06:41.5
#21 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:07:35.1
#20 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:08:02.0
#20 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:08:55.5
#23 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:09:22.6
#18 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:10:16.0
#21 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:10:43.1
#22 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:11:36.5
#22 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:12:03.6
#20 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:12:45.6
#30 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:13:03.8
#21 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:14:26.3
#28 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:14:44.4
#21 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:27:57.6
#31 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:28:15.7
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:29:00.3
#21 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:29:26.9
#21 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:30:20.2
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:30:47.2
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:31:40.3
#21 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:32:07.4
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:33:01.2
#22 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:33:27.7
#21 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:34:21.1
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:34:48.0
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:35:41.4
#24 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:36:08.2
#18 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:36:49.8
#28 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:37:08.0
#22 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:38:33.9
#21 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:39:00.9
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:39:54.9
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:40:21.3
#21 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:41:14.9
#30 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:41:32.5
#21 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:42:57.5
#22 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:43:24.3
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:44:17.7
#21 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:44:44.5
#17 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:45:37.8
#18 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:46:04.8
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:46:58.4
#21 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:47:25.1
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:49:05.8
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:49:33.0
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:50:26.3
#21 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:50:53.3
#17 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:51:46.5
#21 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:52:13.5
#17 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:53:07.1
#21 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:53:33.8
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:54:27.2
#22 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:54:54.0
#18 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:55:47.4
#21 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:56:14.3
#18 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:57:20.4
#29 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:57:38.5
#17 of 43 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:59:44.0
#25 of 43 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:00:02.2
#16 of 43 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:00:46.5
#18 of 43 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:01:13.3
#18 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:02:06.7
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:02:33.6
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:03:26.7
#17 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:03:53.9
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:04:47.1
#17 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:05:14.1
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:06:07.3
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:06:34.4
#18 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:07:27.6
#17 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:07:54.6
#18 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:08:37.4
#37 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:08:55.5
#21 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:10:21.6
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:10:48.8
#18 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:11:42.3
#23 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:12:09.0
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:13:02.3
#25 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:13:20.2
#18 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:14:49.4
#33 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:15:07.2
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:15:51.5
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:16:18.4
#22 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:17:11.8
#21 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:17:38.7
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:18:32.2
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:18:58.9
#18 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:20:42.2
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:21:09.4
#18 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:22:03.0
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:22:29.7
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:23:23.0
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:23:50.0
#17 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:24:43.3
#18 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:25:10.2
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:26:03.6
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:26:30.5
#18 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:27:23.8
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:27:50.7
#22 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:28:38.0
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:29:09.5
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:31:22.7
#35 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:31:40.9
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:32:25.6
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:32:52.1
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:33:45.5
#23 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:34:12.3
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:35:05.7
#21 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:35:32.6
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:36:26.4
#18 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:36:52.9
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:37:46.3
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:38:13.1
#18 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:39:06.4
#23 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:39:33.4
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:40:14.8
#33 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:40:32.9
#22 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:41:49.2
#24 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:29:52.8
#30 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:30:10.9
#23 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:30:56.3
#25 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:31:23.6
#18 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:32:18.0
#25 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:32:45.3
#23 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:33:40.0
#25 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:34:07.0
#21 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:35:01.6
#24 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:35:28.7
#24 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:36:23.1
#19 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:36:50.4
#23 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:37:44.8
#20 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:38:12.1
#28 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:38:47.4
#22 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:40:19.3
#24 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:40:46.6
#27 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:41:40.9
#23 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:42:08.2
#23 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:43:02.6
#27 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:43:29.7
#21 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:45:05.1
#25 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:45:23.3
#25 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:46:08.6
#20 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:46:35.8
#24 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:47:30.3
#25 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:47:57.4
#24 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:48:51.8
#26 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:49:19.0
#20 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:51:02.5
#24 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:51:29.8
#22 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:52:24.1
#23 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:52:51.4
#22 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:53:45.8
#24 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:54:13.0
#20 of 42 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:55:07.3
#21 of 42 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:55:34.6
#24 of 42 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:56:28.9
#21 of 42 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:56:56.2
#22 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:57:50.5
#27 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:58:17.8
#21 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:59:06.1
#25 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:59:38.3
#20 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:01:52.9
#29 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:02:11.0
#23 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:02:56.6
#23 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:03:23.6
#23 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:04:18.0
#24 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:04:45.2
#25 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:05:39.6
#25 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:06:06.8
#19 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:07:01.2
#24 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:07:28.4
#23 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:08:22.8
#19 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:08:50.0
#20 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:09:44.6
#21 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:10:11.6
#28 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:10:47.5
#23 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:12:20.0
#20 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:12:47.2
#22 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:13:41.6
#23 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:14:08.7
#20 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:15:03.1
#33 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:15:21.2
#21 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:16:47.0
#28 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:17:05.2
#22 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:17:50.6
#23 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:18:17.6
#22 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:19:11.9
#19 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:19:39.1
#21 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:20:33.4
#23 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:21:00.6
#22 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:22:44.7
#24 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:23:11.9
#20 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:24:06.2
#24 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:24:33.4
#22 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:25:27.8
#23 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:25:54.9
#23 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:26:49.2
#22 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:27:16.4
#26 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:28:10.7
#21 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:28:37.9
#23 of 41 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:29:32.2
#26 of 41 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:29:59.4
#20 of 38 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:30:47.9
#21 of 38 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:31:18.7
#20 of 34 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:33:31.2
#24 of 34 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:33:49.4
#23 of 35 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:34:34.5
#23 of 35 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:35:01.7
#22 of 35 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:35:55.9
#24 of 35 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:36:23.1
#20 of 36 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:37:17.3
#25 of 36 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:37:44.5
#22 of 36 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:38:38.7
#24 of 36 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:39:05.9
#24 of 35 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:40:00.2
#23 of 35 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:40:27.4
#22 of 35 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:41:21.7
#21 of 35 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:41:48.8
#25 of 35 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:42:25.2
#22 of 35 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:43:55.6
#26 of 35 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:44:22.8
#22 of 34 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:45:10.9
#19 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:27:24.4
#31 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:27:42.6
#20 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:28:27.9
#22 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:28:55.1
#24 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:29:49.5
#22 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:30:16.7
#19 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:31:11.2
#25 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:31:38.3
#25 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:32:32.7
#20 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:32:59.9
#23 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:33:54.3
#22 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:34:21.5
#21 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:35:15.9
#23 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:35:43.1
#21 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:36:18.8
#19 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:37:52.6
#22 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:38:19.8
#24 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:39:14.1
#24 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:39:41.2
#25 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:40:35.7
#23 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:41:03.0
#27 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:42:39.9
#32 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:42:58.1
#25 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:43:43.4
#26 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:44:10.6
#22 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:45:05.1
#25 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:45:32.2
#18 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:46:26.7
#26 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:46:53.8
#21 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:48:36.1
#22 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:49:03.3
#20 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:49:57.8
#26 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:50:24.9
#24 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:51:19.3
#23 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:51:46.5
#21 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:52:40.9
#22 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:53:08.1
#18 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:54:02.5
#23 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:54:29.7
#23 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:55:24.1
#24 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:55:51.3
#23 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:56:39.6
#25 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:57:11.9
#24 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:59:24.3
#37 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:59:42.5
#23 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:00:27.7
#28 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:00:54.8
#21 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:01:49.1
#27 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:02:16.4
#25 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:03:10.6
#28 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:03:37.9
#24 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:04:32.1
#21 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:04:59.4
#21 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:05:53.8
#23 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:06:20.9
#22 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:07:15.2
#23 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:07:42.4
#32 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:08:18.0
#23 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:09:50.1
#29 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:10:17.4
#26 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:11:11.6
#25 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:11:38.8
#28 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:12:33.0
#27 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:13:00.2
#30 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:14:35.2
#30 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:14:53.3
#29 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:15:38.4
#25 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:16:05.7
#19 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:16:59.8
#27 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:17:27.1
#27 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:18:21.3
#33 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:18:48.5
#27 of 44 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:20:33.0
#21 of 44 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:21:00.2
#25 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:21:54.4
#29 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:22:21.6
#26 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:23:15.8
#23 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:23:43.0
#26 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:24:37.2
#33 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:25:04.4
#28 of 44 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:25:58.6
#23 of 44 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:26:25.8
#26 of 44 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:27:20.0
#25 of 44 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:27:47.2
#24 of 44 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:28:35.4
#23 of 44 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:29:07.1
#22 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:31:21.7
#35 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:31:39.9
#25 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:32:24.9
#32 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:32:52.1
#23 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:33:46.2
#25 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:34:13.4
#27 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:35:07.5
#30 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:35:34.7
#29 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:36:28.8
#29 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:36:56.1
#25 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:37:50.2
#32 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:38:17.4
#25 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:39:11.8
#29 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:39:38.7
#35 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:40:14.9
#21 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:41:45.4
#29 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:42:12.6
#21 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:43:00.0
#30 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:43:18.1
#27 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:44:03.0
#23 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:44:30.2
#20 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:45:24.1
#24 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:45:51.3
#18 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:46:45.6
#23 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:47:12.4
#19 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:48:06.3
#23 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:48:33.5
#22 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:49:27.5
#23 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:49:54.6
#19 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:50:48.6
#21 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:51:15.8
#25 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:51:53.9
#23 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:53:30.3
#23 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:53:57.6
#20 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:54:51.6
#24 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:55:18.6
#22 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:56:12.4
#36 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:56:30.5
#23 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:57:58.2
#32 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:58:16.1
#24 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:59:00.8
#23 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:59:28.0
#20 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:00:21.9
#25 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:00:49.1
#24 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:01:42.9
#25 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:02:10.1
#23 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:03:58.6
#23 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:04:25.9
#20 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:05:19.6
#23 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:05:46.8
#23 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:06:40.5
#24 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:07:07.7
#22 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:08:01.5
#19 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:08:28.7
#24 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:09:22.5
#28 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:09:49.4
#21 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:10:43.3
#27 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:11:10.5
#27 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:12:19.2
#40 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:12:37.4
#25 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:14:51.1
#31 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:15:09.3
#21 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:15:53.9
#20 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:16:21.0
#19 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:17:14.9
#26 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:17:41.9
#28 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:18:35.5
#24 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:19:02.7
#23 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:19:56.3
#23 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:20:23.5
#24 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:21:17.1
#25 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:21:44.3
#22 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:22:38.0
#22 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:23:05.2
#31 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:23:42.3
#21 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:25:18.2
#20 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:25:45.4
#23 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:26:39.3
#20 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:27:06.2
#25 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:28:00.1
#25 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:28:17.9
#24 of 43 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:29:45.2
#23 of 43 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:30:03.3
#21 of 41 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:30:47.7
#24 of 41 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:31:14.8
#22 of 41 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:32:08.3
#26 of 41 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:32:35.5
#22 of 40 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:33:29.1
#23 of 40 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:33:56.1
#19 of 38 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:35:44.7
#28 of 38 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:36:12.0
#20 of 39 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:37:05.4
#23 of 39 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:37:32.6
#22 of 39 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:38:26.1
#25 of 39 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:38:53.2
#24 of 39 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:39:46.7
#29 of 39 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:40:13.9
#22 of 37 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:41:07.4
#22 of 37 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:41:34.5
#21 of 37 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:42:28.0
#23 of 37 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:42:55.2
#30 of 40 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:43:42.6
#22 of 40 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:44:15.7
#21 of 40 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:46:35.4
#32 of 40 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:46:53.5
#26 of 40 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:47:37.8
#22 of 40 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:48:04.9
#24 of 39 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:48:58.2
#17 of 39 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:49:25.3
#18 of 38 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:50:18.6
#19 of 38 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:50:45.8
#20 of 37 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:51:39.5
#24 of 37 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:52:06.2
#23 of 37 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:52:59.6
#23 of 37 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:53:26.7
#22 of 38 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:54:20.0
#19 of 38 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:54:47.1
#26 of 38 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:55:25.6
#22 of 37 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:57:00.0
#16 of 37 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:57:27.3
#22 of 34 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:58:20.7
#24 of 34 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:58:38.5
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:39:43.5
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:40:13.8
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:40:44.0
# 1 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:41:14.2
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:41:44.5
# 1 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:42:14.7
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:42:45.0
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:43:15.2
# 1 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:43:45.4
# 1 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:44:15.7
# 1 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:44:47.5
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:45:17.7
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:47:03.8
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:47:34.1
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:48:04.3
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:48:34.6
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:49:04.8
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:49:35.0
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:50:05.3
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:50:35.5
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:51:05.8
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:51:36.0
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:52:07.8
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:52:38.0
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:53:08.3
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:54:47.5
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:55:17.7
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:55:48.0
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:58:09.8
# 6 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:58:40.1
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:59:10.3
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:59:40.6
# 6 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:00:10.8
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:00:41.0
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:01:11.3
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:01:41.5
# 6 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:02:11.8
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:02:42.0
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:03:13.8
# 6 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:03:44.0
# 6 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:04:14.3
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:08:46.0
# 6 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:09:16.3
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:09:46.5
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:10:16.8
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:10:47.0
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:11:17.2
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:11:47.5
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:12:17.7
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:12:48.0
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:13:18.2
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:13:50.0
# 7 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:14:20.2
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:14:50.5
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:17:07.6
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:17:37.8
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:18:08.1
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:18:38.3
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:19:08.5
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:19:38.8
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:20:09.0
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:20:39.3
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:21:09.5
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:21:39.7
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:22:11.5
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:22:41.8
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:23:12.0
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:25:26.5
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:25:56.8
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:27:34.8
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:28:05.1
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:28:35.3
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:29:05.5
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:29:35.8
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:30:06.0
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:30:36.3
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:30:37.0
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:31:07.2
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:32:37.9
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:34:40.4
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:35:40.8
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:37:45.4
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:38:15.6
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:39:16.1
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:39:46.2
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:40:16.6
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:41:17.0
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:42:19.1
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:42:49.1
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:47:38.5
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:48:08.8
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:48:39.0
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:49:39.5
# 3 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:50:09.7
# 3 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:50:39.6
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:51:10.3
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:51:40.5
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:52:12.3
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:52:42.5
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:53:44.4
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:54:14.8
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:54:45.0
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:58:03.6
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:59:04.1
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:59:34.4
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:00:04.6
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:00:34.4
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:01:35.1
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:02:05.6
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:02:35.8
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:03:07.6
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:03:37.8
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:07:13.6
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:07:43.0
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:08:44.6
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:10:15.4
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:11:17.3
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:11:47.4
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:13:32.2
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/03:56:49.6
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/03:57:19.8
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/03:57:50.1
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/03:58:20.3
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/03:58:50.5
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:00:51.5
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:01:21.7
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:01:53.5
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:02:54.0
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:04:29.0
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:04:59.2
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:06:29.9
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:07:00.2
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:07:30.4
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:08:00.7
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:09:01.1
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:09:32.9
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:10:03.2
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:10:33.4
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:11:47.1
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:12:47.6
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:14:22.7
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:14:52.9
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:15:23.1
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:15:53.4
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:16:23.6
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:16:53.9
# 3 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:17:24.0
# 3 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:17:54.3
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:18:24.6
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:18:54.8
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:19:26.6
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:19:56.8
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:20:59.0
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:21:29.1
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:21:59.3
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:25:18.4
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:25:48.6
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:26:18.8
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:27:19.3
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:28:19.8
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:28:50.1
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:29:20.3
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:29:50.5
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:31:22.8
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:32:57.5
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:33:27.7
# 3 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:33:58.0
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:34:28.2
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:34:58.4
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:35:28.7
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:35:58.9
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:36:29.2
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:36:59.3
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:37:29.4
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:38:01.2
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:38:31.7
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:39:01.8
# 3 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:40:14.4
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:40:44.6
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:41:14.9
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:43:19.0
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/03:48:08.0
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/03:52:14.8
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/03:54:15.8
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/03:57:48.9
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/04:13:05.1
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/04:15:06.1
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/04:18:09.4
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/04:30:34.6

#Inspect gain tables interactively and decide whether to manually flag something
#plotms(LB_p5,xaxis='time',yaxis='GainPhase',iteraxis='spw',showgui=True)
#Print calibration png file
#plotms(
#    LB_p5,xaxis='time', yaxis='GainPhase',iteraxis='spw',exprange='all',
#    overwrite=True,showgui=False,
#    plotfile=os.path.join(LB_selfcal_folder,f'{prefix}_LB_p5_gain_phase_vs_time.png')
#)

#Inspect gain tables and decide whether to flag something
plotms(
    LB_p5,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(LB_selfcal_folder,f'{prefix}_LB_p5_gain_phase_vs_time.png'),
)

#Flag problematic antennas
flagdata(vis=LB_p5,mode='manual',spw='0', timerange='2017/10/10/09:05:30~2017/10/10/09:05:40',antenna='DV12')
flagdata(vis=LB_p5,mode='manual',spw='0', timerange='2017/10/10/09:05:55~2017/10/10/09:06:05',antenna='DV24')
flagdata(vis=LB_p5,mode='manual',spw='0', timerange='2017/10/10/09:06:50~2017/10/10/09:07:00',antenna='DV24')
flagdata(vis=LB_p5,mode='manual',spw='0', timerange='2017/10/10/09:08:10~2017/10/10/09:08:20',antenna='DA45')
flagdata(vis=LB_p5,mode='manual',spw='0', timerange='2017/10/10/09:09:40~2017/10/10/09:09:50',antenna='DA58,DV22,DA45,DA65,DV13')
flagdata(vis=LB_p5,mode='manual',spw='0', timerange='2017/10/10/09:11:00~2017/10/10/09:11:10',antenna='DA63')
flagdata(vis=LB_p5,mode='manual',spw='0', timerange='2017/10/10/09:12:30~2017/10/10/09:12:40',antenna='DA54')
flagdata(vis=LB_p5,mode='manual',spw='0', timerange='2017/10/10/09:14:00~2017/10/10/09:14:10',antenna='DA63,PM03')
flagdata(vis=LB_p5,mode='manual',spw='0', timerange='2017/10/10/09:15:30~2017/10/10/09:15:40',antenna='DA54,DA56,DA60,PM04')
flagdata(vis=LB_p5,mode='manual',spw='0', timerange='2017/10/10/09:15:50~2017/10/10/09:16:00',antenna='DA59,DA60')
flagdata(vis=LB_p5,mode='manual',spw='0', timerange='2017/10/10/09:16:30~2017/10/10/09:16:40',antenna='DA52,DV22')
flagdata(vis=LB_p5,mode='manual',spw='0', timerange='2017/10/10/09:22:45~2017/10/10/09:22:55',antenna='DA47,DA63,DV16')
flagdata(vis=LB_p5,mode='manual',spw='0', timerange='2017/10/10/09:23:10~2017/10/10/09:23:20',antenna='DA65')
flagdata(vis=LB_p5,mode='manual',spw='0', timerange='2017/10/10/09:27:15~2017/10/10/09:27:25',antenna='DV05')
flagdata(vis=LB_p5,mode='manual',spw='0', timerange='2017/10/10/09:28:10~2017/10/10/09:28:20',antenna='DA47,DA60')
flagdata(vis=LB_p5,mode='manual',spw='0', timerange='2017/10/10/09:29:20~2017/10/10/09:29:30',antenna='DA65,DV24')
flagdata(vis=LB_p5,mode='manual',spw='0', timerange='2017/10/10/09:32:20~2017/10/10/09:32:30',antenna='DV23')
flagdata(vis=LB_p5,mode='manual',spw='0', timerange='2017/10/10/09:35:50~2017/10/10/09:36:00',antenna='DA57,DA60')
flagdata(vis=LB_p5,mode='manual',spw='0', timerange='2017/10/10/09:38:30~2017/10/10/09:38:40',antenna='DA52')
flagdata(vis=LB_p5,mode='manual',spw='0', timerange='2017/10/10/09:39:00~2017/10/10/09:39:10',antenna='DV05')
flagdata(vis=LB_p5,mode='manual',spw='0', timerange='2017/10/10/09:41:20~2017/10/10/09:41:30',antenna='DA57,DV24')
flagdata(vis=LB_p5,mode='manual',spw='0', timerange='2017/10/10/09:44:10~2017/10/10/09:44:20',antenna='DA52')
flagdata(vis=LB_p5,mode='manual',spw='0', timerange='2017/10/10/09:45:30~2017/10/10/09:45:40',antenna='DV05')
flagdata(vis=LB_p5,mode='manual',spw='0', timerange='2017/10/10/09:45:50~2017/10/10/09:46:00',antenna='DA52,DA64')
flagdata(vis=LB_p5,mode='manual',spw='0', timerange='2017/10/10/09:49:40~2017/10/10/09:49:50',antenna='DA57')
flagdata(vis=LB_p5,mode='manual',spw='0', timerange='2017/10/10/09:56:10~2017/10/10/09:56:20',antenna='DV05')
flagdata(vis=LB_p5,mode='manual',spw='0', timerange='2017/10/10/09:59:50~2017/10/10/10:00:00',antenna='DV04')
flagdata(vis=LB_p5,mode='manual',spw='0', timerange='2017/10/10/10:01:30~2017/10/10/10:01:40',antenna='DV22')
flagdata(vis=LB_p5,mode='manual',spw='0', timerange='2017/10/10/10:03:50~2017/10/10/10:04:00',antenna='DA57')
flagdata(vis=LB_p5,mode='manual',spw='0', timerange='2017/10/10/10:05:20~2017/10/10/10:05:30',antenna='DA57')
flagdata(vis=LB_p5,mode='manual',spw='0', timerange='2017/10/10/10:07:30~2017/10/10/10:07:40',antenna='DA47')
flagdata(vis=LB_p5,mode='manual',spw='0', timerange='2017/10/10/10:08:00~2017/10/10/10:08:10',antenna='DV04')
flagdata(vis=LB_p5,mode='manual',spw='0', timerange='2017/10/10/10:09:20~2017/10/10/10:09:30',antenna='DA54,DV04')
flagdata(vis=LB_p5,mode='manual',spw='0', timerange='2017/10/10/10:10:10~2017/10/10/10:10:20',antenna='DA64,DV12')
flagdata(vis=LB_p5,mode='manual',spw='0', timerange='2017/10/10/10:13:00~2017/10/10/10:13:10',antenna='DA59,DV20')

flagdata(vis=LB_p5,mode='manual',spw='4', timerange='2017/10/10/10:27:50~2017/10/10/10:28:00',antenna='DV22')
flagdata(vis=LB_p5,mode='manual',spw='4', timerange='2017/10/10/10:28:10~2017/10/10/10:28:20',antenna='DA60')
flagdata(vis=LB_p5,mode='manual',spw='4', timerange='2017/10/10/10:30:40~2017/10/10/10:30:50',antenna='DA52')
flagdata(vis=LB_p5,mode='manual',spw='4', timerange='2017/10/10/10:31:40~2017/10/10/10:31:50',antenna='DV24')
flagdata(vis=LB_p5,mode='manual',spw='4', timerange='2017/10/10/10:34:40~2017/10/10/10:34:50',antenna='DV03')
flagdata(vis=LB_p5,mode='manual',spw='4', timerange='2017/10/10/10:36:00~2017/10/10/10:36:10',antenna='DA60')
flagdata(vis=LB_p5,mode='manual',spw='4', timerange='2017/10/10/10:37:00~2017/10/10/10:37:10',antenna='DA45,DV10')
flagdata(vis=LB_p5,mode='manual',spw='4', timerange='2017/10/10/10:39:00~2017/10/10/10:39:10',antenna='DA47')
flagdata(vis=LB_p5,mode='manual',spw='4', timerange='2017/10/10/10:40:20~2017/10/10/10:40:30',antenna='DA58')
flagdata(vis=LB_p5,mode='manual',spw='4', timerange='2017/10/10/10:41:30~2017/10/10/10:41:40',antenna='DA45')
flagdata(vis=LB_p5,mode='manual',spw='4', timerange='2017/10/10/10:46:00~2017/10/10/10:46:10',antenna='DA48,DV08')
flagdata(vis=LB_p5,mode='manual',spw='4', timerange='2017/10/10/10:50:20~2017/10/10/10:50:30',antenna='DA44,DA54,DV15')
flagdata(vis=LB_p5,mode='manual',spw='4', timerange='2017/10/10/10:50:50~2017/10/10/10:51:00',antenna='DA44,DV15')
flagdata(vis=LB_p5,mode='manual',spw='4', timerange='2017/10/10/10:51:40~2017/10/10/10:51:50',antenna='DA44,DA54,DA64')
flagdata(vis=LB_p5,mode='manual',spw='4', timerange='2017/10/10/11:00:40~2017/10/10/11:00:50',antenna='DV08,DV17')
flagdata(vis=LB_p5,mode='manual',spw='4', timerange='2017/10/10/11:03:50~2017/10/10/11:04:00',antenna='DA59,DV03')
flagdata(vis=LB_p5,mode='manual',spw='4', timerange='2017/10/10/11:07:50~2017/10/10/11:08:00',antenna='DV12')
flagdata(vis=LB_p5,mode='manual',spw='4', timerange='2017/10/10/11:08:30~2017/10/10/11:08:40',antenna='DA48,DA65')
flagdata(vis=LB_p5,mode='manual',spw='4', timerange='2017/10/10/11:11:40~2017/10/10/11:11:50',antenna='DA43,DA44')
flagdata(vis=LB_p5,mode='manual',spw='4', timerange='2017/10/10/11:12:00~2017/10/10/11:12:10',antenna='DV13')
flagdata(vis=LB_p5,mode='manual',spw='4', timerange='2017/10/10/11:13:00~2017/10/10/11:13:10',antenna='DA53,DV22')
flagdata(vis=LB_p5,mode='manual',spw='4', timerange='2017/10/10/11:13:15~2017/10/10/11:13:25',antenna='DA53,DA56')
flagdata(vis=LB_p5,mode='manual',spw='4', timerange='2017/10/10/11:16:10~2017/10/10/11:16:20',antenna='DA59')
flagdata(vis=LB_p5,mode='manual',spw='4', timerange='2017/10/10/11:18:30~2017/10/10/11:18:40',antenna='DV16')
flagdata(vis=LB_p5,mode='manual',spw='4', timerange='2017/10/10/11:20:40~2017/10/10/11:20:50',antenna='DA60')
flagdata(vis=LB_p5,mode='manual',spw='4', timerange='2017/10/10/11:22:00~2017/10/10/11:22:10',antenna='DA60')
flagdata(vis=LB_p5,mode='manual',spw='4', timerange='2017/10/10/11:23:45~2017/10/10/11:23:55',antenna='DA52')
flagdata(vis=LB_p5,mode='manual',spw='4', timerange='2017/10/10/11:24:40~2017/10/10/11:24:50',antenna='DA45,DA56')
flagdata(vis=LB_p5,mode='manual',spw='4', timerange='2017/10/10/11:26:00~2017/10/10/11:26:10',antenna='DA48')
flagdata(vis=LB_p5,mode='manual',spw='4', timerange='2017/10/10/11:26:25~2017/10/10/11:26:35',antenna='DV20')
flagdata(vis=LB_p5,mode='manual',spw='4', timerange='2017/10/10/11:27:50~2017/10/10/11:28:00',antenna='DA48')
flagdata(vis=LB_p5,mode='manual',spw='4', timerange='2017/10/10/11:31:20~2017/10/10/11:31:30',antenna='DA63')
flagdata(vis=LB_p5,mode='manual',spw='4', timerange='2017/10/10/11:31:40~2017/10/10/11:31:50',antenna='DA65')
flagdata(vis=LB_p5,mode='manual',spw='4', timerange='2017/10/10/11:36:20~2017/10/10/11:36:30',antenna='DV22')
flagdata(vis=LB_p5,mode='manual',spw='4', timerange='2017/10/10/11:36:50~2017/10/10/11:37:00',antenna='DA48,DA54,DV08,DV22')
flagdata(vis=LB_p5,mode='manual',spw='4', timerange='2017/10/10/11:39:00~2017/10/10/11:39:10',antenna='DA48')

flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/10/11/07:32:10~2017/10/11/07:32:20',antenna='DA64,DV04,DV22')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/10/11/07:32:40~2017/10/11/07:32:50',antenna='DA43')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/10/11/07:33:35~2017/10/11/07:33:45',antenna='DA43,DA44,DV15,DV20')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/10/11/07:34:00~2017/10/11/07:34:10',antenna='DA43,DV11')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/10/11/07:35:00~2017/10/11/07:35:10',antenna='DA44,DV05,DV16')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/10/11/07:35:20~2017/10/11/07:35:30',antenna='DA60,DV15')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/10/11/07:36:20~2017/10/11/07:36:30',antenna='DA58,DA64')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/10/11/07:36:45~2017/10/11/07:36:55',antenna='DA56')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/10/11/07:37:40~2017/10/11/07:37:50',antenna='DA44,DV15')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/10/11/07:38:10~2017/10/11/07:38:20',antenna='DV08')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/10/11/07:38:45~2017/10/11/07:38:50',antenna='DV04')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/10/11/07:40:10~2017/10/11/07:40:20',antenna='DA63')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/10/11/07:43:00~2017/10/11/07:43:10',antenna='DA52,DV15')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/10/11/07:45:20~2017/10/11/07:45:30',antenna='DA45,DA52,DV22')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/10/11/07:47:30~2017/10/11/07:47:35',antenna='DV16')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/10/11/07:53:45~2017/10/11/07:53:50',antenna='DA60')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/10/11/07:55:30~2017/10/11/07:55:40',antenna='DV17')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/10/11/07:56:20~2017/10/11/07:56:30',antenna='DA64,DV20')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/10/11/07:56:50~2017/10/11/07:57:00',antenna='DA52,DA53')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/10/11/08:02:10~2017/10/11/08:02:20',antenna='DA60,DA63,DV01')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/10/11/08:06:00~2017/10/11/08:06:10',antenna='DA45')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/10/11/08:07:00~2017/10/11/08:07:10',antenna='DV15')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/10/11/08:07:20~2017/10/11/08:07:30',antenna='DA63')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/10/11/08:08:20~2017/10/11/08:08:30',antenna='DA54,DV15')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/10/11/08:08:40~2017/10/11/08:08:50',antenna='DA56,DV11,DV17,DV22')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/10/11/08:09:40~2017/10/11/08:09:50',antenna='DV11,DV15')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/10/11/08:10:10~2017/10/11/08:10:20',antenna='DV05,DV15')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/10/11/08:12:15~2017/10/11/08:12:25',antenna='DV16')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/10/11/08:15:20~2017/10/11/08:15:25',antenna='DA46,DA50,DA54,DA63,DV11')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/10/11/08:17:50~2017/10/11/08:18:00',antenna='DA54,DA56,DV17')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/10/11/08:19:10~2017/10/11/08:19:20',antenna='DV20')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/10/11/08:19:30~2017/10/11/08:19:40',antenna='DA45,DA65,DV03')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/10/11/08:21:00~2017/10/11/08:21:10',antenna='DV04')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/10/11/08:24:00~2017/10/11/08:24:10',antenna='DA45')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/10/11/08:25:50~2017/10/11/08:26:00',antenna='DA52,DA53')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/10/11/08:29:30~2017/10/11/08:29:40',antenna='DA44,DV13')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/10/11/08:29:50~2017/10/11/08:30:00',antenna='DV13,DV17')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/10/11/08:29:55~2017/10/11/08:30:00',antenna='DV13,DV17')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/10/11/08:30:45~2017/10/11/08:30:50',antenna='PM04')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/10/11/08:31:15~2017/10/11/08:31:25',antenna='DV05,DV08')

flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/10/11/08:33:00~2017/10/11/08:53:00',antenna='')

flagdata(vis=LB_p5,mode='manual',spw='12',timerange='',antenna='')

flagdata(vis=LB_p5,mode='manual',spw='16',timerange='',antenna='')

flagdata(vis=LB_p5,mode='manual',spw='20',timerange='2017/12/09/04:40:40~2017/12/09/04:40:50',antenna='DV10')
flagdata(vis=LB_p5,mode='manual',spw='20',timerange='2017/12/09/04:41:40~2017/12/09/04:41:50',antenna='DA46')
flagdata(vis=LB_p5,mode='manual',spw='20',timerange='2017/12/09/04:43:40~2017/12/09/04:43:50',antenna='DV15')
flagdata(vis=LB_p5,mode='manual',spw='20',timerange='2017/12/09/04:44:45~2017/12/09/04:44:50',antenna='DV15')
flagdata(vis=LB_p5,mode='manual',spw='20',timerange='2017/12/09/04:47:30~2017/12/09/04:47:40',antenna='DA46,DV11')
flagdata(vis=LB_p5,mode='manual',spw='20',timerange='2017/12/09/04:54:45~2017/12/09/04:54:50',antenna='DV11,DV15')
flagdata(vis=LB_p5,mode='manual',spw='20',timerange='2017/12/09/04:55:45~2017/12/09/04:55:50',antenna='DV11')
flagdata(vis=LB_p5,mode='manual',spw='20',timerange='2017/12/09/04:59:40~2017/12/09/04:59:50',antenna='DV15')
flagdata(vis=LB_p5,mode='manual',spw='20',timerange='2017/12/09/05:03:40~2017/12/09/05:03:50',antenna='DA49,DV13')
flagdata(vis=LB_p5,mode='manual',spw='20',timerange='2017/12/09/05:04:10~2017/12/09/05:04:15',antenna='DA49')
flagdata(vis=LB_p5,mode='manual',spw='20',timerange='2017/12/09/05:10:45~2017/12/09/05:10:50',antenna='DV10')
flagdata(vis=LB_p5,mode='manual',spw='20',timerange='2017/12/09/05:12:15~2017/12/09/05:12:20',antenna='DA46,DV07')
flagdata(vis=LB_p5,mode='manual',spw='20',timerange='2017/12/09/05:13:40~2017/12/09/05:13:50',antenna='DA46,DV13')
flagdata(vis=LB_p5,mode='manual',spw='20',timerange='2017/12/09/05:17:30~2017/12/09/05:17:40',antenna='DV10')
flagdata(vis=LB_p5,mode='manual',spw='20',timerange='2017/12/09/05:18:00~2017/12/09/05:18:10',antenna='DV10')
flagdata(vis=LB_p5,mode='manual',spw='20',timerange='2017/12/09/05:17:00~2017/12/09/05:17:10',antenna='DV11')
flagdata(vis=LB_p5,mode='manual',spw='20',timerange='2017/12/09/05:19:00~2017/12/09/05:19:10',antenna='DA50,DV13')
flagdata(vis=LB_p5,mode='manual',spw='20',timerange='2017/12/09/05:20:00~2017/12/09/05:20:10',antenna='DV15')
flagdata(vis=LB_p5,mode='manual',spw='20',timerange='2017/12/09/05:20:30~2017/12/09/05:20:40',antenna='DA50,DV10,DV15,DA45')
flagdata(vis=LB_p5,mode='manual',spw='20',timerange='2017/12/09/05:21:00~2017/12/09/05:21:10',antenna='DA46')
flagdata(vis=LB_p5,mode='manual',spw='20',timerange='2017/12/09/05:21:20~2017/12/09/05:21:40',antenna='DV11,DA46')
flagdata(vis=LB_p5,mode='manual',spw='20',timerange='2017/12/09/05:22:40~2017/12/09/05:22:50',antenna='DA46,DA55')
flagdata(vis=LB_p5,mode='manual',spw='20',timerange='2017/12/09/05:23:10~2017/12/09/05:23:20',antenna='DA46')
flagdata(vis=LB_p5,mode='manual',spw='20',timerange='2017/12/09/05:25:20~2017/12/09/05:25:30',antenna='DA46,DA51,DA55')
flagdata(vis=LB_p5,mode='manual',spw='20',timerange='2017/12/09/05:25:50~2017/12/09/05:26:00',antenna='DA51')
flagdata(vis=LB_p5,mode='manual',spw='20',timerange='2017/12/09/05:27:30~2017/12/09/05:27:40',antenna='DA46,DA51,DV11,DA50')
flagdata(vis=LB_p5,mode='manual',spw='20',timerange='2017/12/09/05:28:00~2017/12/09/05:28:20',antenna='DA51,DV11')
flagdata(vis=LB_p5,mode='manual',spw='20',timerange='2017/12/09/05:29:00~2017/12/09/05:29:20',antenna='DV15')
flagdata(vis=LB_p5,mode='manual',spw='20',timerange='2017/12/09/05:30:30~2017/12/09/05:30:40',antenna='DA46')

flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/12/17/05:49:05~2017/12/17/05:49:10',antenna='DA55')
flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/12/17/05:51:30~2017/12/17/05:51:50',antenna='DA44,DV07')
flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/12/17/05:53:00~2017/12/17/05:53:20',antenna='DA43,DA55')
flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/12/17/05:53:40~2017/12/17/05:53:50',antenna='DA43,DA44,DV13')
flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/12/17/05:54:10~2017/12/17/05:54:20',antenna='DA43,DA55,DV13')
flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/12/17/05:54:40~2017/12/17/05:54:50',antenna='DA41,DA42,DA43,DA44,DA55,DV13')
flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/12/17/05:58:00~2017/12/17/05:58:10',antenna='DV11,DV13,DV07,DA44')
flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/12/17/05:58:30~2017/12/17/05:58:40',antenna='DV13,DA55')
flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/12/17/05:59:00~2017/12/17/05:59:10',antenna='DA55,DV13,DV07')
flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/12/17/05:59:30~2017/12/17/05:59:40',antenna='DV13,DV07')
flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/12/17/06:00:00~2017/12/17/06:00:10',antenna='DV13')
flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/12/17/06:00:30~2017/12/17/06:00:40',antenna='DV13')
flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/12/17/06:01:00~2017/12/17/06:01:10',antenna='DV13,DV07')
flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/12/17/06:01:30~2017/12/17/06:01:40',antenna='DV13,DV07')
flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/12/17/06:02:00~2017/12/17/06:02:10',antenna='DV13,DV07')
flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/12/17/06:02:30~2017/12/17/06:02:40',antenna='DV13,DV07')
flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/12/17/06:03:00~2017/12/17/06:03:10',antenna='DV13,DV07')
flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/12/17/06:03:30~2017/12/17/06:03:40',antenna='DV13,DV07')
flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/12/17/06:04:00~2017/12/17/06:04:10',antenna='DV13')
flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/12/17/06:08:10~2017/12/17/06:08:20',antenna='DV07')

flagdata(vis=LB_p5,mode='manual',spw='28',timerange='2017/12/27/04:08:00~2017/12/27/04:08:05',antenna='DA42,DV13')
flagdata(vis=LB_p5,mode='manual',spw='28',timerange='2017/12/27/04:08:30~2017/12/27/04:08:35',antenna='DA42,DA55,DV13')
flagdata(vis=LB_p5,mode='manual',spw='28',timerange='2017/12/27/04:09:00~2017/12/27/04:09:05',antenna='DA42,DV13')
flagdata(vis=LB_p5,mode='manual',spw='28',timerange='2017/12/27/04:09:30~2017/12/27/04:09:35',antenna='DA42,DA44,DV13')
flagdata(vis=LB_p5,mode='manual',spw='28',timerange='2017/12/27/04:10:00~2017/12/27/04:10:05',antenna='DA42,DA44,DV11,DV13')
flagdata(vis=LB_p5,mode='manual',spw='28',timerange='2017/12/27/04:10:30~2017/12/27/04:10:35',antenna='DA42,DA44,DV13')
flagdata(vis=LB_p5,mode='manual',spw='28',timerange='2017/12/27/04:11:45~2017/12/27/04:11:50',antenna='DA42,DA44,DV13')
flagdata(vis=LB_p5,mode='manual',spw='28',timerange='2017/12/27/04:25:45~2017/12/27/04:25:50',antenna='DV11')
flagdata(vis=LB_p5,mode='manual',spw='28',timerange='2017/12/27/04:32:50~2017/12/27/04:33:00',antenna='DA44')
flagdata(vis=LB_p5,mode='manual',spw='28',timerange='2017/12/27/04:33:25~2017/12/27/04:33:30',antenna='DA44')
flagdata(vis=LB_p5,mode='manual',spw='28',timerange='2017/12/27/04:34:20~2017/12/27/04:34:30',antenna='DA44')
flagdata(vis=LB_p5,mode='manual',spw='28',timerange='2017/12/27/04:34:50~2017/12/27/04:35:00',antenna='DA44')
flagdata(vis=LB_p5,mode='manual',spw='28',timerange='2017/12/27/04:35:20~2017/12/27/04:35:30',antenna='DA44')
flagdata(vis=LB_p5,mode='manual',spw='28',timerange='2017/12/27/04:35:50~2017/12/27/04:36:00',antenna='DA44')
flagdata(vis=LB_p5,mode='manual',spw='28',timerange='2017/12/27/04:36:20~2017/12/27/04:36:30',antenna='DA44')
flagdata(vis=LB_p5,mode='manual',spw='28',timerange='2017/12/27/04:36:50~2017/12/27/04:37:00',antenna='DA44')
flagdata(vis=LB_p5,mode='manual',spw='28',timerange='2017/12/27/04:37:20~2017/12/27/04:37:30',antenna='DA44')
flagdata(vis=LB_p5,mode='manual',spw='28',timerange='2017/12/27/04:38:00~2017/12/27/04:38:10',antenna='DA44')
flagdata(vis=LB_p5,mode='manual',spw='28',timerange='2017/12/27/04:38:30~2017/12/27/04:38:40',antenna='DA44')
flagdata(vis=LB_p5,mode='manual',spw='28',timerange='2017/12/27/04:39:00~2017/12/27/04:39:10',antenna='DA44')
flagdata(vis=LB_p5,mode='manual',spw='28',timerange='2017/12/27/04:40:40~2017/12/27/04:40:45',antenna='DA55')
flagdata(vis=LB_p5,mode='manual',spw='28',timerange='2017/12/27/04:42:45~2017/12/27/04:42:50',antenna='DV11')
flagdata(vis=LB_p5,mode='manual',spw='28',timerange='2017/12/27/04:43:15~2017/12/27/04:43:20',antenna='DA43')

flagdata(vis=LB_p5,mode='manual',spw='32',timerange='2017/12/28/03:46:00~2017/12/28/03:46:20',antenna='DV11')
flagdata(vis=LB_p5,mode='manual',spw='32',timerange='2017/12/28/03:46:30~2017/12/28/03:46:40',antenna='DA46')
flagdata(vis=LB_p5,mode='manual',spw='32',timerange='2017/12/28/03:48:30~2017/12/28/03:48:40',antenna='DV11')
flagdata(vis=LB_p5,mode='manual',spw='32',timerange='2017/12/28/03:54:40~2017/12/28/03:54:50',antenna='DA55,DV11')
flagdata(vis=LB_p5,mode='manual',spw='32',timerange='2017/12/28/03:55:15~2017/12/28/03:55:20',antenna='DA55')
flagdata(vis=LB_p5,mode='manual',spw='32',timerange='2017/12/28/03:55:45~2017/12/28/03:55:50',antenna='DV11')
flagdata(vis=LB_p5,mode='manual',spw='32',timerange='2017/12/28/03:56:15~2017/12/28/03:56:20',antenna='DA55,DV11')
flagdata(vis=LB_p5,mode='manual',spw='32',timerange='2017/12/28/03:56:45~2017/12/28/03:56:50',antenna='DA55,DV11')
flagdata(vis=LB_p5,mode='manual',spw='32',timerange='2017/12/28/03:57:15~2017/12/28/03:57:20',antenna='DA55')
flagdata(vis=LB_p5,mode='manual',spw='32',timerange='2017/12/28/03:59:00~2017/12/28/03:59:05',antenna='DA55,DV11')
flagdata(vis=LB_p5,mode='manual',spw='32',timerange='2017/12/28/04:00:00~2017/12/28/04:00:10',antenna='DA44,DA55')
flagdata(vis=LB_p5,mode='manual',spw='32',timerange='2017/12/28/04:03:10~2017/12/28/04:03:15',antenna='DA44')
flagdata(vis=LB_p5,mode='manual',spw='32',timerange='2017/12/28/04:03:40~2017/12/28/04:03:45',antenna='DA44,DA55,DA46')
flagdata(vis=LB_p5,mode='manual',spw='32',timerange='2017/12/28/04:04:05~2017/12/28/04:04:15',antenna='DA44,DA55')
flagdata(vis=LB_p5,mode='manual',spw='32',timerange='2017/12/28/04:04:40~2017/12/28/04:04:45',antenna='DA44')
flagdata(vis=LB_p5,mode='manual',spw='32',timerange='2017/12/28/04:05:10~2017/12/28/04:05:15',antenna='DA44')
flagdata(vis=LB_p5,mode='manual',spw='32',timerange='2017/12/28/04:05:30~2017/12/28/04:06:00',antenna='DA44')
flagdata(vis=LB_p5,mode='manual',spw='32',timerange='2017/12/28/04:09:15~2017/12/28/04:09:20',antenna='DA42,DV11,DV13,DA55')
flagdata(vis=LB_p5,mode='manual',spw='32',timerange='2017/12/28/04:12:30~2017/12/28/04:12:40',antenna='DA42,DA55,DV13')
flagdata(vis=LB_p5,mode='manual',spw='32',timerange='2017/12/28/04:13:00~2017/12/28/04:13:10',antenna='DV11')
flagdata(vis=LB_p5,mode='manual',spw='32',timerange='2017/12/28/04:13:30~2017/12/28/04:13:40',antenna='DA55,DV11')
flagdata(vis=LB_p5,mode='manual',spw='32',timerange='2017/12/28/04:17:35~2017/12/28/04:17:40',antenna='DA44,DA55')
flagdata(vis=LB_p5,mode='manual',spw='32',timerange='2017/12/28/04:18:00~2017/12/28/04:18:10',antenna='DA44')
flagdata(vis=LB_p5,mode='manual',spw='32',timerange='2017/12/28/04:18:30~2017/12/28/04:18:40',antenna='DA55')
flagdata(vis=LB_p5,mode='manual',spw='32',timerange='2017/12/28/04:20:00~2017/12/28/04:20:20',antenna='DA44')
flagdata(vis=LB_p5,mode='manual',spw='32',timerange='2017/12/28/04:20:40~2017/12/28/04:20:45',antenna='DA44,DV11,DV13,DA51')
flagdata(vis=LB_p5,mode='manual',spw='32',timerange='2017/12/28/04:21:10~2017/12/28/04:21:15',antenna='DA44,DV11,DV13')
flagdata(vis=LB_p5,mode='manual',spw='32',timerange='2017/12/28/04:21:40~2017/12/28/04:21:45',antenna='DA44,DV13,DV11')
flagdata(vis=LB_p5,mode='manual',spw='32',timerange='2017/12/28/04:22:10~2017/12/28/04:22:15',antenna='DV11,DV13')
flagdata(vis=LB_p5,mode='manual',spw='32',timerange='2017/12/28/04:22:40~2017/12/28/04:22:45',antenna='DV11,DV13')
flagdata(vis=LB_p5,mode='manual',spw='32',timerange='2017/12/28/04:23:10~2017/12/28/04:23:15',antenna='DV11,DV13')
flagdata(vis=LB_p5,mode='manual',spw='32',timerange='2017/12/28/04:28:00~2017/12/28/04:28:05',antenna='DV11')
flagdata(vis=LB_p5,mode='manual',spw='32',timerange='2017/12/28/04:30:00~2017/12/28/04:30:10',antenna='DV13')
flagdata(vis=LB_p5,mode='manual',spw='32',timerange='2017/12/28/04:30:30~2017/12/28/04:30:40',antenna='DV13')

#Inspect gain tables to check if flagging worked
plotms(
    LB_p5,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(LB_selfcal_folder,f'{prefix}_LB_p5_gain_phase_vs_time_flagged.png'),
)

#Apply the calibration gains and split-off the corrected ms
applycal(
    vis=LB_cont_p4+'.ms',spw=LB_contspws,spwmap=LB_spw_mapping,
    gaintable=[LB_p5],interp='linearPD',calwt=True,applymode='calonly'
)
LB_cont_p5 = LB_cont_p4.replace('p4','p5')
os.system('rm -rf %s.ms*' %LB_cont_p5)
split(vis=LB_cont_p4+'.ms',outputvis=LB_cont_p5+'.ms',datacolumn='corrected')

#CLEAN, again with a ~6sigma threshold
tclean_wrapper(
    vis       = LB_cont_p5+'.ms',
    imagename = LB_cont_p5,
    threshold = '0.0392mJy',
    **LB_tclean_wrapper_kwargs
)
estimate_SNR(LB_cont_p5+'.image',disk_mask=LB_mask,noise_mask=noise_annulus_LB)
generate_image_png(
    LB_cont_p5+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_LB,10*rms_LB],
    save_folder=LB_selfcal_folder
)
#MWC_758_SBLB_contp5.image
#Beam 0.065 arcsec x 0.045 arcsec (7.20 deg)
#Flux inside disk mask: 66.30 mJy
#Peak intensity of source: 1.22 mJy/beam
#rms: 6.54e-03 mJy/beam
#Peak SNR: 186.87

#Sixth round of phase-only self-cal
LB_p6 = LB_p5.replace('p5','p6')
os.system('rm -rf '+LB_p6)
gaincal(
    vis=LB_cont_p5+'.ms',caltable=LB_p6,
    gaintype='T',combine='spw',calmode='p',solint='18s',
    spw=LB_contspws,refant=LB_refant,
    minsnr=3.,minblperant=4
)
#19 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:00:19.1
#24 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:00:37.2
#23 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:01:22.1
#24 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:01:40.3
#25 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:01:58.2
#25 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:02:43.5
#20 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:03:01.5
#24 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:03:19.6
#24 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:04:04.6
#23 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:04:22.7
#24 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:04:40.8
#24 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:05:25.9
#22 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:05:43.9
#27 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:06:02.0
#24 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:06:47.0
#20 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:07:05.1
#23 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:07:23.3
#27 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:08:08.3
#24 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:08:26.4
#26 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:08:44.5
#25 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:09:17.2
#26 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:09:35.3
#21 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:11:03.1
#23 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:11:21.2
#28 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:11:39.4
#24 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:12:24.2
#24 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:12:42.3
#26 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:13:00.5
#24 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:13:45.3
#22 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:14:03.3
#27 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:15:27.8
#22 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:15:46.0
#25 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:16:30.8
#25 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:16:49.0
#24 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:17:07.1
#25 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:17:52.0
#26 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:18:10.1
#28 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:18:28.2
#25 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:19:13.3
#22 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:19:31.2
#25 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:19:49.3
#23 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:21:22.7
#24 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:21:40.8
#28 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:21:58.9
#27 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:22:43.9
#22 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:23:01.8
#27 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:23:20.0
#24 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:24:04.7
#26 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:24:22.8
#26 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:24:41.0
#21 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:25:25.7
#26 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:25:43.9
#21 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:26:02.0
#21 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:26:46.8
#25 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:27:04.9
#23 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:27:23.0
#26 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:28:07.8
#24 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:28:25.9
#21 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:28:44.1
#24 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:29:28.8
#27 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:30:00.4
#27 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:32:05.4
#24 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:32:23.6
#24 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:33:08.7
#24 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:33:26.4
#25 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:33:44.5
#24 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:34:29.1
#23 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:34:47.3
#24 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:35:05.4
#24 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:35:50.1
#24 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:36:08.2
#25 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:36:26.4
#22 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:37:11.2
#23 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:37:29.1
#23 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:37:47.3
#27 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:38:32.3
#24 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:38:50.1
#20 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:39:08.2
#25 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:39:52.9
#22 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:40:11.0
#27 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:40:29.1
#25 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:41:02.4
#21 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:41:20.5
#29 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:42:47.3
#23 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:43:05.4
#27 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:43:23.6
#24 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:44:08.3
#22 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:44:26.3
#24 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:44:44.4
#23 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:45:29.3
#24 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:45:47.1
#23 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:47:12.7
#20 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:47:30.5
#26 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:48:15.0
#23 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:48:33.1
#24 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:48:51.3
#25 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:49:35.8
#23 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:49:53.9
#23 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:50:12.0
#27 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:50:56.6
#26 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:51:14.6
#25 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:51:32.7
#27 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:53:05.2
#22 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:53:23.4
#29 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:53:41.1
#26 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:54:26.4
#25 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:54:44.0
#26 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:55:02.2
#23 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:55:47.0
#23 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:56:04.7
#27 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:56:22.8
#24 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:57:07.7
#27 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:57:25.3
#26 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:57:43.4
#25 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:58:28.2
#29 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:58:45.9
#25 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:59:04.1
#25 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:59:49.0
#23 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:00:06.6
#25 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:00:24.7
#24 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:01:09.6
#22 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:01:39.6
#20 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:03:45.3
#24 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:04:03.4
#23 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:04:48.1
#27 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:05:05.8
#26 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:05:23.8
#24 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:06:08.5
#26 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:06:26.4
#24 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:06:44.5
#20 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:07:29.0
#23 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:07:46.9
#23 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:08:05.0
#19 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:08:49.4
#23 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:09:07.4
#21 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:09:25.6
#21 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:10:09.9
#22 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:10:28.0
#24 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:10:46.1
#22 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:11:30.4
#25 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:11:48.5
#23 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:12:06.7
#22 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:12:39.6
#21 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:12:57.7
#20 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:14:20.2
#26 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:14:38.4
#21 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:27:51.6
#21 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:28:09.7
#24 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:28:54.3
#21 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:29:11.8
#22 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:29:30.0
#26 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:30:14.1
#23 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:30:32.1
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:30:50.2
#22 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:31:34.3
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:31:52.3
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:32:10.5
#22 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:32:55.2
#23 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:33:12.6
#25 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:33:30.7
#22 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:34:15.1
#21 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:34:32.8
#23 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:34:51.0
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:35:35.3
#24 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:35:53.1
#21 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:36:11.2
#18 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:36:43.8
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:37:01.9
#24 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:38:27.8
#23 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:38:46.0
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:39:03.9
#23 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:39:48.8
#24 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:40:06.2
#25 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:40:24.3
#24 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:41:08.8
#22 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:41:26.5
#22 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:42:51.4
#22 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:43:09.2
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:43:27.3
#22 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:44:11.6
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:44:29.4
#21 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:44:47.6
#21 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:45:31.8
#17 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:45:49.7
#21 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:46:07.8
#22 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:46:52.3
#22 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:47:09.9
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:47:28.1
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:48:59.8
#26 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:49:17.9
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:49:36.0
#21 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:50:20.2
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:50:38.2
#22 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:50:56.3
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:51:40.5
#18 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:51:58.4
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:52:16.6
#22 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:53:01.1
#24 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:53:18.7
#22 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:53:36.8
#21 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:54:21.1
#21 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:54:38.9
#22 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:54:57.1
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:55:41.4
#22 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:55:59.2
#24 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:56:17.3
#21 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:57:14.4
#23 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:57:32.5
#22 of 43 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:59:38.0
#20 of 43 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:59:56.1
#18 of 43 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:00:40.4
#18 of 43 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:00:58.2
#20 of 43 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:01:16.4
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:02:00.7
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:02:18.5
#22 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:02:36.6
#23 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:03:20.7
#23 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:03:38.7
#22 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:03:56.9
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:04:41.1
#22 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:04:59.0
#18 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:05:17.1
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:06:01.2
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:06:19.3
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:06:37.4
#22 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:07:21.5
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:07:39.5
#21 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:07:57.7
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:08:31.3
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:08:49.4
#21 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:10:15.5
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:10:33.6
#24 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:10:51.8
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:11:36.3
#17 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:11:53.9
#23 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:12:12.1
#23 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:12:56.2
#21 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:13:14.2
#21 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:14:43.3
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:15:01.2
#21 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:15:45.5
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:16:03.3
#23 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:16:21.5
#22 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:17:05.7
#23 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:17:23.6
#24 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:17:41.7
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:18:26.1
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:18:43.8
#21 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:19:02.0
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:20:36.2
#18 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:20:54.3
#22 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:21:12.4
#18 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:21:57.0
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:22:14.6
#25 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:22:32.7
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:23:17.0
#23 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:23:34.8
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:23:53.0
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:24:37.2
#23 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:24:55.1
#23 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:25:13.3
#21 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:25:57.6
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:26:15.4
#22 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:26:33.5
#22 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:27:17.7
#21 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:27:35.6
#21 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:27:53.8
#22 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:28:38.0
#23 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:29:09.5
#22 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:31:16.7
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:31:34.8
#22 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:32:19.5
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:32:37.0
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:32:55.1
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:33:39.5
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:33:57.2
#23 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:34:15.4
#18 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:34:59.6
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:35:17.5
#21 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:35:35.6
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:36:20.3
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:36:37.7
#21 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:36:55.9
#22 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:37:40.2
#22 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:37:58.0
#22 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:38:16.1
#21 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:39:00.3
#22 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:39:18.2
#24 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:39:36.4
#22 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:40:08.7
#21 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:40:26.9
#22 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:41:49.2
#23 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:29:46.8
#25 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:30:04.9
#25 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:30:50.3
#27 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:31:08.4
#26 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:31:26.6
#26 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:32:12.0
#30 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:32:30.1
#25 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:32:48.3
#28 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:33:33.9
#24 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:33:51.8
#32 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:34:10.0
#24 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:34:55.5
#27 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:35:13.5
#26 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:35:31.7
#27 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:36:17.1
#23 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:36:35.2
#23 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:36:53.4
#29 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:37:38.8
#26 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:37:56.9
#21 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:38:15.1
#26 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:38:47.4
#22 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:40:13.3
#25 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:40:31.4
#24 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:40:49.6
#26 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:41:34.9
#29 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:41:53.0
#24 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:42:11.2
#29 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:42:56.5
#25 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:43:14.6
#31 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:43:32.7
#27 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:44:59.1
#23 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:45:17.2
#30 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:46:02.5
#26 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:46:20.7
#21 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:46:38.8
#24 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:47:24.2
#25 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:47:42.3
#26 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:48:00.4
#25 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:48:45.7
#26 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:49:03.9
#26 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:49:22.0
#26 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:50:56.5
#23 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:51:14.6
#22 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:51:32.8
#27 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:52:18.1
#23 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:52:36.2
#25 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:52:54.4
#24 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:53:39.7
#27 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:53:57.8
#26 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:54:16.0
#29 of 42 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:55:01.3
#22 of 42 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:55:19.4
#18 of 42 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:55:37.6
#25 of 42 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:56:22.9
#22 of 42 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:56:41.0
#23 of 42 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:56:59.2
#35 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:57:44.5
#25 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:58:02.6
#27 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:58:20.8
#29 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:59:06.1
#25 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:59:38.3
#24 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:01:46.8
#22 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:02:05.0
#21 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:02:50.5
#25 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:03:08.5
#25 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:03:26.6
#24 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:04:11.9
#24 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:04:30.0
#31 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:04:48.2
#27 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:05:33.5
#26 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:05:51.7
#25 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:06:09.8
#23 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:06:55.2
#22 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:07:13.2
#27 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:07:31.4
#22 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:08:16.7
#31 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:08:34.8
#21 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:08:53.0
#17 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:09:38.6
#24 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:09:56.5
#23 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:10:14.6
#29 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:10:47.5
#25 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:12:14.0
#22 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:12:32.1
#23 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:12:50.3
#25 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:13:35.6
#25 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:13:53.6
#27 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:14:11.8
#24 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:14:57.0
#27 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:15:15.1
#23 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:16:41.0
#22 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:16:59.1
#24 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:17:44.5
#27 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:18:02.5
#23 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:18:20.6
#30 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:19:05.9
#25 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:19:24.0
#27 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:19:42.1
#22 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:20:27.4
#26 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:20:45.5
#25 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:21:03.6
#26 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:22:38.6
#26 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:22:56.8
#26 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:23:14.9
#33 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:24:00.1
#23 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:24:18.3
#27 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:24:36.4
#25 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:25:21.7
#20 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:25:39.8
#29 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:25:57.9
#23 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:26:43.2
#28 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:27:01.3
#33 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:27:19.4
#32 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:28:04.7
#28 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:28:22.8
#25 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:28:40.9
#29 of 41 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:29:26.2
#21 of 41 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:29:44.3
#25 of 41 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:30:02.5
#22 of 38 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:30:47.9
#25 of 38 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:31:18.7
#22 of 34 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:33:25.2
#26 of 34 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:33:43.3
#22 of 35 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:34:28.5
#25 of 35 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:34:46.6
#24 of 35 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:35:04.8
#21 of 35 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:35:49.9
#23 of 35 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:36:08.0
#27 of 35 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:36:26.2
#23 of 36 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:37:11.3
#17 of 36 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:37:29.4
#22 of 36 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:37:47.6
#24 of 36 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:38:32.7
#24 of 36 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:38:50.8
#28 of 36 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:39:09.0
#23 of 35 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:39:54.1
#27 of 35 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:40:12.2
#27 of 35 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:40:30.4
#21 of 35 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:41:15.7
#23 of 35 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:41:33.6
#20 of 35 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:41:51.8
#23 of 35 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:42:25.2
#23 of 35 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:43:49.5
#27 of 35 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:44:07.7
#26 of 35 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:44:25.8
#22 of 34 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:45:10.9
#26 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:27:18.4
#25 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:27:36.5
#22 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:28:21.8
#22 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:28:40.0
#26 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:28:58.1
#25 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:29:43.4
#26 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:30:01.6
#26 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:30:19.7
#26 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:31:05.1
#27 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:31:23.2
#26 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:31:41.3
#28 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:32:26.6
#26 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:32:44.8
#22 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:33:02.9
#25 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:33:48.2
#24 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:34:06.4
#28 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:34:24.5
#26 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:35:09.8
#19 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:35:28.0
#38 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:35:46.1
#20 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:36:18.8
#24 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:37:46.5
#29 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:38:04.7
#27 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:38:22.8
#29 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:39:08.1
#22 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:39:26.2
#27 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:39:44.2
#27 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:40:29.7
#28 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:40:47.9
#23 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:41:06.0
#25 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:42:33.9
#25 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:42:52.0
#22 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:43:37.3
#27 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:43:55.5
#29 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:44:13.6
#27 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:44:59.1
#25 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:45:17.1
#26 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:45:35.2
#22 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:46:20.6
#24 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:46:38.7
#29 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:46:56.8
#23 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:48:30.0
#28 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:48:48.2
#22 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:49:06.3
#25 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:49:51.8
#30 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:50:09.8
#24 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:50:27.9
#23 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:51:13.3
#26 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:51:31.4
#27 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:51:49.5
#20 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:52:34.8
#24 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:52:53.0
#22 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:53:11.1
#24 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:53:56.4
#22 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:54:14.6
#26 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:54:32.7
#26 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:55:18.0
#25 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:55:36.2
#26 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:55:54.3
#25 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:56:39.6
#28 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:57:11.9
#28 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:59:18.3
#26 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:59:36.4
#30 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:00:21.7
#23 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:00:39.8
#28 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:00:57.9
#29 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:01:43.1
#29 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:02:01.3
#28 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:02:19.4
#27 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:03:04.6
#29 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:03:22.8
#29 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:03:40.9
#30 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:04:26.1
#27 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:04:44.3
#27 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:05:02.4
#27 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:05:47.8
#28 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:06:05.8
#26 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:06:23.9
#26 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:07:09.2
#27 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:07:27.3
#25 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:07:45.4
#32 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:08:18.0
#25 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:09:44.1
#28 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:10:02.3
#27 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:10:20.4
#27 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:11:05.5
#30 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:11:23.7
#20 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:11:41.8
#28 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:12:26.9
#26 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:12:45.1
#27 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:13:03.2
#37 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:14:29.1
#32 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:14:47.3
#26 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:15:32.4
#31 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:15:50.5
#28 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:16:08.7
#24 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:16:53.8
#25 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:17:12.0
#30 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:17:30.1
#25 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:18:15.2
#31 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:18:33.4
#27 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:18:51.5
#30 of 44 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:20:26.9
#33 of 44 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:20:45.1
#20 of 44 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:21:03.2
#31 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:21:48.3
#30 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:22:06.5
#34 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:22:24.6
#28 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:23:09.7
#29 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:23:27.9
#28 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:23:46.0
#29 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:24:31.2
#34 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:24:49.3
#29 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:25:07.4
#33 of 44 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:25:52.6
#28 of 44 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:26:10.7
#28 of 44 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:26:28.8
#29 of 44 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:27:14.0
#32 of 44 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:27:32.1
#29 of 44 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:27:50.3
#25 of 44 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:28:35.4
#23 of 44 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:29:07.1
#25 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:31:15.7
#29 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:31:33.8
#26 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:32:18.9
#24 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:32:37.0
#32 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:32:55.1
#34 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:33:40.2
#28 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:33:58.3
#27 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:34:16.5
#32 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:35:01.5
#29 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:35:19.6
#30 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:35:37.8
#30 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:36:22.8
#28 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:36:40.9
#29 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:36:59.1
#29 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:37:44.1
#30 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:38:02.3
#28 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:38:20.4
#27 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:39:05.7
#26 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:39:23.6
#30 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:39:41.7
#33 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:40:14.9
#28 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:41:39.4
#29 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:41:57.5
#26 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:42:15.6
#28 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:42:54.0
#25 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:43:12.1
#24 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:43:56.9
#26 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:44:15.1
#27 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:44:33.2
#22 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:45:18.0
#24 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:45:36.2
#26 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:45:54.3
#24 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:46:39.5
#30 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:46:57.3
#26 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:47:15.5
#25 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:48:00.3
#25 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:48:18.4
#22 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:48:36.6
#23 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:49:21.5
#24 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:49:39.6
#24 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:49:57.6
#30 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:50:42.5
#26 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:51:00.7
#22 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:51:18.8
#26 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:51:53.9
#24 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:53:24.3
#23 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:53:42.4
#28 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:54:00.6
#26 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:54:45.6
#24 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:55:03.5
#24 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:55:21.6
#26 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:56:06.3
#24 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:56:24.5
#26 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:57:52.1
#26 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:58:10.0
#22 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:58:54.8
#26 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:59:12.9
#20 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:59:31.1
#23 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:00:15.9
#25 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:00:33.9
#24 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:00:52.1
#25 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:01:36.9
#29 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:01:55.0
#27 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:02:13.1
#28 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:03:52.6
#25 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:04:10.7
#28 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:04:28.9
#27 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:05:13.5
#24 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:05:31.7
#25 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:05:49.8
#26 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:06:34.5
#28 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:06:52.6
#25 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:07:10.7
#26 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:07:55.5
#26 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:08:13.5
#22 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:08:31.7
#27 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:09:16.5
#30 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:09:34.5
#30 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:09:52.4
#25 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:10:37.2
#28 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:10:55.4
#28 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:11:13.5
#29 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:12:13.2
#27 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:12:31.3
#29 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:14:45.1
#27 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:15:03.2
#21 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:15:47.9
#31 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:16:05.9
#23 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:16:24.0
#23 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:17:08.9
#25 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:17:26.7
#26 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:17:44.9
#23 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:18:29.5
#30 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:18:47.6
#27 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:19:05.7
#28 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:19:50.3
#26 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:20:08.4
#25 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:20:26.5
#26 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:21:11.1
#30 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:21:29.2
#26 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:21:47.4
#29 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:22:31.9
#24 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:22:50.1
#25 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:23:08.2
#26 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:23:42.3
#24 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:25:12.2
#28 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:25:30.3
#26 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:25:48.5
#27 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:26:33.2
#22 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:26:51.1
#23 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:27:09.2
#29 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:27:54.0
#24 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:28:11.8
#24 of 43 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:29:39.1
#30 of 43 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:29:57.2
#27 of 41 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:30:41.6
#27 of 41 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:30:59.7
#25 of 41 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:31:17.8
#23 of 41 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:32:02.3
#25 of 41 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:32:20.4
#27 of 41 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:32:38.5
#24 of 40 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:33:23.0
#28 of 40 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:33:41.0
#23 of 40 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:33:59.1
#22 of 38 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:35:38.7
#30 of 38 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:35:56.8
#25 of 38 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:36:15.0
#25 of 39 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:36:59.3
#30 of 39 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:37:17.5
#24 of 39 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:37:35.6
#28 of 39 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:38:20.0
#18 of 39 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:38:38.1
#25 of 39 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:38:56.3
#28 of 39 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:39:40.6
#29 of 39 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:39:58.8
#31 of 39 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:40:16.9
#27 of 37 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:41:01.3
#27 of 37 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:41:19.4
#20 of 37 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:41:37.5
#27 of 37 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:42:21.9
#22 of 37 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:42:40.0
#26 of 37 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:42:58.2
#28 of 40 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:43:42.6
#27 of 40 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:44:15.7
#28 of 40 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:46:29.3
#27 of 40 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:46:47.5
#23 of 40 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:47:31.8
#31 of 40 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:47:49.8
#28 of 40 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:48:07.9
#23 of 39 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:48:52.2
#24 of 39 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:49:10.2
#23 of 39 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:49:28.4
#24 of 38 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:50:12.6
#24 of 38 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:50:30.7
#27 of 38 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:50:48.8
#30 of 37 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:51:33.5
#24 of 37 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:51:51.1
#26 of 37 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:52:09.3
#27 of 37 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:52:53.5
#24 of 37 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:53:11.6
#24 of 37 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:53:29.7
#23 of 38 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:54:14.0
#26 of 38 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:54:32.0
#23 of 38 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:54:50.2
#26 of 38 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:55:25.6
#22 of 37 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:56:54.0
#23 of 37 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:57:12.1
#21 of 37 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:57:30.3
#23 of 34 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:58:14.7
#30 of 34 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:58:32.5
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:39:37.5
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:39:55.6
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:40:13.8
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:40:31.9
# 6 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:40:50.1
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:41:08.2
# 6 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:41:26.3
# 6 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:41:44.5
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:42:02.6
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:42:20.8
# 7 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:42:38.9
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:42:57.1
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:43:15.2
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:43:33.4
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:43:51.5
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:44:09.6
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:44:28.3
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:44:47.5
# 6 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:45:05.6
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:45:23.8
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:46:57.8
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:47:16.0
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:47:34.1
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:47:52.3
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:48:10.4
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:48:28.5
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:48:46.7
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:49:04.8
# 6 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:49:23.0
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:49:41.1
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:49:59.3
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:50:17.4
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:50:35.6
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:50:53.7
# 7 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:51:11.8
# 6 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:51:30.0
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:51:48.6
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:52:07.8
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:52:26.0
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:52:44.1
# 7 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:53:02.2
# 7 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:53:17.4
# 7 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:54:41.5
# 6 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:54:59.6
# 6 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:55:17.7
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:55:35.9
# 8 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:55:54.0
# 6 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:58:03.8
# 6 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:58:21.9
# 7 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:58:40.1
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:58:58.2
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:59:16.4
# 6 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:59:34.5
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:59:52.7
# 6 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:00:10.8
# 7 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:00:28.9
# 6 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:00:47.1
# 6 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:01:05.2
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:01:23.4
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:01:41.5
# 7 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:01:59.7
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:02:17.8
# 6 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:02:36.0
# 6 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:02:54.6
# 9 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:03:13.8
# 9 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:03:31.9
# 8 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:03:50.1
# 9 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:04:08.2
#10 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:04:23.3
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:08:40.0
# 7 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:08:58.1
# 8 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:09:16.3
# 7 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:09:34.4
# 6 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:09:52.6
# 7 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:10:10.7
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:10:28.8
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:10:47.0
# 6 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:11:05.1
# 6 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:11:23.3
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:11:41.4
# 6 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:11:59.6
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:12:17.7
# 7 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:12:35.9
# 6 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:12:54.0
# 7 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:13:12.1
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:13:30.8
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:13:50.0
# 6 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:14:08.1
# 9 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:14:26.3
# 9 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:14:44.4
#10 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:14:59.5
# 6 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:17:01.5
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:17:19.7
# 7 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:17:37.8
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:17:56.0
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:18:14.1
# 6 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:18:32.3
# 8 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:18:50.4
# 7 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:19:08.5
# 9 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:19:26.7
# 8 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:19:44.8
# 8 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:20:03.0
# 7 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:20:21.1
# 7 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:20:39.3
# 7 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:20:57.4
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:21:15.6
# 6 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:21:33.7
# 6 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:21:52.4
# 6 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:22:11.5
# 7 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:22:29.7
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:22:47.8
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:23:06.0
#10 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:23:21.1
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:25:20.5
# 6 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:25:38.6
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:25:56.8
#12 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:26:08.9
# 8 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:27:28.8
# 7 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:27:46.9
# 7 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:28:05.1
# 8 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:28:23.2
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:28:41.4
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:28:59.5
# 7 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:29:17.7
# 6 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:29:35.8
# 6 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:29:53.9
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:30:12.1
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:30:30.2
# 9 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:30:45.4
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:29:30.4
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:29:48.6
# 4 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:30:24.9
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:30:43.0
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:31:01.2
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:31:19.3
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:31:55.6
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:32:13.7
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:32:31.9
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:32:50.0
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:33:44.5
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:34:02.6
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:34:21.3
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:34:40.3
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:34:58.6
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:35:16.9
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:35:34.7
# 3 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:35:50.1
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:37:27.2
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:37:45.4
# 3 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:38:03.5
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:38:39.8
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:39:16.1
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:39:34.2
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:39:52.2
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:40:10.5
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:40:28.7
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:40:46.7
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:41:04.9
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:41:23.1
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:41:41.2
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:41:59.9
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:42:19.1
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:42:37.2
# 3 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:42:55.7
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:43:13.5
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:43:28.5
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:44:27.9
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:44:46.0
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:45:22.3
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:47:02.3
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:47:38.5
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:47:56.7
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:48:14.8
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:48:33.0
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:48:51.1
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:49:09.3
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:49:27.3
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:49:45.6
# 5 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:50:03.7
# 3 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:50:21.9
# 4 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:50:39.7
# 3 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:50:58.0
# 3 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:51:16.5
# 4 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:51:34.6
# 3 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:51:53.2
# 3 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:52:12.2
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:52:30.4
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:52:48.5
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:53:06.7
# 5 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:53:21.7
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:53:38.4
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:53:56.4
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:54:14.4
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:54:32.4
# 3 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:54:51.0
# 3 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:57:57.6
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:58:15.7
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:58:33.9
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:59:10.2
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:59:28.3
# 3 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:59:46.5
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:00:04.7
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:00:22.6
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:00:40.7
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:00:59.6
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:01:17.3
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:01:35.1
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:01:53.7
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:02:11.6
# 3 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:02:29.8
# 4 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:02:48.5
# 3 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:03:07.6
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:03:25.5
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:03:43.9
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:04:02.2
# 3 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:04:17.3
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:05:55.0
# 3 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:06:13.1
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:06:31.3
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:06:49.4
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:07:07.6
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:07:25.7
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:07:43.5
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:08:02.5
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:08:20.6
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:08:38.2
# 3 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:08:56.4
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:09:14.6
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:09:32.7
# 3 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:09:50.9
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:10:09.4
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:10:46.5
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:11:05.2
# 3 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:11:23.7
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:11:41.1
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:11:56.4
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:12:55.9
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:13:14.0
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:13:32.2
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:13:50.3
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:14:08.4
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:15:31.4
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:15:49.5
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:16:07.7
# 3 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:16:19.8
# 3 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/03:56:43.5
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/03:57:01.7
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/03:57:19.8
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/03:57:38.0
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/03:57:56.1
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/03:58:14.3
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/03:58:32.4
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/03:58:50.5
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/03:59:08.7
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/03:59:26.8
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/03:59:45.0
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:00:03.1
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:00:21.3
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:00:39.4
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:00:57.5
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:01:15.7
# 3 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:01:34.3
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:01:53.5
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:02:11.7
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:02:29.8
# 4 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:02:47.9
# 4 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:03:03.1
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:04:22.9
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:04:41.1
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:04:59.2
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:05:17.4
# 3 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:05:35.5
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:05:53.7
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:06:11.8
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:06:29.9
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:06:48.1
# 5 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:07:06.2
# 3 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:07:24.4
# 3 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:07:42.5
# 5 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:08:00.7
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:08:18.8
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:08:37.0
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:08:55.1
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:09:13.8
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:09:32.9
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:09:51.1
# 3 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:10:09.2
# 4 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:10:27.4
# 4 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:10:42.5
# 3 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:11:41.0
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:11:59.2
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:12:17.3
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:12:35.5
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:12:53.6
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:14:16.6
# 3 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:14:34.8
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:14:52.9
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:15:11.0
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:15:29.2
# 4 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:15:47.3
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:16:05.5
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:16:23.6
# 3 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:16:41.8
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:16:59.9
# 7 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:17:18.0
# 3 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:17:36.2
# 3 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:17:54.3
# 3 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:18:12.5
# 3 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:18:30.6
# 4 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:18:48.8
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:19:07.4
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:19:26.6
# 3 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:19:44.7
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:20:02.9
# 4 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:20:36.1
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:20:53.0
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:21:10.9
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:21:29.1
# 3 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:21:47.1
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:22:05.4
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:25:12.3
# 3 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:25:30.5
# 3 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:25:48.6
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:26:06.8
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:26:24.9
# 3 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:26:43.1
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:27:01.2
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:27:19.3
# 3 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:27:37.5
# 3 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:28:13.8
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:28:31.9
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:28:50.1
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:29:08.2
# 3 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:29:26.4
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:29:44.5
# 3 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:30:03.0
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:30:22.3
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:30:40.5
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:30:58.6
# 3 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:31:16.5
# 4 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:31:31.8
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:32:51.4
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:33:09.6
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:33:27.7
# 3 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:33:45.9
# 3 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:34:04.0
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:34:22.2
# 4 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:34:40.3
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:34:58.5
# 3 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:35:16.6
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:35:34.7
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:35:52.9
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:36:11.0
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:36:29.2
# 3 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:36:47.3
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:37:05.4
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:37:23.6
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:37:42.3
# 5 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:38:01.4
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:38:19.9
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:38:37.9
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:38:55.9
# 5 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:39:11.1
# 3 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:40:08.3
# 4 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:40:26.5
# 3 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:40:44.6
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:41:02.8
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:41:20.9
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:43:00.9
# 3 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:43:19.0
# 7 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:43:31.1
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/03:44:00.1
# 2 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/03:44:18.3
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/03:44:36.3
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/03:44:54.5
# 2 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/03:45:30.8
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/03:46:07.1
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/03:47:19.6
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/03:47:37.8
# 2 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/03:47:55.9
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/03:48:14.1
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/03:48:32.2
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/03:48:50.9
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/03:49:46.3
# 2 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/03:50:04.5
# 2 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/03:51:38.5
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/03:51:56.6
# 2 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/03:52:14.8
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/03:52:32.9
# 2 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/03:54:03.7
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/03:54:21.8
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/03:54:39.9
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/03:55:34.4
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/03:55:52.5
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/03:56:10.7
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/03:57:06.6
# 2 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/03:57:42.9
# 2 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/03:57:58.0
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/03:58:56.9
# 2 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/03:59:15.1
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/04:01:33.2
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/04:01:51.4
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/04:02:09.5
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/04:03:58.4
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/04:05:11.0
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/04:05:47.2
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/04:06:05.4
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/04:06:43.2
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/04:07:19.5
# 2 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/04:07:52.8
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/04:08:09.6
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/04:08:45.5
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/04:09:03.7
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/04:09:22.1
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/04:12:46.9
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/04:13:05.1
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/04:13:59.5
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/04:14:17.7
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/04:14:54.0
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/04:15:48.4
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/04:16:24.7
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/04:17:19.6
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/04:17:38.9
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/04:17:56.5
# 2 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/04:18:15.1
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/04:20:25.2
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/04:21:37.7
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/04:21:55.9
# 2 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/04:24:20.8
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/04:24:57.6
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/04:25:16.9
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/04:25:35.4
# 2 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/04:26:26.5
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/04:28:18.8
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/04:28:37.0
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/04:30:16.4
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/04:30:34.6
# 2 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/04:30:46.7

#Inspect gain tables interactively and decide whether to manually flag something
#plotms(LB_p6,xaxis='time',yaxis='GainPhase',iteraxis='spw',showgui=True)
#Print calibration png file
#plotms(
#    LB_p6,xaxis='time', yaxis='GainPhase',iteraxis='spw',exprange='all',
#    overwrite=True,showgui=False,
#    plotfile=os.path.join(LB_selfcal_folder,f'{prefix}_LB_p6_gain_phase_vs_time.png')
#)

#Inspect gain tables and decide whether to flag something
plotms(
    LB_p6,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(LB_selfcal_folder,f'{prefix}_LB_p6_gain_phase_vs_time.png'),
)

#Flag problematic antennas
flagdata(vis=LB_p6,mode='manual',spw='0', timerange='2017/10/10/09:01:20~2017/10/10/09:01:30',antenna='DV20')
flagdata(vis=LB_p6,mode='manual',spw='0', timerange='2017/10/10/09:01:35~2017/10/10/09:01:45',antenna='DA60')
flagdata(vis=LB_p6,mode='manual',spw='0', timerange='2017/10/10/09:03:00~2017/10/10/09:03:10',antenna='DA52,DA64')
flagdata(vis=LB_p6,mode='manual',spw='0', timerange='2017/10/10/09:03:10~2017/10/10/09:03:20',antenna='DV17')
flagdata(vis=LB_p6,mode='manual',spw='0', timerange='2017/10/10/09:04:00~2017/10/10/09:04:10',antenna='DA52')
flagdata(vis=LB_p6,mode='manual',spw='0', timerange='2017/10/10/09:07:00~2017/10/10/09:07:10',antenna='DA52')
flagdata(vis=LB_p6,mode='manual',spw='0', timerange='2017/10/10/09:07:20~2017/10/10/09:07:30',antenna='DV04')
flagdata(vis=LB_p6,mode='manual',spw='0', timerange='2017/10/10/09:11:00~2017/10/10/09:11:10',antenna='DA63,DV20')
flagdata(vis=LB_p6,mode='manual',spw='0', timerange='2017/10/10/09:15:20~2017/10/10/09:15:30',antenna='DA54')
flagdata(vis=LB_p6,mode='manual',spw='0', timerange='2017/10/10/09:15:40~2017/10/10/09:15:50',antenna='DA56,DV05')
flagdata(vis=LB_p6,mode='manual',spw='0', timerange='2017/10/10/09:18:05~2017/10/10/09:18:15',antenna='DA48')
flagdata(vis=LB_p6,mode='manual',spw='0', timerange='2017/10/10/09:21:20~2017/10/10/09:21:30',antenna='DV12')
flagdata(vis=LB_p6,mode='manual',spw='0', timerange='2017/10/10/09:22:40~2017/10/10/09:22:50',antenna='DA47')
flagdata(vis=LB_p6,mode='manual',spw='0', timerange='2017/10/10/09:23:00~2017/10/10/09:23:10',antenna='DV12,DV16')
flagdata(vis=LB_p6,mode='manual',spw='0', timerange='2017/10/10/09:23:15~2017/10/10/09:23:25',antenna='DA65')
flagdata(vis=LB_p6,mode='manual',spw='0', timerange='2017/10/10/09:26:40~2017/10/10/09:26:50',antenna='DV20')
flagdata(vis=LB_p6,mode='manual',spw='0', timerange='2017/10/10/09:28:00~2017/10/10/09:28:10',antenna='DA60')
flagdata(vis=LB_p6,mode='manual',spw='0', timerange='2017/10/10/09:28:40~2017/10/10/09:28:50',antenna='DA45,DV04')
flagdata(vis=LB_p6,mode='manual',spw='0', timerange='2017/10/10/09:32:20~2017/10/10/09:32:30',antenna='DA48,DA60,PM04')
flagdata(vis=LB_p6,mode='manual',spw='0', timerange='2017/10/10/09:34:20~2017/10/10/09:34:30',antenna='DA54,DA58')
flagdata(vis=LB_p6,mode='manual',spw='0', timerange='2017/10/10/09:34:40~2017/10/10/09:34:50',antenna='DV05,DV08')
flagdata(vis=LB_p6,mode='manual',spw='0', timerange='2017/10/10/09:35:00~2017/10/10/09:35:10',antenna='DA56')
flagdata(vis=LB_p6,mode='manual',spw='0', timerange='2017/10/10/09:35:45~2017/10/10/09:35:55',antenna='DA57,DV03')
flagdata(vis=LB_p6,mode='manual',spw='0', timerange='2017/10/10/09:38:45~2017/10/10/09:38:55',antenna='DA58,DV12')
flagdata(vis=LB_p6,mode='manual',spw='0', timerange='2017/10/10/09:41:20~2017/10/10/09:41:30',antenna='DA52,DV20')
flagdata(vis=LB_p6,mode='manual',spw='0', timerange='2017/10/10/09:49:50~2017/10/10/09:50:00',antenna='DA65')
flagdata(vis=LB_p6,mode='manual',spw='0', timerange='2017/10/10/09:53:20~2017/10/10/09:53:25',antenna='DA64')
flagdata(vis=LB_p6,mode='manual',spw='0', timerange='2017/10/10/10:01:30~2017/10/10/10:01:40',antenna='DV22')
flagdata(vis=LB_p6,mode='manual',spw='0', timerange='2017/10/10/10:03:40~2017/10/10/10:03:50',antenna='DA57,DV08,DV17')
flagdata(vis=LB_p6,mode='manual',spw='0', timerange='2017/10/10/10:05:20~2017/10/10/10:05:30',antenna='DA57')
flagdata(vis=LB_p6,mode='manual',spw='0', timerange='2017/10/10/10:09:20~2017/10/10/10:09:30',antenna='DV04,DV24')
flagdata(vis=LB_p6,mode='manual',spw='0', timerange='2017/10/10/10:10:00~2017/10/10/10:10:10',antenna='DA47,DA52')
flagdata(vis=LB_p6,mode='manual',spw='0', timerange='2017/10/10/10:12:50~2017/10/10/10:13:00',antenna='DV20')

flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/10:27:50~2017/10/10/10:28:00',antenna='DV22')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/10:28:50~2017/10/10/10:29:00',antenna='DA54')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/10:30:50~2017/10/10/10:31:00',antenna='DA47,DV08')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/10:31:30~2017/10/10/10:31:40',antenna='DA54,DV24')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/10:32:10~2017/10/10/10:32:20',antenna='DA58,DA63')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/10:32:50~2017/10/10/10:33:00',antenna='DV24')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/10:34:10~2017/10/10/10:34:20',antenna='DV17,DV24')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/10:34:30~2017/10/10/10:34:40',antenna='DV04')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/10:35:30~2017/10/10/10:35:40',antenna='DA65,DV24')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/10:37:00~2017/10/10/10:37:10',antenna='DA56,DA60')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/10:41:20~2017/10/10/10:41:30',antenna='DA63')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/10:44:40~2017/10/10/10:44:50',antenna='DA54')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/10:45:45~2017/10/10/10:45:55',antenna='DA52')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/10:46:00~2017/10/10/10:46:10',antenna='DV08')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/10:47:20~2017/10/10/10:47:30',antenna='DV16')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/10:50:15~2017/10/10/10:50:25',antenna='DV15')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/10:50:20~2017/10/10/10:50:30',antenna='DA44')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/10:50:30~2017/10/10/10:50:40',antenna='DA44,DV15')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/10:50:50~2017/10/10/10:51:00',antenna='DA44,DV15')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/10:51:40~2017/10/10/10:51:50',antenna='DA44')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/11:03:50~2017/10/10/11:04:00',antenna='DA59')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/11:04:35~2017/10/10/11:04:45',antenna='DA48')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/11:08:30~2017/10/10/11:08:40',antenna='DA48')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/11:10:10~2017/10/10/11:10:20',antenna='DA63')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/11:10:30~2017/10/10/11:10:40',antenna='DV04')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/11:11:30~2017/10/10/11:11:40',antenna='DA43,DA44,DA60')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/11:11:50~2017/10/10/11:12:00',antenna='DA43,DA44,DA48,DV13,DV17,PM04')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/11:12:10~2017/10/10/11:12:20',antenna='DA53,DV13')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/11:12:50~2017/10/10/11:13:00',antenna='DA44,DV22')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/11:13:10~2017/10/10/11:13:20',antenna='DA53')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/11:17:00~2017/10/10/11:17:10',antenna='DA60')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/11:17:20~2017/10/10/11:17:30',antenna='DA48,DA60')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/11:18:40~2017/10/10/11:18:50',antenna='DA63')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/11:20:30~2017/10/10/11:20:40',antenna='DA60,DA65,DV12')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/11:20:50~2017/10/10/11:21:00',antenna='DV04,DV24')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/11:21:50~2017/10/10/11:22:00',antenna='DA60,DA65')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/11:23:10~2017/10/10/11:23:20',antenna='DA64,DV05')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/11:28:30~2017/10/10/11:28:40',antenna='DA54,DA65')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/11:31:30~2017/10/10/11:31:40',antenna='DA63,DA65')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/11:32:10~2017/10/10/11:32:20',antenna='DA47')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/11:32:30~2017/10/10/11:32:40',antenna='DA63')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/11:32:50~2017/10/10/11:33:30',antenna='DV20')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/11:33:50~2017/10/10/11:34:00',antenna='DA45,DV17')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/11:34:50~2017/10/10/11:35:00',antenna='DA45,DV17')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/11:35:30~2017/10/10/11:35:40',antenna='DV05')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/11:36:15~2017/10/10/11:36:25',antenna='DV08,DV22')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/11:36:30~2017/10/10/11:36:40',antenna='DA47,DA64,DV22')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/11:36:50~2017/10/10/11:37:00',antenna='DV16,DV17,DV22')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/11:37:35~2017/10/10/11:37:45',antenna='DA48')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/11:39:10~2017/10/10/11:39:20',antenna='DA48')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/11:39:30~2017/10/10/11:39:40',antenna='DV17')
flagdata(vis=LB_p6,mode='manual',spw='4', timerange='2017/10/10/11:40:00~2017/10/10/11:40:10',antenna='DA63,DV05')

flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/07:30:00~2017/10/11/07:30:10',antenna='DV04')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/07:30:45~2017/10/11/07:30:55',antenna='DA52,DA64')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/07:31:00~2017/10/11/07:31:10',antenna='DA47')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/07:31:20~2017/10/11/07:31:30',antenna='DV16,DV20')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/07:32:45~2017/10/11/07:32:55',antenna='DA43')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/07:33:30~2017/10/11/07:33:40',antenna='DA43,DA44')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/07:33:50~2017/10/11/07:34:00',antenna='DA43,DV04,DV15')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/07:34:50~2017/10/11/07:35:00',antenna='DA44,DA65,DV16')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/07:35:10~2017/10/11/07:35:20',antenna='DA60,DV16')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/07:35:30~2017/10/11/07:35:40',antenna='PM04')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/07:36:10~2017/10/11/07:36:20',antenna='DA64,DV23')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/07:36:30~2017/10/11/07:36:40',antenna='DA64,DV15')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/07:36:50~2017/10/11/07:37:00',antenna='DA56')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/07:37:30~2017/10/11/07:37:40',antenna='DV15')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/07:37:50~2017/10/11/07:38:00',antenna='DA44')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/07:38:10~2017/10/11/07:38:20',antenna='DV08')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/07:38:40~2017/10/11/07:38:50',antenna='DV04')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/07:40:10~2017/10/11/07:40:20',antenna='DV03')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/07:40:30~2017/10/11/07:40:40',antenna='DA63,DV08')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/07:41:50~2017/10/11/07:42:00',antenna='DA54,DA60')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/07:42:50~2017/10/11/07:43:00',antenna='DA58,DV04,PM04')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/07:46:30~2017/10/11/07:46:40',antenna='DA63,DV01')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/07:47:20~2017/10/11/07:47:30',antenna='DV12')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/07:47:40~2017/10/11/07:47:50',antenna='DV16')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/07:49:00~2017/10/11/07:49:10',antenna='DA65')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/07:52:30~2017/10/11/07:52:40',antenna='DA59')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/07:53:50~2017/10/11/07:54:00',antenna='DV03')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/07:54:10~2017/10/11/07:54:20',antenna='DA63,DV10')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/07:55:00~2017/10/11/07:55:10',antenna='DV08,DV22')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/07:55:30~2017/10/11/07:55:40',antenna='DA47,DV12')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/07:56:40~2017/10/11/07:56:50',antenna='DA47,DA64')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/07:57:40~2017/10/11/07:57:50',antenna='DV11,PM04')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/07:59:00~2017/10/11/07:59:10',antenna='DA58')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/07:59:30~2017/10/11/07:59:40',antenna='DA64,DV22')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/08:04:10~2017/10/11/08:04:20',antenna='DA58')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/08:06:50~2017/10/11/08:07:00',antenna='DV20')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/08:08:30~2017/10/11/08:08:40',antenna='DA45,DA63')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/08:08:50~2017/10/11/08:09:00',antenna='DA56,DV04,DV11')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/08:09:30~2017/10/11/08:09:40',antenna='DA54,DV11,DV15')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/08:09:50~2017/10/11/08:10:00',antenna='DV15')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/08:10:10~2017/10/11/08:10:20',antenna='DV15')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/08:12:10~2017/10/11/08:12:20',antenna='DA64')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/08:13:30~2017/10/11/08:13:40',antenna='DA59,DV12')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/08:13:50~2017/10/11/08:14:00',antenna='DV04,DV20')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/08:14:50~2017/10/11/08:15:00',antenna='DV05')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/08:15:10~2017/10/11/08:15:20',antenna='DA47')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/08:16:40~2017/10/11/08:16:50',antenna='DA64,DA65')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/08:17:40~2017/10/11/08:17:50',antenna='DA54')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/08:18:00~2017/10/11/08:18:10',antenna='DA63')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/08:18:15~2017/10/11/08:18:25',antenna='DA45,DV20')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/08:19:20~2017/10/11/08:19:30',antenna='DA54,DV15')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/08:19:40~2017/10/11/08:19:50',antenna='DA45,DA65,DV03')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/08:24:10~2017/10/11/08:24:20',antenna='DA45,DA52')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/08:25:35~2017/10/11/08:25:45',antenna='DA52,DA63')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/08:27:15~2017/10/11/08:27:25',antenna='DA60')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/08:28:20~2017/10/11/08:28:30',antenna='DV05,DV13')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/08:28:40~2017/10/11/08:28:50',antenna='DA44,DA52,DV13')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/08:29:20~2017/10/11/08:29:30',antenna='DV13')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/08:29:40~2017/10/11/08:29:50',antenna='DA44,DV08,DV13,DV15,DV16,PM04')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/08:30:00~2017/10/11/08:30:10',antenna='DV08,DV12,DV13,DV20,PM04')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/08:31:10~2017/10/11/08:31:20',antenna='DA58,DA61,DA63,DV03,DV10')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/10/11/08:33:00~2017/10/11/08:53:00',antenna='')

flagdata(vis=LB_p6,mode='manual',spw='12',timerange='2017/10/15/07:32:40~2017/10/15/07:32:50',antenna='DV17')
flagdata(vis=LB_p6,mode='manual',spw='12',timerange='2017/10/15/07:33:00~2017/10/15/07:33:10',antenna='DV15,DV22')
flagdata(vis=LB_p6,mode='manual',spw='12',timerange='2017/10/15/07:34:00~2017/10/15/07:34:10',antenna='DA58')
flagdata(vis=LB_p6,mode='manual',spw='12',timerange='2017/10/15/07:35:20~2017/10/15/07:35:30',antenna='DA57')
flagdata(vis=LB_p6,mode='manual',spw='12',timerange='2017/10/15/07:37:00~2017/10/15/08:47:00',antenna='')

flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/10/16/08:43:10~2017/10/16/08:43:20',antenna='DA60,DV04,DV12')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/10/16/08:49:20~2017/10/16/08:49:30',antenna='DA51,DV25')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/10/16/08:45:50~2017/10/16/08:46:00',antenna='DA60')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/10/16/08:46:35~2017/10/16/08:46:45',antenna='PM03')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/10/16/08:48:30~2017/10/16/08:48:40',antenna='DA42')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/10/16/08:49:20~2017/10/16/08:49:30',antenna='DV05')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/10/16/08:49:30~2017/10/16/08:49:40',antenna='DV24')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/10/16/08:50:55~2017/10/16/08:51:05',antenna='DV24')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/10/16/08:51:50~2017/10/16/08:52:00',antenna='DV10')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/10/16/08:53:20~2017/10/16/08:53:30',antenna='DA52')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/10/16/08:54:00~2017/10/16/08:54:10',antenna='DV11')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/10/16/08:58:05~2017/10/16/08:58:15',antenna='DA63')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/10/16/08:59:30~2017/10/16/08:59:40',antenna='DA59,DV10')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/10/16/09:00:10~2017/10/16/09:00:20',antenna='DA44,DV10')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/10/16/09:00:50~2017/10/16/09:01:00',antenna='DA52,DA53')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/10/16/09:03:00~2017/10/16/10:03:00',antenna='')

flagdata(vis=LB_p6,mode='manual',spw='20',timerange='2017/12/09/04:46:00~2017/12/09/05:46:00',antenna='')

flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/12/17/05:50:20~2017/12/17/05:50:25',antenna='DA44')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/12/17/05:53:35~2017/12/17/05:53:45',antenna='DA44,DV13')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/12/17/05:53:50~2017/12/17/05:54:00',antenna='DA43,DA44')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/12/17/05:54:10~2017/12/17/05:54:20',antenna='DA43,DV13')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/12/17/05:54:30~2017/12/17/05:54:40',antenna='DA41,DA42,DA43,DA44,DV13')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/12/17/05:54:45~2017/12/17/05:54:55',antenna='DA41,DA42,DA43,DA44')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/12/17/05:57:55~2017/12/17/05:58:00',antenna='DV13,DV07')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/12/17/05:58:15~2017/12/17/05:58:20',antenna='DV11,DV13')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/12/17/05:58:30~2017/12/17/05:58:35',antenna='DV13')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/12/17/05:58:50~2017/12/17/05:58:55',antenna='DV13')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/12/17/05:59:10~2017/12/17/05:59:15',antenna='DA55,DV13,DV07')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/12/17/05:59:20~2017/12/17/05:59:40',antenna='DV13,DV07')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/12/17/05:59:40~2017/12/17/06:00:00',antenna='DA44,DA45,DV07')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/12/17/06:00:00~2017/12/17/06:00:05',antenna='DV13')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/12/17/06:00:20~2017/12/17/06:00:25',antenna='DV13')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/12/17/06:00:35~2017/12/17/06:00:45',antenna='DV13')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/12/17/06:00:55~2017/12/17/06:01:05',antenna='DV13')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/12/17/06:01:10~2017/12/17/06:01:20',antenna='DV13,DV07')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/12/17/06:01:30~2017/12/17/06:01:40',antenna='DV13,DV07')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/12/17/06:01:50~2017/12/17/06:01:55',antenna='DV13,DV07')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/12/17/06:02:10~2017/12/17/06:02:15',antenna='DV13,DV07')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/12/17/06:02:20~2017/12/17/06:02:30',antenna='DV07')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/12/17/06:02:45~2017/12/17/06:02:50',antenna='DV13')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/12/17/06:03:05~2017/12/17/06:03:10',antenna='DV13,DV07')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/12/17/06:03:25~2017/12/17/06:03:30',antenna='DV13')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/12/17/06:03:40~2017/12/17/06:03:45',antenna='DV13,DV07')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/12/17/06:04:00~2017/12/17/06:04:05',antenna='DV13')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/12/17/06:11:50~2017/12/17/06:12:00',antenna='DV11')

flagdata(vis=LB_p6,mode='manual',spw='28',timerange='2017/12/27/04:08:00~2017/12/27/04:08:05',antenna='DA42')
flagdata(vis=LB_p6,mode='manual',spw='28',timerange='2017/12/27/04:08:15~2017/12/27/04:08:20',antenna='DA42,DV13')
flagdata(vis=LB_p6,mode='manual',spw='28',timerange='2017/12/27/04:08:35~2017/12/27/04:08:40',antenna='DA42,DV13')
flagdata(vis=LB_p6,mode='manual',spw='28',timerange='2017/12/27/04:08:50~2017/12/27/04:09:00',antenna='DA42,DV13')
flagdata(vis=LB_p6,mode='manual',spw='28',timerange='2017/12/27/04:09:10~2017/12/27/04:09:20',antenna='DA42,DV13')
flagdata(vis=LB_p6,mode='manual',spw='28',timerange='2017/12/27/04:09:30~2017/12/27/04:09:35',antenna='DA42,DA44,DV13')
flagdata(vis=LB_p6,mode='manual',spw='28',timerange='2017/12/27/04:09:50~2017/12/27/04:09:55',antenna='DA42,DA44,DV13')
flagdata(vis=LB_p6,mode='manual',spw='28',timerange='2017/12/27/04:10:05~2017/12/27/04:10:10',antenna='DA42,DA44,DV13')
flagdata(vis=LB_p6,mode='manual',spw='28',timerange='2017/12/27/04:10:25~2017/12/27/04:10:30',antenna='DA42,DA44,DV13')
flagdata(vis=LB_p6,mode='manual',spw='28',timerange='2017/12/27/04:10:40~2017/12/27/04:10:45',antenna='DA42,DA44,DA50,DV13')
flagdata(vis=LB_p6,mode='manual',spw='28',timerange='2017/12/27/04:11:40~2017/12/27/04:11:45',antenna='DA42,DA44,DV13')
flagdata(vis=LB_p6,mode='manual',spw='28',timerange='2017/12/27/04:11:55~2017/12/27/04:12:00',antenna='DA44')
flagdata(vis=LB_p6,mode='manual',spw='28',timerange='2017/12/27/04:21:25~2017/12/27/04:21:30',antenna='DA55')
flagdata(vis=LB_p6,mode='manual',spw='28',timerange='2017/12/27/04:34:00~2017/12/27/04:34:10',antenna='DA44')
flagdata(vis=LB_p6,mode='manual',spw='28',timerange='2017/12/27/04:34:20~2017/12/27/04:34:30',antenna='DA44')
flagdata(vis=LB_p6,mode='manual',spw='28',timerange='2017/12/27/04:34:50~2017/12/27/04:35:00',antenna='DA44')
flagdata(vis=LB_p6,mode='manual',spw='28',timerange='2017/12/27/04:35:30~2017/12/27/04:35:40',antenna='DA44')
flagdata(vis=LB_p6,mode='manual',spw='28',timerange='2017/12/27/04:35:50~2017/12/27/04:36:00',antenna='DA44')
flagdata(vis=LB_p6,mode='manual',spw='28',timerange='2017/12/27/04:36:10~2017/12/27/04:36:15',antenna='DA44')
flagdata(vis=LB_p6,mode='manual',spw='28',timerange='2017/12/27/04:36:25~2017/12/27/04:36:30',antenna='DA44')
flagdata(vis=LB_p6,mode='manual',spw='28',timerange='2017/12/27/04:36:45~2017/12/27/04:36:50',antenna='DA44')
flagdata(vis=LB_p6,mode='manual',spw='28',timerange='2017/12/27/04:37:05~2017/12/27/04:37:10',antenna='DA44,DA50')
flagdata(vis=LB_p6,mode='manual',spw='28',timerange='2017/12/27/04:37:20~2017/12/27/04:37:25',antenna='DA44')
flagdata(vis=LB_p6,mode='manual',spw='28',timerange='2017/12/27/04:37:40~2017/12/27/04:37:45',antenna='DA44')
flagdata(vis=LB_p6,mode='manual',spw='28',timerange='2017/12/27/04:38:15~2017/12/27/04:38:25',antenna='DA44')
flagdata(vis=LB_p6,mode='manual',spw='28',timerange='2017/12/27/04:38:35~2017/12/27/04:38:40',antenna='DA44')
flagdata(vis=LB_p6,mode='manual',spw='28',timerange='2017/12/27/04:38:50~2017/12/27/04:39:00',antenna='DA44')
flagdata(vis=LB_p6,mode='manual',spw='28',timerange='2017/12/27/04:39:10~2017/12/27/04:39:15',antenna='DA44')
flagdata(vis=LB_p6,mode='manual',spw='28',timerange='2017/12/27/04:40:40~2017/12/27/04:40:50',antenna='DA50,DA55')

flagdata(vis=LB_p6,mode='manual',spw='32',timerange='2017/12/28/03:55:50~2017/12/28/03:55:55',antenna='DV11')
flagdata(vis=LB_p6,mode='manual',spw='32',timerange='2017/12/28/03:56:25~2017/12/28/03:56:30',antenna='DA55,DV11')
flagdata(vis=LB_p6,mode='manual',spw='32',timerange='2017/12/28/03:58:55~2017/12/28/03:59:00',antenna='DA55')
flagdata(vis=LB_p6,mode='manual',spw='32',timerange='2017/12/28/04:03:00~2017/12/28/04:03:10',antenna='DA44')
flagdata(vis=LB_p6,mode='manual',spw='32',timerange='2017/12/28/04:03:20~2017/12/28/04:03:25',antenna='DA44')
flagdata(vis=LB_p6,mode='manual',spw='32',timerange='2017/12/28/04:03:40~2017/12/28/04:03:45',antenna='DA44,DA55,DA46')
flagdata(vis=LB_p6,mode='manual',spw='32',timerange='2017/12/28/04:03:55~2017/12/28/04:04:05',antenna='DA44,DA46')
flagdata(vis=LB_p6,mode='manual',spw='32',timerange='2017/12/28/04:04:15~2017/12/28/04:04:20',antenna='DA44')
flagdata(vis=LB_p6,mode='manual',spw='32',timerange='2017/12/28/04:04:30~2017/12/28/04:04:35',antenna='DA44')
flagdata(vis=LB_p6,mode='manual',spw='32',timerange='2017/12/28/04:04:50~2017/12/28/04:04:55',antenna='DA44')
flagdata(vis=LB_p6,mode='manual',spw='32',timerange='2017/12/28/04:05:10~2017/12/28/04:05:15',antenna='DA44')
flagdata(vis=LB_p6,mode='manual',spw='32',timerange='2017/12/28/04:05:25~2017/12/28/04:05:30',antenna='DA44')
flagdata(vis=LB_p6,mode='manual',spw='32',timerange='2017/12/28/04:09:20~2017/12/28/04:09:30',antenna='DA42,DV13')
flagdata(vis=LB_p6,mode='manual',spw='32',timerange='2017/12/28/04:20:40~2017/12/28/04:20:45',antenna='DA44,DV11,DV13')
flagdata(vis=LB_p6,mode='manual',spw='32',timerange='2017/12/28/04:21:00~2017/12/28/04:21:05',antenna='DA44,DV11,DV13')
flagdata(vis=LB_p6,mode='manual',spw='32',timerange='2017/12/28/04:21:15~2017/12/28/04:21:20',antenna='DA44,DV11,DV13')
flagdata(vis=LB_p6,mode='manual',spw='32',timerange='2017/12/28/04:21:35~2017/12/28/04:21:40',antenna='DV13')
flagdata(vis=LB_p6,mode='manual',spw='32',timerange='2017/12/28/04:21:55~2017/12/28/04:22:00',antenna='DV13')
flagdata(vis=LB_p6,mode='manual',spw='32',timerange='2017/12/28/04:22:10~2017/12/28/04:22:20',antenna='DV13,DV11')
flagdata(vis=LB_p6,mode='manual',spw='32',timerange='2017/12/28/04:22:45~2017/12/28/04:22:55',antenna='DV11')
flagdata(vis=LB_p6,mode='manual',spw='32',timerange='2017/12/28/04:30:10~2017/12/28/04:30:20',antenna='DV13')
flagdata(vis=LB_p6,mode='manual',spw='32',timerange='2017/12/28/04:30:30~2017/12/28/04:30:40',antenna='DV13')
flagdata(vis=LB_p6,mode='manual',spw='32',timerange='2017/12/28/04:30:40~2017/12/28/04:30:50',antenna='DV13')

#Inspect gain tables to check if flagging worked
plotms(
    LB_p6,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(LB_selfcal_folder,f'{prefix}_LB_p6_gain_phase_vs_time_flagged.png'),
)

#Apply the calibration gains and split-off the corrected ms
applycal(
    vis=LB_cont_p5+'.ms',spw=LB_contspws,spwmap=LB_spw_mapping,
    gaintable=[LB_p6],interp='linearPD',calwt=True,applymode='calonly'
)
LB_cont_p6 = LB_cont_p5.replace('p5','p6')
os.system('rm -rf %s.ms*' %LB_cont_p6)
split(vis=LB_cont_p5+'.ms',outputvis=LB_cont_p6+'.ms',datacolumn='corrected')

#CLEAN, again with a ~6sigma threshold
tclean_wrapper(
    vis       = LB_cont_p6+'.ms',
    imagename = LB_cont_p6,
    threshold = '0.0392mJy',
    **LB_tclean_wrapper_kwargs
)
estimate_SNR(LB_cont_p6+'.image',disk_mask=LB_mask,noise_mask=noise_annulus_LB)
generate_image_png(
    LB_cont_p6+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_LB,10*rms_LB],
    save_folder=LB_selfcal_folder
)
#MWC_758_SBLB_contp6.image
#Beam 0.065 arcsec x 0.045 arcsec (7.20 deg)
#Flux inside disk mask: 65.99 mJy
#Peak intensity of source: 1.23 mJy/beam
#rms: 6.57e-03 mJy/beam
#Peak SNR: 186.88

#Maybe worth stopping after p4
#Check that the applycal p4 -> p5 didn't change the .ms structure (i.e., CLEANing to 6sigma, do I recover the same p5 gain solutions?)

split(vis=LB_cont_p3+'.ms',outputvis=LB_cont_p4+'.ms',datacolumn='corrected')

#CLEAN, again with a ~6sigma threshold
tclean_wrapper(
    vis       = LB_cont_p4+'.ms',
    imagename = LB_cont_p4,
    threshold = '0.0393mJy',
    **LB_tclean_wrapper_kwargs
)
estimate_SNR(LB_cont_p4+'.image',disk_mask=LB_mask,noise_mask=noise_annulus_LB)
generate_image_png(
    LB_cont_p4+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_LB,10*rms_LB],
    save_folder=LB_selfcal_folder
)
#MWC_758_SBLB_contp4.image
#Beam 0.065 arcsec x 0.045 arcsec (7.20 deg)
#Flux inside disk mask: 66.13 mJy
#Peak intensity of source: 1.21 mJy/beam
#rms: 6.54e-03 mJy/beam
#Peak SNR: 185.23

#MWC_758_SBLB_contp4.image
#Beam 0.065 arcsec x 0.045 arcsec (7.20 deg)
#Flux inside disk mask: 66.13 mJy
#Peak intensity of source: 1.21 mJy/beam
#rms: 6.54e-03 mJy/beam
#Peak SNR: 185.23

#Check same phase gain solutions as fifth round of phase-only self-cal
LB_p5_check = LB_p4.replace('p4','p5_check')
os.system('rm -rf '+LB_p5_check)
gaincal(
    vis=LB_cont_p4+'.ms',caltable=LB_p5_check,
    gaintype='T',combine='spw',calmode='p',solint='30s',
    spw=LB_contspws,refant=LB_refant,
    minsnr=3.,minblperant=4
)

#Inspect gain tables interactively and decide whether to manually flag something
#plotms(LB_p5_check,xaxis='time',yaxis='GainPhase',iteraxis='spw',showgui=True)
#Print calibration png file
#plotms(
#    LB_p5_check,xaxis='time', yaxis='GainPhase',iteraxis='spw',exprange='all',
#    overwrite=True,showgui=False,
#    plotfile=os.path.join(LB_selfcal_folder,f'{prefix}_LB_p5_check_gain_phase_vs_time.png')
#)

#Inspect gain tables and decide whether to flag something
plotms(
    LB_p5_check,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(LB_selfcal_folder,f'{prefix}_LB_p5_check_gain_phase_vs_time.png'),
)
"""
#CLEAN, again with a ~1sigma threshold before amplitude selfcal
tclean_wrapper(
    vis       = LB_cont_p4+'.ms',
    imagename = LB_cont_p4,
    threshold = '6.49e-03mJy',
    **LB_tclean_wrapper_kwargs
)
estimate_SNR(LB_cont_p4+'.image',disk_mask=LB_mask,noise_mask=noise_annulus_LB)
generate_image_png(
    LB_cont_p4+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_LB,10*rms_LB],
    save_folder=LB_selfcal_folder
)
#MWC_758_SBLB_contp4.image
#Beam 0.065 arcsec x 0.045 arcsec (7.20 deg)
#Flux inside disk mask: 66.22 mJy
#Peak intensity of source: 1.22 mJy/beam
#rms: 6.52e-03 mJy/beam
#Peak SNR: 186.57

#MWC_758_SBLB_contp4.image
#Beam 0.065 arcsec x 0.045 arcsec (7.20 deg)
#Flux inside disk mask: 57.75 mJy
#Peak intensity of source: 1.15 mJy/beam
#rms: 6.47e-03 mJy/beam
#Peak SNR: 177.82

#First round of amplitude self-cal
LB_ap0 = prefix+'_SBLB.ap0'
os.system('rm -rf '+LB_ap0)
gaincal(
    vis=LB_cont_p4+'.ms',caltable=LB_ap0,
    gaintype='T',combine='scan,spw',calmode='ap',solint='inf',
    spw=LB_contspws,refant=LB_refant,
    minsnr=5.0,minblperant=4,solnorm=False
)
# 3 of 43 solutions flagged due to SNR < 5 in spw=8  at 2017/10/11/08:03:57.8
# 6 of 46 solutions flagged due to SNR < 5 in spw=12 at 2017/10/15/07:57:58.7
# 9 of 44 solutions flagged due to SNR < 5 in spw=16 at 2017/10/16/09:13:45.9

#Inspect gain tables interactively and decide whether to manually flag something
#plotms(LB_ap0,xaxis='time', yaxis='GainPhase',iteraxis='spw')
#plotms(LB_ap0,xaxis='time', yaxis='GainAmp',iteraxis='spw')
#
#Print calibration png file
#plotms(LB_ap0,xaxis='time', yaxis='GainPhase',iteraxis='spw',exprange='all',
#       overwrite=True,showgui=False,
#       plotfile=os.path.join(LB_selfcal_folder,f'{prefix}_LB_gain_ap0_phase_vs_time.png'))
#plotms(LB_ap0,xaxis='time', yaxis='GainAmp',iteraxis='spw',exprange='all',
#       overwrite=True,showgui=False,
#       plotfile=os.path.join(LB_selfcal_folder,f'{prefix}_LB_gain_ap0_amp_vs_time.png'))

#Inspect gain tables and decide whether to flag something
plotms(
    LB_ap0,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(LB_selfcal_folder,f'{prefix}_LB_ap0_gain_phase_vs_time.png'),
)
plotms(
    LB_ap0,xaxis='time',yaxis='GainAmp',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(LB_selfcal_folder,f'{prefix}_LB_ap0_gain_amp_vs_time.png'),
)

#Flag problematic antennas
flagdata(vis=LB_ap0,mode='clip',clipminmax=[0.8,1.2],clipoutside=True,datacolumn='CPARAM')
flagdata(vis=LB_ap0,mode='manual',spw='0', antenna='DA54,DA56,DV24')
flagdata(vis=LB_ap0,mode='manual',spw='4', antenna='DA54,DV05,DV24')
flagdata(vis=LB_ap0,mode='manual',spw='8', antenna='DA45')
flagdata(vis=LB_ap0,mode='manual',spw='16',antenna='DV12')
flagdata(vis=LB_ap0,mode='manual',spw='24',antenna='DV07')

#Inspect gain tables to check if flagging worked
plotms(
    LB_ap0,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(LB_selfcal_folder,f'{prefix}_LB_ap0_gain_phase_vs_time_flagged.png'),
)
plotms(
    LB_ap0,xaxis='time',yaxis='GainAmp',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(LB_selfcal_folder,f'{prefix}_LB_ap0_gain_amp_vs_time_flagged.png'),
)

#Apply the calibration gains and split-off the corrected ms
applycal(
    vis=LB_cont_p4+'.ms',spw=LB_contspws,spwmap=LB_spw_mapping,
    gaintable=[LB_ap0],interp='linearPD',calwt=True,applymode='calonly'
)
LB_cont_ap0 = prefix+'_SBLB_contap0'
os.system('rm -rf %s.ms*' %LB_cont_ap0)
split(vis=LB_cont_p4+'.ms',outputvis=LB_cont_ap0+'.ms',datacolumn='corrected')

#CLEAN with a ~1sigma threshold before amplitude self-cal
tclean_wrapper(
    vis=LB_cont_ap0+'.ms',
    imagename = LB_cont_ap0,
    threshold = '6.47e-03mJy',
    **LB_tclean_wrapper_kwargs
)
estimate_SNR(LB_cont_ap0+'.image',disk_mask=LB_mask,noise_mask=noise_annulus_LB)
generate_image_png(
    LB_cont_ap0+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_LB,10*rms_LB],
    save_folder=LB_selfcal_folder
)
#MWC_758_SBLB_contp4.image
#Beam 0.065 arcsec x 0.045 arcsec (7.20 deg)
#Flux inside disk mask: 57.75 mJy
#Peak intensity of source: 1.15 mJy/beam
#rms: 6.47e-03 mJy/beam
#Peak SNR: 177.82

#MWC_758_SBLB_contap0.image
#Beam 0.068 arcsec x 0.045 arcsec (9.81 deg)
#Flux inside disk mask: 57.68 mJy
#Peak intensity of source: 1.21 mJy/beam
#rms: 6.49e-03 mJy/beam
#Peak SNR: 185.99

#For each step of the self-cal, check how it improved things
non_self_caled_LB_vis = LB_cont_p0
self_caled_LB_visibilities = {
    'p1' :LB_cont_p1,
    'p2' :LB_cont_p2,
    'p3' :LB_cont_p3,
    'p4' :LB_cont_p4,
    #'p5' :LB_cont_p5,
    #'p6' :LB_cont_p6,
    'ap0':LB_cont_ap0,
}

#LB_EBs = ('EB0','EB1','EB2','EB3','EB4','EB5','EB6','EB7','EB8','EB9')
#LB_EB_spws = ('0,1,2,3', '4,5,6,7', '8,9,10,11', '12,13,14,15', '16,17,18,19', '20,21,22,23', '24,25,26,27', '28,29,30,31', '32,33,34,35') #fill out by referring to listobs output
#
#for self_cal_step,self_caled_vis in self_caled_LB_visibilities.items():
#    for EB_key,spw in zip(LB_EBs,LB_EB_spws):
#
#        print('selfcal_step: {}'.format(self_cal_step)+', EB: {}'.format(EB_key))
#
#        nametemplate = f'{prefix}_LB_{EB_key}_{self_cal_step}_compare_amp_vs_time'
#        visibilities = [self_caled_vis+'.ms',non_self_caled_LB_vis+'.ms']
#
#        plot_amp_vs_time_comparison(
#            nametemplate=nametemplate,visibilities=visibilities,spw=spw,
#            uvrange=uv_ranges['LB'],output_folder=LB_selfcal_folder
#        )

#Set to the EB of the combined SBLB data that corresponds to flux_ref_EB
SBLB_flux_ref_EB = 8 #this is SB_EB3

all_LB_visibilities = self_caled_LB_visibilities.copy()
all_LB_visibilities['p0'] = LB_cont_p0

total_number_of_EBs = number_of_EBs['SB'] + number_of_EBs['LB']
for self_cal_step,vis_name in all_LB_visibilities.items():
    #Split out SB EBs
    vis_ms = vis_name+'.ms'
    nametemplate = vis_ms.replace('.ms','_EB')
    split_all_obs(msfile=vis_ms,nametemplate=nametemplate)

    exported_ms = []
    for i in range(total_number_of_EBs):
        EB_vis = f'{nametemplate}{i}.ms'
        listobs(vis=EB_vis,listfile=EB_vis+'.listobs.txt',overwrite=True)
        #Export MS contents into numpy save files
        export_MS(EB_vis)
        exported_ms.append(EB_vis.replace('.ms','.vis.npz'))

    for i,exp_ms in enumerate(exported_ms):
        png_filename = f'flux_comparison_EB{i}_to_EB{SBLB_flux_ref_EB}'+f'_SBLB_{self_cal_step}.png'
        plot_label = os.path.join(LB_selfcal_folder,png_filename)
        estimate_flux_scale(
            reference=f'{nametemplate}{SBLB_flux_ref_EB}.vis.npz',
            comparison=exp_ms,incl=incl,PA=PA,plot_label=plot_label,#uvbins=np.arange(40.,300.,20.),
        )

    fluxscale = [1.,]*total_number_of_EBs
    plot_label = os.path.join(LB_selfcal_folder,f'deprojected_vis_profiles_SBLB_{self_cal_step}.png')
    plot_deprojected(filelist=exported_ms,fluxscale=fluxscale,PA=PA,incl=incl,show_err=True,plot_label=plot_label)

#Redo the flux comparison images without the uvbins parameters, to have clearer plots
#for self_cal_step,vis in self_caled_LB_visibilities.items():
#    nametemplate = f'{prefix}_SBLB_cont{self_cal_step}_EB'
#    for i in range(total_number_of_EBs):
#        reference = f'{nametemplate}{SBLB_flux_ref_EB}.vis.npz'
#        output = f'flux_comparison_2_EB{i}_to_EB{SBLB_flux_ref_EB}'+f'_SBLB_{self_cal_step}.png'
#        plot_label = os.path.join(LB_selfcal_folder,output)
#        estimate_flux_scale(
#           reference=reference,comparison=f'{nametemplate}{i}.vis.npz',
#           incl=incl PA=PA,plot_label=plot_label
#        )

# In the concatenated SBLB file:
# EB0 = LB EB0
# EB1 = LB EB1
# EB2 = LB EB2
# EB3 = LB EB3
# EB4 = LB EB4
# EB5 = SB EB0
# EB6 = SB EB1
# EB7 = SB EB2
# EB8 = SB EB3

#ratio          = [0.88988,0.89305,0.89593,0.90156,0.90895,0.92577] #MWC_758_SB_contp0...p5_EB0.vis.npz vs MWC_758_SB_EB3_initcont.vis.npz
#scaling_factor = [0.943  ,0.945  ,0.947  ,0.950  ,0.953  ,0.962  ]
#ratio          = [0.98340,0.97806,0.98762,0.99364,0.99671,0.99951] #MWC_758_SB_contp0...p5_EB1.vis.npz vs MWC_758_SB_EB3_initcont.vis.npz
#scaling_factor = [0.992  ,0.989  ,0.994  ,0.997  ,0.998  ,1.000  ]
#ratio          = [0.93677,0.96069,0.95945,0.96593,0.96917,0.97148] #MWC_758_SB_contp0...p5_EB2.vis.npz vs MWC_758_SB_EB3_initcont.vis.npz
#scaling_factor = [0.968  ,0.980  ,0.980  ,0.983  ,0.984  ,0.986  ]

#ratio          = [0.82477,0.91568,0.98312,1.00422,1.03841,1.03459] #MWC_758_SBLB_contp0...p5_EB0.vis.npz vs MWC_758_SBLB_contp0...p5_EB8.vis.npz
#scaling_factor = [0.908  ,0.957  ,0.992  ,1.002  ,1.019  ,1.017  ]
#ratio          = [0.82298,0.91481,0.94454,0.96959,1.02252,1.06883] #MWC_758_SBLB_contp0...p5_EB1.vis.npz vs MWC_758_SBLB_contp0...p5_EB8.vis.npz
#scaling_factor = [0.907  ,0.956  ,0.972  ,0.985  ,1.011  ,1.034  ]
#ratio          = [0.80492,0.92314,0.90848,0.91003,0.93038,0.94416] #MWC_758_SBLB_contp0...p5_EB2.vis.npz vs MWC_758_SBLB_contp0...p5_EB8.vis.npz
#scaling_factor = [0.897  ,0.961  ,0.953  ,0.954  ,0.965  ,0.972  ]
#ratio          = [0.94428,0.99588,1.00153,0.98510,1.01808,1.05846] #MWC_758_SBLB_contp0...p5_EB3.vis.npz vs MWC_758_SBLB_contp0...p5_EB8.vis.npz
#scaling_factor = [0.972  ,0.998  ,1.001  ,0.993  ,1.009  ,1.029  ]
#ratio          = [1.03013,1.06972,1.04402,1.04893,1.07551,1.07234] #MWC_758_SBLB_contp0...p5_EB4.vis.npz vs MWC_758_SBLB_contp0...p5_EB8.vis.npz
#scaling_factor = [1.015  ,1.034  ,1.022  ,1.024  ,1.037  ,1.036  ]
#ratio          = [0.92268,0.92291,0.92254,0.92251,0.92266,0.99209] #MWC_758_SBLB_contp0...p5_EB5.vis.npz vs MWC_758_SBLB_contp0...p5_EB8.vis.npz
#scaling_factor = [0.961  ,0.961  ,0.960  ,0.960  ,0.961  ,0.996  ]
#ratio          = [1.00050,1.00059,1.00074,0.99988,0.99991,0.98880] #MWC_758_SBLB_contp0...p5_EB6.vis.npz vs MWC_758_SBLB_contp0...p5_EB8.vis.npz
#scaling_factor = [1.000  ,1.000  ,1.000  ,1.000  ,1.000  ,0.994  ]
#ratio          = [0.97132,0.97124,0.97121,0.97124,0.97127,0.99629] #MWC_758_SBLB_contp0...p5_EB7.vis.npz vs MWC_758_SBLB_contp0...p5_EB8.vis.npz
#scaling_factor = [0.986  ,0.986  ,0.986  ,0.986  ,0.986  ,0.998  ]

# It might be worth stopping at p4 for the SBLBs because the gain solutions become consistently strange for contp5,6 and sloint <30s 
# (maybe due to the very very long baselines?)
# SB_EB0,LB_EB2,LB4 substantial improvement but flux offsets >4%, all the other EBs show flux offsets <4%
# Amplitude selfcal after p4 affects the flux ratios substantially:
#   - All the SB_EBs are essentially in line with each other (rescale below just to be as consistent as possible)
#   - All the LB_EBs except LB_EB0 show flux offsets >4% compared to SB_EB3. Rescale (also LB_EB0 just to be as consistent as possible)
# Stress the difference from the flux ratios you got CLEANing with robust=0.5 (does a higher SNR really matter, then?)
# Correct all EBs with flux ratios after ap0_afterp4

#END of COMBINED SB+LB phase-only self-cal iteration 1

#Rescaling flux of the non-self-caled SB and LB EBs

#First check that you get the same flux ratios as above
for baseline_key,n_EB in number_of_EBs.items():
    for i in range(n_EB):
        vis = f'{prefix}_{baseline_key}_EB{i}_initcont_selfcal.ms'#vis = f'{prefix}_{baseline_key}_EB{i}_initcont_shift.ms'
        listobs(vis=vis,listfile=f'{vis}.listobs.txt',overwrite=True)

flux_ref_EB = 'SB_EB3' 
for params in data_params.values():
    plot_label = os.path.join(LB_selfcal_folder,'check_flux_comparison_'+params['name']+f'_to_{flux_ref_EB}.png')
    estimate_flux_scale(
        reference=f'{prefix}_{flux_ref_EB}_initcont_selfcal.vis.npz',#reference=f'{prefix}_{flux_ref_EB}_initcont_shift.vis.npz',
        comparison=prefix+'_'+params['name']+'_initcont_selfcal.vis.npz',#comparison=prefix+'_'+params['name']+'_initcont_shift.vis.npz',
        incl=incl,PA=PA,plot_label=plot_label
    )

#The ratio of the fluxes of MWC_758_LB_EB0_initcont_selfcal.vis.npz to
#MWC_758_SB_EB3_initcont_selfcal.vis.npz is 0.95742
#The scaling factor for gencal is 0.978 for your comparison measurement
#The error on the weighted mean ratio is 4.902e-03, although it's likely that
#the weights in the measurement sets are too off by some constant factor

#The ratio of the fluxes of MWC_758_LB_EB1_initcont_selfcal.vis.npz to
#MWC_758_SB_EB3_initcont_selfcal.vis.npz is 0.93723
#The scaling factor for gencal is 0.968 for your comparison measurement
#The error on the weighted mean ratio is 4.833e-03, although it's likely that
#the weights in the measurement sets are too off by some constant factor

#The ratio of the fluxes of MWC_758_LB_EB2_initcont_selfcal.vis.npz to
#MWC_758_SB_EB3_initcont_selfcal.vis.npz is 0.90976
#The scaling factor for gencal is 0.954 for your comparison measurement
#The error on the weighted mean ratio is 6.064e-03, although it's likely that
#the weights in the measurement sets are too off by some constant factor

#The ratio of the fluxes of MWC_758_LB_EB3_initcont_selfcal.vis.npz to
#MWC_758_SB_EB3_initcont_selfcal.vis.npz is 1.01642
#The scaling factor for gencal is 1.008 for your comparison measurement
#The error on the weighted mean ratio is 4.064e-03, although it's likely that
#the weights in the measurement sets are too off by some constant factor

#The ratio of the fluxes of MWC_758_LB_EB4_initcont_selfcal.vis.npz to
#MWC_758_SB_EB3_initcont_selfcal.vis.npz is 1.11297
#The scaling factor for gencal is 1.055 for your comparison measurement
#The error on the weighted mean ratio is 5.086e-03, although it's likely that
#the weights in the measurement sets are too off by some constant factor

#The ratio of the fluxes of MWC_758_SB_EB0_initcont_selfcal.vis.npz to
#MWC_758_SB_EB3_initcont_selfcal.vis.npz is 0.88668
#The scaling factor for gencal is 0.942 for your comparison measurement
#The error on the weighted mean ratio is 1.119e-03, although it's likely that
#the weights in the measurement sets are too off by some constant factor

#The ratio of the fluxes of MWC_758_SB_EB1_initcont_selfcal.vis.npz to
#MWC_758_SB_EB3_initcont_selfcal.vis.npz is 0.98104
#The scaling factor for gencal is 0.990 for your comparison measurement
#The error on the weighted mean ratio is 7.696e-04, although it's likely that
#the weights in the measurement sets are too off by some constant factor

#The ratio of the fluxes of MWC_758_SB_EB2_initcont_selfcal.vis.npz to
#MWC_758_SB_EB3_initcont_selfcal.vis.npz is 0.93757
#The scaling factor for gencal is 0.968 for your comparison measurement
#The error on the weighted mean ratio is 8.283e-04, although it's likely that
#the weights in the measurement sets are too off by some constant factor

#The ratio of the fluxes of MWC_758_SB_EB3_initcont_selfcal.vis.npz to
#MWC_758_SB_EB3_initcont_selfcal.vis.npz is 1.00000
#The scaling factor for gencal is 1.000 for your comparison measurement
#The error on the weighted mean ratio is 6.559e-04, although it's likely that
#the weights in the measurement sets are too off by some constant factor

# os.system('rm -rf '+prefix+'_SB_EB0_initcont_shift_rescaled.ms')
# rescale_flux(vis=prefix+'_SB_EB0_initcont_shift.ms', gencalparameter=[1.150])
# Splitting out rescaled values into new MS: CQ_Tau_SB_EB0_initcont_shift_rescaled.ms
# listobs(vis=prefix+'_SB_EB0_initcont_shift_rescaled.ms',listfile=prefix+'_SB_EB0_initcont_shift_rescaled.ms.listobs.txt',overwrite=True)
#
# export_MS(prefix+'_SB_EB0_initcont_shift_rescaled.ms') 
# Measurement set exported to CQ_Tau_SB_EB0_initcont_shift_rescaled.vis.npz

# Rescale the LBs
for i,gencalpar in enumerate([1.017,1.034,0.972,1.029,1.036]): #[0.995,1.029,0.970,1.030,1.031]
    os.system('rm -rf '+prefix+f'_LB_EB{i}_initcont_selfcal_rescaled.ms')
    rescale_flux(vis=prefix+f'_LB_EB{i}_initcont_selfcal.ms', gencalparameter=[gencalpar])
    #Splitting out rescaled values into new MS: MWC_758_SB_LB0_initcont_selfcal_rescaled.ms
    #Splitting out rescaled values into new MS: MWC_758_SB_LB1_initcont_selfcal_rescaled.ms
    #Splitting out rescaled values into new MS: MWC_758_SB_LB2_initcont_selfcal_rescaled.ms
    #Splitting out rescaled values into new MS: MWC_758_SB_LB3_initcont_selfcal_rescaled.ms
    #Splitting out rescaled values into new MS: MWC_758_SB_LB4_initcont_selfcal_rescaled.ms
    listobs(vis=prefix+f'_LB_EB{i}_initcont_selfcal_rescaled.ms',listfile=prefix+f'_LB_EB{i}_initcont_selfcal_rescaled.ms.listobs.txt',overwrite=True)

    export_MS(prefix+f'_LB_EB{i}_initcont_selfcal_rescaled.ms')
    #Measurement set exported to MWC_758_SB_EB0_initcont_selfcal_rescaled.vis.npz
    #Measurement set exported to MWC_758_SB_EB1_initcont_selfcal_rescaled.vis.npz
    #Measurement set exported to MWC_758_SB_EB2_initcont_selfcal_rescaled.vis.npz
    #Measurement set exported to MWC_758_SB_EB3_initcont_selfcal_rescaled.vis.npz
    #Measurement set exported to MWC_758_SB_EB4_initcont_selfcal_rescaled.vis.npz

# Rescale the SBs
for i,gencalpar in enumerate([0.996,0.994,0.998]):
    os.system('rm -rf '+prefix+f'_SB_EB{i}_initcont_selfcal_rescaled.ms')
    rescale_flux(vis=prefix+f'_SB_EB{i}_initcont_selfcal.ms', gencalparameter=[gencalpar])
    #Splitting out rescaled values into new MS: MWC_758_SB_EB0_initcont_selfcal_rescaled.ms
    #Splitting out rescaled values into new MS: MWC_758_SB_EB1_initcont_selfcal_rescaled.ms
    #Splitting out rescaled values into new MS: MWC_758_SB_EB2_initcont_selfcal_rescaled.ms
    listobs(vis=prefix+f'_SB_EB{i}_initcont_selfcal_rescaled.ms',listfile=prefix+f'_SB_EB{i}_initcont_selfcal_rescaled.ms.listobs.txt',overwrite=True)

    export_MS(prefix+f'_SB_EB{i}_initcont_selfcal_rescaled.ms')
    #Measurement set exported to MWC_758_SB_EB0_initcont_selfcal_rescaled.vis.npz
    #Measurement set exported to MWC_758_SB_EB1_initcont_selfcal_rescaled.vis.npz
    #Measurement set exported to MWC_758_SB_EB2_initcont_selfcal_rescaled.vis.npz

# output = f'flux_comparison_LB_EB{i}_rescaled_to_SB_EB3.png'
# plot_label = os.path.join(LB_selfcal_folder,output)
# estimate_flux_scale(
#     reference  = prefix+f'_SB_EB3_initcont_selfcal.vis.npz',#prefix+'_SB_EB3_initcont_shift.vis.npz',
#     comparison = prefix+f'_LB_EB{i}_initcont_selfcal_rescaled.vis.npz',#prefix+'_LB_EB0_initcont_shift_rescaled.vis.npz',
#     incl=incl,PA=PA,plot_label=plot_label
# )
# #The ratio of the fluxes of MWC_758_LB_EB0_initcont_selfcal_rescaled.vis.npz to
# #MWC_758_SB_EB3_initcont_selfcal.vis.npz is 0.91307
# #The scaling factor for gencal is 0.956 for your comparison measurement
# #The error on the weighted mean ratio is 4.675e-03, although it's likely that
# #the weights in the measurement sets are too off by some constant factor

flux_ref_EB = 'SB_EB3' 
for params in data_params.values():
    plot_label = os.path.join(LB_selfcal_folder,'flux_comparison_'+params['name']+f'_rescaled_to_{flux_ref_EB}.png')
    if params['name'] == flux_ref_EB:
        estimate_flux_scale(
            reference=f'{prefix}_{flux_ref_EB}_initcont_selfcal.vis.npz',
            comparison=prefix+'_'+params['name']+'_initcont_selfcal.vis.npz',
            incl=incl,PA=PA,plot_label=plot_label
        )
    else:
        estimate_flux_scale(
            reference=f'{prefix}_{flux_ref_EB}_initcont_selfcal.vis.npz',
            comparison=prefix+'_'+params['name']+'_initcont_selfcal_rescaled.vis.npz',
            incl=incl,PA=PA,plot_label=plot_label
        )

#The ratio of the fluxes of MWC_758_LB_EB0_initcont_selfcal_rescaled.vis.npz to
#MWC_758_SB_EB3_initcont_selfcal.vis.npz is 0.92568
#The scaling factor for gencal is 0.962 for your comparison measurement
#The error on the weighted mean ratio is 4.740e-03, although it's likely that
#the weights in the measurement sets are too off by some constant factor

#The ratio of the fluxes of MWC_758_LB_EB1_initcont_selfcal_rescaled.vis.npz to
#MWC_758_SB_EB3_initcont_selfcal.vis.npz is 0.87661
#The scaling factor for gencal is 0.936 for your comparison measurement
#The error on the weighted mean ratio is 4.520e-03, although it's likely that
#the weights in the measurement sets are too off by some constant factor

#The ratio of the fluxes of MWC_758_LB_EB2_initcont_selfcal_rescaled.vis.npz to
#MWC_758_SB_EB3_initcont_selfcal.vis.npz is 0.96292
#The scaling factor for gencal is 0.981 for your comparison measurement
#The error on the weighted mean ratio is 6.418e-03, although it's likely that
#the weights in the measurement sets are too off by some constant factor

#The ratio of the fluxes of MWC_758_LB_EB3_initcont_selfcal_rescaled.vis.npz to
#MWC_758_SB_EB3_initcont_selfcal.vis.npz is 0.95994
#The scaling factor for gencal is 0.980 for your comparison measurement
#The error on the weighted mean ratio is 3.838e-03, although it's likely that
#the weights in the measurement sets are too off by some constant factor

#The ratio of the fluxes of MWC_758_LB_EB4_initcont_selfcal_rescaled.vis.npz to
#MWC_758_SB_EB3_initcont_selfcal.vis.npz is 1.03697
#The scaling factor for gencal is 1.018 for your comparison measurement
#The error on the weighted mean ratio is 4.739e-03, although it's likely that
#the weights in the measurement sets are too off by some constant factor

#The ratio of the fluxes of MWC_758_SB_EB0_initcont_selfcal_rescaled.vis.npz to
#MWC_758_SB_EB3_initcont_selfcal.vis.npz is 0.89382
#The scaling factor for gencal is 0.945 for your comparison measurement
#The error on the weighted mean ratio is 1.128e-03, although it's likely that
#the weights in the measurement sets are too off by some constant factor

#The ratio of the fluxes of MWC_758_SB_EB1_initcont_selfcal_rescaled.vis.npz to
#MWC_758_SB_EB3_initcont_selfcal.vis.npz is 0.99292
#The scaling factor for gencal is 0.996 for your comparison measurement
#The error on the weighted mean ratio is 7.790e-04, although it's likely that
#the weights in the measurement sets are too off by some constant factor

#The ratio of the fluxes of MWC_758_SB_EB2_initcont_selfcal_rescaled.vis.npz to
#MWC_758_SB_EB3_initcont_selfcal.vis.npz is 0.94133
#The scaling factor for gencal is 0.970 for your comparison measurement
#The error on the weighted mean ratio is 8.316e-04, although it's likely that
#the weights in the measurement sets are too off by some constant factor

#The ratio of the fluxes of MWC_758_SB_EB3_initcont_selfcal.vis.npz to
#MWC_758_SB_EB3_initcont_selfcal.vis.npz is 1.00000
#The scaling factor for gencal is 1.000 for your comparison measurement
#The error on the weighted mean ratio is 6.559e-04, although it's likely that
#the weights in the measurement sets are too off by some constant factor

#Begin of SB self-cal - iteration 2
#For phase self-cal, clean down to ~6sigma
SB_selfcal_iteration2_folder = get_figures_folderpath('7.1_selfcal_SB_iteration2_figures')
make_figures_folder(SB_selfcal_iteration2_folder)

SB_iteration2_cont_p0 = prefix+'_SB_iteration2_contp0'
os.system('rm -rf %s.ms*' %SB_iteration2_cont_p0)
concat(
    vis=[
        f'{prefix}_SB_EB0_initcont_selfcal_rescaled.ms',#f'{prefix}_SB_EB0_initcont_shift_rescaled.ms',
        f'{prefix}_SB_EB1_initcont_selfcal_rescaled.ms',#f'{prefix}_SB_EB1_initcont_shift_rescaled.ms',
        f'{prefix}_SB_EB2_initcont_selfcal_rescaled.ms',#f'{prefix}_SB_EB2_initcont_shift_rescaled.ms',
        f'{prefix}_SB_EB3_initcont_selfcal.ms',#f'{prefix}_SB_EB3_initcont_shift.ms',
    ],
    concatvis=SB_iteration2_cont_p0+'.ms',dirtol='0.1arcsec',copypointing=False
)
listobs(vis=SB_iteration2_cont_p0+'.ms',listfile=SB_iteration2_cont_p0+'.ms.listobs.txt',overwrite=True)
#2025-01-29 06:42:49     WARN    concat::::casa  The setup of the input MSs is not fully consistent. The concatenation may fail
#2025-01-29 06:42:49     WARN    concat::::casa  and/or the affected columns may contain partially only default data.
#2025-01-29 06:42:49     WARN    concat::::casa  
#    {'MWC_758_SB_EB3_initcont_selfcal.ms': {'Main': {'present_a': True, 'present_b': True, 'missingcol_a': ['MODEL_DATA'], 'missingcol_b': []}}}

#Define new SB mask using the same centre as before (checked and agrees with the listobs one)
mask_pa        = PA
mask_semimajor = 0.8
mask_semiminor = 0.8
mask_ra  = '05h30m27.535799s' 
mask_dec = '25d19m56.611160s'

SB_mask= f"ellipse[[{mask_ra},{mask_dec}], [{mask_semimajor:.3f}arcsec, {mask_semiminor:.3f}arcsec], {mask_pa:.1f}deg]"

noise_annulus_SB = f"annulus[[{mask_ra}, {mask_dec}],['4.arcsec', '6.arcsec']]"

SB_tclean_wrapper_kwargs = {
    'deconvolver':'multiscale','scales':[0,2,4,8,12],
    'smallscalebias':0.6,'gain':0.3,'cycleniter':300,
    'niter':1000000,'mask':SB_mask,
    'cellsize':'0.014arcsec','imsize':3600,
    'parallel':use_parallel,'savemodel':'modelcolumn',
    'robust':0.5,'interactive':False,
    'gridder':'standard',
}

tclean_wrapper(
    vis       = SB_iteration2_cont_p0+'.ms',
    imagename = SB_iteration2_cont_p0,
    threshold = '0.0960mJy', 
    **SB_tclean_wrapper_kwargs
)
estimate_SNR(SB_iteration2_cont_p0+'.image',disk_mask=SB_mask,noise_mask=noise_annulus_SB)
rms_iteration2_SB = imstat(imagename=SB_iteration2_cont_p0+'.image',region=noise_annulus_SB)['rms'][0]
generate_image_png(
    SB_iteration2_cont_p0+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_iteration2_SB,10*rms_iteration2_SB],
    save_folder=SB_selfcal_iteration2_folder
)
#MWC_758_SB_contp0.image
#Beam 0.193 arcsec x 0.146 arcsec (-28.31 deg)
#Flux inside disk mask: 56.86 mJy
#Peak intensity of source: 5.03 mJy/beam
#rms: 1.58e-02 mJy/beam
#Peak SNR: 317.43

#MWC_758_SB_iteration2_contp0.image
#Beam 0.193 arcsec x 0.146 arcsec (-28.29 deg)
#Flux inside disk mask: 57.17 mJy
#Peak intensity of source: 5.06 mJy/beam
#rms: 1.59e-02 mJy/beam
#Peak SNR: 318.54

#First round of phase-only self-cal
#NOTE: you need .p1 instead of _p1 in the caltable name if you want flagdata to work (i.e., flagging problematic antennas non interactively)...
SB_iteration2_p1 = prefix+'_SB_iteration2.p1'
os.system('rm -rf '+SB_iteration2_p1)
gaincal(
    vis=SB_iteration2_cont_p0+'.ms',caltable=SB_iteration2_p1,
    gaintype='G',combine='scan,spw',calmode='p',solint='inf',
    spw=SB_contspws,refant=SB_refant,
    minsnr=3.,minblperant=4
)

#Inspect gain tables interactively and decide whether to manually flag something
#plotms(SB_iteration2_p1,xaxis='time',yaxis='GainPhase',iteraxis='spw',showgui=True)
#Print calibration png file
#plotms(
#    SB_iteration2_p1,xaxis='time',yaxis='GainPhase',iteraxis='spw',exprange='all',
#    overwrite=True,showgui=False,
#    plotfile=os.path.join(SB_selfcal_iteration2_folder,f'{prefix}_SB_iteration2_p1_gain_phase_vs_time.png')
#)

#Inspect gain tables and decide whether to flag something
plotms(
    SB_iteration2_p1,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(SB_selfcal_iteration2_folder,f'{prefix}_SB_iteration2_p1_gain_phase_vs_time.png'),
)

#Apply the calibration gains and split-off the corrected ms
applycal(
    vis=SB_iteration2_cont_p0+'.ms',spw=SB_contspws,spwmap=SB_spw_mapping,
    gaintable=[SB_iteration2_p1],interp='linearPD',calwt=True,applymode='calonly'
)
SB_iteration2_cont_p1 = SB_iteration2_cont_p0.replace('p0','p1')
os.system('rm -rf %s.ms*' %SB_iteration2_cont_p1)
split(vis=SB_iteration2_cont_p0+'.ms',outputvis=SB_iteration2_cont_p1+'.ms',datacolumn='corrected')

#CLEAN, again with a ~6sigma threshold
tclean_wrapper(
    vis       = SB_iteration2_cont_p1+'.ms',
    imagename = SB_iteration2_cont_p1,
    threshold = '0.0924mJy',
    **SB_tclean_wrapper_kwargs
)
estimate_SNR(SB_iteration2_cont_p1+'.image',disk_mask=SB_mask,noise_mask=noise_annulus_SB)
generate_image_png(
    SB_iteration2_cont_p1+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_iteration2_SB,10*rms_iteration2_SB],
    save_folder=SB_selfcal_iteration2_folder
)
#MWC_758_SB_contp1.image
#Beam 0.193 arcsec x 0.146 arcsec (-28.31 deg)
#Flux inside disk mask: 56.81 mJy
#Peak intensity of source: 5.06 mJy/beam
#rms: 1.53e-02 mJy/beam
#Peak SNR: 330.85

#MWC_758_SB_iteration2_contp1.image
#Beam 0.193 arcsec x 0.146 arcsec (-28.29 deg)
#Flux inside disk mask: 57.16 mJy
#Peak intensity of source: 5.09 mJy/beam
#rms: 1.53e-02 mJy/beam
#Peak SNR: 331.53

#Re-run the first step of phase-only selfcal with infinite integration interval to improve the model and get a better EB alignment
#One step is enough, there is not much improvement later on.
#
#Step .p1_bis
SB_iteration2_p1_bis = SB_iteration2_p1.replace('p1','p1_bis')
os.system('rm -rf '+SB_iteration2_p1_bis)
gaincal(
    vis=SB_iteration2_cont_p1+'.ms',caltable=SB_iteration2_p1_bis,
    gaintype='T',combine='scan,spw',calmode='p',solint='inf',
    spw=SB_contspws,refant=SB_refant,
    minsnr=3.,minblperant=4
)

#Inspect gain tables interactively and decide whether to manually flag something
#plotms(SB_iteration2_p1_bis,xaxis='time',yaxis='GainPhase',iteraxis='spw',showgui=True)
#Print calibration png file
#plotms(
#    SB_iteration2_p1_bis,xaxis='time',yaxis='GainPhase',iteraxis='spw',exprange='all',
#    overwrite=True,showgui=False,
#    plotfile=os.path.join(SB_selfcal_iteration2_folder,f'{prefix}_SB_iteration2_p1_bis_gain_phase_vs_time.png')
#)

#Inspect gain tables and decide whether to flag something
plotms(
    SB_iteration2_p1_bis,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(SB_selfcal_iteration2_folder,f'{prefix}_SB_iteration2_p1_bis_gain_phase_vs_time.png'),
)

#Apply the calibration gains and split-off the corrected ms
applycal(
    vis=SB_iteration2_cont_p1+'.ms',spw=SB_contspws,spwmap=SB_spw_mapping,
    gaintable=[SB_iteration2_p1_bis],interp='linearPD',calwt=True,applymode='calonly'
)
SB_iteration2_cont_p1_bis = SB_iteration2_cont_p1.replace('p1','p1_bis')
os.system('rm -rf %s.ms*' %SB_iteration2_cont_p1_bis)
split(vis=SB_iteration2_cont_p1+'.ms',outputvis=SB_iteration2_cont_p1_bis+'.ms',datacolumn='corrected')

#CLEAN, again with a ~6sigma threshold
tclean_wrapper(
    vis       = SB_iteration2_cont_p1_bis+'.ms',
    imagename = SB_iteration2_cont_p1_bis,
    threshold = '0.0918mJy',
    **SB_tclean_wrapper_kwargs
)
estimate_SNR(SB_iteration2_cont_p1_bis+'.image',disk_mask=SB_mask,noise_mask=noise_annulus_SB)
generate_image_png(
    SB_iteration2_cont_p1_bis+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_iteration2_SB,10*rms_iteration2_SB],
    save_folder=SB_selfcal_iteration2_folder
)
#MWC_758_SB_contp1_bis.image
#Beam 0.193 arcsec x 0.146 arcsec (-28.31 deg)
#Flux inside disk mask: 56.82 mJy
#Peak intensity of source: 5.06 mJy/beam
#rms: 1.53e-02 mJy/beam
#Peak SNR: 330.75

#MWC_758_SB_iteration2_contp1_bis.image
#Beam 0.193 arcsec x 0.146 arcsec (-28.29 deg)
#Flux inside disk mask: 57.16 mJy
#Peak intensity of source: 5.09 mJy/beam
#rms: 1.54e-02 mJy/beam
#Peak SNR: 331.05

#Second round of phase-only self-cal
SB_iteration2_p2 = SB_iteration2_p1_bis.replace('p1_bis','p2')
os.system('rm -rf '+SB_iteration2_p2)
gaincal(
    vis=SB_iteration2_cont_p1_bis+'.ms',caltable=SB_iteration2_p2,
    gaintype='T',combine='scan,spw',calmode='p',solint='360s',
    spw=SB_contspws,refant=SB_refant,
    minsnr=3.,minblperant=4
)
#4 of 43 solutions flagged due to SNR < 3 in spw=0 at 2017/12/09/04:53:08.3
#1 of 43 solutions flagged due to SNR < 3 in spw=0 at 2017/12/09/05:02:34.0
#4 of 43 solutions flagged due to SNR < 3 in spw=0 at 2017/12/09/05:14:50.5
#3 of 43 solutions flagged due to SNR < 3 in spw=0 at 2017/12/09/05:23:12.0
#1 of 45 solutions flagged due to SNR < 3 in spw=4 at 2017/12/17/06:11:47.4
#2 of 46 solutions flagged due to SNR < 3 in spw=8 at 2017/12/27/04:10:33.4
#1 of 46 solutions flagged due to SNR < 3 in spw=8 at 2017/12/27/04:31:22.8
#1 of 46 solutions flagged due to SNR < 3 in spw=8 at 2017/12/27/04:35:39.6
#1 of 46 solutions flagged due to SNR < 3 in spw=8 at 2017/12/27/04:39:01.8

#Inspect gain tables interactively and decide whether to manually flag something
#plotms(SB_iteration2_p2,xaxis='time',yaxis='GainPhase',iteraxis='spw',showgui=True)
#Print calibration png file
#plotms(
#    SB_iteration2_p2,xaxis='time',yaxis='GainPhase',iteraxis='spw',exprange='all',
#    overwrite=True,showgui=False,
#    plotfile=os.path.join(SB_selfcal_iteration2_folder,f'{prefix}_SB_iteration2_p2_gain_phase_vs_time.png')
#)

#Inspect gain tables and decide whether to flag something
plotms(
    SB_iteration2_p2,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(SB_selfcal_iteration2_folder,f'{prefix}_SB_iteration2_p2_gain_phase_vs_time.png'),
)

#Flag problematic antennas
flagdata(vis=SB_iteration2_p2,mode='manual',spw='4', timerange='2017/12/17/06:00:40~2017/12/17/06:01:00',antenna='DV13')
flagdata(vis=SB_iteration2_p2,mode='manual',spw='8', timerange='2017/12/27/04:10:00~2017/12/27/04:11:00',antenna='DA42,DA44,DV13')
flagdata(vis=SB_iteration2_p2,mode='manual',spw='8', timerange='2017/12/27/04:35:00~2017/12/27/04:36:00',antenna='DA44')
flagdata(vis=SB_iteration2_p2,mode='manual',spw='8', timerange='2017/12/27/04:38:00~2017/12/27/04:40:00',antenna='DA44')

#Inspect gain tables to check if flagging worked
plotms(
    SB_iteration2_p2,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(SB_selfcal_iteration2_folder,f'{prefix}_SB_iteration2_p2_gain_phase_vs_time_flagged.png'),
)

#Apply the calibration gains and split-off the corrected ms
applycal(
    vis=SB_iteration2_cont_p1_bis+'.ms',spw=SB_contspws,spwmap=SB_spw_mapping,
    gaintable=[SB_iteration2_p2],interp='linearPD',calwt=True,applymode='calonly'
)
SB_iteration2_cont_p2 = SB_iteration2_cont_p1_bis.replace('p1_bis','p2')
os.system('rm -rf %s.ms*' %SB_iteration2_cont_p2)
split(vis=SB_iteration2_cont_p1_bis+'.ms',outputvis=SB_iteration2_cont_p2+'.ms',datacolumn='corrected')

#CLEAN, again with a ~6sigma threshold
tclean_wrapper(
    vis       = SB_iteration2_cont_p2+'.ms',
    imagename = SB_iteration2_cont_p2,
    threshold = '0.0792mJy',
    **SB_tclean_wrapper_kwargs
)
estimate_SNR(SB_iteration2_cont_p2+'.image',disk_mask=SB_mask,noise_mask=noise_annulus_SB)
generate_image_png(
    SB_iteration2_cont_p2+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_iteration2_SB,10*rms_iteration2_SB],
    save_folder=SB_selfcal_iteration2_folder
)
#MWC_758_SB_contp2.image
#Beam 0.193 arcsec x 0.146 arcsec (-28.31 deg)
#Flux inside disk mask: 57.35 mJy
#Peak intensity of source: 5.34 mJy/beam
#rms: 1.31e-02 mJy/beam
#Peak SNR: 406.88

#MWC_758_SB_iteration2_contp2.image
#Beam 0.193 arcsec x 0.146 arcsec (-28.29 deg)
#Flux inside disk mask: 57.68 mJy
#Peak intensity of source: 5.37 mJy/beam
#rms: 1.32e-02 mJy/beam
#Peak SNR: 407.73

#Third round of phase-only self-cal
SB_iteration2_p3 = SB_iteration2_p2.replace('p2','p3')
os.system('rm -rf '+SB_iteration2_p3)
gaincal(
    vis=SB_iteration2_cont_p2+'.ms',caltable=SB_iteration2_p3,
    gaintype='T',combine='scan,spw',calmode='p',solint='120s', #diff from exoALMA, kept combine='scan' because some scans for the SB_EBs~60sec
    spw=SB_contspws,refant=SB_refant,
    minsnr=3.,minblperant=4
)
#1 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:40:28.5
#1 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:42:29.5
#2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:51:51.1
#4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:53:08.3
#1 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:55:17.4
#1 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:58:55.3
#2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:00:56.3
#2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:02:58.0
#6 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:04:14.3
#1 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:09:31.4
#3 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:11:32.4
#2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:13:34.1
#5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:14:50.5
#1 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:17:53.0
#2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:21:55.7
#2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:23:12.0
#2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:25:41.7
#2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:30:05.7
#1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:49:53.8
#1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:51:57.5
#3 of 43 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:54:57.1
#1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/03:56:58.6
#1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:10:33.4
#8 of 45 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:22:11.4
#1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:31:22.8
#1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:39:01.8
#1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:40:44.5
#1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:57:48.9
#3 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:09:28.0

#Inspect gain tables interactively and decide whether to manually flag something
#plotms(SB_iteration2_p3,xaxis='time',yaxis='GainPhase',iteraxis='spw',showgui=True)
#Print calibration png file
#plotms(
#    SB_iteration2_p3,xaxis='time',yaxis='GainPhase',iteraxis='spw',exprange='all',
#    overwrite=True,showgui=False,
#    plotfile=os.path.join(SB_selfcal_iteration2_folder,f'{prefix}_SB_iteration2_p3_gain_phase_vs_time.png')
#)

#Inspect gain tables and decide whether to flag something
plotms(
    SB_iteration2_p3,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(SB_selfcal_iteration2_folder,f'{prefix}_SB_iteration2_p3_gain_phase_vs_time.png'),
)

#Flag problematic antennas
flagdata(vis=SB_iteration2_p3,mode='manual',spw='0', timerange='2017/12/09/04:58:30~2017/12/09/04:59:00',antenna='DV15')
flagdata(vis=SB_iteration2_p3,mode='manual',spw='0', timerange='2017/12/09/05:04:00~2017/12/09/05:04:30',antenna='DA45,DA49,PM02,DV10,DA50')
flagdata(vis=SB_iteration2_p3,mode='manual',spw='0', timerange='2017/12/09/05:11:00~2017/12/09/05:12:00',antenna='DA46')
flagdata(vis=SB_iteration2_p3,mode='manual',spw='0', timerange='2017/12/09/05:13:00~2017/12/09/05:14:00',antenna='DA55')
flagdata(vis=SB_iteration2_p3,mode='manual',spw='0', timerange='2017/12/09/05:17:30~2017/12/09/05:18:00',antenna='DV10')
flagdata(vis=SB_iteration2_p3,mode='manual',spw='0', timerange='2017/12/09/05:19:30~2017/12/09/05:20:00',antenna='DA46')
flagdata(vis=SB_iteration2_p3,mode='manual',spw='0', timerange='2017/12/09/05:23:00~2017/12/09/05:23:20',antenna='DA55')
flagdata(vis=SB_iteration2_p3,mode='manual',spw='0', timerange='2017/12/09/05:25:30~2017/12/09/05:26:00',antenna='DA51')

flagdata(vis=SB_iteration2_p3,mode='manual',spw='4', timerange='2017/12/17/05:53:30~2017/12/17/05:54:00',antenna='DA43,DA44')
flagdata(vis=SB_iteration2_p3,mode='manual',spw='4', timerange='2017/12/17/05:54:30~2017/12/17/05:55:00',antenna='DA41,DA42,DA43,DA44')
flagdata(vis=SB_iteration2_p3,mode='manual',spw='4', timerange='2017/12/17/05:58:30~2017/12/17/05:59:00',antenna='DA44,DA45,DA50,DA55,DV07,DV13')
flagdata(vis=SB_iteration2_p3,mode='manual',spw='4', timerange='2017/12/17/06:00:30~2017/12/17/06:01:00',antenna='DV13')
flagdata(vis=SB_iteration2_p3,mode='manual',spw='4', timerange='2017/12/17/06:02:30~2017/12/17/06:03:00',antenna='DV07,DV13')
flagdata(vis=SB_iteration2_p3,mode='manual',spw='4', timerange='2017/12/17/06:04:00~2017/12/17/06:04:20',antenna='DV13')
flagdata(vis=SB_iteration2_p3,mode='manual',spw='4', timerange='2017/12/17/06:08:20~2017/12/17/06:08:40',antenna='DA55,DV07,DV11,DA50')
flagdata(vis=SB_iteration2_p3,mode='manual',spw='4', timerange='2017/12/17/06:11:40~2017/12/17/06:12:00',antenna='DV11')

flagdata(vis=SB_iteration2_p3,mode='manual',spw='8', timerange='2017/12/27/04:09:00~2017/12/27/04:09:30',antenna='DA42,DV13')
flagdata(vis=SB_iteration2_p3,mode='manual',spw='8', timerange='2017/12/27/04:10:30~2017/12/27/04:11:00',antenna='DA42,DA44,DV13')
flagdata(vis=SB_iteration2_p3,mode='manual',spw='8', timerange='2017/12/27/04:33:30~2017/12/27/04:34:00',antenna='DA44')
flagdata(vis=SB_iteration2_p3,mode='manual',spw='8', timerange='2017/12/27/04:35:00~2017/12/27/04:36:00',antenna='DA44,DV11')
flagdata(vis=SB_iteration2_p3,mode='manual',spw='8', timerange='2017/12/27/04:37:00~2017/12/27/04:38:00',antenna='DA44')
flagdata(vis=SB_iteration2_p3,mode='manual',spw='8', timerange='2017/12/27/04:39:00~2017/12/27/04:40:00',antenna='DA44')

flagdata(vis=SB_iteration2_p3,mode='manual',spw='12',timerange='2017/12/28/03:56:20~2017/12/28/03:56:40',antenna='DA55')
flagdata(vis=SB_iteration2_p3,mode='manual',spw='12',timerange='2017/12/28/03:57:30~2017/12/28/03:58:00',antenna='DA55')
flagdata(vis=SB_iteration2_p3,mode='manual',spw='12',timerange='2017/12/28/03:59:20~2017/12/28/03:59:40',antenna='DA55')
flagdata(vis=SB_iteration2_p3,mode='manual',spw='12',timerange='2017/12/28/04:04:00~2017/12/28/04:04:30',antenna='DA44,DA55,DA46,DV07')
flagdata(vis=SB_iteration2_p3,mode='manual',spw='12',timerange='2017/12/28/04:09:00~2017/12/28/04:09:30',antenna='DA42,DV13')
flagdata(vis=SB_iteration2_p3,mode='manual',spw='12',timerange='2017/12/28/04:15:00~2017/12/28/04:15:30',antenna='DA44,DV11,DV13')
flagdata(vis=SB_iteration2_p3,mode='manual',spw='12',timerange='2017/12/28/04:20:30~2017/12/28/04:21:00',antenna='DA44,DV11,DV13')
flagdata(vis=SB_iteration2_p3,mode='manual',spw='12',timerange='2017/12/28/04:30:00~2017/12/28/04:30:30',antenna='DV13')

#Inspect gain tables to check if flagging worked
plotms(
    SB_iteration2_p3,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(SB_selfcal_iteration2_folder,f'{prefix}_SB_iteration2_p3_gain_phase_vs_time_flagged.png'),
)

#Apply the calibration gains and split-off the corrected ms
applycal(
    vis=SB_iteration2_cont_p2+'.ms',spw=SB_contspws,spwmap=SB_spw_mapping,
    gaintable=[SB_iteration2_p3],interp='linearPD',calwt=True,applymode='calonly'
)
SB_iteration2_cont_p3 = SB_iteration2_cont_p2.replace('p2','p3')
os.system('rm -rf %s.ms*' %SB_iteration2_cont_p3)
split(vis=SB_iteration2_cont_p2+'.ms',outputvis=SB_iteration2_cont_p3+'.ms',datacolumn='corrected')

#CLEAN, again with a ~6sigma threshold
tclean_wrapper(
    vis       = SB_iteration2_cont_p3+'.ms',
    imagename = SB_iteration2_cont_p3,
    threshold = '0.0768mJy',
    **SB_tclean_wrapper_kwargs
)
estimate_SNR(SB_iteration2_cont_p3+'.image',disk_mask=SB_mask,noise_mask=noise_annulus_SB)
generate_image_png(
    SB_iteration2_cont_p3+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_iteration2_SB,10*rms_iteration2_SB],
    save_folder=SB_selfcal_iteration2_folder
)
#MWC_758_SB_contp3.image
#Beam 0.193 arcsec x 0.146 arcsec (-28.31 deg)
#Flux inside disk mask: 57.63 mJy
#Peak intensity of source: 5.52 mJy/beam
#rms: 1.28e-02 mJy/beam
#Peak SNR: 432.31

#MWC_758_SB_iteration2_contp3.image
#Beam 0.193 arcsec x 0.146 arcsec (-28.29 deg)
#Flux inside disk mask: 57.95 mJy
#Peak intensity of source: 5.55 mJy/beam
#rms: 1.28e-02 mJy/beam
#Peak SNR: 432.78

#Fourth round of phase-only self-cal
SB_iteration2_p4 = SB_iteration2_p3.replace('p3','p4')
os.system('rm -rf '+SB_iteration2_p4)
gaincal(
    vis=SB_iteration2_cont_p3+'.ms',caltable=SB_iteration2_p4,
    gaintype='T',combine='spw',calmode='p',solint='60s',
    spw=SB_contspws,refant=SB_refant,
    minsnr=3.,minblperant=4
)
#2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:39:58.6
#2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:40:59.0
#2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:41:59.5
#2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:43:00.0
#1 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:44:00.5
#1 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:45:02.5
#2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:48:19.3
#2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:50:20.3
#3 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:51:20.7
#3 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:52:22.8
#4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:53:08.3
#2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:55:02.5
#5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:55:48.0
#4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:58:25.0
#2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:59:25.5
#3 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:00:25.9
#3 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:01:26.4
#2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:02:26.9
#4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:03:28.9
#4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:04:14.3
#4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:09:01.2
#2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:10:01.6
#5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:11:02.1
#1 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:12:02.6
#3 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:13:03.1
#2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:14:05.1
#5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:14:50.5
#1 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:17:22.7
#3 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:18:23.2
#2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:19:23.7
#2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:21:24.6
#2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:22:26.6
#3 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:23:12.0
#3 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:25:41.7
#3 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:28:50.3
#5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:30:36.2
#1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:35:40.8
#2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:50:24.3
#1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:51:25.3
#1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:52:27.4
#2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:02:54.0
#1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:07:45.5
#1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:12:47.6
#2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:17:39.2
#1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:21:59.3
#1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:31:22.8
#1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:34:13.0
#1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:35:13.5
#1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:36:14.0
#1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:37:14.3
#1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:38:16.3
#1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:39:01.8
#1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:40:29.5
#1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:41:14.9
#1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:57:48.9
#1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:09:16.6

#Inspect gain tables interactively and decide whether to manually flag something
#plotms(SB_iteration2_p4,xaxis='time',yaxis='GainPhase',iteraxis='spw',showgui=True)
#Print calibration png file
#plotms(
#    SB_iteration2_p4,xaxis='time',yaxis='GainPhase',iteraxis='spw',exprange='all',
#    overwrite=True,showgui=False,
#    plotfile=os.path.join(SB_selfcal_iteration2_folder,f'{prefix}_SB_iteration2_p4_gain_phase_vs_time.png')
#)

#Inspect gain tables and decide whether to flag something
plotms(
    SB_iteration2_p4,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(SB_selfcal_iteration2_folder,f'{prefix}_SB_iteration2_p4_gain_phase_vs_time.png'),
)

#Flag problematic antennas
flagdata(vis=SB_iteration2_p4,mode='manual',spw='0', timerange='2017/12/09/04:51:20~2017/12/09/04:51:25',antenna='DV15')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='0', timerange='2017/12/09/04:52:20~2017/12/09/04:52:25',antenna='DV10')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='0', timerange='2017/12/09/04:54:30~2017/12/09/04:55:30',antenna='DV15')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='0', timerange='2017/12/09/04:55:30~2017/12/09/04:56:00',antenna='DV11')
## flagdata(vis=SB_iteration2_p4,mode='manual',spw='0', timerange='2017/12/09/05:00:20~2017/12/09/05:00:30',antenna='DV10')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='0', timerange='2017/12/09/05:04:00~2017/12/09/05:04:20',antenna='DA45,DA50,PM02,DA49,DV10')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='0', timerange='2017/12/09/05:14:00~2017/12/09/05:14:20',antenna='DA55')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='0', timerange='2017/12/09/05:17:20~2017/12/09/05:17:40',antenna='DV10')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='0', timerange='2017/12/09/05:19:20~2017/12/09/05:19:30',antenna='DV15')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='0', timerange='2017/12/09/05:20:20~2017/12/09/05:20:40',antenna='DA46')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='0', timerange='2017/12/09/05:21:20~2017/12/09/05:21:40',antenna='DA46')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='0', timerange='2017/12/09/05:22:00~2017/12/09/05:23:00',antenna='DA46')
## flagdata(vis=SB_iteration2_p4,mode='manual',spw='0', timerange='2017/12/09/05:23:10~2017/12/09/05:23:20',antenna='DA55,DV10')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='0', timerange='2017/12/09/05:25:00~2017/12/09/05:26:00',antenna='DA51,DA46')

flagdata(vis=SB_iteration2_p4,mode='manual',spw='4', timerange='2017/12/17/05:49:20~2017/12/17/05:49:30',antenna='DA55')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='4', timerange='2017/12/17/05:53:10~2017/12/17/05:53:30',antenna='DA55')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='4', timerange='2017/12/17/05:53:30~2017/12/17/05:54:00',antenna='DA43,DA44')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='4', timerange='2017/12/17/05:54:30~2017/12/17/05:55:00',antenna='DA41,DA42,DA43,DA44,DA55,DV13')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='4', timerange='2017/12/17/05:58:10~2017/12/17/05:58:20',antenna='DA44,DA45,DA55,DA50,DA55,DV07,DV13')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='4', timerange='2017/12/17/05:59:10~2017/12/17/05:59:20',antenna='DA44,DA45,DA55,DA50,DA55,DV07,DV13')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='4', timerange='2017/12/17/06:00:10~2017/12/17/06:00:20',antenna='DV13')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='4', timerange='2017/12/17/06:01:20~2017/12/17/06:01:30',antenna='DV07,DV13')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='4', timerange='2017/12/17/06:02:20~2017/12/17/06:02:30',antenna='DV07,DV13')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='4', timerange='2017/12/17/06:03:20~2017/12/17/06:03:30',antenna='DV07,DV13')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='4', timerange='2017/12/17/06:04:00~2017/12/17/06:04:10',antenna='DV13')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='4', timerange='2017/12/17/06:08:40~2017/12/17/06:09:00',antenna='DA55,DV07,DA50')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='4', timerange='2017/12/17/06:09:40~2017/12/17/06:10:00',antenna='DA55,DV07')

flagdata(vis=SB_iteration2_p4,mode='manual',spw='8', timerange='2017/12/27/04:08:30~2017/12/27/04:09:00',antenna='DA42,DV13')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='8', timerange='2017/12/27/04:09:30~2017/12/27/04:10:00',antenna='DA42,DA44,DV11,DV13')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='8', timerange='2017/12/27/04:10:30~2017/12/27/04:11:00',antenna='DA42,DA44,DV13')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='8', timerange='2017/12/27/04:25:30~2017/12/27/04:25:40',antenna='DV11')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='8', timerange='2017/12/27/04:33:00~2017/12/27/04:33:20',antenna='DA44')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='8', timerange='2017/12/27/04:34:00~2017/12/27/04:34:20',antenna='DA44,DV11')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='8', timerange='2017/12/27/04:35:00~2017/12/27/04:35:20',antenna='DA44')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='8', timerange='2017/12/27/04:36:00~2017/12/27/04:36:20',antenna='DA44')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='8', timerange='2017/12/27/04:37:00~2017/12/27/04:37:20',antenna='DA44')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='8', timerange='2017/12/27/04:38:00~2017/12/27/04:38:20',antenna='DA44')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='8', timerange='2017/12/27/04:39:00~2017/12/27/04:39:20',antenna='DA44')

flagdata(vis=SB_iteration2_p4,mode='manual',spw='12',timerange='2017/12/28/03:46:20~2017/12/28/03:46:30',antenna='DA46')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='12',timerange='2017/12/28/03:55:00~2017/12/28/03:55:10',antenna='DA55')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='12',timerange='2017/12/28/03:56:00~2017/12/28/03:56:20',antenna='DV11,DA55')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='12',timerange='2017/12/28/03:57:00~2017/12/28/03:57:10',antenna='DV11,DA55')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='12',timerange='2017/12/28/03:59:00~2017/12/28/03:59:30',antenna='DA55')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='12',timerange='2017/12/28/04:00:00~2017/12/28/04:00:20',antenna='DA44,DA55,DV07')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='12',timerange='2017/12/28/04:02:50~2017/12/28/04:03:00',antenna='DA44,DV07,DA46')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='12',timerange='2017/12/28/04:03:40~2017/12/28/04:04:00',antenna='DA44,DA55,DA46,DV07')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='12',timerange='2017/12/28/04:04:40~2017/12/28/04:05:00',antenna='DA44')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='12',timerange='2017/12/28/04:09:10~2017/12/28/04:09:20',antenna='DA42,DV13')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='12',timerange='2017/12/28/04:21:20~2017/12/28/04:21:40',antenna='DA44,DV11,DV13')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='12',timerange='2017/12/28/04:22:20~2017/12/28/04:22:40',antenna='DV11,DV13')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='12',timerange='2017/12/28/04:30:10~2017/12/28/04:30:30',antenna='DV13')

#Inspect gain tables to check if flagging worked
plotms(
    SB_iteration2_p4,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(SB_selfcal_iteration2_folder,f'{prefix}_SB_iteration2_p4_gain_phase_vs_time_flagged.png'),
)

#Apply the calibration gains and split-off the corrected ms
applycal(
    vis=SB_iteration2_cont_p3+'.ms',spw=SB_contspws,spwmap=SB_spw_mapping,
    gaintable=[SB_iteration2_p4],interp='linearPD',calwt=True,applymode='calonly'
)
SB_iteration2_cont_p4 = SB_iteration2_cont_p3.replace('p3','p4')
os.system('rm -rf %s.ms*' %SB_iteration2_cont_p4)
split(vis=SB_iteration2_cont_p3+'.ms',outputvis=SB_iteration2_cont_p4+'.ms',datacolumn='corrected')

#CLEAN, again with a ~6sigma threshold
tclean_wrapper(
    vis       = SB_iteration2_cont_p4+'.ms',
    imagename = SB_iteration2_cont_p4,
    threshold = '0.0768mJy',
    **SB_tclean_wrapper_kwargs
)
estimate_SNR(SB_iteration2_cont_p4+'.image',disk_mask=SB_mask,noise_mask=noise_annulus_SB)
generate_image_png(
    SB_iteration2_cont_p4+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_iteration2_SB,10*rms_iteration2_SB],
    save_folder=SB_selfcal_iteration2_folder
)
#MWC_758_SB_contp4.image
#Beam 0.193 arcsec x 0.146 arcsec (-28.31 deg)
#Flux inside disk mask: 57.80 mJy
#Peak intensity of source: 5.60 mJy/beam
#rms: 1.28e-02 mJy/beam
#Peak SNR: 437.52

#MWC_758_SB_iteration2_contp4.image
#Beam 0.193 arcsec x 0.146 arcsec (-28.29 deg)
#Flux inside disk mask: 58.17 mJy
#Peak intensity of source: 5.63 mJy/beam
#rms: 1.29e-02 mJy/beam
#Peak SNR: 437.79

#Fifth round of phase-only selfcal
SB_iteration2_p5 = SB_iteration2_p4.replace('p4','p5')
os.system('rm -rf '+SB_iteration2_p5)
gaincal(
    vis=SB_iteration2_cont_p4+'.ms',caltable=SB_iteration2_p5,
    gaintype='T',combine='spw',calmode='p',solint='20s',
    spw=SB_contspws,refant=SB_refant,
    minsnr=3.,minblperant=4
)
# 4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:39:40.5
# 4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:40:01.7
# 5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:40:19.8
# 3 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:40:41.0
# 5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:41:02.2
# 4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:41:20.3
# 5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:41:41.5
# 4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:42:02.6
# 3 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:42:20.8
# 5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:42:41.9
# 4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:43:03.1
# 4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:43:21.3
# 3 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:43:39.4
# 2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:44:00.6
# 2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:44:21.7
# 3 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:44:41.4
# 4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:44:59.6
# 2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:45:20.7
# 2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:47:00.8
# 4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:47:22.0
# 5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:47:40.2
# 5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:48:01.3
# 4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:48:22.5
# 6 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:48:40.6
# 3 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:49:01.8
# 6 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:49:23.0
# 5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:49:41.1
# 2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:50:02.3
# 5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:50:23.5
# 6 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:50:41.6
# 2 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:50:59.7
# 4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:51:20.9
# 9 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:51:42.1
# 5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:52:01.8
# 5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:52:19.9
# 5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:52:41.1
# 7 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:53:02.2
# 8 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:53:17.4
# 4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:54:44.5
# 6 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:55:05.7
# 7 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:55:23.8
# 6 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:55:45.0
#13 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:56:00.1
# 5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:58:06.8
# 7 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:58:28.0
# 7 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:58:46.1
# 4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:59:07.3
# 6 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:59:28.5
# 7 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/04:59:46.6
# 6 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:00:07.8
# 8 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:00:28.9
# 6 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:00:47.1
# 4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:01:08.3
# 7 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:01:29.4
# 3 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:01:47.6
# 9 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:02:05.7
# 3 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:02:26.9
# 8 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:02:48.0
# 9 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:03:07.7
# 9 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:03:25.9
# 7 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:03:47.0
#10 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:04:08.2
#10 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:04:23.3
# 5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:08:43.0
# 7 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:09:04.2
# 6 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:09:22.3
# 4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:09:43.5
# 6 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:10:04.7
# 4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:10:22.8
# 5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:10:44.0
# 6 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:11:05.1
# 6 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:11:23.3
# 5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:11:44.5
# 9 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:12:05.6
# 6 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:12:23.8
# 9 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:12:41.9
# 4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:13:03.1
# 6 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:13:24.2
# 7 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:13:43.9
# 4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:14:02.1
# 6 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:14:23.2
# 8 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:14:44.4
# 9 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:14:59.5
# 4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:17:04.6
# 8 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:17:25.7
# 5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:17:43.9
# 4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:18:05.0
# 7 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:18:26.2
# 7 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:18:44.4
# 5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:19:05.5
# 9 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:19:26.7
# 6 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:19:44.8
# 6 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:20:06.0
# 3 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:20:27.2
# 8 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:20:45.3
# 6 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:21:03.5
# 5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:21:24.6
# 3 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:21:45.8
# 6 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:22:05.5
# 4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:22:23.6
# 3 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:22:44.8
# 3 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:23:06.0
#10 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:23:21.1
# 3 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:25:23.5
# 6 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:25:44.7
# 7 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:26:02.8
# 7 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:27:31.8
# 5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:27:53.0
# 9 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:28:11.1
# 5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:28:32.3
# 4 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:28:53.5
# 7 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:29:11.6
# 5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:29:32.8
# 6 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:29:53.9
# 5 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:30:12.1
# 3 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:30:33.2
#11 of 43 solutions flagged due to SNR < 3 in spw=0  at 2017/12/09/05:30:48.4
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:29:33.5
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:29:54.6
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:30:12.8
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:30:33.9
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:30:55.1
# 4 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:31:13.2
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:31:34.4
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:31:55.6
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:32:13.7
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:32:34.9
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:32:56.1
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:33:32.3
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:33:53.5
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:34:14.7
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:34:34.4
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:34:52.6
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:35:13.7
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:35:34.7
# 3 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:35:50.1
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:37:12.1
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:37:33.3
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:37:51.4
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:38:12.6
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:38:33.7
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:38:51.9
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:39:13.1
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:39:34.2
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:39:52.2
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:40:34.6
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:40:52.8
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:41:11.0
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:41:53.3
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:42:13.0
# 3 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:42:31.2
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:42:52.3
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:43:13.5
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:43:28.5
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:44:30.9
# 6 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:45:46.5
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:47:26.4
# 3 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:47:44.6
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:48:05.8
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:48:26.9
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:48:45.1
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:49:27.3
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:49:45.6
# 4 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:50:06.6
# 4 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:50:27.9
# 4 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:50:45.9
# 4 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:51:04.3
# 3 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:51:25.2
# 4 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:51:46.3
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:52:06.3
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:52:24.3
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:52:45.5
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:53:06.7
# 5 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:53:21.7
# 3 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:54:02.8
# 3 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:54:21.0
# 3 of 43 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:54:57.1
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:58:00.6
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:58:21.8
# 3 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:58:39.9
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:59:01.1
# 4 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:59:22.3
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/05:59:40.4
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:00:01.6
# 3 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:00:22.6
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:00:40.7
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:01:23.2
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:01:41.0
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:01:59.6
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:02:20.6
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:02:41.9
# 3 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:03:01.4
# 3 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:03:19.5
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:03:40.8
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:04:02.2
# 3 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:04:17.3
# 3 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:06:01.1
# 3 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:06:19.2
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:06:40.4
# 4 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:07:01.5
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:07:19.7
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:07:40.5
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:08:02.5
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:08:20.6
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:08:41.4
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:09:02.5
# 3 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:09:20.6
# 3 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:09:59.6
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:10:21.1
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:10:40.8
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:10:58.9
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:11:20.4
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:11:41.1
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:11:56.4
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:12:58.9
# 2 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:13:20.1
# 3 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:13:38.2
# 6 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:14:14.5
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:15:34.4
# 3 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:15:55.6
# 1 of 45 solutions flagged due to SNR < 3 in spw=4  at 2017/12/17/06:16:13.7
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/03:56:46.6
# 4 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/03:57:07.7
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/03:57:25.9
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/03:57:47.0
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/03:58:08.2
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/03:58:26.3
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/03:58:47.5
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/03:59:08.7
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:00:27.3
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:00:45.5
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:01:06.6
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:01:27.8
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:01:47.5
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:02:05.6
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:02:26.8
# 4 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:02:47.9
# 4 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:03:03.1
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:04:26.0
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:04:47.1
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:05:05.3
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:05:47.6
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:06:05.8
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:06:26.9
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:06:48.1
# 5 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:07:06.2
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:07:27.4
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:07:48.6
# 4 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:08:06.7
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:08:24.9
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:09:07.2
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:09:26.9
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:09:45.0
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:10:06.2
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:10:27.4
# 4 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:10:42.5
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:11:44.1
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:12:05.2
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:12:23.4
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:12:44.5
# 5 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:12:59.7
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:14:19.6
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:14:40.8
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:14:58.9
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:15:20.1
# 4 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:15:41.3
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:15:59.4
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:16:20.6
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:16:41.8
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:16:59.9
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:17:21.1
# 4 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:17:42.3
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:18:00.4
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:18:18.5
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:18:39.7
# 4 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:19:00.9
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:19:20.5
# 4 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:19:38.7
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:19:59.9
# 4 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:20:36.1
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:20:56.0
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:21:17.0
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:21:35.1
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:21:56.3
# 8 of 45 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:22:11.4
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:25:15.4
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:25:36.5
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:25:54.7
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:26:15.8
# 4 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:26:37.0
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:26:55.1
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:27:16.3
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:27:37.5
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:28:16.8
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:28:38.0
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:28:56.1
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:29:35.4
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:29:56.6
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:30:16.5
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:30:34.4
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:30:55.6
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:31:16.5
# 4 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:31:31.8
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:32:54.5
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:33:15.6
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:33:33.8
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:33:54.9
# 4 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:34:16.1
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:34:34.3
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:34:55.4
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:35:16.6
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:35:34.7
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:35:55.9
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:36:17.1
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:36:35.2
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:36:53.4
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:37:14.6
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:37:35.4
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:37:55.5
# 4 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:38:13.4
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:38:34.7
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:38:55.9
# 4 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:39:11.1
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:40:11.4
# 4 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:40:32.5
# 5 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:40:50.7
# 1 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:41:11.8
# 6 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:41:27.0
# 2 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:43:06.9
# 3 of 46 solutions flagged due to SNR < 3 in spw=8  at 2017/12/27/04:43:25.1
# 2 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:44:03.0
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:45:03.6
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:45:24.7
# 2 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:45:42.9
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:47:25.7
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:47:43.8
# 2 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:48:02.0
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:48:23.1
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:48:44.3
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:49:22.1
# 2 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:50:04.5
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:52:02.7
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:52:20.8
# 2 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:54:03.7
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:54:21.8
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:55:40.4
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:56:01.6
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:57:00.6
# 2 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:57:42.9
# 2 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:57:58.0
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:59:00.0
# 2 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/03:59:21.1
# 3 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:00:15.6
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:01:36.3
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:03:58.4
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:05:35.1
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:06:55.3
# 2 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:07:52.8
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:09:13.5
# 3 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:09:28.0
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:13:11.1
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:14:54.0
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:16:12.6
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:17:13.1
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:18:48.3
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:20:31.2
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:20:49.3
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:21:31.7
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:22:11.0
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:24:51.5
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:25:11.3
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:25:28.9
# 2 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:26:26.5
# 3 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:28:43.0
# 1 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/12/28/04:30:22.5

#Inspect gain tables interactively and decide whether to manually flag something
#plotms(SB_iteration2_p5,xaxis='time',yaxis='GainPhase',iteraxis='spw',showgui=True)
#Print calibration png file
#plotms(
#    SB_iteration2_p5,xaxis='time',yaxis='GainPhase',iteraxis='spw',exprange='all',
#    overwrite=True,showgui=False,
#    plotfile=os.path.join(SB_selfcal_iteration2_folder,f'{prefix}_SB_iteration2_p5_gain_phase_vs_time.png')
#)

#Inspect gain tables and decide whether to flag something
plotms(
    SB_iteration2_p5,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(SB_selfcal_iteration2_folder,f'{prefix}_SB_iteration2_p5_gain_phase_vs_time.png'),
)

#Flag problematic antennas
flagdata(vis=SB_iteration2_p5,mode='manual',spw='0', timerange='2017/12/09/04:40:35~2017/12/09/04:40:45',antenna='DV10')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='0', timerange='2017/12/09/04:41:00~2017/12/09/04:41:05',antenna='DV11')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='0', timerange='2017/12/09/04:47:39~2017/12/09/04:47:41',antenna='DA46')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='0', timerange='2017/12/09/04:54:40~2017/12/09/04:54:50',antenna='DV12,DV11,DV15')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='0', timerange='2017/12/09/05:02:05~2017/12/09/05:02:10',antenna='DA45,DA50,DA51')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='0', timerange='2017/12/09/05:02:45~2017/12/09/05:02:50',antenna='DA45')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='0', timerange='2017/12/09/05:03:20~2017/12/09/05:03:30',antenna='DA49')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='0', timerange='2017/12/09/05:04:20~2017/12/09/05:04:25',antenna='DA45,DA49,DA46')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='0', timerange='2017/12/09/05:08:40~2017/12/09/05:08:45',antenna='DA45,DA49,DA50,DV07,PM02')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='0', timerange='2017/12/09/05:09:00~2017/12/09/05:09:05',antenna='PM02')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='0', timerange='2017/12/09/05:17:40~2017/12/09/05:17:50',antenna='DV10')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='0', timerange='2017/12/09/05:18:00~2017/12/09/05:18:10',antenna='DV10')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='0', timerange='2017/12/09/05:20:20~2017/12/09/05:20:40',antenna='DA46')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='0', timerange='2017/12/09/05:20:40~2017/12/09/05:20:50',antenna='DA50')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='0', timerange='2017/12/09/05:22:20~2017/12/09/05:22:25',antenna='DA46')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='0', timerange='2017/12/09/05:22:40~2017/12/09/05:22:50',antenna='DA55')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='0', timerange='2017/12/09/05:25:20~2017/12/09/05:25:40',antenna='DA51')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='0', timerange='2017/12/09/05:25:40~2017/12/09/05:25:50',antenna='DA51')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='0', timerange='2017/12/09/05:26:00~2017/12/09/05:26:10',antenna='DA51')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='0', timerange='2017/12/09/05:30:30~2017/12/09/05:30:35',antenna='DA46')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='0', timerange='2017/12/09/05:30:40~2017/12/09/05:30:50',antenna='DA55')

flagdata(vis=SB_iteration2_p5,mode='manual',spw='4', timerange='2017/12/17/05:49:05~2017/12/17/05:49:10',antenna='DA55')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='4', timerange='2017/12/17/05:53:40~2017/12/17/05:53:45',antenna='DA43,DA44')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='4', timerange='2017/12/17/05:54:00~2017/12/17/05:54:05',antenna='DA43,DA44')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='4', timerange='2017/12/17/05:54:20~2017/12/17/05:54:25',antenna='DA43')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='4', timerange='2017/12/17/05:54:40~2017/12/17/05:54:45',antenna='DA41,DA42,DA43,DA44,DA55,DV13')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='4', timerange='2017/12/17/05:54:55~2017/12/17/05:55:00',antenna='DA41,DA42,DA43,DA44')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='4', timerange='2017/12/17/05:58:00~2017/12/17/05:58:05',antenna='DA44,DV07,DV13')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='4', timerange='2017/12/17/05:58:20~2017/12/17/05:58:25',antenna='DA44,DA45,DA50,DA55,DV07,DV13')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='4', timerange='2017/12/17/05:58:20~2017/12/17/05:58:25',antenna='DA44,DA45,DA50,DA55,DV07,DV13')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='4', timerange='2017/12/17/05:58:35~2017/12/17/05:58:40',antenna='DA45,DA50,DV07,DV13')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='4', timerange='2017/12/17/05:59:00~2017/12/17/05:59:05',antenna='DA44,DA45,DA50,DV07,DV13')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='4', timerange='2017/12/17/05:59:20~2017/12/17/05:59:30',antenna='DA45,DA50,DV07')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='4', timerange='2017/12/17/05:59:40~2017/12/17/05:59:45',antenna='DV07,DV13')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='4', timerange='2017/12/17/06:00:00~2017/12/17/06:00:10',antenna='DV13')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='4', timerange='2017/12/17/06:00:20~2017/12/17/06:00:25',antenna='DV13')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='4', timerange='2017/12/17/06:00:30~2017/12/17/06:00:45',antenna='DV13')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='4', timerange='2017/12/17/06:01:00~2017/12/17/06:01:05',antenna='DV13')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='4', timerange='2017/12/17/06:01:20~2017/12/17/06:01:25',antenna='DV07,DV13')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='4', timerange='2017/12/17/06:01:40~2017/12/17/06:01:45',antenna='DV07,DV13')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='4', timerange='2017/12/17/06:01:55~2017/12/17/06:02:00',antenna='DV07,DV13')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='4', timerange='2017/12/17/06:02:10~2017/12/17/06:02:30',antenna='DV07,DV13')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='4', timerange='2017/12/17/06:02:40~2017/12/17/06:02:50',antenna='DV07,DV13')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='4', timerange='2017/12/17/06:03:00~2017/12/17/06:03:10',antenna='DV07,DV13')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='4', timerange='2017/12/17/06:03:10~2017/12/17/06:03:20',antenna='DV13')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='4', timerange='2017/12/17/06:03:40~2017/12/17/06:03:50',antenna='DV07,DV13')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='4', timerange='2017/12/17/06:04:00~2017/12/17/06:04:10',antenna='DV13')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='4', timerange='2017/12/17/06:08:40~2017/12/17/06:08:50',antenna='DA55,DV07')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='4', timerange='2017/12/17/06:09:00~2017/12/17/06:09:10',antenna='DA55,DV07')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='4', timerange='2017/12/17/06:09:20~2017/12/17/06:09:30',antenna='DV07')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='4', timerange='2017/12/17/06:09:30~2017/12/17/06:09:40',antenna='DA45,DV07,DA50')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='4', timerange='2017/12/17/06:04:00~2017/12/17/06:04:10',antenna='DV13')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='4', timerange='2017/12/17/06:09:50~2017/12/17/06:10:00',antenna='DV07')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='4', timerange='2017/12/17/06:10:20~2017/12/17/06:10:25',antenna='DV07')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='4', timerange='2017/12/17/06:11:50~2017/12/17/06:12:00',antenna='DV11')

flagdata(vis=SB_iteration2_p5,mode='manual',spw='8', timerange='2017/12/27/04:07:40~2017/12/27/04:07:50',antenna='DA42')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8', timerange='2017/12/27/04:08:00~2017/12/27/04:08:10',antenna='DA42')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8', timerange='2017/12/27/04:08:20~2017/12/27/04:08:30',antenna='DA42,DV13')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8', timerange='2017/12/27/04:08:40~2017/12/27/04:08:50',antenna='DA42,DV13')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8', timerange='2017/12/27/04:09:00~2017/12/27/04:09:10',antenna='DA42,DV13')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8', timerange='2017/12/27/04:09:20~2017/12/27/04:09:30',antenna='DA42,DV13')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8', timerange='2017/12/27/04:09:40~2017/12/27/04:09:50',antenna='DA42,DA44,DV13')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8', timerange='2017/12/27/04:10:00~2017/12/27/04:10:10',antenna='DA42,DA44,DV13')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8', timerange='2017/12/27/04:10:20~2017/12/27/04:10:30',antenna='DA42,DA44,DA55,DV13')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8', timerange='2017/12/27/04:10:40~2017/12/27/04:10:50',antenna='DA42,DA44,DV13')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8', timerange='2017/12/27/04:11:40~2017/12/27/04:11:50',antenna='DA42,DA44,DV13')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8', timerange='2017/12/27/04:25:30~2017/12/27/04:25:40',antenna='DV11')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8', timerange='2017/12/27/04:25:50~2017/12/27/04:26:00',antenna='DA55')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8', timerange='2017/12/27/04:29:50~2017/12/27/04:30:00',antenna='DV07,DV11')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8', timerange='2017/12/27/04:33:10~2017/12/27/04:33:20',antenna='DA44')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8', timerange='2017/12/27/04:33:30~2017/12/27/04:33:40',antenna='DA44,DA45')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8', timerange='2017/12/27/04:34:10~2017/12/27/04:34:20',antenna='DA44')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8', timerange='2017/12/27/04:34:50~2017/12/27/04:35:00',antenna='DA44')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8', timerange='2017/12/27/04:35:30~2017/12/27/04:35:40',antenna='DA44')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8', timerange='2017/12/27/04:35:50~2017/12/27/04:36:00',antenna='DA44')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8', timerange='2017/12/27/04:36:10~2017/12/27/04:36:20',antenna='DA44')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8', timerange='2017/12/27/04:36:30~2017/12/27/04:36:40',antenna='DA44')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8', timerange='2017/12/27/04:37:10~2017/12/27/04:37:20',antenna='DA44')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8', timerange='2017/12/27/04:37:30~2017/12/27/04:37:40',antenna='DA44')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8', timerange='2017/12/27/04:37:30~2017/12/27/04:37:40',antenna='DA44')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8', timerange='2017/12/27/04:37:50~2017/12/27/04:38:00',antenna='DA44')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8', timerange='2017/12/27/04:38:30~2017/12/27/04:38:40',antenna='DA44,DV13')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8', timerange='2017/12/27/04:38:50~2017/12/27/04:39:00',antenna='DA44')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8', timerange='2017/12/27/04:39:10~2017/12/27/04:39:20',antenna='DA44')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8', timerange='2017/12/27/04:41:10~2017/12/27/04:41:20',antenna='DV11')

flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/03:45:40~2017/12/28/03:45:45',antenna='DA55')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/03:46:40~2017/12/28/03:46:50',antenna='DA55,DA46')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/03:48:00~2017/12/28/03:48:10',antenna='DA46')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/03:49:00~2017/12/28/03:49:10',antenna='DA44')
## flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/03:54:40~2017/12/28/03:54:45',antenna='DA55')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/03:55:00~2017/12/28/03:55:10',antenna='DA55')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/03:55:20~2017/12/28/03:55:30',antenna='DA55,DV11')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/03:55:40~2017/12/28/03:55:50',antenna='DA55')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/03:56:00~2017/12/28/03:56:10',antenna='DA55,DV07')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/03:56:20~2017/12/28/03:56:30',antenna='DA55,DV11')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/03:56:40~2017/12/28/03:56:50',antenna='DA55,DV11')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/03:57:20~2017/12/28/03:57:25',antenna='DA55')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/03:58:50~2017/12/28/03:59:00',antenna='DA45,DA50,DA55')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/03:59:10~2017/12/28/03:59:30',antenna='DA50')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/03:59:30~2017/12/28/03:59:40',antenna='DA55')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/03:59:50~2017/12/28/04:00:10',antenna='DA44,DA55,DV07')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/04:00:10~2017/12/28/04:00:20',antenna='DA42,DA44,DV07')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/04:01:30~2017/12/28/04:01:40',antenna='DA44,DA56,DV07,DA46')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/04:01:50~2017/12/28/04:02:00',antenna='DV07')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/04:02:10~2017/12/28/04:02:20',antenna='DA55,DA48,DA51')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/04:03:10~2017/12/28/04:03:20',antenna='DA44')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/04:03:20~2017/12/28/04:03:40',antenna='DA44,DA55,DA46')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/04:03:50~2017/12/28/04:04:00',antenna='DA44,DA46')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/04:04:10~2017/12/28/04:04:20',antenna='DA44')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/04:04:30~2017/12/28/04:04:40',antenna='DA44')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/04:04:50~2017/12/28/04:05:00',antenna='DA44')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/04:05:10~2017/12/28/04:05:20',antenna='DA44')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/04:09:10~2017/12/28/04:09:20',antenna='DA42,DV11,DV13')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/04:09:20~2017/12/28/04:09:30',antenna='DA42,DV13')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/04:12:20~2017/12/28/04:12:40',antenna='DA42,DA55,DV13')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/04:12:50~2017/12/28/04:13:00',antenna='DA55')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/04:13:10~2017/12/28/04:13:20',antenna='DV11')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/04:13:30~2017/12/28/04:13:35',antenna='DV11')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/04:13:50~2017/12/28/04:14:00',antenna='DV11')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/04:14:10~2017/12/28/04:14:20',antenna='DV11')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/04:17:10~2017/12/28/04:17:15',antenna='DA44,DV13')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/04:17:30~2017/12/28/04:17:40',antenna='DA44,DA55,DV02')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/04:18:10~2017/12/28/04:18:20',antenna='DA44,DA55')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/04:20:40~2017/12/28/04:20:50',antenna='DA44,DV13,DA51')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/04:21:10~2017/12/28/04:21:20',antenna='DA44,DV13')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/04:21:30~2017/12/28/04:21:40',antenna='DA44,DV13')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/04:21:45~2017/12/28/04:21:55',antenna='DV11,DV13')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/04:22:05~2017/12/28/04:22:15',antenna='DV11,DV13')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/04:22:30~2017/12/28/04:22:40',antenna='DA55,DV13')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/04:22:45~2017/12/28/04:22:55',antenna='DV11')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/04:28:00~2017/12/28/04:28:10',antenna='DV11')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/04:28:40~2017/12/28/04:28:50',antenna='DA44')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/04:30:00~2017/12/28/04:30:10',antenna='DV13')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/04:30:20~2017/12/28/04:30:30',antenna='DV13')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='12',timerange='2017/12/28/04:30:35~2017/12/28/04:30:45',antenna='DV13')

#Inspect gain tables to check if flagging worked
plotms(
    SB_iteration2_p5,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(SB_selfcal_iteration2_folder,f'{prefix}_SB_iteration2_p5_gain_phase_vs_time_flagged.png'),
)

#Apply the calibration gains and split-off the corrected ms
applycal(
    vis=SB_iteration2_cont_p4+'.ms',spw=SB_contspws,spwmap=SB_spw_mapping,
    gaintable=[SB_iteration2_p5],interp='linearPD',calwt=True,applymode='calonly'
)

SB_iteration2_cont_p5 = SB_iteration2_cont_p4.replace('p4','p5')
os.system('rm -rf %s.ms*' %SB_iteration2_cont_p5)
split(vis=SB_iteration2_cont_p4+'.ms',outputvis=SB_iteration2_cont_p5+'.ms',datacolumn='corrected')

#CLEAN, again with a ~6sigma threshold
tclean_wrapper(
    vis       = SB_iteration2_cont_p5+'.ms',
    imagename = SB_iteration2_cont_p5,
    threshold = '0.0768mJy',
    **SB_tclean_wrapper_kwargs
)
estimate_SNR(SB_iteration2_cont_p5+'.image',disk_mask=SB_mask,noise_mask=noise_annulus_SB)
generate_image_png(
    SB_iteration2_cont_p5+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_iteration2_SB,10*rms_iteration2_SB],
    save_folder=SB_selfcal_iteration2_folder
)
#MWC_758_SB_contp5.image
#Beam 0.193 arcsec x 0.146 arcsec (-28.31 deg)
#Flux inside disk mask: 58.12 mJy
#Peak intensity of source: 5.67 mJy/beam
#rms: 1.28e-02 mJy/beam
#Peak SNR: 442.55

#MWC_758_SB_iteration2_contp5.image
#Beam 0.193 arcsec x 0.146 arcsec (-28.29 deg)
#Flux inside disk mask: 58.45 mJy
#Peak intensity of source: 5.71 mJy/beam
#rms: 1.29e-02 mJy/beam
#Peak SNR: 443.37

#Check again how SB phase-only selfcal improved things at each step
non_self_caled_SB_iteration2_vis = SB_iteration2_cont_p0
self_caled_SB_iteration2_visibilities = {
    'p1':SB_iteration2_cont_p1,
    'p2':SB_iteration2_cont_p2,
    'p3':SB_iteration2_cont_p3,
    'p4':SB_iteration2_cont_p4,
    'p5':SB_iteration2_cont_p5
}

#SB_EBs = ('EB0','EB1','EB2','EB3')
#SB_EB_spws = ('0,1,2,3','4,5,6,7','8,9,10,11','12,13,14,15') #fill out by referring to listobs output
#
#for self_cal_step,self_caled_vis in self_caled_SB_iteration2_visibilities.items():
#    for EB_key,spw in zip(SB_EBs,SB_EB_spws):
#
#        print('selfcal_step: {}'.format(self_cal_step)+', EB: {}'.format(EB_key))
#
#        nametemplate = f'{prefix}_SB_iteration2_{EB_key}_{self_cal_step}_compare_amp_vs_time'
#        visibilities = [self_caled_vis+'.ms',non_self_caled_SB_iteration2_vis+'.ms']
#
#        plot_amp_vs_time_comparison(
#            nametemplate=nametemplate,visibilities=visibilities,spw=spw,
#            uvrange=uv_ranges['SB'],output_folder=SB_selfcal_iteration2_folder
#        )

all_SB_iteration2_visibilities = self_caled_SB_iteration2_visibilities.copy()
all_SB_iteration2_visibilities['p0'] = SB_iteration2_cont_p0

SB_flux_ref_EB = 3 #this is SB_EB3

for self_cal_step,vis_name in all_SB_iteration2_visibilities.items():
    #Split out SB EBs
    vis_ms = vis_name+'.ms'
    nametemplate = vis_ms.replace('.ms','_EB')
    split_all_obs(msfile=vis_ms,nametemplate=nametemplate)

    exported_ms = []
    for i in range(number_of_EBs['SB']):
        EB_vis = f'{nametemplate}{i}.ms'
        listobs(vis=EB_vis,listfile=EB_vis+'.listobs.txt',overwrite=True)
        #Export MS contents into numpy save files
        export_MS(EB_vis)
        exported_ms.append(EB_vis.replace('.ms','.vis.npz'))

    for i,exp_ms in enumerate(exported_ms):
        # png_filename = f'iteration2_flux_comparison_SB_EB{i}_{self_cal_step}_to_{flux_ref_EB}.png'
        png_filename = f'iteration2_flux_comparison_SB_EB{i}_{self_cal_step}_to_EB{SB_flux_ref_EB}.png'
        plot_label = os.path.join(SB_selfcal_iteration2_folder,png_filename)
        estimate_flux_scale(
            #reference=f'{prefix}_{flux_ref_EB}_initcont_shift.vis.npz',
            reference=f'{nametemplate}{SB_flux_ref_EB}.vis.npz',
            comparison=exp_ms,incl=incl,PA=PA,plot_label=plot_label
        )

    fluxscale = [1.,]*number_of_EBs['SB']
    plot_label = os.path.join(SB_selfcal_iteration2_folder,f'deprojected_vis_profiles_SB_iteration2_{self_cal_step}.png')
    plot_deprojected(filelist=exported_ms,fluxscale=fluxscale,PA=PA,incl=incl,show_err=True,plot_label=plot_label)

#iteration_1                
#ratio          = [0.88988,0.89305,0.89593,0.90156,0.90895,0.92577] #MWC_758_SB_contp0...p5_EB0.vis.npz vs MWC_758_SB_EB3_initcont.vis.npz
#scaling_factor = [0.943  ,0.945  ,0.947  ,0.950  ,0.953  ,0.962  ]
#ratio          = [0.98340,0.97806,0.98762,0.99364,0.99671,0.99951] #MWC_758_SB_contp0...p5_EB1.vis.npz vs MWC_758_SB_EB3_initcont.vis.npz
#scaling_factor = [0.992  ,0.989  ,0.994  ,0.997  ,0.998  ,1.000  ]
#ratio          = [0.93677,0.96069,0.95945,0.96593,0.96917,0.97148] #MWC_758_SB_contp0...p5_EB2.vis.npz vs MWC_758_SB_EB3_initcont.vis.npz
#scaling_factor = [0.968  ,0.980  ,0.980  ,0.983  ,0.984  ,0.986  ]

#iteration_2
#ratio          = [0.89704,0.90017,0.90316,0.90890,0.91631,0.93330] #MWC_758_SB_contp0...p5_EB0.vis.npz vs MWC_758_SB_EB3_initcont.vis.npz
#scaling_factor = [0.947  ,0.949  ,0.950  ,0.953  ,0.957  ,0.966  ]
#ratio          = [0.99531,0.98994,0.99961,1.00568,1.00880,1.01165] #MWC_758_SB_contp0...p5_EB1.vis.npz vs MWC_758_SB_EB3_initcont.vis.npz
#scaling_factor = [0.998  ,0.995  ,1.000  ,1.003  ,1.004  ,1.006  ]
#ratio          = [0.94053,0.96452,0.96327,0.96981,0.97305,0.97537] #MWC_758_SB_contp0...p5_EB2.vis.npz vs MWC_758_SB_EB3_initcont.vis.npz
#scaling_factor = [0.970  ,0.982  ,0.981  ,0.985  ,0.986  ,0.988  ]

#Begin of SB+LB self-cal - iteration 2
#For phase self-cal, clean down to ~6sigma; for amplitude self-cal, clean down to ~1 sigma
LB_selfcal_iteration2_folder = get_figures_folderpath('8.1_selfcal_SBLB_iteration2_figures')
make_figures_folder(LB_selfcal_iteration2_folder)

LB_iteration2_cont_p0 = prefix+'_SBLB_iteration2_contp0'
os.system('rm -rf %s.ms*' %LB_iteration2_cont_p0)
concat(
    vis=[
        SB_iteration2_cont_p5+'.ms',
        f'{prefix}_LB_EB0_initcont_selfcal_rescaled.ms',#f'{prefix}_LB_EB0_initcont_shift_rescaled.ms',
        f'{prefix}_LB_EB1_initcont_selfcal_rescaled.ms',#f'{prefix}_LB_EB1_initcont_shift_rescaled.ms',
        f'{prefix}_LB_EB2_initcont_selfcal_rescaled.ms',#f'{prefix}_LB_EB2_initcont_shift_rescaled.ms',
        f'{prefix}_LB_EB3_initcont_selfcal_rescaled.ms',#f'{prefix}_LB_EB3_initcont_shift_rescaled.ms',
        f'{prefix}_LB_EB4_initcont_selfcal_rescaled.ms',#f'{prefix}_LB_EB4_initcont_shift_rescaled.ms',
    ],
    concatvis=LB_iteration2_cont_p0+'.ms',dirtol='0.1arcsec',copypointing=False
)
listobs(vis=LB_iteration2_cont_p0+'.ms',listfile=LB_iteration2_cont_p0+'.ms.listobs.txt',overwrite=True)
#2025-01-29 14:41:16     SEVERE  getcell::TIME   Exception Reported: TableProxy::getCell: no such row
#2025-01-29 14:41:20     WARN    concat::::casa  Some but not all of the input MSs are lacking a populated POINTING table:
#2025-01-29 14:41:20     WARN    concat::::casa     0: MWC_758_SB_iteration2_contp5.ms
#2025-01-29 14:41:20     WARN    concat::::casa  The joint dataset will not have a valid POINTING table.
#2025-01-29 14:41:20     WARN    concat::::casa  The setup of the input MSs is not fully consistent. The concatenation may fail
#2025-01-29 14:41:20     WARN    concat::::casa  and/or the affected columns may contain partially only default data.
#2025-01-29 14:41:20     WARN    concat::::casa  {
#    'MWC_758_LB_EB0_initcont_selfcal_rescaled.ms': {'Main': {'present_a': True, 'present_b': True, 'missingcol_a': [], 'missingcol_b': ['MODEL_DATA']}},
#    'MWC_758_LB_EB1_initcont_selfcal_rescaled.ms': {'Main': {'present_a': True, 'present_b': True, 'missingcol_a': [], 'missingcol_b': ['MODEL_DATA']}},
#    'MWC_758_LB_EB2_initcont_selfcal_rescaled.ms': {'Main': {'present_a': True, 'present_b': True, 'missingcol_a': [], 'missingcol_b': ['MODEL_DATA']}},
#    'MWC_758_LB_EB3_initcont_selfcal_rescaled.ms': {'Main': {'present_a': True, 'present_b': True, 'missingcol_a': [], 'missingcol_b': ['MODEL_DATA']}},
#    'MWC_758_LB_EB4_initcont_selfcal_rescaled.ms': {'Main': {'present_a': True, 'present_b': True, 'missingcol_a': [], 'missingcol_b': ['MODEL_DATA']}}
#}

#Define new SB mask using the same centre as before (checked and agrees with the listobs one)
mask_pa        = PA
mask_semimajor = 0.8
mask_semiminor = 0.8
mask_ra  = '05h30m27.535947s' 
mask_dec = '25d19m56.615500s'

LB_mask = f'ellipse[[{mask_ra},{mask_dec}], [{mask_semimajor:.3f}arcsec, {mask_semiminor:.3f}arcsec], {mask_pa:.1f}deg]'

noise_annulus_LB = f"annulus[[{mask_ra}, {mask_dec}],['4.arcsec', '6.arcsec']]"

LB_tclean_wrapper_kwargs = {
    'deconvolver':'multiscale','scales':[0,2,4,8,16,24,32],
    'smallscalebias':0.6,'gain':0.3,'cycleniter':300,
    'niter':1000000,'mask':LB_mask,
    'cellsize':'0.005arcsec','imsize':8000,
    'parallel':use_parallel,'savemodel':'modelcolumn',
    'robust':1.0,'interactive':False, #diff from exoALMA, chose robust=1.0 to increase a bit the peak SNR and aid self-cal
    'gridder':'standard',
}

tclean_wrapper(
    vis       = LB_iteration2_cont_p0+'.ms', 
    imagename = LB_iteration2_cont_p0,
    threshold = '0.0408mJy',
    **LB_tclean_wrapper_kwargs
)
estimate_SNR(LB_iteration2_cont_p0+'.image',disk_mask=LB_mask,noise_mask=noise_annulus_LB)
rms_iteration2_LB = imstat(imagename=LB_iteration2_cont_p0+'.image',region=noise_annulus_LB)['rms'][0]
generate_image_png(
    LB_iteration2_cont_p0+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_iteration2_LB,10*rms_iteration2_LB],
    save_folder=LB_selfcal_iteration2_folder
)
#MWC_758_SBLB_contp0.image
#Beam 0.064 arcsec x 0.045 arcsec (7.62 deg)
#Flux inside disk mask: 67.18 mJy
#Peak intensity of source: 1.14 mJy/beam
#rms: 6.80e-03 mJy/beam
#Peak SNR: 168.11

#MWC_758_SBLB_iteration2_contp0.image
#Beam 0.063 arcsec x 0.044 arcsec (7.28 deg)
#Flux inside disk mask: 67.77 mJy
#Peak intensity of source: 1.12 mJy/beam
#rms: 6.68e-03 mJy/beam
#Peak SNR: 167.33

#First round of phase-only self-cal
#NOTE: you need .p1 instead of _p1 in the caltable name if you want flagdata to work (i.e., flagging problematic antennas non interactively)...
LB_iteration2_p1 = prefix+'_SBLB_iteration2.p1'
os.system('rm -rf '+LB_iteration2_p1)
#If you get many flagged solutions, change gaintype to 'T'
gaincal(
    vis=LB_iteration2_cont_p0+'.ms',caltable=LB_iteration2_p1,
    gaintype='G',combine='scan,spw',calmode='p',solint='inf',
    spw=LB_contspws,refant=LB_refant,
    minsnr=3.,minblperant=4
)
# 9 of 92 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:36:59.3
#19 of 88 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:03:31.9
#20 of 86 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:04:22.6
#25 of 92 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:58:16.0
#26 of 88 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:14:38.9

#Inspect gain tables interactively and decide whether to manually flag something
#plotms(LB_iteration2_p1,xaxis='time',yaxis='GainPhase',iteraxis='spw',showgui=True)
#Print calibration png file
#plotms(
#    LB_iteration2_p1,xaxis='time', yaxis='GainPhase',iteraxis='spw',exprange='all',
#    overwrite=True,showgui=False,
#    plotfile=os.path.join(LB_selfcal_iteration2_folder,f'{prefix}_LB_iteration2_p1_gain_phase_vs_time.png')
#)

#Inspect gain tables and decide whether to flag something
plotms(
    LB_iteration2_p1,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(LB_selfcal_iteration2_folder,f'{prefix}_LB_iteration2_p1_gain_phase_vs_time.png'),
)

#Flag problematic antennas
flagdata(vis=LB_iteration2_p1,mode='manual',spw='0', antenna='DV24')
flagdata(vis=LB_iteration2_p1,mode='manual',spw='4', antenna='DA45,DA54,DV24')
flagdata(vis=LB_iteration2_p1,mode='manual',spw='8', antenna='DA45,DA47')
flagdata(vis=LB_iteration2_p1,mode='manual',spw='20',antenna='DV15')
flagdata(vis=LB_iteration2_p1,mode='manual',spw='32',antenna='DA44,DV11')

#Inspect gain tables to check if flagging worked
plotms(
    LB_iteration2_p1,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(LB_selfcal_iteration2_folder,f'{prefix}_LB_iteration2_p1_gain_phase_vs_time_flagged.png'),
)

#Apply the calibration gains and split-off the corrected ms
applycal(
    vis=LB_iteration2_cont_p0+'.ms',spw=LB_contspws,spwmap=LB_spw_mapping,
    gaintable=[LB_iteration2_p1],interp='linearPD',calwt=True,applymode='calonly'
)
LB_iteration2_cont_p1 = LB_iteration2_cont_p0.replace('p0','p1')
os.system('rm -rf %s.ms*' %LB_iteration2_cont_p1)
split(vis=LB_iteration2_cont_p0+'.ms',outputvis=LB_iteration2_cont_p1+'.ms',datacolumn='corrected')

#CLEAN, again with a ~6sigma threshold
tclean_wrapper(
    vis       = LB_iteration2_cont_p1+'.ms',
    imagename = LB_iteration2_cont_p1,
    threshold = '0.0399mJy',
    **LB_tclean_wrapper_kwargs
)
estimate_SNR(LB_iteration2_cont_p1+'.image',disk_mask=LB_mask,noise_mask=noise_annulus_LB)
generate_image_png(
    LB_iteration2_cont_p1+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_iteration2_LB,10*rms_iteration2_LB],
    save_folder=LB_selfcal_iteration2_folder
)
#MWC_758_SBLB_contp1.image
#Beam 0.065 arcsec x 0.045 arcsec (7.20 deg)
#Flux inside disk mask: 66.27 mJy
#Peak intensity of source: 1.16 mJy/beam
#rms: 6.66e-03 mJy/beam
#Peak SNR: 174.00

#MWC_758_SBLB_iteration2_contp1.image
#Beam 0.063 arcsec x 0.044 arcsec (7.28 deg)
#Flux inside disk mask: 67.19 mJy
#Peak intensity of source: 1.08 mJy/beam
#rms: 6.53e-03 mJy/beam
#Peak SNR: 165.95

#Re-run the first step of phase-only selfcal with infinite integration interval to improve the model and get a better EB alignment
#One step is enough, there is not much improvement later on.
#
#Step .p1_bis
LB_iteration2_p1_bis = LB_iteration2_p1.replace('p1','p1_bis')
os.system('rm -rf '+LB_iteration2_p1_bis)
gaincal(
    vis=LB_iteration2_cont_p1+'.ms',caltable=LB_iteration2_p1_bis,
    gaintype='T',combine='scan,spw',calmode='p',solint='inf',
    spw=LB_contspws,refant=LB_refant,
    minsnr=3.,minblperant=4
)
# 5 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:36:58.3
# 4 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:03:34.4
#10 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:04:20.9
# 9 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:58:14.3
#12 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:14:40.3

#Inspect gain tables interactively and decide whether to manually flag something
#plotms(LB_iteration2_p1_bis,xaxis='time',yaxis='GainPhase',iteraxis='spw',showgui=True)
#Print calibration png file
#plotms(
#    LB_iteration2_p1_bis,xaxis='time',yaxis='GainPhase',iteraxis='spw',exprange='all',
#    overwrite=True,showgui=False,
#    plotfile=os.path.join(LB_selfcal_iteration2_folder,f'{prefix}_LB_iteration2_p1_bis_gain_phase_vs_time.png')
#)

#Inspect gain tables and decide whether to flag something
plotms(
    LB_iteration2_p1_bis,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(LB_selfcal_iteration2_folder,f'{prefix}_LB_iteration2_p1_bis_gain_phase_vs_time.png'),
)

#Flag problematic antennas
flagdata(vis=LB_iteration2_p1_bis,mode='manual',spw='4', antenna='DA54,DV24')
flagdata(vis=LB_iteration2_p1_bis,mode='manual',spw='12',antenna='DA63,DV17')
flagdata(vis=LB_iteration2_p1_bis,mode='manual',spw='16',antenna='DV04,DV16')
flagdata(vis=LB_iteration2_p1_bis,mode='manual',spw='20',antenna='DV15')
flagdata(vis=LB_iteration2_p1_bis,mode='manual',spw='32',antenna='DA44,DV11')

#Inspect gain tables to check if flagging worked
plotms(
    LB_iteration2_p1_bis,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(LB_selfcal_iteration2_folder,f'{prefix}_LB_iteration2_p1_bis_gain_phase_vs_time_flagged.png'),
)

#Apply the calibration gains and split-off the corrected ms
applycal(
    vis=LB_iteration2_cont_p1+'.ms',spw=LB_contspws,spwmap=LB_spw_mapping,
    gaintable=[LB_iteration2_p1_bis],interp='linearPD',calwt=True,applymode='calonly'
)
LB_iteration2_cont_p1_bis = LB_iteration2_cont_p1.replace('p1','p1_bis')
os.system('rm -rf %s.ms*' %LB_iteration2_cont_p1_bis)
split(vis=LB_iteration2_cont_p1+'.ms',outputvis=LB_iteration2_cont_p1_bis+'.ms',datacolumn='corrected')

#CLEAN, again with a ~6sigma threshold
tclean_wrapper(
    vis       = LB_iteration2_cont_p1_bis+'.ms',
    imagename = LB_iteration2_cont_p1_bis,
    threshold = '0.0392mJy',
    **LB_tclean_wrapper_kwargs
)
estimate_SNR(LB_iteration2_cont_p1_bis+'.image',disk_mask=LB_mask,noise_mask=noise_annulus_LB)
generate_image_png(
    LB_iteration2_cont_p1_bis+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_iteration2_LB,10*rms_iteration2_LB],
    save_folder=LB_selfcal_iteration2_folder
)
#MWC_758_SBLB_contp1_bis.image
#Beam 0.065 arcsec x 0.045 arcsec (7.20 deg)
#Flux inside disk mask: 66.43 mJy
#Peak intensity of source: 1.15 mJy/beam
#rms: 6.67e-03 mJy/beam
#Peak SNR: 173.01

#MWC_758_SBLB_iteration2_contp1_bis.image
#Beam 0.063 arcsec x 0.044 arcsec (7.28 deg)
#Flux inside disk mask: 66.89 mJy
#Peak intensity of source: 1.11 mJy/beam
#rms: 6.52e-03 mJy/beam
#Peak SNR: 169.84

#Second round of phase-only self-cal
LB_iteration2_p2 = LB_iteration2_p1_bis.replace('p1_bis','p2')
os.system('rm -rf '+LB_iteration2_p2)
gaincal(
    vis=LB_iteration2_cont_p1_bis+'.ms',caltable=LB_iteration2_p2,
    gaintype='T',combine='scan,spw',calmode='p',solint='360s',
    spw=LB_contspws,refant=LB_refant,
    minsnr=3.,minblperant=4
)
#19 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:03:16.6
#19 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:09:25.3
#22 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:12:54.4
#20 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:16:49.9
#23 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:19:47.3
#22 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:24:11.2
#25 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:27:23.0
#20 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:30:56.4
#17 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:37:16.9
#21 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:40:57.4
#20 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:46:41.8
#24 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:48:50.3
#17 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:52:18.5
#16 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:58:26.4
#27 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:01:43.6
#19 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:06:38.2
#17 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:12:14.4
#21 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:27:52.6
#22 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:28:10.7
#18 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:31:42.3
#22 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:34:53.0
#13 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:38:23.5
#28 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:41:31.5
#15 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:45:10.8
#14 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:51:49.3
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:54:59.1
#16 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:58:26.5
#18 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:04:48.8
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:08:29.0
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:14:20.6
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:16:19.4
#17 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:19:49.8
#16 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:25:55.0
#21 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:29:13.5
#16 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:34:09.4
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:39:23.1
#24 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:29:55.8
#20 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:33:37.6
#25 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:36:52.4
#19 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:40:34.9
#32 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:43:35.8
#18 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:47:18.7
#19 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:51:16.6
#18 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:55:02.4
#30 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:58:19.8
#16 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:02:31.1
#18 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:08:07.6
#20 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:14:49.6
#25 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:18:17.6
#17 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:21:51.6
#18 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:27:52.7
#24 of 38 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:31:20.7
#21 of 37 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:37:30.9
#19 of 35 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:42:17.6
#20 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:27:27.5
#18 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:31:13.1
#20 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:34:23.5
#16 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:38:11.2
#29 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:41:08.0
#16 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:44:51.8
#20 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:48:50.1
#17 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:52:35.3
#23 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:55:53.3
#19 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:59:57.3
#21 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:05:35.3
#20 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:12:10.9
#24 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:15:56.6
#17 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:19:32.4
#18 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:25:45.9
#21 of 44 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:29:09.1
#16 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:34:11.2
#18 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:39:48.6
#20 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:43:03.0
#13 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:46:45.6
#25 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:49:57.7
#15 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:53:38.8
#15 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:00:06.3
#13 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:06:38.3
#28 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:09:52.4
#16 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:13:29.5
#28 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:16:31.1
#15 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:19:57.7
#23 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:23:26.8
#17 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:27:41.8
#24 of 41 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:31:15.8
#16 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:34:37.3
#16 of 43 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:42:14.6
#34 of 40 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:44:18.8
#13 of 40 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:49:09.6
#17 of 39 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:55:27.9
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:40:07.6
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:53:08.3
# 1 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:57:48.9
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:02:34.0
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:14:50.5
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:23:12.0
# 1 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:28:15.7
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:11:47.4
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:10:33.4
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:31:22.8
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:35:39.6
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:39:01.8
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/03:57:48.9

#Inspect gain tables interactively and decide whether to manually flag something
#plotms(LB_iteration2_p2,xaxis='time',yaxis='GainPhase',iteraxis='spw',showgui=True)
#Print calibration png file
#plotms(
#    LB_iteration2_p2,xaxis='time', yaxis='GainPhase',iteraxis='spw',exprange='all',
#    overwrite=True,showgui=False,
#    plotfile=os.path.join(LB_selfcal_iteration2_folder,f'{prefix}_LB_iteration2_p2_gain_phase_vs_time.png')
#)

#Inspect gain tables and decide whether to flag something
plotms(
    LB_iteration2_p2,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(LB_selfcal_iteration2_folder,f'{prefix}_LB_iteration2_p2_gain_phase_vs_time.png'),
)

#Flag problematic antennas
flagdata(vis=LB_iteration2_p2,mode='manual',spw='0', timerange='2017/10/10/09:09:20~2017/10/10/09:09:30',antenna='DA47')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='0', timerange='2017/10/10/09:12:50~2017/10/10/09:13:00',antenna='DA45,DA60')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='0', timerange='2017/10/10/09:30:50~2017/10/10/09:31:00',antenna='DV20,DV22')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='0', timerange='2017/10/10/09:48:50~2017/10/10/09:49:00',antenna='DA47,DA64')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='0', timerange='2017/10/10/09:58:20~2017/10/10/09:58:30',antenna='DA45,DA58')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='0', timerange='2017/10/10/10:01:40~2017/10/10/10:01:50',antenna='DA53')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='0', timerange='2017/10/10/10:06:30~2017/10/10/10:06:40',antenna='DA64')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='0', timerange='2017/10/10/10:12:10~2017/10/10/10:12:20',antenna='DV17')

flagdata(vis=LB_iteration2_p2,mode='manual',spw='4', timerange='2017/10/10/10:27:50~2017/10/10/10:28:00',antenna='DV22')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='4', timerange='2017/10/10/10:38:20~2017/10/10/10:38:30',antenna='DA54,DV24')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='4', timerange='2017/10/10/10:51:30~2017/10/10/10:52:00',antenna='DA44,DA53,DA54,DV15,DV24,PM04')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='4', timerange='2017/10/10/10:54:40~2017/10/10/10:55:00',antenna='DA54')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='4', timerange='2017/10/10/11:16:10~2017/10/10/11:16:20',antenna='DA59')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='4', timerange='2017/10/10/11:34:00~2017/10/10/11:34:10',antenna='DA45')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='4', timerange='2017/10/10/11:39:20~2017/10/10/11:39:30',antenna='DV17')

flagdata(vis=LB_iteration2_p2,mode='manual',spw='8', timerange='2017/10/11/07:36:50~2017/10/11/07:37:00',antenna='DA44')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='8', timerange='2017/10/11/07:43:30~2017/10/11/07:43:40',antenna='DV04')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='8', timerange='2017/10/11/08:02:30~2017/10/11/08:02:40',antenna='DA59')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='8', timerange='2017/10/11/08:27:50~2017/10/11/08:28:00',antenna='DA52')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='8', timerange='2017/10/11/08:31:20~2017/10/11/08:31:30',antenna='DA63')

flagdata(vis=LB_iteration2_p2,mode='manual',spw='12',timerange='2017/10/15/07:34:20~2017/10/15/07:34:30',antenna='DA58,DV12,DV17,DV13')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='12',timerange='2017/10/15/07:41:00~2017/10/15/07:41:10',antenna='DA63')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='12',timerange='2017/10/15/07:55:50~2017/10/15/07:56:00',antenna='DA52,DA58')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='12',timerange='2017/10/15/08:29:00~2017/10/15/08:29:10',antenna='DV03')

flagdata(vis=LB_iteration2_p2,mode='manual',spw='16',timerange='2017/10/16/09:00:00~2017/10/16/09:00:30',antenna='DV15')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='16',timerange='2017/10/16/09:06:30~2017/10/16/09:06:40',antenna='DA56')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='16',timerange='2017/10/16/09:09:50~2017/10/16/09:10:00',antenna='DV15,PM01')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='16',timerange='2017/10/16/09:16:25~2017/10/16/09:16:35',antenna='DV12,DV15,DV25')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='16',timerange='2017/10/16/09:31:00~2017/10/16/09:31:20',antenna='DA51,DA55,DV23')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='16',timerange='2017/10/16/09:34:30~2017/10/16/09:34:40',antenna='DA46,DA52')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='16',timerange='2017/10/16/09:42:10~2017/10/16/09:42:20',antenna='DV10')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='16',timerange='2017/10/16/09:44:10~2017/10/16/09:44:20',antenna='DA50')

flagdata(vis=LB_iteration2_p2,mode='manual',spw='20',timerange='2017/12/09/04:49:30~2017/12/09/04:50:00',antenna='DV15')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='20',timerange='2017/12/09/04:57:30~2017/12/09/04:58:00',antenna='DV11')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='20',timerange='2017/12/09/05:11:30~2017/12/09/05:12:00',antenna='DA46,DA55,DV10,DV15')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='20',timerange='2017/12/09/05:19:40~2017/12/09/05:20:00',antenna='DV10')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='20',timerange='2017/12/09/05:23:10~2017/12/09/05:23:20',antenna='DA46')

flagdata(vis=LB_iteration2_p2,mode='manual',spw='24',timerange='2017/12/17/06:00:30~2017/12/17/06:01:00',antenna='DV07,DV13')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='24',timerange='2017/12/17/06:04:00~2017/12/17/06:04:20',antenna='DV13')

flagdata(vis=LB_iteration2_p2,mode='manual',spw='28',timerange='2017/12/27/04:10:30~2017/12/27/04:10:40',antenna='DA42,DA44,DV13')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='28',timerange='2017/12/27/04:35:30~2017/12/27/04:35:40',antenna='DA44')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='28',timerange='2017/12/27/04:39:00~2017/12/27/04:39:10',antenna='DA44')

flagdata(vis=LB_iteration2_p2,mode='manual',spw='32',timerange='2017/12/28/03:44:50~2017/12/28/03:45:00',antenna='DA44')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='32',timerange='2017/12/28/03:54:30~2017/12/28/03:54:40',antenna='DA44,DA55')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='32',timerange='2017/12/28/03:57:30~2017/12/28/03:58:00',antenna='DA44,DA55,DV11')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='32',timerange='2017/12/28/04:02:00~2017/12/28/04:02:20',antenna='DA44,DA46,DA55,DV11,DV13,DV07')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='32',timerange='2017/12/28/04:07:00~2017/12/28/04:07:20',antenna='DA44,DV11')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='32',timerange='2017/12/28/04:15:00~2017/12/28/04:15:40',antenna='DA44,DV11,DV13')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='32',timerange='2017/12/28/04:22:50~2017/12/28/04:23:00',antenna='DV11,DV13')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='32',timerange='2017/12/28/04:26:10~2017/12/28/04:26:20',antenna='DA44,DA55,DV11')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='32',timerange='2017/12/28/04:28:40~2017/12/28/04:29:00',antenna='DA44,DA55,DV11,DV13')

#Inspect gain tables to check if flagging worked
plotms(
    LB_iteration2_p2,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(LB_selfcal_iteration2_folder,f'{prefix}_LB_iteration2_p2_gain_phase_vs_time_flagged.png'),
)

#Apply the calibration gains and split-off the corrected ms
applycal(
    vis=LB_iteration2_cont_p1_bis+'.ms',spw=LB_contspws,spwmap=LB_spw_mapping,
    gaintable=[LB_iteration2_p2],interp='linearPD',calwt=True,applymode='calonly'
)
LB_iteration2_cont_p2 = LB_iteration2_cont_p1_bis.replace('p1_bis','p2')
os.system('rm -rf %s.ms*' %LB_iteration2_cont_p2)
split(vis=LB_iteration2_cont_p1_bis+'.ms',outputvis=LB_iteration2_cont_p2+'.ms',datacolumn='corrected')

#CLEAN, again with a ~6sigma threshold
tclean_wrapper(
    vis       = LB_iteration2_cont_p2+'.ms',
    imagename = LB_iteration2_cont_p2,
    threshold = '0.0391mJy',
    **LB_tclean_wrapper_kwargs
)
estimate_SNR(LB_iteration2_cont_p2+'.image',disk_mask=LB_mask,noise_mask=noise_annulus_LB)
generate_image_png(
    LB_iteration2_cont_p2+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_iteration2_LB,10*rms_iteration2_LB],
    save_folder=LB_selfcal_iteration2_folder
)
#MWC_758_SBLB_contp2.image
#Beam 0.065 arcsec x 0.045 arcsec (7.20 deg)
#Flux inside disk mask: 65.99 mJy
#Peak intensity of source: 1.16 mJy/beam
#rms: 6.61e-03 mJy/beam
#Peak SNR: 175.02

#MWC_758_SBLB_iteration2_contp2.image
#Beam 0.063 arcsec x 0.044 arcsec (7.28 deg)
#Flux inside disk mask: 66.78 mJy
#Peak intensity of source: 1.08 mJy/beam
#rms: 6.49e-03 mJy/beam
#Peak SNR: 166.96

#Third round of phase-only self-cal
LB_iteration2_p3 = LB_iteration2_p2.replace('p2','p3')
os.system('rm -rf '+LB_iteration2_p3)
gaincal(
    vis=LB_iteration2_cont_p2+'.ms',caltable=LB_iteration2_p3,
    gaintype='T',combine='scan,spw',calmode='p',solint='120s',
    spw=LB_contspws,refant=LB_refant,
    minsnr=3.,minblperant=4
)

#Inspect gain tables interactively and decide whether to manually flag something
#plotms(LB_iteration2_p3,xaxis='time',yaxis='GainPhase',iteraxis='spw',showgui=True)
#Print calibration png file
#plotms(
#    LB_iteration2_p3,xaxis='time', yaxis='GainPhase',iteraxis='spw',exprange='all',
#    overwrite=True,showgui=False,
#    plotfile=os.path.join(LB_selfcal_iteration2_folder,f'{prefix}_LB_iteration2_p3_gain_phase_vs_time.png')
#)

#Inspect gain tables and decide whether to flag something
plotms(
    LB_iteration2_p3,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(LB_selfcal_iteration2_folder,f'{prefix}_LB_iteration2_p3_gain_phase_vs_time.png'),
)

#Flag problematic antennas
flagdata(vis=LB_iteration2_p3,mode='manual',spw='0', timerange='2017/10/10/09:03:30~2017/10/10/09:03:40',antenna='DA52')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='0', timerange='2017/10/10/09:07:20~2017/10/10/09:07:30',antenna='DA65')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='0', timerange='2017/10/10/09:08:50~2017/10/10/09:09:00',antenna='DA47')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='0', timerange='2017/10/10/09:13:00~2017/10/10/09:13:10',antenna='DA47,DA52')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='0', timerange='2017/10/10/09:14:20~2017/10/10/09:14:30',antenna='DV03')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='0', timerange='2017/10/10/09:15:40~2017/10/10/09:15:50',antenna='DA56')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='0', timerange='2017/10/10/09:18:25~2017/10/10/09:18:35',antenna='DA41,DA42,DV17')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='0', timerange='2017/10/10/09:29:40~2017/10/10/09:29:50',antenna='DA47')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='0', timerange='2017/10/10/09:40:30~2017/10/10/09:40:40',antenna='DA59')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='0', timerange='2017/10/10/09:43:30~2017/10/10/09:43:40',antenna='DV05,DV15,DV16')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='0', timerange='2017/10/10/09:47:25~2017/10/10/09:47:35',antenna='DA65,DV15,DV22')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='0', timerange='2017/10/10/09:57:40~2017/10/10/09:58:00',antenna='DA63')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='0', timerange='2017/10/10/10:00:20~2017/10/10/10:00:30',antenna='DV15')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='0', timerange='2017/10/10/10:01:20~2017/10/10/10:01:40',antenna='DV22')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='0', timerange='2017/10/10/10:08:00~2017/10/10/10:08:10',antenna='DA65,DV15,DV22')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='0', timerange='2017/10/10/10:09:30~2017/10/10/10:09:40',antenna='DA54')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='0', timerange='2017/10/10/10:10:40~2017/10/10/10:10:50',antenna='DA64')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='0', timerange='2017/10/10/10:12:10~2017/10/10/10:12:20',antenna='DV17')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='0', timerange='2017/10/10/10:14:25~2017/10/10/10:14:35',antenna='DV04')

flagdata(vis=LB_iteration2_p3,mode='manual',spw='4', timerange='2017/10/10/10:27:50~2017/10/10/10:28:00',antenna='DV22')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='4', timerange='2017/10/10/10:28:10~2017/10/10/10:28:20',antenna='DV05')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='4', timerange='2017/10/10/10:30:50~2017/10/10/10:31:00',antenna='DA60')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='4', timerange='2017/10/10/10:32:20~2017/10/10/10:32:30',antenna='DV24')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='4', timerange='2017/10/10/10:35:00~2017/10/10/10:35:10',antenna='DV16')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='4', timerange='2017/10/10/10:36:40~2017/10/10/10:36:50',antenna='DA60')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='4', timerange='2017/10/10/10:39:10~2017/10/10/10:39:20',antenna='DA52')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='4', timerange='2017/10/10/10:46:05~2017/10/10/10:46:15',antenna='DA56')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='4', timerange='2017/10/10/10:50:50~2017/10/10/10:51:00',antenna='DA44,DA53,PM04')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='4', timerange='2017/10/10/10:56:15~2017/10/10/10:56:25',antenna='DA58')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='4', timerange='2017/10/10/10:57:20~2017/10/10/10:57:30',antenna='DV13')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='4', timerange='2017/10/10/11:12:10~2017/10/10/11:12:20',antenna='DA42,DV01')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='4', timerange='2017/10/10/11:17:40~2017/10/10/11:17:50',antenna='DA41,DA44,DA53,DA60')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='4', timerange='2017/10/10/11:22:30~2017/10/10/11:22:40',antenna='DA48')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='4', timerange='2017/10/10/11:25:10~2017/10/10/11:25:20',antenna='DA47')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='4', timerange='2017/10/10/11:27:50~2017/10/10/11:28:00',antenna='DA47')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='4', timerange='2017/10/10/11:28:50~2017/10/10/11:29:00',antenna='DV24')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='4', timerange='2017/10/10/11:35:30~2017/10/10/11:35:40',antenna='DV05')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='4', timerange='2017/10/10/11:38:10~2017/10/10/11:38:20',antenna='DV22')

flagdata(vis=LB_iteration2_p3,mode='manual',spw='8', timerange='2017/10/11/07:29:35~2017/10/11/07:29:45',antenna='DV08')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='8', timerange='2017/10/11/07:31:30~2017/10/11/07:31:40',antenna='DV15,DA53')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='8', timerange='2017/10/11/07:32:40~2017/10/11/07:32:50',antenna='DV22,DA41,DA42')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='8', timerange='2017/10/11/07:34:20~2017/10/11/07:34:30',antenna='DV13,DA53')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='8', timerange='2017/10/11/07:37:00~2017/10/11/07:37:10',antenna='DA44,DV15,PM04')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='8', timerange='2017/10/11/07:45:50~2017/10/11/07:46:00',antenna='DA63')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='8', timerange='2017/10/11/07:52:30~2017/10/11/07:53:00',antenna='DA49')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='8', timerange='2017/10/11/07:54:50~2017/10/11/07:55:00',antenna='DA59')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='8', timerange='2017/10/11/07:57:40~2017/10/11/07:57:50',antenna='DA47')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='8', timerange='2017/10/11/07:58:20~2017/10/11/07:58:30',antenna='DA54')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='8', timerange='2017/10/11/08:14:10~2017/10/11/08:14:20',antenna='DV03')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='8', timerange='2017/10/11/08:19:40~2017/10/11/08:20:00',antenna='DA53,DA54')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='8', timerange='2017/10/11/08:20:40~2017/10/11/08:20:50',antenna='DA46')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='8', timerange='2017/10/11/08:23:20~2017/10/11/08:23:30',antenna='DA45')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='8', timerange='2017/10/11/08:24:30~2017/10/11/08:24:40',antenna='DA58')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='8', timerange='2017/10/11/08:27:20~2017/10/11/08:27:30',antenna='DA60')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='8', timerange='2017/10/11/08:28:50~2017/10/11/08:29:00',antenna='DA63')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='8', timerange='2017/10/11/08:30:00~2017/10/11/08:30:10',antenna='DV05,DV12')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='8', timerange='2017/10/11/08:36:40~2017/10/11/08:36:50',antenna='DV20,PM01')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='8', timerange='2017/10/11/08:37:40~2017/10/11/08:37:50',antenna='')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='8', timerange='2017/10/11/08:40:30~2017/10/11/08:40:40',antenna='DA59,DA64,DV16,DA65')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='8', timerange='2017/10/11/08:41:40~2017/10/11/08:41:50',antenna='DV04')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='8', timerange='2017/10/11/08:44:20~2017/10/11/08:44:30',antenna='DA59')

flagdata(vis=LB_iteration2_p3,mode='manual',spw='12',timerange='2017/10/15/07:27:20~2017/10/15/07:27:30',antenna='DA60')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='12',timerange='2017/10/15/07:29:35~2017/10/15/07:29:45',antenna='DV03')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='12',timerange='2017/10/15/07:33:00~2017/10/15/07:33:10',antenna='DA52')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='12',timerange='2017/10/15/07:39:00~2017/10/15/07:39:10',antenna='DV03')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='12',timerange='2017/10/15/07:39:40~2017/10/15/07:39:50',antenna='DA54,DA57')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='12',timerange='2017/10/15/07:50:25~2017/10/15/07:50:35',antenna='DA46,DV04')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='12',timerange='2017/10/15/07:53:10~2017/10/15/07:53:20',antenna='DV04')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='12',timerange='2017/10/15/07:55:50~2017/10/15/07:56:00',antenna='DA58,DV04,DV17')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='12',timerange='2017/10/15/08:03:40~2017/10/15/08:03:50',antenna='DA60,DV12,DV24')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='12',timerange='2017/10/15/08:11:40~2017/10/15/08:11:50',antenna='DV03')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='12',timerange='2017/10/15/08:15:20~2017/10/15/08:15:30',antenna='DA58')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='12',timerange='2017/10/15/08:22:20~2017/10/15/08:22:30',antenna='PM03')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='12',timerange='2017/10/15/08:32:00~2017/10/15/08:32:10',antenna='DV03')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='12',timerange='2017/10/15/08:35:35~2017/10/15/08:35:45',antenna='DA44,DA48,DA52,DV23')

flagdata(vis=LB_iteration2_p3,mode='manual',spw='16',timerange='2017/10/16/08:43:00~2017/10/16/08:43:10',antenna='DA60')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='16',timerange='2017/10/16/08:48:30~2017/10/16/08:48:40',antenna='DA42,DA45')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='16',timerange='2017/10/16/08:50:10~2017/10/16/08:50:20',antenna='DA55')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='16',timerange='2017/10/16/08:58:00~2017/10/16/08:58:10',antenna='DV12')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='16',timerange='2017/10/16/09:00:50~2017/10/16/09:01:00',antenna='DV08,DV20')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='16',timerange='2017/10/16/09:21:40~2017/10/16/09:21:50',antenna='DA60,DA64')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='16',timerange='2017/10/16/09:23:00~2017/10/16/09:23:10',antenna='DV10')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='16',timerange='2017/10/16/09:25:00~2017/10/16/10:05:00',antenna='')

flagdata(vis=LB_iteration2_p3,mode='manual',spw='20',timerange='2017/12/09/04:40:00~2017/12/09/04:40:10',antenna='DV15')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='20',timerange='2017/12/09/04:51:40~2017/12/09/04:52:00',antenna='DV10,DV11,DV15')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='20',timerange='2017/12/09/04:55:00~2017/12/09/04:55:20',antenna='DV10,DV15')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='20',timerange='2017/12/09/04:58:30~2017/12/09/04:59:00',antenna='DV11,DV15')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='20',timerange='2017/12/09/05:00:30~2017/12/09/05:01:00',antenna='DV10,DV11,DV15')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='20',timerange='2017/12/09/05:04:00~2017/12/09/05:04:30',antenna='DA49')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='20',timerange='2017/12/09/05:09:30~2017/12/09/05:09:40',antenna='DV11')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='20',timerange='2017/12/09/05:11:30~2017/12/09/05:11:40',antenna='DA46,DV15')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='20',timerange='2017/12/09/05:13:30~2017/12/09/05:13:40',antenna='DA46,DA55,DV10')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='20',timerange='2017/12/09/05:17:40~2017/12/09/05:18:00',antenna='DV10')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='20',timerange='2017/12/09/05:23:00~2017/12/09/05:23:20',antenna='DA46,DV10')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='20',timerange='2017/12/09/05:25:40~2017/12/09/05:26:00',antenna='DA46,DA51,DA55,DV15')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='20',timerange='2017/12/09/05:28:10~2017/12/09/05:28:20',antenna='DV10')

flagdata(vis=LB_iteration2_p3,mode='manual',spw='24',timerange='2017/12/17/05:53:30~2017/12/17/05:54:00',antenna='DA43,DA44')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='24',timerange='2017/12/17/05:54:30~2017/12/17/05:55:00',antenna='DA41,DA42,DA43,DA44,DV11,DV13')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='24',timerange='2017/12/17/05:58:40~2017/12/17/05:59:00',antenna='DA55,DV13,DV07')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='24',timerange='2017/12/17/06:00:40~2017/12/17/06:01:00',antenna='DV13')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='24',timerange='2017/12/17/06:02:40~2017/12/17/06:03:00',antenna='DV13,DV07')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='24',timerange='2017/12/17/06:04:00~2017/12/17/06:04:20',antenna='DV13')

flagdata(vis=LB_iteration2_p3,mode='manual',spw='28',timerange='2017/12/27/04:09:00~2017/12/27/04:09:30',antenna='DA42,DV13')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='28',timerange='2017/12/27/04:10:30~2017/12/27/04:11:00',antenna='DA42,DA44,DV13')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='28',timerange='2017/12/27/04:33:30~2017/12/27/04:34:00',antenna='DA44,DV11')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='28',timerange='2017/12/27/04:35:30~2017/12/27/04:36:00',antenna='DA44,DV11')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='28',timerange='2017/12/27/04:37:30~2017/12/27/04:38:00',antenna='DA44,DV11')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='28',timerange='2017/12/27/04:38:50~2017/12/27/04:39:10',antenna='DA44,DV11')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='28',timerange='2017/12/27/04:40:40~2017/12/27/04:40:50',antenna='DV11')

flagdata(vis=LB_iteration2_p3,mode='manual',spw='32',timerange='2017/12/28/03:56:30~2017/12/28/03:56:40',antenna='DA55')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='32',timerange='2017/12/28/03:57:40~2017/12/28/03:57:50',antenna='DA55')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='32',timerange='2017/12/28/03:59:30~2017/12/28/03:59:40',antenna='DA55')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='32',timerange='2017/12/28/04:04:20~2017/12/28/04:04:30',antenna='DA44,DA46,DA55')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='32',timerange='2017/12/28/04:09:20~2017/12/28/04:09:30',antenna='DA42,DV11,DV13')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='32',timerange='2017/12/28/04:18:30~2017/12/28/04:18:40',antenna='DA55')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='32',timerange='2017/12/28/04:20:50~2017/12/28/04:21:00',antenna='DA44,DV11,DV13')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='32',timerange='2017/12/28/04:30:10~2017/12/28/04:30:20',antenna='DV13')

#Inspect gain tables to check if flagging worked
plotms(
    LB_iteration2_p3,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(LB_selfcal_iteration2_folder,f'{prefix}_LB_iteration2_p3_gain_phase_vs_time_flagged.png'),
)

#Apply the calibration gains and split-off the corrected ms
applycal(
    vis=LB_iteration2_cont_p2+'.ms',spw=LB_contspws,spwmap=LB_spw_mapping,
    gaintable=[LB_iteration2_p3],interp='linearPD',calwt=True,applymode='calonly'
)
LB_iteration2_cont_p3 = LB_iteration2_cont_p2.replace('p2','p3')
os.system('rm -rf %s.ms*' %LB_iteration2_cont_p3)
split(vis=LB_iteration2_cont_p2+'.ms',outputvis=LB_iteration2_cont_p3+'.ms',datacolumn='corrected')

#CLEAN, again with a ~6sigma threshold
tclean_wrapper(
    vis       = LB_iteration2_cont_p3+'.ms',
    imagename = LB_iteration2_cont_p3,
    threshold = '0.0389mJy',
    **LB_tclean_wrapper_kwargs
)
estimate_SNR(LB_iteration2_cont_p3+'.image',disk_mask=LB_mask,noise_mask=noise_annulus_LB)
generate_image_png(
    LB_iteration2_cont_p3+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_iteration2_LB,10*rms_iteration2_LB],
    save_folder=LB_selfcal_iteration2_folder
)
#MWC_758_SBLB_contp3.image
#Beam 0.065 arcsec x 0.045 arcsec (7.20 deg)
#Flux inside disk mask: 65.99 mJy
#Peak intensity of source: 1.20 mJy/beam
#rms: 6.56e-03 mJy/beam
#Peak SNR: 182.50

#MWC_758_SBLB_iteration2_contp3.image
#Beam 0.063 arcsec x 0.044 arcsec (7.28 deg)
#Flux inside disk mask: 66.89 mJy
#Peak intensity of source: 1.14 mJy/beam
#rms: 6.44e-03 mJy/beam
#Peak SNR: 176.43

#Fourth round of phase-only self-cal
LB_iteration2_p4 = LB_iteration2_p3.replace('p3','p4')
os.system('rm -rf '+LB_iteration2_p4)
gaincal(
    vis=LB_iteration2_cont_p3+'.ms',caltable=LB_iteration2_p4,
    gaintype='T',combine='scan,spw',calmode='p',solint='60s',
    spw=LB_contspws,refant=LB_refant,
    minsnr=3.,minblperant=4
)
#22 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:00:28.1
#21 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:01:40.1
#22 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:03:01.7
#20 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:04:22.7
#23 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:05:44.0
#21 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:07:05.1
#21 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:08:26.4
#21 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:09:26.3
#20 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:11:21.2
#19 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:12:42.4
#22 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:13:54.2
#24 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:15:36.9
#23 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:16:48.9
#21 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:18:10.2
#21 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:19:31.5
#19 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:21:40.8
#23 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:23:02.0
#23 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:24:22.9
#19 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:25:43.9
#19 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:27:04.9
#20 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:28:25.9
#24 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:29:44.6
#24 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:32:14.5
#21 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:33:26.9
#20 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:34:47.3
#21 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:36:08.3
#19 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:37:29.3
#18 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:38:50.5
#22 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:40:11.0
#22 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:41:11.4
#23 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:43:05.4
#20 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:44:26.4
#21 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:45:38.4
#21 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:47:21.8
#19 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:48:33.1
#21 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:49:54.0
#22 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:51:14.8
#21 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:53:23.0
#20 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:54:44.5
#22 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:56:05.2
#20 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:57:25.8
#22 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/09:58:46.4
#21 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:00:07.2
#22 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:01:25.1
#21 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:03:54.3
#22 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:05:06.1
#19 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:06:26.7
#19 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:07:47.2
#19 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:09:07.6
#18 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:10:28.0
#18 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:11:48.6
#19 of 45 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:12:48.7
#20 of 46 solutions flagged due to SNR < 3 in spw=0  at 2017/10/10/10:14:29.3
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:27:57.6
#29 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:28:15.7
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:29:12.4
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:30:32.3
#18 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:31:52.4
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:33:13.3
#17 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:34:33.2
#18 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:35:53.5
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:36:52.9
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:38:45.8
#18 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:40:07.0
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:41:17.9
#17 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:43:09.6
#18 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:44:29.8
#17 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:45:50.0
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:47:10.5
#18 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:49:17.9
#17 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:50:38.4
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:51:58.7
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:53:19.2
#18 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:54:39.3
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:55:59.5
#21 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:57:23.4
#17 of 43 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/10:59:47.0
#16 of 43 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:00:58.5
#18 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:02:18.8
#17 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:03:38.8
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:04:59.2
#18 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:06:19.4
#17 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:07:39.7
#18 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:08:40.4
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:10:33.6
#26 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:11:54.5
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:13:05.3
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:14:52.4
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:16:03.6
#18 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:17:23.9
#20 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:18:44.3
#17 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:20:54.3
#18 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:22:15.1
#18 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:23:35.1
#15 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:24:55.4
#15 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:26:15.7
#17 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:27:35.9
#18 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:28:54.0
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:31:25.8
#18 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:32:37.7
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:33:57.6
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:35:17.8
#17 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:36:38.4
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:37:58.3
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:39:18.5
#18 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:40:17.8
#19 of 44 solutions flagged due to SNR < 3 in spw=4  at 2017/10/10/11:41:49.3
#30 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:29:39.7
#22 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:29:57.8
#19 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:31:08.4
#18 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:32:30.1
#20 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:33:52.1
#25 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:35:13.6
#19 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:36:35.2
#20 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:37:56.9
#24 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:38:47.4
#20 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:40:31.4
#22 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:41:53.0
#24 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:43:14.6
#21 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:45:08.1
#18 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:46:20.7
#19 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:47:42.4
#20 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:49:03.9
#20 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:51:14.6
#21 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:52:36.2
#18 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:53:57.9
#19 of 42 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:55:19.4
#19 of 42 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:56:41.0
#20 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:58:02.6
#20 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/07:59:22.2
#21 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:01:55.9
#22 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:03:08.7
#20 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:04:30.0
#20 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:05:51.7
#20 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:07:13.3
#20 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:08:34.9
#18 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:09:56.7
#23 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:10:47.5
#21 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:12:32.1
#19 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:13:53.7
#19 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:15:06.1
#19 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:16:50.1
#22 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:18:02.6
#20 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:19:24.0
#22 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:20:45.5
#21 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:22:56.8
#21 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:24:18.3
#19 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:25:39.9
#19 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:27:01.3
#20 of 43 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:28:22.8
#19 of 41 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:29:44.4
#17 of 38 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:31:03.5
#23 of 34 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:33:34.3
#23 of 35 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:34:46.6
#21 of 35 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:36:08.0
#19 of 36 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:37:29.4
#20 of 36 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:38:50.8
#21 of 35 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:40:12.3
#21 of 35 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:41:33.8
#24 of 35 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:42:25.2
#18 of 35 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:44:07.7
#23 of 34 solutions flagged due to SNR < 3 in spw=8  at 2017/10/11/08:45:10.9
#21 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:27:27.5
#20 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:28:39.9
#20 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:30:01.6
#18 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:31:23.2
#21 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:32:44.7
#19 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:34:06.4
#19 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:35:28.0
#20 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:36:18.8
#17 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:38:04.6
#21 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:39:26.0
#24 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:40:47.8
#23 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:42:42.9
#21 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:43:55.5
#20 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:45:17.2
#21 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:46:38.8
#19 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:48:48.2
#16 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:50:09.9
#17 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:51:31.4
#17 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:52:52.9
#19 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:54:14.5
#22 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:55:36.1
#26 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:56:55.8
#23 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/07:59:27.3
#18 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:00:39.6
#24 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:02:01.1
#23 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:03:22.5
#21 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:04:44.1
#22 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:06:05.8
#22 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:07:27.2
#33 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:08:18.0
#22 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:10:02.2
#20 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:11:23.7
#24 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:12:45.1
#27 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:14:38.2
#24 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:15:50.5
#18 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:17:11.9
#26 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:18:33.4
#21 of 44 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:20:45.0
#24 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:22:06.5
#23 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:23:27.8
#24 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:24:49.3
#20 of 44 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:26:10.7
#24 of 44 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:27:32.1
#19 of 44 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:28:51.2
#24 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:31:24.8
#23 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:32:37.0
#23 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:33:58.3
#23 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:35:19.6
#25 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:36:40.9
#25 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:38:02.3
#19 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:39:23.9
#27 of 46 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:40:14.9
#17 of 45 solutions flagged due to SNR < 3 in spw=12 at 2017/10/15/08:41:57.5
#20 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:43:03.0
#19 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:44:15.1
#18 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:45:36.2
#19 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:46:57.7
#17 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:48:18.4
#19 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:49:39.6
#18 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:51:00.7
#26 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:51:53.9
#19 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:53:42.4
#20 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:55:03.7
#21 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:56:15.4
#20 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:58:01.2
#16 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/08:59:12.9
#21 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:00:34.1
#24 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:01:55.1
#20 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:04:10.7
#15 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:05:31.7
#21 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:06:52.5
#19 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:08:13.6
#20 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:09:34.3
#19 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:10:55.3
#23 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:12:22.3
#22 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:14:54.2
#21 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:16:06.0
#18 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:17:27.0
#18 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:18:47.6
#21 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:20:08.4
#21 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:21:29.2
#20 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:22:50.1
#31 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:23:42.3
#18 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:25:30.3
#21 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:26:51.4
#22 of 44 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:28:03.1
#22 of 43 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:29:48.2
#18 of 41 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:30:59.7
#20 of 41 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:32:20.4
#20 of 40 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:33:41.2
#21 of 38 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:35:56.8
#16 of 39 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:37:17.4
#18 of 39 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:38:38.2
#19 of 39 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:39:58.8
#14 of 37 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:41:19.5
#17 of 37 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:42:40.1
#18 of 40 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:43:59.2
#19 of 40 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:46:38.4
#22 of 40 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:47:49.9
#18 of 39 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:49:10.3
#16 of 38 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:50:30.7
#17 of 37 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:51:51.6
#20 of 37 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:53:11.7
#21 of 38 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:54:32.2
#31 of 38 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:55:25.6
#17 of 37 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:57:12.1
#22 of 34 solutions flagged due to SNR < 3 in spw=16 at 2017/10/16/09:58:23.8
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:39:46.5
# 1 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:40:34.9
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:41:35.3
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:42:35.8
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:43:36.3
# 1 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:44:37.7
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:45:20.7
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:47:18.8
# 1 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:48:19.3
# 1 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:49:19.8
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:50:20.3
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:51:20.7
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:52:22.8
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:53:08.3
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:55:02.5
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:55:48.0
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/04:58:25.0
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:00:25.9
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:01:26.4
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:02:26.9
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:03:28.9
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:04:14.3
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:09:01.2
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:10:01.6
# 4 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:11:02.1
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:12:02.6
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:13:03.1
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:14:05.1
# 5 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:14:50.5
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:17:22.7
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:18:23.2
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:19:23.7
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:21:24.6
# 1 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:22:26.7
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:23:12.0
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:25:41.7
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:27:49.8
# 2 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:28:50.3
# 1 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:29:50.8
# 3 of 43 solutions flagged due to SNR < 3 in spw=20 at 2017/12/09/05:30:36.3
# 2 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:50:24.3
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:51:25.3
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:52:27.4
# 3 of 43 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/05:54:57.1
# 1 of 45 solutions flagged due to SNR < 3 in spw=24 at 2017/12/17/06:11:47.4
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/03:56:58.6
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:02:47.9
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:07:45.5
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:10:33.4
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:12:47.6
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:14:37.8
# 2 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:17:39.2
# 8 of 45 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:22:11.4
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:31:22.8
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:35:13.5
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:36:14.0
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:37:14.3
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:38:16.3
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:39:01.8
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:40:29.5
# 1 of 46 solutions flagged due to SNR < 3 in spw=28 at 2017/12/27/04:41:14.9
# 1 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/03:57:48.9
# 3 of 45 solutions flagged due to SNR < 3 in spw=32 at 2017/12/28/04:09:28.0

#Inspect gain tables interactively and decide whether to manually flag something
#plotms(LB_iteration2_p4,xaxis='time',yaxis='GainPhase',iteraxis='spw',showgui=True)
#Print calibration png file
#plotms(
#    LB_iteration2_p4,xaxis='time', yaxis='GainPhase',iteraxis='spw',exprange='all',
#    overwrite=True,showgui=False,
#    plotfile=os.path.join(LB_selfcal_iteration2_folder,f'{prefix}_LB_iteration2_p4_gain_phase_vs_time.png')
#)

#Inspect gain tables and decide whether to flag something
plotms(
    LB_iteration2_p4,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(LB_selfcal_iteration2_folder,f'{prefix}_LB_iteration2_p4_gain_phase_vs_time.png'),
)

#Flag problematic antennas
flagdata(vis=LB_iteration2_p4,mode='manual',spw='0', timerange='2017/10/10/09:04:20~2017/10/10/09:04:30',antenna='DA64,DV24')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='0', timerange='2017/10/10/09:07:00~2017/10/10/09:07:10',antenna='DA60,DA64')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='0', timerange='2017/10/10/09:11:20~2017/10/10/09:11:30',antenna='DA63')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='0', timerange='2017/10/10/09:12:40~2017/10/10/09:12:50',antenna='DA65')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='0', timerange='2017/10/10/09:19:30~2017/10/10/09:19:40',antenna='DA47')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='0', timerange='2017/10/10/09:23:00~2017/10/10/09:23:10',antenna='DA65')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='0', timerange='2017/10/10/09:25:40~2017/10/10/09:25:50',antenna='DV04')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='0', timerange='2017/10/10/09:33:20~2017/10/10/09:33:30',antenna='DA64')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='0', timerange='2017/10/10/09:38:50~2017/10/10/09:39:00',antenna='DA59,DV13')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='0', timerange='2017/10/10/09:51:10~2017/10/10/09:51:20',antenna='DV22')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='0', timerange='2017/10/10/10:01:20~2017/10/10/10:01:30',antenna='DV22')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='0', timerange='2017/10/10/10:09:00~2017/10/10/10:09:10',antenna='DA64')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='0', timerange='2017/10/10/10:14:20~2017/10/10/10:14:30',antenna='DV04')

flagdata(vis=LB_iteration2_p4,mode='manual',spw='4', timerange='2017/10/10/10:27:50~2017/10/10/10:28:00',antenna='DV22')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='4', timerange='2017/10/10/10:28:10~2017/10/10/10:28:20',antenna='DA52,DA60')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='4', timerange='2017/10/10/10:31:50~2017/10/10/10:32:00',antenna='DV24')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='4', timerange='2017/10/10/10:34:30~2017/10/10/10:34:40',antenna='DV16')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='4', timerange='2017/10/10/10:35:50~2017/10/10/10:36:00',antenna='DV12')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='4', timerange='2017/10/10/10:36:50~2017/10/10/10:37:00',antenna='DA60')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='4', timerange='2017/10/10/10:38:40~2017/10/10/10:38:50',antenna='DV22')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='4', timerange='2017/10/10/10:50:30~2017/10/10/10:50:40',antenna='DA44,DA53,DA54,DV12,DV15,PM04')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='4', timerange='2017/10/10/10:57:20~2017/10/10/10:57:30',antenna='DV13')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='4', timerange='2017/10/10/11:03:30~2017/10/10/11:03:40',antenna='DA54')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='4', timerange='2017/10/10/11:06:10~2017/10/10/11:06:20',antenna='DV15')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='4', timerange='2017/10/10/11:11:50~2017/10/10/11:12:00',antenna='DA41,DA42,DV01,DV13')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='4', timerange='2017/10/10/11:13:00~2017/10/10/11:13:10',antenna='DA42,DA48,DA53,DV01,DV17,DV22')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='4', timerange='2017/10/10/11:17:20~2017/10/10/11:17:30',antenna='DA41,DA44,DA48,DA53,DA60')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='4', timerange='2017/10/10/11:27:30~2017/10/10/11:27:40',antenna='DA59')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='4', timerange='2017/10/10/11:28:50~2017/10/10/11:29:00',antenna='DA65,DV17')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='4', timerange='2017/10/10/11:31:20~2017/10/10/11:31:30',antenna='DA47')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='4', timerange='2017/10/10/11:32:30~2017/10/10/11:32:40',antenna='DA59')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='4', timerange='2017/10/10/11:37:50~2017/10/10/11:38:00',antenna='DA58')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='4', timerange='2017/10/10/11:39:10~2017/10/10/11:39:20',antenna='DA48')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='4', timerange='2017/10/10/11:40:10~2017/10/10/11:40:20',antenna='DA60,DV02')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='4', timerange='2017/10/10/11:41:40~2017/10/10/11:41:50',antenna='DA52,DA60,DV22')

flagdata(vis=LB_iteration2_p4,mode='manual',spw='8', timerange='2017/10/11/07:29:30~2017/10/11/07:29:40',antenna='DV08')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='8', timerange='2017/10/11/07:29:50~2017/10/11/07:30:00',antenna='DA53')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='8', timerange='2017/10/11/07:31:00~2017/10/11/07:31:10',antenna='DA53,DV15,DV22')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='8', timerange='2017/10/11/07:32:30~2017/10/11/07:32:40',antenna='DA64,DV04')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='8', timerange='2017/10/11/07:33:50~2017/10/11/07:34:00',antenna='DA53,DV13')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='8', timerange='2017/10/11/07:36:30~2017/10/11/07:36:40',antenna='DA44,DV15')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='8', timerange='2017/10/11/07:38:40~2017/10/11/07:38:50',antenna='DV08')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='8', timerange='2017/10/11/07:45:00~2017/10/11/07:45:10',antenna='DA53')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='8', timerange='2017/10/11/07:47:40~2017/10/11/07:47:50',antenna='DV12')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='8', timerange='2017/10/11/07:49:00~2017/10/11/07:49:10',antenna='DA46')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='8', timerange='2017/10/11/07:53:50~2017/10/11/07:54:00',antenna='DA45,DA60')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='8', timerange='2017/10/11/08:10:40~2017/10/11/08:10:50',antenna='DV05,DV16')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='8', timerange='2017/10/11/08:15:00~2017/10/11/08:15:10',antenna='DV20')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='8', timerange='2017/10/11/08:16:45~2017/10/11/08:16:55',antenna='DA65')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='8', timerange='2017/10/11/08:19:20~2017/10/11/08:19:30',antenna='DA58')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='8', timerange='2017/10/11/08:28:20~2017/10/11/08:28:30',antenna='DV05,DA52')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='8', timerange='2017/10/11/08:29:40~2017/10/11/08:29:50',antenna='DA42,DA44,DA53,DV15,DV20')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='8', timerange='2017/10/11/08:31:00~2017/10/11/08:31:10',antenna='DA65')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='8', timerange='2017/10/11/08:36:00~2017/10/11/08:36:10',antenna='DV20,PM01')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='8', timerange='2017/10/11/08:37:20~2017/10/11/08:37:30',antenna='DV03')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='8', timerange='2017/10/11/08:38:45~2017/10/11/08:38:55',antenna='DV08')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='8', timerange='2017/10/11/08:42:20~2017/10/11/08:42:30',antenna='DV16')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='8', timerange='2017/10/11/08:45:10~2017/10/11/08:45:15',antenna='DA47')

flagdata(vis=LB_iteration2_p4,mode='manual',spw='12',timerange='2017/10/15/07:27:20~2017/10/15/07:27:30',antenna='DA58,DV05')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='12',timerange='2017/10/15/07:28:30~2017/10/15/07:28:40',antenna='DV02')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='12',timerange='2017/10/15/07:31:20~2017/10/15/07:31:30',antenna='DA65')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='12',timerange='2017/10/15/07:32:40~2017/10/15/07:32:50',antenna='DV13,DV15')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='12',timerange='2017/10/15/07:39:20~2017/10/15/07:39:30',antenna='DV03')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='12',timerange='2017/10/15/07:59:20~2017/10/15/07:59:30',antenna='DA44,DA53')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='12',timerange='2017/10/15/08:03:20~2017/10/15/08:03:30',antenna='DA44,PM04')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='12',timerange='2017/10/15/08:04:40~2017/10/15/08:04:50',antenna='DA53,DV13')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='12',timerange='2017/10/15/08:06:00~2017/10/15/08:06:10',antenna='DA42,DA63')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='12',timerange='2017/10/15/08:11:20~2017/10/15/08:11:30',antenna='DV05')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='12',timerange='2017/10/15/08:14:30~2017/10/15/08:14:40',antenna='DA65,DV10,PM03')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='12',timerange='2017/10/15/08:22:00~2017/10/15/08:22:10',antenna='DV12')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='12',timerange='2017/10/15/08:26:10~2017/10/15/08:26:20',antenna='DA48,PM04')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='12',timerange='2017/10/15/08:28:50~2017/10/15/08:29:00',antenna='DA58')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='12',timerange='2017/10/15/08:31:20~2017/10/15/08:31:30',antenna='DV24')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='12',timerange='2017/10/15/08:38:00~2017/10/15/08:38:10',antenna='DA63')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='12',timerange='2017/10/15/08:40:10~2017/10/15/08:40:20',antenna='DA51,DV05,DV12')

flagdata(vis=LB_iteration2_p4,mode='manual',spw='16',timerange='2017/10/16/08:48:10~2017/10/16/08:48:20',antenna='DA42')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='16',timerange='2017/10/16/08:49:30~2017/10/16/08:49:40',antenna='DA42,DA51,DA55,PM01,PM03')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='16',timerange='2017/10/16/08:51:50~2017/10/16/08:52:00',antenna='DV10')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='16',timerange='2017/10/16/08:53:40~2017/10/16/08:53:50',antenna='DV24')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='16',timerange='2017/10/16/08:55:00~2017/10/16/08:55:10',antenna='DV16')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='16',timerange='2017/10/16/08:58:00~2017/10/16/08:58:10',antenna='DV15')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='16',timerange='2017/10/16/09:00:30~2017/10/16/09:00:40',antenna='DA44,DA53')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='16',timerange='2017/10/16/09:04:10~2017/10/16/09:04:20',antenna='DV23')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='16',timerange='2017/10/16/09:09:30~2017/10/16/09:09:40',antenna='DV11')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='16',timerange='2017/10/16/09:17:20~2017/10/16/09:17:30',antenna='DA60')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='16',timerange='2017/10/16/09:18:40~2017/10/16/09:18:50',antenna='DV08')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='16',timerange='2017/10/16/09:21:20~2017/10/16/09:21:30',antenna='DA53')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='16',timerange='2017/10/16/09:23:40~2017/10/16/09:23:50',antenna='DA44,DA53')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='16',timerange='2017/10/16/09:25:00~2017/10/16/10:05:00',antenna='')

flagdata(vis=LB_iteration2_p4,mode='manual',spw='20',timerange='2017/12/09/04:48:10~2017/12/09/04:48:20',antenna='DV15')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='20',timerange='2017/12/09/04:49:15~2017/12/09/04:49:25',antenna='DV15')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='20',timerange='2017/12/09/04:50:15~2017/12/09/04:50:25',antenna='DV15')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='20',timerange='2017/12/09/04:51:15~2017/12/09/04:51:25',antenna='DV15')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='20',timerange='2017/12/09/04:52:15~2017/12/09/04:52:25',antenna='DV11')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='20',timerange='2017/12/09/04:55:00~2017/12/09/04:55:10',antenna='DV15')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='20',timerange='2017/12/09/04:55:45~2017/12/09/04:55:50',antenna='DV11')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='20',timerange='2017/12/09/04:59:20~2017/12/09/04:59:30',antenna='DV11,DV15')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='20',timerange='2017/12/09/05:00:20~2017/12/09/05:00:30',antenna='DV10')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='20',timerange='2017/12/09/05:01:25~2017/12/09/05:01:30',antenna='DV10')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='20',timerange='2017/12/09/05:02:20~2017/12/09/05:02:30',antenna='DV10')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='20',timerange='2017/12/09/05:04:10~2017/12/09/05:04:20',antenna='DA49,DA47')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='20',timerange='2017/12/09/05:10:00~2017/12/09/05:10:05',antenna='DV10')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='20',timerange='2017/12/09/05:11:00~2017/12/09/05:11:05',antenna='DV10')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='20',timerange='2017/12/09/05:14:00~2017/12/09/05:14:10',antenna='DA46,DA55,DV11')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='20',timerange='2017/12/09/05:17:20~2017/12/09/05:17:25',antenna='DV10')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='20',timerange='2017/12/09/05:20:20~2017/12/09/05:20:25',antenna='DA46,DV15')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='20',timerange='2017/12/09/05:21:20~2017/12/09/05:21:25',antenna='DA46')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='20',timerange='2017/12/09/05:22:25~2017/12/09/05:22:30',antenna='DA46,DA55')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='20',timerange='2017/12/09/05:23:10~2017/12/09/05:23:25',antenna='DA46')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='20',timerange='2017/12/09/05:25:40~2017/12/09/05:25:45',antenna='DA46,DA51,DA55')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='20',timerange='2017/12/09/05:27:45~2017/12/09/05:27:50',antenna='DA46,DA51,DV10')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='20',timerange='2017/12/09/05:30:35~2017/12/09/05:40:00',antenna='DA46')

flagdata(vis=LB_iteration2_p4,mode='manual',spw='24',timerange='2017/12/17/05:49:20~2017/12/17/05:49:30',antenna='DA55')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='24',timerange='2017/12/17/05:54:20~2017/12/17/05:54:30',antenna='DA41,DA42,DA43,DA44,DA45@W205,DA55,DV13')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='24',timerange='2017/12/17/05:54:50~2017/12/17/05:55:00',antenna='DA41,DA42,DA43,DA44')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='24',timerange='2017/12/17/05:58:15~2017/12/17/05:58:20',antenna='DV13,DA44,DA55,DV07')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='24',timerange='2017/12/17/05:59:15~2017/12/17/05:59:20',antenna='DA55,DV13,DV07')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='24',timerange='2017/12/17/06:00:15~2017/12/17/06:00:20',antenna='DV13')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='24',timerange='2017/12/17/06:01:20~2017/12/17/06:01:25',antenna='DV13,DV07')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='24',timerange='2017/12/17/06:02:20~2017/12/17/06:02:25',antenna='DV13,DV07')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='24',timerange='2017/12/17/06:03:20~2017/12/17/06:03:25',antenna='DV13,DV07')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='24',timerange='2017/12/17/06:04:05~2017/12/17/06:04:10',antenna='DV13')

flagdata(vis=LB_iteration2_p4,mode='manual',spw='28',timerange='2017/12/27/04:08:45~2017/12/27/04:08:50',antenna='DA42,DV13')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='28',timerange='2017/12/27/04:09:45~2017/12/27/04:09:50',antenna='DA42,DA44,DV11,DV13')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='28',timerange='2017/12/27/04:10:30~2017/12/27/04:10:35',antenna='DA42,DA44,DV13')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='28',timerange='2017/12/27/04:12:00~2017/12/27/04:12:10',antenna='DA44')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='28',timerange='2017/12/27/04:21:30~2017/12/27/04:21:40',antenna='DA55')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='28',timerange='2017/12/27/04:25:25~2017/12/27/04:25:35',antenna='DV11')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='28',timerange='2017/12/27/04:26:30~2017/12/27/04:26:40',antenna='DV11')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='28',timerange='2017/12/27/04:33:10~2017/12/27/04:33:20',antenna='DA44')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='28',timerange='2017/12/27/04:34:10~2017/12/27/04:34:20',antenna='DA44,DV11')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='28',timerange='2017/12/27/04:35:00~2017/12/27/04:35:20',antenna='DA45,DA44')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='28',timerange='2017/12/27/04:36:00~2017/12/27/04:36:20',antenna='DA44')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='28',timerange='2017/12/27/04:37:00~2017/12/27/04:37:20',antenna='DA44')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='28',timerange='2017/12/27/04:38:00~2017/12/27/04:38:20',antenna='DA44')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='28',timerange='2017/12/27/04:39:00~2017/12/27/04:39:20',antenna='DA44')

flagdata(vis=LB_iteration2_p4,mode='manual',spw='32',timerange='2017/12/28/03:45:00~2017/12/28/03:45:10',antenna='DA55')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='32',timerange='2017/12/28/03:55:00~2017/12/28/03:55:10',antenna='DA55')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='32',timerange='2017/12/28/03:56:00~2017/12/28/03:56:10',antenna='DV11,DA55')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='32',timerange='2017/12/28/03:57:00~2017/12/28/03:57:10',antenna='DA55')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='32',timerange='2017/12/28/03:59:10~2017/12/28/03:59:20',antenna='DA55')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='32',timerange='2017/12/28/04:00:00~2017/12/28/04:00:10',antenna='DA55,DA44')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='32',timerange='2017/12/28/04:02:50~2017/12/28/04:03:00',antenna='DA44')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='32',timerange='2017/12/28/04:03:50~2017/12/28/04:04:00',antenna='DA44,DA55,DA46')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='32',timerange='2017/12/28/04:04:50~2017/12/28/04:05:00',antenna='DA44')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='32',timerange='2017/12/28/04:09:20~2017/12/28/04:09:30',antenna='DA42,DV13')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='32',timerange='2017/12/28/04:12:45~2017/12/28/04:12:55',antenna='DA42,DA55,DV11,DV13')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='32',timerange='2017/12/28/04:17:50~2017/12/28/04:18:00',antenna='DA55')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='32',timerange='2017/12/28/04:18:35~2017/12/28/04:18:45',antenna='DA55')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='32',timerange='2017/12/28/04:20:25~2017/12/28/04:20:35',antenna='DV11')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='32',timerange='2017/12/28/04:21:25~2017/12/28/04:21:35',antenna='DV11,DV13,DA44')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='32',timerange='2017/12/28/04:22:25~2017/12/28/04:22:35',antenna='DV11,DV13')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='32',timerange='2017/12/28/04:30:10~2017/12/28/04:30:30',antenna='DV13')

#Inspect gain tables to check if flagging worked
plotms(
    LB_iteration2_p4,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(LB_selfcal_iteration2_folder,f'{prefix}_LB_iteration2_p4_gain_phase_vs_time_flagged.png'),
)

#Apply the calibration gains and split-off the corrected ms
applycal(
    vis=LB_iteration2_cont_p3+'.ms',spw=LB_contspws,spwmap=LB_spw_mapping,
    gaintable=[LB_iteration2_p4],interp='linearPD',calwt=True,applymode='calonly'
)

LB_iteration2_cont_p4 = LB_iteration2_cont_p3.replace('p3','p4')
os.system('rm -rf %s.ms*' %LB_iteration2_cont_p4)
split(vis=LB_iteration2_cont_p3+'.ms',outputvis=LB_iteration2_cont_p4+'.ms',datacolumn='corrected')

#CLEAN, again with a ~6sigma threshold
tclean_wrapper(
    vis       = LB_iteration2_cont_p4+'.ms',
    imagename = LB_iteration2_cont_p4,
    threshold = '0.0386mJy',
    **LB_tclean_wrapper_kwargs
)
estimate_SNR(LB_iteration2_cont_p4+'.image',disk_mask=LB_mask,noise_mask=noise_annulus_LB)
generate_image_png(
    LB_iteration2_cont_p4+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_iteration2_LB,10*rms_iteration2_LB],
    save_folder=LB_selfcal_iteration2_folder
)
#MWC_758_SBLB_contp4.image
#Beam 0.065 arcsec x 0.045 arcsec (7.20 deg)
#Flux inside disk mask: 66.13 mJy
#Peak intensity of source: 1.21 mJy/beam
#rms: 6.54e-03 mJy/beam
#Peak SNR: 185.23

#MWC_758_SBLB_iteration2_contp4.image
#Beam 0.063 arcsec x 0.044 arcsec (7.28 deg)
#Flux inside disk mask: 66.74 mJy
#Peak intensity of source: 1.13 mJy/beam
#rms: 6.41e-03 mJy/beam
#Peak SNR: 176.44

#CLEAN with a ~1sigma threshold before amplitude self-cal
tclean_wrapper(
    vis       = LB_iteration2_cont_p4+'.ms',
    imagename = LB_iteration2_cont_p4,
    threshold = '6.41e-03mJy',
    **LB_tclean_wrapper_kwargs
)
estimate_SNR(LB_iteration2_cont_p4+'.image',disk_mask=LB_mask,noise_mask=noise_annulus_LB)
generate_image_png(
    LB_iteration2_cont_p4+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_iteration2_LB,10*rms_iteration2_LB],
    save_folder=LB_selfcal_iteration2_folder
)
#MWC_758_SBLB_contp4.image
#Beam 0.065 arcsec x 0.045 arcsec (7.20 deg)
#Flux inside disk mask: 57.78 mJy
#Peak intensity of source: 1.14 mJy/beam
#rms: 6.49e-03 mJy/beam
#Peak SNR: 175.43

#MWC_758_SBLB_iteration2_contp4.image
#Beam 0.063 arcsec x 0.044 arcsec (7.28 deg)
#Flux inside disk mask: 58.05 mJy
#Peak intensity of source: 1.07 mJy/beam
#rms: 6.36e-03 mJy/beam
#Peak SNR: 167.98

#First round of amplitude self-cal
LB_iteration2_ap0 = prefix+'_SBLB_iteration2.ap0'
os.system('rm -rf '+LB_iteration2_ap0)
gaincal(
    vis=LB_iteration2_cont_p4+'.ms',caltable=LB_iteration2_ap0,
    gaintype='T',combine='scan,spw',calmode='ap',solint='inf',
    spw=LB_contspws,refant=LB_refant,
    minsnr=5.0,minblperant=4,solnorm=False
)
#4 of 43 solutions flagged due to SNR < 5 in spw=8  at 2017/10/11/08:03:57.6
#7 of 46 solutions flagged due to SNR < 5 in spw=12 at 2017/10/15/07:57:59.8
#9 of 44 solutions flagged due to SNR < 5 in spw=16 at 2017/10/16/09:13:43.8

#Inspect gain tables interactively and decide whether to manually flag something
#plotms(LB_iteration2_ap0,xaxis='time', yaxis='GainPhase',iteraxis='spw')
#plotms(LB_iteration2_ap0,xaxis='time', yaxis='GainAmp',iteraxis='spw')
#
#Print calibration png file
#plotms(LB_iteration2_ap0,xaxis='time', yaxis='GainPhase',iteraxis='spw',exprange='all',
#       overwrite=True,showgui=False,
#       plotfile=os.path.join(LB_selfcal_iteration2_folder,f'{prefix}_LB_iteration2_gain_ap0_phase_vs_time.png'))
#plotms(LB_iteration2_ap0,xaxis='time', yaxis='GainAmp',iteraxis='spw',exprange='all',
#       overwrite=True,showgui=False,
#       plotfile=os.path.join(LB_selfcal_iteration2_folder,f'{prefix}_LB_iteration2_gain_ap0_amp_vs_time.png'))

#Inspect gain tables and decide whether to flag something
plotms(
    LB_iteration2_ap0,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(LB_selfcal_iteration2_folder,f'{prefix}_LB_iteration2_ap0_gain_phase_vs_time.png'),
)
plotms(
    LB_iteration2_ap0,xaxis='time',yaxis='GainAmp',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(LB_selfcal_iteration2_folder,f'{prefix}_LB_iteration2_ap0_gain_amp_vs_time.png'),
)

#Flag problematic antennas
flagdata(vis=LB_iteration2_ap0,mode='clip',clipminmax=[0.8,1.2],clipoutside=True,datacolumn='CPARAM')
flagdata(vis=LB_iteration2_ap0,mode='manual',spw='4', antenna='DV24')
flagdata(vis=LB_iteration2_ap0,mode='manual',spw='8', antenna='DA63,DA64,DV04,DV17')
flagdata(vis=LB_iteration2_ap0,mode='manual',spw='12',antenna='DA59')
flagdata(vis=LB_iteration2_ap0,mode='manual',spw='24',antenna='DV07')
flagdata(vis=LB_iteration2_ap0,mode='manual',spw='32',antenna='DA44,DV11')

#Inspect gain tables to check if flagging worked
plotms(
    LB_iteration2_ap0,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(LB_selfcal_iteration2_folder,f'{prefix}_LB_iteration2_ap0_gain_phase_vs_time_flagged.png'),
)
plotms(
    LB_iteration2_ap0,xaxis='time',yaxis='GainAmp',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(LB_selfcal_iteration2_folder,f'{prefix}_LB_iteration2_ap0_gain_amp_vs_time_flagged.png'),
)

#Apply the calibration gains and split-off the corrected ms
applycal(
    vis=LB_iteration2_cont_p4+'.ms',spw=LB_contspws,spwmap=LB_spw_mapping,
    gaintable=[LB_iteration2_ap0],interp='linearPD',calwt=True,applymode='calonly'
)
LB_iteration2_cont_ap0 = prefix+'_SBLB_iteration2_contap0'
os.system('rm -rf %s.ms*' %LB_iteration2_cont_ap0)
split(vis=LB_iteration2_cont_p4+'.ms',outputvis=LB_iteration2_cont_ap0+'.ms',datacolumn='corrected')

#CLEAN with a ~1sigma threshold before amplitude self-cal
tclean_wrapper(
    vis=LB_iteration2_cont_ap0+'.ms',
    imagename = LB_iteration2_cont_ap0,
    threshold = '6.36e-03mJy',
    **LB_tclean_wrapper_kwargs
)
estimate_SNR(LB_iteration2_cont_ap0+'.image',disk_mask=LB_mask,noise_mask=noise_annulus_LB)
generate_image_png(
    LB_iteration2_cont_ap0+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_iteration2_LB,10*rms_iteration2_LB],
    save_folder=LB_selfcal_iteration2_folder
)
#MWC_758_SBLB_contap0.image
#Beam 0.067 arcsec x 0.046 arcsec (7.20 deg)
#Flux inside disk mask: 57.75 mJy
#Peak intensity of source: 1.19 mJy/beam
#rms: 6.50e-03 mJy/beam
#Peak SNR: 183.45

#MWC_758_SBLB_iteration2_contap0.image
#Beam 0.066 arcsec x 0.045 arcsec (8.24 deg)
#Flux inside disk mask: 57.93 mJy
#Peak intensity of source: 1.16 mJy/beam
#rms: 6.44e-03 mJy/beam
#Peak SNR: 179.98
"""
#Try ampl self-cal on scan length intervals
LB_iteration2_ap1 = prefix+'_SBLB_iteration2.ap1'
os.system('rm -rf '+LB_iteration2_ap1)
gaincal(
    vis=LB_iteration2_cont_ap0+'.ms',caltable=LB_iteration2_ap1,
    gaintype='T',combine='spw',calmode='ap',solint='inf',
    spw=LB_contspws,refant=LB_refant,
    minsnr=5.0,minblperant=4,solnorm=False
)

#Inspect gain tables interactively and decide whether to manually flag something
#plotms(LB_iteration2_ap1,xaxis='time', yaxis='GainPhase',iteraxis='spw')
#plotms(LB_iteration2_ap1,xaxis='time', yaxis='GainAmp',iteraxis='spw')
#
#Print calibration png file
#plotms(LB_iteration2_ap1,xaxis='time', yaxis='GainPhase',iteraxis='spw',exprange='all',
#       overwrite=True,showgui=False,
#       plotfile=os.path.join(LB_selfcal_iteration2_folder,f'{prefix}_LB_iteration2_gain_ap1_phase_vs_time.png'))
#plotms(LB_iteration2_ap1,xaxis='time', yaxis='GainAmp',iteraxis='spw',exprange='all',
#       overwrite=True,showgui=False,
#       plotfile=os.path.join(LB_selfcal_iteration2_folder,f'{prefix}_LB_iteration2_gain_ap1_amp_vs_time.png'))

#Inspect gain tables and decide whether to flag something
plotms(
    LB_iteration2_ap1,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(LB_selfcal_iteration2_folder,f'{prefix}_LB_iteration2_ap1_gain_phase_vs_time.png'),
)
plotms(
    LB_iteration2_ap1,xaxis='time',yaxis='GainAmp',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(LB_selfcal_iteration2_folder,f'{prefix}_LB_iteration2_ap1_gain_amp_vs_time.png'),
)

#Flag problematic antennas
flagdata(vis=LB_iteration2_ap1,mode='clip',clipminmax=[0.8,1.2],clipoutside=True,datacolumn='CPARAM')

flagdata(vis=LB_iteration2_ap1,mode='manual',spw='0', timerange='2017/10/10/09:07:00~2017/10/10/09:07:10',antenna='DA44,DV11,PM04')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='0', timerange='2017/10/10/09:08:25~2017/10/10/09:08:30',antenna='DA45')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='0', timerange='2017/10/10/09:11:20~2017/10/10/09:11:25',antenna='DA54')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='0', timerange='2017/10/10/09:13:50~2017/10/10/09:14:00',antenna='DA53')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='0', timerange='2017/10/10/09:24:20~2017/10/10/09:24:30',antenna='DA44,DV22,PM03')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='0', timerange='2017/10/10/09:27:00~2017/10/10/09:27:10',antenna='DA57')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='0', timerange='2017/10/10/09:32:10~2017/10/10/09:32:20',antenna='DA53,PM04')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='0', timerange='2017/10/10/09:36:00~2017/10/10/09:36:10',antenna='DA47,DV15')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='0', timerange='2017/10/10/09:37:20~2017/10/10/09:37:30',antenna='DV15')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='0', timerange='2017/10/10/09:38:50~2017/10/10/09:39:00',antenna='PM04')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='0', timerange='2017/10/10/09:44:20~2017/10/10/09:44:30',antenna='DV15')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='0', timerange='2017/10/10/09:49:50~2017/10/10/09:50:00',antenna='DA42,DA44,DA53,DV15,PM04')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='0', timerange='2017/10/10/09:54:40~2017/10/10/09:54:50',antenna='DV24')

flagdata(vis=LB_iteration2_ap1,mode='manual',spw='4', timerange='2017/10/10/10:29:10~2017/10/10/10:29:20',antenna='DA46')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='4', timerange='2017/10/10/10:34:30~2017/10/10/10:34:40',antenna='DA53,DA54,DV24')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='4', timerange='2017/10/10/10:35:50~2017/10/10/10:36:00',antenna='DA53,DA54,DV15,DV24,PM04')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='4', timerange='2017/10/10/10:36:50~2017/10/10/10:37:00',antenna='DV15')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='4', timerange='2017/10/10/10:40:00~2017/10/10/10:40:10',antenna='DA53,PM04')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='4', timerange='2017/10/10/10:44:20~2017/10/10/10:44:40',antenna='DV24')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='4', timerange='2017/10/10/10:47:00~2017/10/10/10:47:20',antenna='DA54,DV24,DA44')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='4', timerange='2017/10/10/10:49:10~2017/10/10/10:49:20',antenna='DA53,DV15,PM04')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='4', timerange='2017/10/10/10:50:30~2017/10/10/10:50:40',antenna='DA54')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='4', timerange='2017/10/10/10:51:50~2017/10/10/10:52:00',antenna='DA54,DV15,DV24')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='4', timerange='2017/10/10/10:53:10~2017/10/10/10:53:20',antenna='DA44,DV15,PM04,DA53')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='4', timerange='2017/10/10/10:54:30~2017/10/10/10:54:40',antenna='DA54,DV24')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='4', timerange='2017/10/10/10:55:50~2017/10/10/10:56:00',antenna='DA44,DA53,DV15,PM04')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='4', timerange='2017/10/10/10:57:20~2017/10/10/10:57:30',antenna='PM04')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='4', timerange='2017/10/10/10:59:40~2017/10/10/10:59:50',antenna='DA54,PM04')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='4', timerange='2017/10/10/11:02:10~2017/10/10/11:02:20',antenna='DV24')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='4', timerange='2017/10/10/11:06:10~2017/10/10/11:06:20',antenna='DV22')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='4', timerange='2017/10/10/11:07:30~2017/10/10/11:07:40',antenna='DA53,DV15,DA44,PM04')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='4', timerange='2017/10/10/11:10:30~2017/10/10/11:10:40',antenna='DA53')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='4', timerange='2017/10/10/11:11:50~2017/10/10/11:12:00',antenna='DA44,DA53')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='4', timerange='2017/10/10/11:13:00~2017/10/10/11:13:10',antenna='DA48,DA53')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='4', timerange='2017/10/10/11:14:50~2017/10/10/11:15:00',antenna='DA44')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='4', timerange='2017/10/10/11:24:50~2017/10/10/11:25:00',antenna='DV24,DA53,DV15')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='4', timerange='2017/10/10/11:27:30~2017/10/10/11:27:40',antenna='DV15')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='4', timerange='2017/10/10/11:37:50~2017/10/10/11:38:00',antenna='DA54,DV24')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='4', timerange='2017/10/10/11:40:10~2017/10/10/11:40:20',antenna='DV02')

flagdata(vis=LB_iteration2_ap1,mode='manual',spw='8', timerange='2017/10/11/07:31:00~2017/10/11/07:31:10',antenna='DA53,DA56,DA60,DV15,PM04')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='8', timerange='2017/10/11/07:33:50~2017/10/11/07:34:00',antenna='DA44,DA47,DV11,DV15')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='8', timerange='2017/10/11/07:46:20~2017/10/11/07:46:30',antenna='DA53')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='8', timerange='2017/10/11/07:53:50~2017/10/11/07:54:00',antenna='DV08,PM03')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='8', timerange='2017/10/11/07:56:40~2017/10/11/07:56:50',antenna='DV04,DA53')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='8', timerange='2017/10/11/07:58:00~2017/10/11/07:58:10',antenna='PM04')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='8', timerange='2017/10/11/08:03:00~2017/10/11/08:03:10',antenna='DA49')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='8', timerange='2017/10/11/08:04:20~2017/10/11/08:04:40',antenna='DV11')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='8', timerange='2017/10/11/08:07:10~2017/10/11/08:07:20',antenna='DA49,DV04,DV22,PM03')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='8', timerange='2017/10/11/08:09:50~2017/10/11/08:10:00',antenna='DA54,DV12,DA53')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='8', timerange='2017/10/11/08:36:00~2017/10/11/08:36:10',antenna='DA58,DV05,DA42')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='8', timerange='2017/10/11/08:37:20~2017/10/11/08:37:30',antenna='DA56,DV04,DA50')

flagdata(vis=LB_iteration2_ap1,mode='manual',spw='12',timerange='2017/10/15/07:30:00~2017/10/15/07:30:10',antenna='DA53')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='12',timerange='2017/10/15/07:32:40~2017/10/15/07:32:50',antenna='PM04')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='12',timerange='2017/10/15/07:38:00~2017/10/15/07:38:10',antenna='PM04')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='12',timerange='2017/10/15/08:12:40~2017/10/15/08:12:50',antenna='DV25')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='12',timerange='2017/10/15/08:15:50~2017/10/15/08:16:00',antenna='DA62')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='12',timerange='2017/10/15/08:38:00~2017/10/15/08:38:10',antenna='DA41,DV25,DV14')

flagdata(vis=LB_iteration2_ap1,mode='manual',spw='20',timerange='2017/12/09/05:25:30~2017/12/09/05:26:00',antenna='DA51')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='24',timerange='2017/12/17/05:54:00~2017/12/17/05:54:30',antenna='DA43,DA44')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='24',timerange='2017/12/17/06:00:55~2017/12/17/06:01:05',antenna='DV13')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='32',timerange='2017/12/28/03:59:30~2017/12/28/03:59:40',antenna='DA55')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='32',timerange='2017/12/28/04:30:00~2017/12/28/04:30:30',antenna='DV13')

#Inspect gain tables to check if flagging worked
plotms(
    LB_iteration2_ap1,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(LB_selfcal_iteration2_folder,f'{prefix}_LB_iteration2_ap1_gain_phase_vs_time_flagged.png'),
)
plotms(
    LB_iteration2_ap1,xaxis='time',yaxis='GainAmp',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(LB_selfcal_iteration2_folder,f'{prefix}_LB_iteration2_ap1_gain_amp_vs_time_flagged.png'),
)

#Apply the calibration gains and split-off the corrected ms
applycal(
    vis=LB_iteration2_cont_ap0+'.ms',spw=LB_contspws,spwmap=LB_spw_mapping,
    gaintable=[LB_iteration2_ap1],interp='linearPD',calwt=True,applymode='calonly'
)
LB_iteration2_cont_ap1 = prefix+'_SBLB_iteration2_contap1'
os.system('rm -rf %s.ms*' %LB_iteration2_cont_ap1)
split(vis=LB_iteration2_cont_ap0+'.ms',outputvis=LB_iteration2_cont_ap1+'.ms',datacolumn='corrected')

#CLEAN with a ~1sigma threshold before amplitude self-cal
tclean_wrapper(
    vis       = LB_iteration2_cont_ap1+'.ms',
    imagename = LB_iteration2_cont_ap1,
    threshold = '0.0060mJy',
    **LB_tclean_wrapper_kwargs
)
estimate_SNR(LB_iteration2_cont_ap1+'.image',disk_mask=LB_mask,noise_mask=noise_annulus_LB)
generate_image_png(
    LB_iteration2_cont_ap1+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_iteration2_LB,10*rms_iteration2_LB],
    save_folder=LB_selfcal_iteration2_folder
)
#CQ_Tau_SBLB_iteration2_contap1.image
#Beam 0.081 arcsec x 0.059 arcsec (11.53 deg)
#Flux inside disk mask: 147.09 mJy
#Peak intensity of source: 2.65 mJy/beam
#rms: 1.18e-02 mJy/beam
#Peak SNR: 224.17
"""
#Check again how LB phase-only selfcal improved things at each step
non_self_caled_LB_iteration2_vis = LB_iteration2_cont_p0
self_caled_LB_iteration2_visibilities = {
    'p1' :LB_iteration2_cont_p1,
    'p2' :LB_iteration2_cont_p2,
    'p3' :LB_iteration2_cont_p3,
    'p4' :LB_iteration2_cont_p4,
    'ap0':LB_iteration2_cont_ap0,
}

#for vis in self_caled_LB_iteration2_visibilities.values(): 
#    listobs(vis=vis+'.ms',listfile=vis+'.ms.listobs.txt',overwrite=True)
#
#LB_EBs = ('EB0','EB1','EB2','EB3','EB4','EB5','EB6','EB7','EB8','EB9')
#LB_EB_spws = ('0,1,2,3', '4,5,6,7', '8,9,10,11', '12,13,14,15', '16,17,18,19', '20,21,22,23', '24,25,26,27', '28,29,30,31', '32,33,34,35') #fill out by referring to listobs output
#
#for self_cal_step,self_caled_vis in self_caled_LB_iteration2_visibilities.items():
#    for EB_key,spw in zip(LB_EBs,LB_EB_spws):
#
#        print('selfcal_step: {}'.format(self_cal_step)+', EB: {}'.format(EB_key))
#
#        nametemplate = f'{prefix}_LB_iteration2_{EB_key}_{self_cal_step}_compare_amp_vs_time'
#
#        visibilities = [self_caled_vis+'.ms',non_self_caled_LB_iteration2_vis+'.ms']
#
#        plot_amp_vs_time_comparison(
#            nametemplate=nametemplate,visibilities=visibilities,spw=spw,
#            uvrange=uv_ranges['LB'],output_folder=LB_selfcal_iteration2_folder
#        )

#Set to the EB of the combined SBLB data that corresponds to flux_ref_EB
SBLB_flux_ref_EB = 8 #this is SB_EB3

all_LB_iteration2_visibilities = self_caled_LB_iteration2_visibilities.copy()
all_LB_iteration2_visibilities['p0'] = LB_iteration2_cont_p0

total_number_of_EBs = number_of_EBs['SB'] + number_of_EBs['LB']
for self_cal_step,vis_name in all_LB_iteration2_visibilities.items():
    #Split out SB EBs
    vis_ms = vis_name+'.ms'
    nametemplate = vis_ms.replace('.ms','_EB')
    split_all_obs(msfile=vis_ms,nametemplate=nametemplate)

    exported_ms = []
    for i in range(total_number_of_EBs):
        EB_vis = f'{nametemplate}{i}.ms'
        listobs(vis=EB_vis,listfile=EB_vis+'.listobs.txt',overwrite=True)
        #Export MS contents into numpy save files
        export_MS(EB_vis)
        exported_ms.append(EB_vis.replace('.ms','.vis.npz'))

    for i,exp_ms in enumerate(exported_ms):
        png_filename = f'iteration2_flux_comparison_EB{i}_to_EB{SBLB_flux_ref_EB}'+f'_SBLB_{self_cal_step}.png'
        plot_label = os.path.join(LB_selfcal_iteration2_folder,png_filename)

        estimate_flux_scale(
            reference=f'{nametemplate}{SBLB_flux_ref_EB}.vis.npz',
            comparison=exp_ms,incl=incl,PA=PA,plot_label=plot_label,#uvbins=np.arange(40.,300.,20.),
        )

    fluxscale = [1.,]*total_number_of_EBs
    plot_label = os.path.join(LB_selfcal_iteration2_folder,f'deprojected_vis_profiles_SBLB_iteration2_{self_cal_step}.png')
    plot_deprojected(filelist=exported_ms,fluxscale=fluxscale,PA=PA,incl=incl,show_err=True,plot_label=plot_label)

# In the concatenated SBLB file:
# EB0 = LB EB0
# EB1 = LB EB1
# EB2 = LB EB2
# EB3 = LB EB3
# EB4 = LB EB4
# EB5 = SB EB0
# EB6 = SB EB1
# EB7 = SB EB2
# EB8 = SB EB3

#iteration_1
#ratio          = [0.88988,0.89305,0.89593,0.90156,0.90895,0.92577] #MWC_758_SB_contp0...p5_EB0.vis.npz vs MWC_758_SB_EB3_initcont.vis.npz
#scaling_factor = [0.943  ,0.945  ,0.947  ,0.950  ,0.953  ,0.962  ]
#ratio          = [0.98340,0.97806,0.98762,0.99364,0.99671,0.99951] #MWC_758_SB_contp0...p5_EB1.vis.npz vs MWC_758_SB_EB3_initcont.vis.npz
#scaling_factor = [0.992  ,0.989  ,0.994  ,0.997  ,0.998  ,1.000  ]
#ratio          = [0.93677,0.96069,0.95945,0.96593,0.96917,0.97148] #MWC_758_SB_contp0...p5_EB2.vis.npz vs MWC_758_SB_EB3_initcont.vis.npz
#scaling_factor = [0.968  ,0.980  ,0.980  ,0.983  ,0.984  ,0.986  ]

#CLEANing with 1.0 robust (p0 -> ap0_afterp4)
#ratio          = [0.82477,0.91568,0.98312,1.00422,1.03841,1.03459] #MWC_758_SBLB_contp0...p5_EB0.vis.npz vs MWC_758_SBLB_contp0...p5_EB8.vis.npz
#scaling_factor = [0.908  ,0.957  ,0.992  ,1.002  ,1.019  ,1.017  ]
#ratio          = [0.82298,0.91481,0.94454,0.96959,1.02252,1.06883] #MWC_758_SBLB_contp0...p5_EB1.vis.npz vs MWC_758_SBLB_contp0...p5_EB8.vis.npz
#scaling_factor = [0.907  ,0.956  ,0.972  ,0.985  ,1.011  ,1.034  ]
#ratio          = [0.80492,0.92314,0.90848,0.91003,0.93038,0.94416] #MWC_758_SBLB_contp0...p5_EB2.vis.npz vs MWC_758_SBLB_contp0...p5_EB8.vis.npz
#scaling_factor = [0.897  ,0.961  ,0.953  ,0.954  ,0.965  ,0.972  ]
#ratio          = [0.94428,0.99588,1.00153,0.98510,1.01808,1.05846] #MWC_758_SBLB_contp0...p5_EB3.vis.npz vs MWC_758_SBLB_contp0...p5_EB8.vis.npz
#scaling_factor = [0.972  ,0.998  ,1.001  ,0.993  ,1.009  ,1.029  ]
#ratio          = [1.03013,1.06972,1.04402,1.04893,1.07551,1.07234] #MWC_758_SBLB_contp0...p5_EB4.vis.npz vs MWC_758_SBLB_contp0...p5_EB8.vis.npz
#scaling_factor = [1.015  ,1.034  ,1.022  ,1.024  ,1.037  ,1.036  ]
#ratio          = [0.92268,0.92291,0.92254,0.92251,0.92266,0.99209] #MWC_758_SBLB_contp0...p5_EB5.vis.npz vs MWC_758_SBLB_contp0...p5_EB8.vis.npz
#scaling_factor = [0.961  ,0.961  ,0.960  ,0.960  ,0.961  ,0.996  ]
#ratio          = [1.00050,1.00059,1.00074,0.99988,0.99991,0.98880] #MWC_758_SBLB_contp0...p5_EB6.vis.npz vs MWC_758_SBLB_contp0...p5_EB8.vis.npz
#scaling_factor = [1.000  ,1.000  ,1.000  ,1.000  ,1.000  ,0.994  ]
#ratio          = [0.97132,0.97124,0.97121,0.97124,0.97127,0.99629] #MWC_758_SBLB_contp0...p5_EB7.vis.npz vs MWC_758_SBLB_contp0...p5_EB8.vis.npz
#scaling_factor = [0.986  ,0.986  ,0.986  ,0.986  ,0.986  ,0.998  ]

#iteration_2
#ratio          = [0.89704,0.90017,0.90316,0.90890,0.91631,0.93330] #MWC_758_SB_contp0...p5_EB0.vis.npz vs MWC_758_SB_EB3_initcont.vis.npz
#scaling_factor = [0.947  ,0.949  ,0.950  ,0.953  ,0.957  ,0.966  ]
#ratio          = [0.99531,0.98994,0.99961,1.00568,1.00880,1.01165] #MWC_758_SB_contp0...p5_EB1.vis.npz vs MWC_758_SB_EB3_initcont.vis.npz
#scaling_factor = [0.998  ,0.995  ,1.000  ,1.003  ,1.004  ,1.006  ]
#ratio          = [0.94053,0.96452,0.96327,0.96981,0.97305,0.97537] #MWC_758_SB_contp0...p5_EB2.vis.npz vs MWC_758_SB_EB3_initcont.vis.npz
#scaling_factor = [0.970  ,0.982  ,0.981  ,0.985  ,0.986  ,0.988  ]

#CLEANing with 1.0 robust (p0 -> ap0_afterp4)
#ratio          = [0.79717,0.92643,0.93885,0.96156,0.99294,1.01935] #MWC_758_SBLB_contp0...p5_EB0.vis.npz vs MWC_758_SBLB_contp0...p5_EB8.vis.npz
#scaling_factor = [0.893  ,0.963  ,0.969  ,0.981  ,0.996  ,1.010  ]
#ratio          = [0.76958,0.87085,0.87812,0.91807,0.94547,1.06423] #MWC_758_SBLB_contp0...p5_EB1.vis.npz vs MWC_758_SBLB_contp0...p5_EB8.vis.npz
#scaling_factor = [0.877  ,0.933  ,0.937  ,0.958  ,0.972  ,1.032  ]
#ratio          = [0.85173,0.97668,0.99459,1.00013,1.04123,1.01879] #MWC_758_SBLB_contp0...p5_EB2.vis.npz vs MWC_758_SBLB_contp0...p5_EB8.vis.npz
#scaling_factor = [0.923  ,0.988  ,0.997  ,1.000  ,1.020  ,1.009  ]
#ratio          = [0.89173,0.94016,0.94277,0.94134,0.96509,1.04963] #MWC_758_SBLB_contp0...p5_EB3.vis.npz vs MWC_758_SBLB_contp0...p5_EB8.vis.npz
#scaling_factor = [0.944  ,0.970  ,0.971  ,0.970  ,0.982  ,1.025  ]
#ratio          = [0.95971,0.99677,0.97350,0.98486,0.99846,1.06931] #MWC_758_SBLB_contp0...p5_EB4.vis.npz vs MWC_758_SBLB_contp0...p5_EB8.vis.npz
#scaling_factor = [0.980  ,0.998  ,0.987  ,0.992  ,0.999  ,1.034  ]
#ratio          = [0.93017,0.93058,0.93013,0.93011,0.93027,0.99197] #MWC_758_SBLB_contp0...p5_EB5.vis.npz vs MWC_758_SBLB_contp0...p5_EB8.vis.npz
#scaling_factor = [0.964  ,0.965  ,0.964  ,0.964  ,0.965  ,0.996  ]
#ratio          = [1.01265,1.01273,1.01289,1.01206,1.01208,0.98880] #MWC_758_SBLB_contp0...p5_EB6.vis.npz vs MWC_758_SBLB_contp0...p5_EB8.vis.npz
#scaling_factor = [1.006  ,1.006  ,1.006  ,1.006  ,1.006  ,0.994  ]
#ratio          = [0.97521,0.97514,0.97513,0.97515,0.97518,0.99614] #MWC_758_SBLB_contp0...p5_EB7.vis.npz vs MWC_758_SBLB_contp0...p5_EB8.vis.npz
#scaling_factor = [0.988  ,0.987  ,0.987  ,0.987  ,0.988  ,0.998  ]

#END of COMBINED SB+LB phase-only self-cal iteration 2

#Split out final continuum ms table, with a 30s timebin
LB_iteration2_cont_averaged = f'{prefix}_time_ave_continuum'
os.system(f'rm -rf {LB_iteration2_cont_averaged}.ms*')
split(vis=LB_iteration2_cont_ap0+'.ms',outputvis=LB_iteration2_cont_averaged+'.ms',datacolumn='data',keepflags=False,timebin='30s')

#Now apply these solutions to the line data
calibrate_linedata_folder = get_figures_folderpath('9_apply_cal_to_lines')
make_figures_folder(calibrate_linedata_folder)

#Check that lines are not flagged in the non-averaged data
for params in data_params.values():
    plotfile = os.path.join(calibrate_linedata_folder,f'{prefix}_{params["name"]}_chan-v-amp_preselfcal_after_flagging.png')
    plotms(
        vis         = params['vis'],
        xaxis       = 'channel',
        yaxis       = 'amplitude',
        field       = params['field'],
        ydatacolumn = 'data',
        avgtime     = '1e8',
        avgscan     = True,
        avgbaseline = True,
        coloraxis   = 'corr',
        iteraxis    = 'spw',
        showgui     = False,
        exprange    = 'all',
        plotfile    = plotfile
    )

    plotfile = os.path.join(calibrate_linedata_folder,f'{prefix}_{params["name"]}_freq-v-amp_preselfcal_after_flagging.png')
    plotms(
        vis         = params['vis'],
        xaxis       = 'frequency',
        yaxis       = 'amplitude',
        field       = params['field'],
        ydatacolumn = 'data',
        avgtime     = '1e8',
        avgscan     = True,
        avgbaseline = False,
        coloraxis   = 'corr',
        iteraxis    = 'spw',
        showgui     = False,
        exprange    = 'all',
        plotfile    = plotfile
    )

#Apply the gaintables of individual EBs
for params in data_params.values():
    single_EB_p1 = prefix+'_'+params['name']+'_initcont.p1'
    vis          = prefix+'_'+params['name']+'.ms'
    applycal(vis=vis,spw=single_EB_contspws,spwmap=single_EB_spw_mapping,gaintable=[single_EB_p1],interp='linearPD',applymode='calonly',calwt=True)
    split(vis=vis,outputvis=prefix+'_'+params['name']+'_no_ave_selfcal.ms',datacolumn='corrected')
"""
#Align the data (skip)
#We re-align the non-averaged data, as we have done for the "initcont" .ms tables (from *_no_ave_selfcal.ms to *_no_ave_shift.ms)
print('Align no_ave data')
reference_ms = {
    'LB':reference_for_LB_alignment,
    'SB':reference_for_SB_alignment,
}
for params in data_params.values():
    unshifted_ms = prefix+'_'+params['name']+'_no_ave_selfcal.ms'
    array_key, _ = params['name'].split('_') #LB or SB
    offset       = alignment_offsets[params['name']]
    #npix and cell_size are not needed because we do not fit any offset
    alignment.align_measurement_sets(
        reference_ms  = reference_ms[array_key],
        align_ms      = [unshifted_ms],
        align_offsets = [offset],
        npix          = None,
        cell_size     = None,
    )
"""
#If you have re-scaled fluxes, you need to re-scale the shifted *no_ave* EBs as well
#rescale_flux(vis=prefix+'_LB_EB0_no_ave_selfcal.ms', gencalparameter=[1.024])#vis=prefix+'_LB_EB0_no_ave_selfcal_shift.ms'
#listobs(vis=prefix+'_LB_EB0_no_ave_selfcal_shift_rescaled.ms',listfile=prefix+'_LB_EB0_no_ave_selfcal_shift_rescaled.ms.listobs.txt',overwrite=True)
#listobs(vis=prefix+'_LB_EB0_no_ave_selfcal_rescaled.ms',listfile=prefix+'_LB_EB0_no_ave_selfcal_rescaled.ms.listobs.txt',overwrite=True)

#rescale_flux(vis=prefix+'_SB_EB0_no_ave_selfcal.ms', gencalparameter=[0.961])#vis=prefix+'_SB_EB0_no_ave_selfcal_shift.ms'
#listobs(vis=prefix+'_SB_EB0_no_ave_selfcal_shift_rescaled.ms',listfile=prefix+'_SB_EB0_no_ave_selfcal_shift_rescaled.ms.listobs.txt',overwrite=True)
#listobs(vis=prefix+'_SB_EB0_no_ave_selfcal_rescaled.ms',listfile=prefix+'_SB_EB0_no_ave_selfcal_rescaled.ms.listobs.txt',overwrite=True)

# Rescale the LBs
for i,gencalpar in enumerate([1.017,1.034,0.972,1.029,1.036]):
    # os.system('rm -rf '+prefix+f'_LB_EB{i}_no_ave_selfcal_rescaled.ms')
    rescale_flux(vis=prefix+f'_LB_EB{i}_no_ave_selfcal.ms', gencalparameter=[gencalpar])
    
    listobs(vis=prefix+f'_LB_EB{i}_no_ave_selfcal_rescaled.ms',listfile=prefix+f'_LB_EB{i}_no_ave_selfcal_rescaled.ms.listobs.txt',overwrite=True)

for i,gencalpar in enumerate([0.996,0.994,0.998]):
    # os.system('rm -rf '+prefix+f'_SB_EB{i}_no_ave_selfcal_rescaled.ms')
    rescale_flux(vis=prefix+f'_SB_EB{i}_no_ave_selfcal.ms', gencalparameter=[gencalpar])
    
    listobs(vis=prefix+f'_SB_EB{i}_no_ave_selfcal_rescaled.ms',listfile=prefix+f'_SB_EB{i}_no_ave_selfcal_rescaled.ms.listobs.txt',overwrite=True)

#Concat the non-averaged SB data
SB_combined = f'{prefix}_SB_no_ave_concat'
os.system('rm -rf %s.ms*' %SB_combined)
concat(
    vis=[
        f'{prefix}_SB_EB0_no_ave_selfcal_rescaled.ms',
        f'{prefix}_SB_EB1_no_ave_selfcal_rescaled.ms',
        f'{prefix}_SB_EB2_no_ave_selfcal_rescaled.ms',
        f'{prefix}_SB_EB3_no_ave_selfcal.ms',
    ],
    concatvis=SB_combined+'.ms',dirtol='0.1arcsec',copypointing=False
)
listobs(vis=SB_combined+'.ms',listfile=SB_combined+'.ms.listobs.txt',overwrite=True)
#2024-12-27 08:02:02     WARN    MSConcat::copySysCal    /data/beegfs/astro-storage/groups/benisty/frzagaria/SO_detections/CQTau/selfcal_products/CQ_Tau_SB_no_ave_concat.ms does not have a valid syscal table,
#    the MS to be appended, however, has one. Result won't have one.
#2024-12-27 08:02:02     WARN    MSConcat::concatenate (file /source/casa6/casatools/casacore/ms/MSOper/MSConcat.cc, line 1000)     Could not merge SysCal subtables 

#BE CAREFUL HERE
#Using gaintables from iteration2
applycal(
    vis        = SB_combined+'.ms',
    gaintable  = [SB_iteration2_p1,SB_iteration2_p1_bis,SB_iteration2_p2,SB_iteration2_p3,SB_iteration2_p4,SB_iteration2_p5],
    spw        = SB_contspws,
    spwmap     = [SB_spw_mapping]*6,
    interp     = ['linearPD']*6, 
    calwt      = True,
    applymode  = 'calonly',
    flagbackup = False,
)
SB_no_ave_selfcal = f'{prefix}_SB_no_ave_selfcal.ms'
os.system(f'rm -rf {SB_no_ave_selfcal}*')
split(vis=SB_combined+'.ms',outputvis=SB_no_ave_selfcal,datacolumn='corrected')
listobs(vis=SB_no_ave_selfcal,listfile=SB_no_ave_selfcal+'.listobs.txt',overwrite=True)

#Concat the non-averaged LB data
LB_combined = f'{prefix}_SBLB_no_ave_concat'
os.system('rm -rf %s.ms*' %LB_combined)
concat(
    vis=[SB_no_ave_selfcal]+[
        f'{prefix}_LB_EB0_no_ave_selfcal_rescaled.ms',
        f'{prefix}_LB_EB1_no_ave_selfcal_rescaled.ms',
        f'{prefix}_LB_EB2_no_ave_selfcal_rescaled.ms',
        f'{prefix}_LB_EB3_no_ave_selfcal_rescaled.ms',
        f'{prefix}_LB_EB4_no_ave_selfcal_rescaled.ms',
    ],
    concatvis=LB_combined+'.ms',dirtol='0.1arcsec',copypointing=False
)
listobs(vis=LB_combined+'.ms',listfile=LB_combined+'.ms.listobs.txt',overwrite=True)
#2024-12-27 08:23:00     SEVERE  getcell::TIME   Exception Reported: TableProxy::getCell: no such row
#2024-12-27 08:23:01     WARN    concat::::casa  Some but not all of the input MSs are lacking a populated POINTING table:
#2024-12-27 08:23:01     WARN    concat::::casa     0: CQ_Tau_SB_no_ave_selfcal.ms
#2024-12-27 08:23:01     WARN    concat::::casa  The joint dataset will not have a valid POINTING table.
#2024-12-27 08:25:23     WARN    MSConcat::copySysCal    /data/beegfs/astro-storage/groups/benisty/frzagaria/SO_detections/CQTau/selfcal_products/CQ_Tau_SBLB_no_ave_concat.ms does not have a valid syscal table,
#    the MS to be appended, however, has one. Result won't have one.
#2024-12-27 08:25:23     WARN    MSConcat::concatenate (file /source/casa6/casatools/casacore/ms/MSOper/MSConcat.cc, line 1000)     Could not merge SysCal subtables 
#2024-12-27 08:31:57     WARN    MSConcat::copySysCal    /data/beegfs/astro-storage/groups/benisty/frzagaria/SO_detections/CQTau/selfcal_products/CQ_Tau_SBLB_no_ave_concat.ms does not have a valid syscal table,
#    the MS to be appended, however, has one. Result won't have one.
#2024-12-27 08:31:57     WARN    MSConcat::concatenate (file /source/casa6/casatools/casacore/ms/MSOper/MSConcat.cc, line 1000)     Could not merge SysCal subtables

#BE CAREFUL HERE
#Using all gaintables from iteration2 even for LB
applycal(
    vis        = LB_combined+'.ms',
    gaintable  = [LB_iteration2_p1,LB_iteration2_p1_bis,LB_iteration2_p2,LB_iteration2_p3,LB_iteration2_p4,LB_iteration2_ap0],
    spw        = LB_contspws,
    spwmap     = [LB_spw_mapping]*6,
    interp     = ['linearPD']*6,
    calwt      = True,
    applymode  = 'calonly',
    flagbackup = False
)
SBLB_no_ave_selfcal = f'{prefix}_SBLB_no_ave_selfcal_time_ave.ms'
os.system(f'rm -rf {SBLB_no_ave_selfcal}*')
split(vis=LB_combined+'.ms',outputvis=SBLB_no_ave_selfcal,datacolumn='corrected',keepflags=False,timebin='30s') #Time average of 30s, tests show there is no difference with data without time average
listobs(vis=SBLB_no_ave_selfcal,listfile=SBLB_no_ave_selfcal+'.listobs.txt',overwrite=True)

#Check that the solutions have been applied correctly by flagging the line data, averaging and imaging continuum
#Continuum has to be the same imaged in the last step of the self-cal
SBLB_no_ave_selfcal = f'{prefix}_SBLB_no_ave_selfcal_time_ave.ms'
complete_dataset_dict = {
    'vis':       SBLB_no_ave_selfcal,
    'name':      'SBLB_concat',
    'field':     'MWC758',
    'line_spws': np.array([
        1,1,1,1,1, 2, 3,3,3,3,3,
        5,5,5,5,5, 6, 7,7,7,7,7,
        9,9,9,9,9, 10, 11,11,11,11,11,
        13,13,13,13,13, 14, 15,15,15,15,15,
        17,17,17,17,17, 18, 19,19,19,19,19,
        21,21,21,21,21, 22, 23,23,23,23,23,
        25,25,25,25,25, 26, 27,27,27,27,27,
        29,29,29,29,29, 30, 31,31,31,31,31,
        33,33,33,33,33, 34, 35,35,35,35,35
    ]), # list of spws containing lines
    'line_freqs':  np.array([
        rest_freq_H2CO_918_919,rest_freq_DCN_F21,rest_freq_DCN_F22,rest_freq_H2CO_303_202,rest_freq_H2CO_322_221,
        rest_freq_12CO,
        rest_freq_H2CO_322_221,rest_freq_H2CO_321_220,rest_freq_C18O,rest_freq_SO,rest_freq_13CO
    ]*9), #frequencies (Hz) corresponding to line_spws
    'spwcont_forplot': ['0']*9,
    'cont_spws':       '0~35',
    'width_array':     [2,2,32,32]*5+[8,8,192,192]*4,
}

#Use the output of get_flagchannels at the beginning of the script to define fitspw
# Flagchannels input string for LB_EB0: '1:123~125, 1:80~82, 1:80~82, 1:17~19, 1:1~3, 2:0~17, 3:1919~1919, 3:1676~1698, 3:856~879, 3:458~480, 3:0~20'

#Remember how the spws are organized in the SBLB files (seen from listobs)
# [0,1,2,3]     -> LB EB0
# [4,5,6,7]     -> LB EB1
# [8,9,10,11]   -> LB EB2
# [12,13,14,15] -> LB EB3
# [16,17,18,19] -> LB EB4
# [20,21,22,23] -> SB EB0
# [24,25,26,27] -> SB EB1
# [28,29,30,31] -> SB EB2
# [32,33,34,35] -> SB EB3
fitspw =  '0:0, 1:123~125, 1:80~82, 1:80~82, 1:17~19, 1:1~3, 2:0~17, 3:1919~1919, 3:1676~1698, 3:856~879, 3:458~480, 3:0~20, '\
         +'4:0, 5:123~125, 5:80~82, 5:80~82, 5:17~19, 5:1~3, 6:0~17, 7:1919~1919, 7:1676~1698, 7:856~879, 7:458~480, 7:0~20, '\
         +'8:0, 9:123~125, 9:80~82, 9:80~82, 9:17~19, 9:1~3, 10:0~17, 11:1919~1919, 11:1676~1698, 11:856~879, 11:458~480, 11:0~20, '\
         +'12:0, 13:123~125, 13:80~82, 13:80~82, 13:17~19, 13:1~3, 14:0~17, 15:1919~1919, 15:1676~1698, 15:856~879, 15:458~480, 15:0~20, '\
         +'16:0, 17:123~125, 17:80~82, 17:80~82, 17:17~19, 17:1~3, 18:0~17, 19:1919~1919, 19:1676~1698, 19:856~879, 19:458~480, 19:0~20, '\
         +'20:0, 21:123~125, 21:80~82, 21:80~82, 21:17~19, 21:1~3, 22:0~17, 23:1919~1919, 23:1676~1698, 23:856~879, 23:458~480, 23:0~20, '\
         +'24:0, 25:123~125, 25:80~82, 25:80~82, 25:17~19, 25:1~3, 26:0~17, 27:1919~1919, 27:1676~1698, 27:856~879, 27:458~480, 27:0~20, '\
         +'28:0, 29:123~125, 29:80~82, 29:80~82, 29:17~19, 29:1~3, 30:0~17, 31:1919~1919, 31:1676~1698, 31:856~879, 31:458~480, 31:0~20, '\
         +'32:0, 33:123~125, 33:80~82, 33:80~82, 33:17~19, 33:1~3, 34:0~17, 35:1919~1919, 35:1676~1698, 35:856~879, 35:458~480, 35:0~20'
     
avg_cont(
    ms_dict=complete_dataset_dict,output_prefix=prefix,flagchannels=fitspw,
    contspws=complete_dataset_dict['cont_spws'],width_array=complete_dataset_dict['width_array']
)

LB_tclean_wrapper_kwargs = {
    'deconvolver':'multiscale','scales':[0,2,4,8,16,24,32],
    'smallscalebias':0.6,'gain':0.3,'cycleniter':300,
    'niter':1000000,'mask':LB_mask,
    'cellsize':'0.005arcsec','imsize':2000,
    'parallel':use_parallel,'savemodel':'modelcolumn',
    'robust':1.0,'interactive':False,
    'gridder':'standard',
}

#Image the avg then cal with same parameters as in last step of self-cal
tclean_wrapper(
    vis       = LB_iteration2_cont_averaged+'.ms',
    imagename = LB_iteration2_cont_averaged+'_image', 
    threshold = '6.44e-03mJy',
    **LB_tclean_wrapper_kwargs
)

estimate_SNR(LB_iteration2_cont_averaged+'_image'+'.image',disk_mask=LB_mask,noise_mask=noise_annulus_LB)
generate_image_png(
    LB_iteration2_cont_averaged+'_image'+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_iteration2_LB,10*rms_iteration2_LB],save_folder=calibrate_linedata_folder
)
#MWC_758_SBLB_iteration2_contap0.image
#Beam 0.066 arcsec x 0.045 arcsec (7.24 deg)
#Flux inside disk mask: 57.93 mJy
#Peak intensity of source: 1.15 mJy/beam
#rms: 6.44e-03 mJy/beam
#Peak SNR: 179.33

#MWC_758_time_ave_continuum_image.image
#Beam 0.066 arcsec x 0.044 arcsec (7.68 deg)
#Flux inside disk mask: 57.87 mJy
#Peak intensity of source: 1.13 mJy/beam
#rms: 6.38e-03 mJy/beam
#Peak SNR: 176.98

#Image the cal then avg with same parameters as in last step of self-cal
complete_dataset_image = prefix+'_'+complete_dataset_dict['name']+'_initcont_image'
tclean_wrapper(
    vis       = prefix+'_'+complete_dataset_dict['name']+'_initcont.ms',
    imagename = complete_dataset_image, 
    threshold = '6.44e-03mJy',
    **LB_tclean_wrapper_kwargs
)
estimate_SNR(complete_dataset_image+'.image',disk_mask=LB_mask,noise_mask=noise_annulus_LB)
generate_image_png(
    complete_dataset_image+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_iteration2_LB,10*rms_iteration2_LB],save_folder=calibrate_linedata_folder
)
#MWC_758_SBLB_iteration2_contap0.image
#Beam 0.066 arcsec x 0.045 arcsec (7.24 deg)
#Flux inside disk mask: 57.93 mJy
#Peak intensity of source: 1.15 mJy/beam
#rms: 6.44e-03 mJy/beam
#Peak SNR: 179.33

#MWC_758_SBLB_concat_initcont_image.image
#Beam 0.066 arcsec x 0.044 arcsec (7.74 deg)
#Flux inside disk mask: 57.84 mJy
#Peak intensity of source: 1.13 mJy/beam
#rms: 6.38e-03 mJy/beam
#Peak SNR: 176.52

#Plot ratio of cal then avg, and avg then cal. It should be equal to ~one (only difference being the time average)
ref_image = LB_iteration2_cont_averaged+'_image'+'.image'
os.system('rm -rf '+complete_dataset_image+'.ratio')
immath(
    imagename=[ref_image,complete_dataset_image+'.image'],mode='evalexpr',
    outfile=complete_dataset_image+'.ratio',
    expr='iif(IM0 > 3*'+str(rms_iteration2_LB)+', IM1/IM0, 0)'
)
generate_image_png(
    f'{complete_dataset_image}.ratio',plot_sizes=[2*mask_semimajor,2*mask_semimajor],
    color_scale_limits=[0.5,1.5],image_units='ratio',
    save_folder=calibrate_linedata_folder
)
#Check OK: the vortices are doing quite well, the ring is too low SNR and the model components change quite drastically, giving strong difference in ratio
"""
#Do continuum subtraction (need CASA 6.2.1-7)
SBLB_no_ave_selfcal = f'{prefix}_SBLB_no_ave_selfcal_time_ave.ms'
contsub_vis = f'{SBLB_no_ave_selfcal}.contsub'
os.system(f'rm -rf {contsub_vis}*')
uvcontsub(
    vis=SBLB_no_ave_selfcal,spw=complete_dataset_dict['cont_spws'],fitspw=fitspw,
    excludechans=True,solint='int',fitorder=1,want_cont=False
)
listobs(vis=SBLB_no_ave_selfcal+'.contsub',listfile=SBLB_no_ave_selfcal+'.contsub.listobs.txt',overwrite=True)
#2025-01-05 07:33:16     WARN    calibrater::setvi(bool,bool)    Forcing use of OLD VisibilityIterator.
#2025-01-05 09:04:50     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [233.22, 235.205] (GHz).
#2025-01-05 09:04:50     WARN    VBContinuumSubtractor::apply+   The frequency range used for the continuum fit was [233.236, 235.205] (GHz).
#2025-01-05 09:04:50     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [216.528, 218.512] (GHz).
#2025-01-05 09:04:50     WARN    VBContinuumSubtractor::apply+   The frequency range used for the continuum fit was [216.528, 218.449] (GHz).
#2025-01-05 09:04:50     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [230.541, 232.415] (GHz).
#2025-01-05 09:04:50     WARN    VBContinuumSubtractor::apply+   The frequency range used for the continuum fit was [230.558, 232.415] (GHz).
#2025-01-05 09:04:51     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [218.541, 220.415] (GHz).
#2025-01-05 09:04:51     WARN    VBContinuumSubtractor::apply+   The frequency range used for the continuum fit was [218.542, 220.394] (GHz).
#2025-01-05 09:06:24     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [233.22, 235.204] (GHz).
#2025-01-05 09:06:24     WARN    VBContinuumSubtractor::apply+   The frequency range used for the continuum fit was [233.236, 235.204] (GHz).
#2025-01-05 09:06:24     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [216.527, 218.512] (GHz).
#2025-01-05 09:06:24     WARN    VBContinuumSubtractor::apply+   The frequency range used for the continuum fit was [216.527, 218.449] (GHz).
#2025-01-05 09:06:24     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [230.541, 232.415] (GHz).
#2025-01-05 09:06:24     WARN    VBContinuumSubtractor::apply+   The frequency range used for the continuum fit was [230.558, 232.415] (GHz).
#2025-01-05 09:06:25     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [218.541, 220.415] (GHz).
#2025-01-05 09:06:25     WARN    VBContinuumSubtractor::apply+   The frequency range used for the continuum fit was [218.542, 220.394] (GHz).
#2025-01-05 09:08:06     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [233.22, 235.205] (GHz).
#2025-01-05 09:08:06     WARN    VBContinuumSubtractor::apply+   The frequency range used for the continuum fit was [233.236, 235.205] (GHz).
#2025-01-05 09:08:06     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [216.528, 218.512] (GHz).
#2025-01-05 09:08:06     WARN    VBContinuumSubtractor::apply+   The frequency range used for the continuum fit was [216.528, 218.449] (GHz).
#2025-01-05 09:08:07     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [230.541, 232.415] (GHz).
#2025-01-05 09:08:07     WARN    VBContinuumSubtractor::apply+   The frequency range used for the continuum fit was [230.558, 232.415] (GHz).
#2025-01-05 09:08:07     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [218.541, 220.415] (GHz).
#2025-01-05 09:08:07     WARN    VBContinuumSubtractor::apply+   The frequency range used for the continuum fit was [218.542, 220.394] (GHz).
#2025-01-05 09:09:52     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [233.219, 235.204] (GHz).
#2025-01-05 09:09:52     WARN    VBContinuumSubtractor::apply+   The frequency range used for the continuum fit was [233.235, 235.204] (GHz).
#2025-01-05 09:09:52     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [216.527, 218.511] (GHz).
#2025-01-05 09:09:52     WARN    VBContinuumSubtractor::apply+   The frequency range used for the continuum fit was [216.527, 218.449] (GHz).
#2025-01-05 09:09:52     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [230.54, 232.414] (GHz).
#2025-01-05 09:09:52     WARN    VBContinuumSubtractor::apply+   The frequency range used for the continuum fit was [230.558, 232.414] (GHz).
#2025-01-05 09:09:53     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [218.54, 220.414] (GHz).
#2025-01-05 09:09:53     WARN    VBContinuumSubtractor::apply+   The frequency range used for the continuum fit was [218.541, 220.394] (GHz).
#2025-01-05 09:11:36     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [233.219, 235.204] (GHz).
#2025-01-05 09:11:36     WARN    VBContinuumSubtractor::apply+   The frequency range used for the continuum fit was [233.235, 235.204] (GHz).
#2025-01-05 09:11:36     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [216.527, 218.511] (GHz).
#2025-01-05 09:11:36     WARN    VBContinuumSubtractor::apply+   The frequency range used for the continuum fit was [216.527, 218.448] (GHz).
#2025-01-05 09:11:37     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [230.54, 232.414] (GHz).
#2025-01-05 09:11:37     WARN    VBContinuumSubtractor::apply+   The frequency range used for the continuum fit was [230.557, 232.414] (GHz).
#2025-01-05 09:11:38     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [218.54, 220.414] (GHz).
#2025-01-05 09:11:38     WARN    VBContinuumSubtractor::apply+   The frequency range used for the continuum fit was [218.541, 220.393] (GHz).
#2025-01-05 09:13:30     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [233.202, 235.186] (GHz).
#2025-01-05 09:13:30     WARN    VBContinuumSubtractor::apply+   The frequency range used for the continuum fit was [233.217, 235.186] (GHz).
#2025-01-05 09:13:31     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [216.51, 218.494] (GHz).
#2025-01-05 09:13:31     WARN    VBContinuumSubtractor::apply+   The frequency range used for the continuum fit was [216.51, 218.431] (GHz).
#2025-01-05 09:13:32     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [230.523, 232.397] (GHz).
#2025-01-05 09:13:32     WARN    VBContinuumSubtractor::apply+   The frequency range used for the continuum fit was [230.54, 232.397] (GHz).
#2025-01-05 09:13:37     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [218.523, 220.397] (GHz).
#2025-01-05 09:13:37     WARN    VBContinuumSubtractor::apply+   The frequency range used for the continuum fit was [218.524, 220.376] (GHz).
#2025-01-05 09:14:47     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [233.198, 235.182] (GHz).
#2025-01-05 09:14:47     WARN    VBContinuumSubtractor::apply+   The frequency range used for the continuum fit was [233.214, 235.182] (GHz).
#2025-01-05 09:14:48     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [216.506, 218.491] (GHz).
#2025-01-05 09:14:48     WARN    VBContinuumSubtractor::apply+   The frequency range used for the continuum fit was [216.506, 218.428] (GHz).
#2025-01-05 09:14:49     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [230.519, 232.393] (GHz).
#2025-01-05 09:14:49     WARN    VBContinuumSubtractor::apply+   The frequency range used for the continuum fit was [230.537, 232.393] (GHz).
#2025-01-05 09:14:54     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [218.519, 220.393] (GHz).
#2025-01-05 09:14:54     WARN    VBContinuumSubtractor::apply+   The frequency range used for the continuum fit was [218.52, 220.373] (GHz).
#2025-01-05 09:16:15     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [233.194, 235.178] (GHz).
#2025-01-05 09:16:15     WARN    VBContinuumSubtractor::apply+   The frequency range used for the continuum fit was [233.21, 235.178] (GHz).
#2025-01-05 09:16:15     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [216.502, 218.487] (GHz).
#2025-01-05 09:16:15     WARN    VBContinuumSubtractor::apply+   The frequency range used for the continuum fit was [216.502, 218.424] (GHz).
#2025-01-05 09:16:16     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [230.516, 232.39] (GHz).
#2025-01-05 09:16:16     WARN    VBContinuumSubtractor::apply+   The frequency range used for the continuum fit was [230.533, 232.39] (GHz).
#2025-01-05 09:16:23     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [218.516, 220.39] (GHz).
#2025-01-05 09:16:23     WARN    VBContinuumSubtractor::apply+   The frequency range used for the continuum fit was [218.517, 220.369] (GHz).
#2025-01-05 09:17:27     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [233.194, 235.178] (GHz).
#2025-01-05 09:17:27     WARN    VBContinuumSubtractor::apply+   The frequency range used for the continuum fit was [233.209, 235.178] (GHz).
#2025-01-05 09:17:28     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [216.502, 218.486] (GHz).
#2025-01-05 09:17:28     WARN    VBContinuumSubtractor::apply+   The frequency range used for the continuum fit was [216.502, 218.424] (GHz).
#2025-01-05 09:17:29     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [230.515, 232.389] (GHz).
#2025-01-05 09:17:29     WARN    VBContinuumSubtractor::apply+   The frequency range used for the continuum fit was [230.533, 232.389] (GHz).
#2025-01-05 09:17:34     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [218.515, 220.389] (GHz).
#2025-01-05 09:17:34     WARN    VBContinuumSubtractor::apply+   The frequency range used for the continuum fit was [218.516, 220.369] (GHz).
#....10....20....30....40....50....60....70....80....90....100%
"""
#Split final ms table into separate spws for each targeted line (added SiS, SiO, double velocity range coverage, just in case, but check that no overlap occurs)
rest_freq_12CO         = 230.5380000e9 #J=2-1
    #2:0~29,6:0~29,10:0~29,14:0~29,18:0~29,22:0~29,26:0~29,30:0~29,34:0~29
rest_freq_13CO         = 220.3986842e9 #J=2-1 (ignoring splitting)
    #3:0~32,7:0~32,11:0~32,15:0~32,19:0~32,23:0~32,27:0~32,31:0~32,35:0~32
rest_freq_C18O         = 219.5603541e9 #J=2-1
    #3:845~890,7:845~890,11:845~890,15:845~890,19:845~890,23:845~890,27:845~890,31:845~890,35:845~890
rest_freq_H2CO_303_202 = 218.2221920e9
    #1:17~19,5:17~19,9:17~19,13:17~19,17:17~19,21:17~19,25:17~19,29:17~19,33:17~19
rest_freq_H2CO_321_220 = 218.7600660e9
    #3:1665~1709,7:1665~1709,11:1665~1709,15:1665~1709,19:1665~1709,23:1665~1709,27:1665~1709,31:1665~1709,35:1665~1709
rest_freq_H2CO_322_221 = 218.4756320e9
    #flagged
rest_freq_H2CO_918_919 = 216.5686510e9
    #flagged
rest_freq_H2S_220_211  = 216.7104365e9
    #1:113~116,5:113~116,9:113~116,13:113~116,17:113~116,21:113~116,25:113~116,29:113~116,33:113~116
rest_freq_SiO          = 217.1049800e9 #J=5-4
    #1:88~91,5:88~91,9:88~91,13:88~91,17:88~91,21:88~91,25:88~91,29:88~91,33:88~91
rest_freq_DCN_F21      = 217.2384000e9 #J=3-2
rest_freq_DCN_F22      = 217.2386307e9 #J=3-2
    #1:80~82,5:80~82,9:80~82,13:80~82,17:80~82,21:80~82,25:80~82,29:80~82,33:80~82
rest_freq_SiS_1211     = 217.8176630e9 #J=12-11
    #1:43~45,5:43~45,9:43~45,13:43~45,17:43~45,21:43~45,25:43~45,29:43~45,33:43~45
rest_freq_SO           = 219.9494420e9 #3Sigma 6(5)-5(4)
    #3:446~492,7:446~492,11:446~492,15:446~492,19:446~492,23:446~492,27:446~492,31:446~492,35:446~492

#SBs and LBs
#0: 233.35, 235.10 continuum
#1: 216.60, 218.40 rest_freq_H2CO_918_919,rest_freq_DCN_F21,rest_freq_DCN_F22,rest_freq_SiS_1211,rest_freq_H2CO_303_202,rest_freq_H2CO_322_221
#2: 230.50, 232.40 rest_freq_12CO
#3: 218.50, 220.40 rest_freq_H2CO_322_221,rest_freq_H2CO_321_220,rest_freq_C18O,rest_freq_SO,rest_freq_13CO

data_params_LB = {
    f'LB{i}': {
        'vis':          f'{prefix}_LB_EB{i}.ms',
        'name':         f'LB_EB{i}',
        'field':        'MWC758',
        'line_spws':    np.array([1,1,1,1,1,1,1,1, 2, 3,3,3,3,3]), #list of spws containing lines
        'line_freqs':   np.array([
            rest_freq_H2CO_918_919,rest_freq_H2S_220_211,rest_freq_SiO,rest_freq_DCN_F21,rest_freq_DCN_F22,rest_freq_SiS_1211,rest_freq_H2CO_303_202,rest_freq_H2CO_322_221,
            rest_freq_12CO,
            rest_freq_H2CO_322_221,rest_freq_H2CO_321_220,rest_freq_C18O,rest_freq_SO,rest_freq_13CO
        ]), #frequencies (Hz) corresponding to line_spws
        'spwcont_forplot': '0',
        'cont_spws':       '0,1,2,3',
        'width_array':     [2,2,32,32],
    } for i in range(number_of_EBs['LB']) #Phase RMS 13.935, 15.003, 21.429, 32.285, 34.154
}

data_params_SB = {
    f'SB{i}': {
        'vis':          f'{prefix}_SB_EB{i}.ms',
        'name':         f'SB_EB{i}',
        'field':        'MWC758',
        'line_spws':    np.array([1,1,1,1,1,1,1,1, 2, 3,3,3,3,3]), #list of spws containing lines
        'line_freqs':   np.array([
            rest_freq_H2CO_918_919,rest_freq_H2S_220_211,rest_freq_SiO,rest_freq_DCN_F21,rest_freq_DCN_F22,rest_freq_SiS_1211,rest_freq_H2CO_303_202,rest_freq_H2CO_322_221,
            rest_freq_12CO,
            rest_freq_H2CO_322_221,rest_freq_H2CO_321_220,rest_freq_C18O,rest_freq_SO,rest_freq_13CO
        ]), #frequencies (Hz) corresponding to line_spws
        'spwcont_forplot': '0',
        'cont_spws':       '0,1,2,3',
        'width_array':     [8,8,192,192],
    } for i in range(number_of_EBs['SB']) #Phase RMS 38.807, 20.869, 24.399, 16.723
}

data_params = data_params_LB.copy()
data_params.update(data_params_SB)

for params in data_params.values():
    flagchannels_string = get_flagchannels(ms_dict=params,output_prefix=prefix,velocity_range=np.array([-30.,30.]) + v_sys)

# Flagchannels input string for LB_EB0: '1:123~125, 1:80~82, 1:80~82, 1:43~45, 1:17~19, 1:0~3, 2:0~29, 3:1919~1919, 3:1665~1709, 3:845~890, 3:446~492, 3:0~32'
# Flagchannels input string for LB_EB1: '1:123~125, 1:80~82, 1:80~82, 1:43~45, 1:17~19, 1:0~3, 2:0~29, 3:1919~1919, 3:1665~1709, 3:845~890, 3:446~492, 3:0~32'
# Flagchannels input string for LB_EB2: '1:123~125, 1:80~82, 1:80~82, 1:43~45, 1:17~19, 1:0~3, 2:0~29, 3:1919~1919, 3:1665~1709, 3:845~890, 3:446~492, 3:0~32'
# Flagchannels input string for LB_EB3: '1:123~125, 1:80~82, 1:80~82, 1:43~45, 1:17~19, 1:0~3, 2:0~29, 3:1919~1919, 3:1665~1709, 3:845~890, 3:446~492, 3:0~32'
# Flagchannels input string for LB_EB4: '1:123~125, 1:80~82, 1:80~82, 1:43~45, 1:17~19, 1:0~3, 2:0~29, 3:1919~1919, 3:1665~1709, 3:845~890, 3:446~492, 3:0~32'
# Flagchannels input string for SB_EB0: '1:122~125, 1:80~82, 1:80~82, 1:43~45, 1:17~19, 1:0~3, 2:0~28, 3:1919~1919, 3:1664~1709, 3:844~889, 3:446~491, 3:0~31'
# Flagchannels input string for SB_EB1: '1:122~125, 1:80~82, 1:80~82, 1:43~45, 1:17~19, 1:0~3, 2:0~28, 3:1919~1919, 3:1664~1709, 3:844~889, 3:446~491, 3:0~31'
# Flagchannels input string for SB_EB2: '1:122~125, 1:80~82, 1:80~82, 1:43~45, 1:17~19, 1:0~3, 2:0~28, 3:1919~1919, 3:1664~1709, 3:844~889, 3:446~491, 3:0~31'
# Flagchannels input string for SB_EB3: '1:122~125, 1:80~82, 1:80~82, 1:43~45, 1:17~19, 1:0~3, 2:0~28, 3:1919~1919, 3:1664~1709, 3:844~889, 3:446~491, 3:0~31'

# rest_freq_H2CO_918_919 1:123~125,5:123~125,9:123~125,13:123~125,17:123~125,21:123~125,25:123~125,29:123~125,33:123~125
# rest_freq_DCN_F21      1:80~82,5:80~82,9:80~82,13:80~82,17:80~82,21:80~82,25:80~82,29:80~82,33:80~82
# rest_freq_DCN_F22      N/A
# rest_freq_SiS_1211     1:43~45,5:43~45,9:43~45,13:43~45,17:43~45,21:43~45,25:43~45,29:43~45,33:43~45
# rest_freq_H2CO_303_202 1:17~19,5:17~19,9:17~19,13:17~19,17:17~19,21:17~19,25:17~19,29:17~19,33:17~19
# rest_freq_H2CO_322_221 1:0~3,5:0~3,9:0~3,13:0~3,17:0~3,21:0~3,25:0~3,29:0~3,33:0~3

# rest_freq_12CO         2:0~29,6:0~29,10:0~29,14:0~29,18:0~29,22:0~29,26:0~29,30:0~29,34:0~29
            
# rest_freq_H2CO_322_221 3:1919~1919,7:1919~1919,11:1919~1919,15:1919~1919,19:1919~1919,23:1919~1919,27:1919~1919,31:1919~1919,35:1919~1919
# rest_freq_H2CO_321_220 3:1665~1709,7:1665~1709,11:1665~1709,15:1665~1709,19:1665~1709,23:1665~1709,27:1665~1709,31:1665~1709,35:1665~1709
# rest_freq_C18O         3:845~890,7:845~890,11:845~890,15:845~890,19:845~890,23:845~890,27:845~890,31:845~890,35:845~890
# rest_freq_SO           3:446~492,7:446~492,11:446~492,15:446~492,19:446~492,23:446~492,27:446~492,31:446~492,35:446~492
# rest_freq_13CO         3:0~32,7:0~32,11:0~32,15:0~32,19:0~32,23:0~32,27:0~32,31:0~32,35:0~32

SBLB_no_ave_selfcal = f'{prefix}_SBLB_no_ave_selfcal_time_ave.ms'
contsub_vis = f'{SBLB_no_ave_selfcal}.contsub'

#12CO
vis_12CO = SBLB_no_ave_selfcal[:-3]+'_12CO.ms'
os.system(f'rm -rf {vis_12CO}*')
spw_12CO = '2:0~29,6:0~29,10:0~29,14:0~29,18:0~29,22:0~29,26:0~29,30:0~29,34:0~29'

split(vis=SBLB_no_ave_selfcal,outputvis=vis_12CO,spw=spw_12CO,datacolumn='data',keepflags=False)
listobs(vis=vis_12CO,listfile=vis_12CO+'.listobs.txt',overwrite=True)

split(vis=contsub_vis,outputvis=f'{vis_12CO}.contsub',spw=spw_12CO,datacolumn='data',keepflags=False)
listobs(vis=vis_12CO+'.contsub',listfile=vis_12CO+'.contsub.listobs.txt',overwrite=True)

#13CO
vis_13CO = SBLB_no_ave_selfcal[:-3]+'_13CO.ms'
os.system(f'rm -rf {vis_13CO}*')
spw_13CO = '3:0~32,7:0~32,11:0~32,15:0~32,19:0~32,23:0~32,27:0~32,31:0~32,35:0~32'

split(vis=SBLB_no_ave_selfcal,outputvis=vis_13CO,spw=spw_13CO,datacolumn='data',keepflags=False)
listobs(vis=vis_13CO,listfile=vis_13CO+'.listobs.txt',overwrite=True)

split(vis=contsub_vis,outputvis=f'{vis_13CO}.contsub',spw=spw_13CO,datacolumn='data',keepflags=False)
listobs(vis=vis_13CO+'.contsub',listfile=vis_13CO+'.contsub.listobs.txt',overwrite=True)

#C18O
vis_C18O = SBLB_no_ave_selfcal[:-3]+'_C18O.ms'
os.system(f'rm -rf {vis_C18O}*')
spw_C18O = '3:845~890,7:845~890,11:845~890,15:845~890,19:845~890,23:845~890,27:845~890,31:845~890,35:845~890'

split(vis=SBLB_no_ave_selfcal,outputvis=vis_C18O,spw=spw_C18O,datacolumn='data',keepflags=False)
listobs(vis=vis_C18O,listfile=vis_C18O+'.listobs.txt',overwrite=True)

split(vis=contsub_vis,outputvis=f'{vis_C18O}.contsub',spw=spw_C18O,datacolumn='data',keepflags=False)
listobs(vis=vis_C18O+'.contsub',listfile=vis_C18O+'.contsub.listobs.txt',overwrite=True)

#H2CO_303_202
vis_H2CO_303_202 = SBLB_no_ave_selfcal[:-3]+'_H2CO_303_202.ms'
os.system(f'rm -rf {vis_H2CO_303_202}*')
spw_H2CO_303_202 = '1:17~19,5:17~19,9:17~19,13:17~19,17:17~19,21:17~19,25:17~19,29:17~19,33:17~19'

split(vis=SBLB_no_ave_selfcal,outputvis=vis_H2CO_303_202,spw=spw_H2CO_303_202,datacolumn='data',keepflags=False)
listobs(vis=vis_H2CO_303_202,listfile=vis_H2CO_303_202+'.listobs.txt',overwrite=True)

split(vis=contsub_vis,outputvis=f'{vis_H2CO_303_202}.contsub',spw=spw_H2CO_303_202,datacolumn='data',keepflags=False)
listobs(vis=vis_H2CO_303_202+'.contsub',listfile=vis_H2CO_303_202+'.contsub.listobs.txt',overwrite=True)

#H2CO_321_220
vis_H2CO_321_220 = SBLB_no_ave_selfcal[:-3]+'_H2CO_321_220.ms'
os.system(f'rm -rf {vis_H2CO_321_220}*')
spw_H2CO_321_220 = '3:1665~1709,7:1665~1709,11:1665~1709,15:1665~1709,19:1665~1709,23:1665~1709,27:1665~1709,31:1665~1709,35:1665~1709'

split(vis=SBLB_no_ave_selfcal,outputvis=vis_H2CO_321_220,spw=spw_H2CO_321_220,datacolumn='data',keepflags=False)
listobs(vis=vis_H2CO_321_220,listfile=vis_H2CO_321_220+'.listobs.txt',overwrite=True)

split(vis=contsub_vis,outputvis=f'{vis_H2CO_321_220}.contsub',spw=spw_H2CO_321_220,datacolumn='data',keepflags=False)
listobs(vis=vis_H2CO_321_220+'.contsub',listfile=vis_H2CO_321_220+'.contsub.listobs.txt',overwrite=True)
"""
#H2CO_322_221 essentially empty because channels were flagged by pipeline crosscal
vis_H2CO_322_221 = SBLB_no_ave_selfcal[:-3]+'_H2CO_322_221.ms'
os.system(f'rm -rf {vis_H2CO_322_221}*')
spw_H2CO_322_221 = '1:0~3,5:0~3,9:0~3,13:0~3,17:0~3,21:0~3,25:0~3,29:0~3,33:0~3,3:1919~1919,7:1919~1919,11:1919~1919,15:1919~1919,19:1919~1919,23:1919~1919,27:1919~1919,31:1919~1919,35:1919~1919'

split(vis=SBLB_no_ave_selfcal,outputvis=vis_H2CO_322_221,spw=spw_H2CO_322_221,datacolumn='data',keepflags=False)
listobs(vis=vis_H2CO_322_221,listfile=vis_H2CO_322_221+'.listobs.txt',overwrite=True)

split(vis=contsub_vis,outputvis=f'{vis_H2CO_322_221}.contsub',spw=spw_H2CO_322_221,datacolumn='data',keepflags=False)
listobs(vis=vis_H2CO_322_221+'.contsub',listfile=vis_H2CO_322_221+'.contsub.listobs.txt',overwrite=True)

#H2CO_918_919 essentially empty because channels were flagged by pipeline crosscal
vis_H2CO_918_919 = SBLB_no_ave_selfcal[:-3]+'_H2CO_918_919.ms'
os.system(f'rm -rf {vis_H2CO_918_919}*')
spw_H2CO_918_919 = '1:123~125,5:123~125,9:123~125,13:123~125,17:123~125,21:123~125,25:123~125,29:123~125,33:123~125'

split(vis=SBLB_no_ave_selfcal,outputvis=vis_H2CO_918_919,spw=spw_H2CO_918_919,datacolumn='data',keepflags=False)
listobs(vis=vis_H2CO_918_919,listfile=vis_H2CO_918_919+'.listobs.txt',overwrite=True)

split(vis=contsub_vis,outputvis=f'{vis_H2CO_918_919}.contsub',spw=spw_H2CO_918_919,datacolumn='data',keepflags=False)
listobs(vis=vis_H2CO_918_919+'.contsub',listfile=vis_H2CO_918_919+'.contsub.listobs.txt',overwrite=True)
"""
#H2S
vis_H2S = SBLB_no_ave_selfcal[:-3]+'_H2S.ms'
os.system(f'rm -rf {vis_H2S}*')
spw_H2S = '1:113~116,5:113~116,9:113~116,13:113~116,17:113~116,21:113~116,25:113~116,29:113~116,33:113~116'

split(vis=SBLB_no_ave_selfcal,outputvis=vis_H2S,spw=spw_H2S,datacolumn='data',keepflags=False)
listobs(vis=vis_H2S,listfile=vis_H2S+'.listobs.txt',overwrite=True)

split(vis=contsub_vis,outputvis=f'{vis_H2S}.contsub',spw=spw_H2S,datacolumn='data',keepflags=False)
listobs(vis=vis_H2S+'.contsub',listfile=vis_H2S+'.contsub.listobs.txt',overwrite=True)

#SiO
vis_SiO = SBLB_no_ave_selfcal[:-3]+'_SiO.ms'
os.system(f'rm -rf {vis_SiO}*')
spw_SiO = '1:88~91,5:88~91,9:88~91,13:88~91,17:88~91,21:88~91,25:88~91,29:88~91,33:88~91'

split(vis=SBLB_no_ave_selfcal,outputvis=vis_SiO,spw=spw_SiO,datacolumn='data',keepflags=False)
listobs(vis=vis_SiO,listfile=vis_SiO+'.listobs.txt',overwrite=True)

split(vis=contsub_vis,outputvis=f'{vis_SiO}.contsub',spw=spw_SiO,datacolumn='data',keepflags=False)
listobs(vis=vis_SiO+'.contsub',listfile=vis_SiO+'.contsub.listobs.txt',overwrite=True)

#DCN
vis_DCN = SBLB_no_ave_selfcal[:-3]+'_DCN.ms'
os.system(f'rm -rf {vis_DCN}*')
spw_DCN = '1:80~82,5:80~82,9:80~82,13:80~82,17:80~82,21:80~82,25:80~82,29:80~82,33:80~82'

split(vis=SBLB_no_ave_selfcal,outputvis=vis_DCN,spw=spw_DCN,datacolumn='data',keepflags=False)
listobs(vis=vis_DCN,listfile=vis_DCN+'.listobs.txt',overwrite=True)

split(vis=contsub_vis,outputvis=f'{vis_DCN}.contsub',spw=spw_DCN,datacolumn='data',keepflags=False)
listobs(vis=vis_DCN+'.contsub',listfile=vis_DCN+'.contsub.listobs.txt',overwrite=True)

#SiS
vis_SiS = SBLB_no_ave_selfcal[:-3]+'_SiS.ms'
os.system(f'rm -rf {vis_SiS}*')
spw_SiS = '1:43~45,5:43~45,9:43~45,13:43~45,17:43~45,21:43~45,25:43~45,29:43~45,33:43~45'

split(vis=SBLB_no_ave_selfcal,outputvis=vis_SiS,spw=spw_SiS,datacolumn='data',keepflags=False)
listobs(vis=vis_SiS,listfile=vis_SiS+'.listobs.txt',overwrite=True)

split(vis=contsub_vis,outputvis=f'{vis_SiS}.contsub',spw=spw_SiS,datacolumn='data',keepflags=False)
listobs(vis=vis_SiS+'.contsub',listfile=vis_SiS+'.contsub.listobs.txt',overwrite=True)

#SO
vis_SO = SBLB_no_ave_selfcal[:-3]+'_SO.ms'
os.system(f'rm -rf {vis_SO}*')
spw_SO = '3:446~492,7:446~492,11:446~492,15:446~492,19:446~492,23:446~492,27:446~492,31:446~492,35:446~492'

split(vis=SBLB_no_ave_selfcal,outputvis=vis_SO,spw=spw_SO,datacolumn='data',keepflags=False)
listobs(vis=vis_SO,listfile=vis_SO+'.listobs.txt',overwrite=True)

split(vis=contsub_vis,outputvis=f'{vis_SO}.contsub',spw=spw_SO,datacolumn='data',keepflags=False)
listobs(vis=vis_SO+'.contsub',listfile=vis_SO+'.contsub.listobs.txt',overwrite=True)
"""
for _vis,_freq in zip(
    [vis_12CO,vis_13CO,vis_C18O,vis_H2CO_303_202,vis_H2CO_321_220,vis_DCN,vis_SiS,vis_SO],#vis_H2CO_322_221,vis_H2CO_918_919,
    [rest_freq_12CO,rest_freq_13CO,rest_freq_C18O,rest_freq_H2CO_303_202,rest_freq_H2CO_321_220,rest_freq_DCN_F21,rest_freq_SiS_1211,rest_freq_SO]#rest_freq_H2CO_322_221,rest_freq_H2CO_918_919,
):
    for _suffix in ['','.contsub']:
        plotfile = os.path.join(calibrate_linedata_folder,_vis[:-3]+'_freq-v-amp_cal'+_suffix+'.png')
        plotms(
            vis         = _vis + _suffix,
            xaxis       = 'frequency',
            restfreq    = '{}Hz'.format(_freq),
            freqframe   = 'LSRK',
            yaxis       = 'amplitude',
            field       = 'CQ_Tau',
            ydatacolumn = 'data',
            avgtime     = '1e8',
            avgscan     = True,
            avgspw      = True,
            avgbaseline = False,
            showgui     = False,
            exprange    = 'all',
            plotfile    = plotfile,
            overwrite   = True,
        )

        plotfile = os.path.join(calibrate_linedata_folder,_vis[:-3]+'_vel-v-amp_cal'+_suffix+'.png')
        plotms(
            vis         = _vis + _suffix,
            xaxis       = 'velocity',
            restfreq    = '{}Hz'.format(_freq),
            freqframe   = 'LSRK',
            yaxis       = 'amplitude',
            field       = 'CQ_Tau',
            ydatacolumn = 'data',
            avgtime     = '1e8',
            avgscan     = True,
            iteraxis    = 'spw',
                #spw         = '0',
                #avgspw      = True,
            avgbaseline = False,
            showgui     = False,
            exprange    = 'all',
            plotfile    = plotfile,
            overwrite   = True,
        )
"""
#Define line mask and noise annulus
mask_pa        = PA
mask_semimajor = 2.0
mask_semiminor = 2.0
mask_ra  = '05h30m27.535947s'
mask_dec = '25d19m56.615500s'

line_mask = f'ellipse[[{mask_ra},{mask_dec}], [{mask_semimajor:.3f}arcsec, {mask_semiminor:.3f}arcsec], {mask_pa:.1f}deg]'

noise_annulus_line = f"annulus[[{mask_ra}, {mask_dec}],['4.arcsec', '6.arcsec']]"

chanstart = '-6.8km/s'
chanwidth = ' 1.3km/s'
nchan     = 15

imagename = vis_12CO[:-3]+'.contsub_image_1.0robust_4.0sigma'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.image.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis=vis_12CO+'.contsub',imagename=imagename,
    specmode='cube',restoringbeam='common',nterms=1,
    deconvolver='multiscale',scales=[0,2,4,8],
    cell='0.006arcsec',imsize=2400,gain=0.1,niter=50000,
    weighting='briggs',robust=1.0,threshold='1.6440mJy',interactive=False, #~4sigma
    start=chanstart,width=chanwidth,nchan=nchan,
    outframe='LSRK',veltype='radio',restfreq='{}Hz'.format(rest_freq_12CO),
    usemask='user',mask=line_mask, 
)
estimate_SNR(imagename+'.image',disk_mask=line_mask,noise_mask=noise_annulus_line)
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#MWC_758_SBLB_no_ave_selfcal_time_ave_12CO.contsub_image_1.0robust_4.0sigma.image
#Beam 0.064 arcsec x 0.044 arcsec (0.13 deg)
#Flux inside disk mask: 13180.53 mJy
#Peak intensity of source: 12.23 mJy/beam
#rms: 4.11e-01 mJy/beam
#Peak SNR: 29.79

chanstart = '-6.0km/s'
chanwidth = ' 1.4km/s'
nchan     = 15

imagename = vis_13CO[:-3]+'.contsub_image_1.0robust_4.0sigma'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.image.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis=vis_13CO+'.contsub',imagename=imagename,
    specmode='cube',restoringbeam='common',nterms=1,
    deconvolver='multiscale',scales=[0,2,4,8],
    cell='0.006arcsec',imsize=1200,gain=0.1,niter=50000,
    weighting='briggs',robust=1.0,threshold='1.9040mJy',interactive=False, #~4sigma
    start=chanstart,width=chanwidth,nchan=nchan,
    outframe='LSRK',veltype='radio',restfreq='{}Hz'.format(rest_freq_13CO),
    usemask='user',mask=line_mask,
    #uvtaper='0.1arcsec',#interactive=True, 
)
estimate_SNR(imagename+'.image',disk_mask=line_mask,noise_mask=noise_annulus_line)
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#MWC_758_SBLB_no_ave_selfcal_time_ave_13CO.contsub_image_1.0robust_4.0sigma.image
#Beam 0.068 arcsec x 0.044 arcsec (0.15 deg)
#Flux inside disk mask: 4100.14 mJy
#Peak intensity of source: 6.93 mJy/beam
#rms: 4.76e-01 mJy/beam
#Peak SNR: 14.57

imagename = vis_C18O[:-3]+'.contsub_image_1.0robust_4.0sigma'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.image.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis=vis_C18O+'.contsub',imagename=imagename,
    specmode='cube',restoringbeam='common',nterms=1,
    deconvolver='multiscale',scales=[0,2,4,8,16],
    cell='0.006arcsec',imsize=1200,gain=0.1,niter=50000,
    weighting='briggs',robust=1.0,threshold='1.2800mJy',interactive=False, #~4sigma
    start=chanstart,width=chanwidth,nchan=nchan,
    outframe='LSRK',veltype='radio',restfreq='{}Hz'.format(rest_freq_C18O),
    usemask='user',mask=line_mask,
    #uvtaper='0.1arcsec',#interactive=True, 
)
estimate_SNR(imagename+'.image',disk_mask=line_mask,noise_mask=noise_annulus_line)
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#MWC_758_SBLB_no_ave_selfcal_time_ave_C18O.contsub_image_1.0robust_4.0sigma.image
#Beam 0.069 arcsec x 0.045 arcsec (0.15 deg)
#Flux inside disk mask: 1216.29 mJy
#Peak intensity of source: 4.25 mJy/beam
#rms: 3.20e-01 mJy/beam
#Peak SNR: 13.31

imagename = vis_SO[:-3]+'.contsub_image_natural_4.0sigma'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.image.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis=vis_SO+'.contsub',imagename=imagename,
    specmode='cube',restoringbeam='common',nterms=1,
    deconvolver='multiscale',scales=[0,2,4],
    cell='0.009arcsec',imsize=1200,gain=0.1,niter=50000,
    weighting='natural',threshold='1.4840mJy',interactive=False, #~4sigma 
    start=chanstart,width=chanwidth,nchan=nchan,
    outframe='LSRK',veltype='radio',restfreq='{}Hz'.format(rest_freq_SO),
    usemask='user',mask=line_mask,
    #uvtaper='0.1arcsec',#interactive=True, 
)
estimate_SNR(imagename+'.image',disk_mask=line_mask,noise_mask=noise_annulus_line)
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#MWC_758_SBLB_no_ave_selfcal_time_ave_SO.contsub_image_natural_4.0sigma.image
#Beam 0.091 arcsec x 0.072 arcsec (-18.61 deg)
#Flux inside disk mask: -29.52 mJy
#Peak intensity of source: 2.35 mJy/beam
#rms: 3.71e-01 mJy/beam
#Peak SNR: 6.32

imagename = vis_SO[:-3]+'.contsub_image_natural_4.0sigma_uvtaper0.1arcsec'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.image.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis=vis_SO+'.contsub',imagename=imagename,
    specmode='cube',restoringbeam='common',nterms=1,
    deconvolver='multiscale',scales=[0,2,4],
    cell='0.025arcsec',imsize=1200,gain=0.1,niter=50000,
    weighting='natural',threshold='2.0960mJy',interactive=False, #~4sigma 
    start=chanstart,width=chanwidth,nchan=nchan,
    outframe='LSRK',veltype='radio',restfreq='{}Hz'.format(rest_freq_SO),
    usemask='user',mask=line_mask,
    uvtaper='0.1arcsec',#interactive=True, 
)
estimate_SNR(imagename+'.image',disk_mask=line_mask,noise_mask=noise_annulus_line)
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#MWC_758_SBLB_no_ave_selfcal_time_ave_SO.contsub_image_natural_4.0sigma_uvtaper0.1arcsec.image
#Beam 0.247 arcsec x 0.201 arcsec (-32.34 deg)
#Flux inside disk mask: 0.82 mJy
#Peak intensity of source: 2.63 mJy/beam
#rms: 5.24e-01 mJy/beam
#Peak SNR: 5.01

imagename = vis_SO[:-3]+'.contsub_image_0.5robust_4.0sigma'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.image.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis=vis_SO+'.contsub',imagename=imagename,
    specmode='cube',restoringbeam='common',nterms=1,
    deconvolver='multiscale',scales=[0,2,4],
    cell='0.004arcsec',imsize=2400,gain=0.1,niter=50000,
    weighting='briggs',robust=0.5,threshold='1.6840mJy',interactive=False, #~4sigma 
    start=chanstart,width=chanwidth,nchan=nchan,
    outframe='LSRK',veltype='radio',restfreq='{}Hz'.format(rest_freq_SO),
    usemask='user',mask=line_mask,
    #uvtaper='0.1arcsec',#interactive=True, 
)
estimate_SNR(imagename+'.image',disk_mask=line_mask,noise_mask=noise_annulus_line)
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#MWC_758_SBLB_no_ave_selfcal_time_ave_SO.contsub_image_0.5robust_4.0sigma.image
#Beam 0.047 arcsec x 0.028 arcsec (0.16 deg)
#Flux inside disk mask: -0.80 mJy
#Peak intensity of source: 2.61 mJy/beam
#rms: 4.21e-01 mJy/beam
#Peak SNR: 6.20

imagename = vis_SO[:-3]+'.contsub_image_0.75robust_4.0sigma'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.image.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis=vis_SO+'.contsub',imagename=imagename,
    specmode='cube',restoringbeam='common',nterms=1,
    deconvolver='multiscale',scales=[0,2,4],
    cell='0.004arcsec',imsize=2400,gain=0.1,niter=50000,
    weighting='briggs',robust=0.75,threshold='1.5680mJy',interactive=False, #~4sigma 
    start=chanstart,width=chanwidth,nchan=nchan,
    outframe='LSRK',veltype='radio',restfreq='{}Hz'.format(rest_freq_SO),
    usemask='user',mask=line_mask,
    #uvtaper='0.1arcsec',#interactive=True, 
)
estimate_SNR(imagename+'.image',disk_mask=line_mask,noise_mask=noise_annulus_line)
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#MWC_758_SBLB_no_ave_selfcal_time_ave_SO.contsub_image_0.75robust_4.0sigma.image
#Beam 0.059 arcsec x 0.033 arcsec (0.22 deg)
#Flux inside disk mask: 9.64 mJy
#Peak intensity of source: 2.37 mJy/beam
#rms: 3.92e-01 mJy/beam
#Peak SNR: 6.04

imagename = vis_SO[:-3]+'.contsub_image_1.0robust_4.0sigma'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.image.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis=vis_SO+'.contsub',imagename=imagename,
    specmode='cube',restoringbeam='common',nterms=1,
    deconvolver='multiscale',scales=[0,2,4],
    cell='0.006arcsec',imsize=1200,gain=0.1,niter=50000,
    weighting='briggs',robust=1.0,threshold='1.5040mJy',interactive=False, #~4sigma 
    start=chanstart,width=chanwidth,nchan=nchan,
    outframe='LSRK',veltype='radio',restfreq='{}Hz'.format(rest_freq_SO),
    usemask='user',mask=line_mask,
    #uvtaper='0.1arcsec',#interactive=True, 
)
estimate_SNR(imagename+'.image',disk_mask=line_mask,noise_mask=noise_annulus_line)
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#MWC_758_SBLB_no_ave_selfcal_time_ave_SO.contsub_image_1.0robust_4.0sigma.image
#Beam 0.068 arcsec x 0.044 arcsec (8.38 deg)
#Flux inside disk mask: 24.47 mJy
#Peak intensity of source: 2.32 mJy/beam
#rms: 3.76e-01 mJy/beam
#Peak SNR: 6.17

imagename = vis_H2CO_321_220[:-3]+'.contsub_image_natural_4.0sigma'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.image.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis=vis_H2CO_321_220+'.contsub',imagename=imagename,
    specmode='cube',restoringbeam='common',nterms=1,
    deconvolver='multiscale',scales=[0,2,4],
    cell='0.009arcsec',imsize=1200,gain=0.1,niter=50000,
    weighting='natural',threshold='1.1200mJy',interactive=False, #~4sigma 
    start=chanstart,width=chanwidth,nchan=nchan,
    outframe='LSRK',veltype='radio',restfreq='{}Hz'.format(rest_freq_H2CO_321_220),
    usemask='user',mask=line_mask,
    #uvtaper='0.1arcsec',#interactive=True, 
)
estimate_SNR(imagename+'.image',disk_mask=line_mask,noise_mask=noise_annulus_line)
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#MWC_758_SBLB_no_ave_selfcal_time_ave_H2CO_321_220.contsub_image_natural_4.0sigma.image
#Beam 0.093 arcsec x 0.072 arcsec (-22.00 deg)
#Flux inside disk mask: -158.63 mJy
#Peak intensity of source: 1.57 mJy/beam
#rms: 2.80e-01 mJy/beam
#Peak SNR: 5.62

chanstart = '-15.00km/s'
chanwidth = ' 22.00km/s'
nchan     = 3

imagename = vis_SiS[:-3]+'.contsub_image_natural_4.0sigma'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.image.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis=vis_SiS+'.contsub',imagename=imagename,
    specmode='cube',restoringbeam='common',nterms=1,
    deconvolver='multiscale',scales=[0,2,4],
    cell='0.016arcsec',imsize=1200,gain=0.1,niter=50000,
    weighting='natural',threshold='1.5640mJy',interactive=False, #~4sigma
    start=chanstart,width=chanwidth,nchan=nchan,
    outframe='LSRK',veltype='radio',restfreq='{}Hz'.format(rest_freq_SiS_1211),
    usemask='user',mask=line_mask,
    #uvtaper='0.1arcsec',#interactive=True, 
)
estimate_SNR(imagename+'.image',disk_mask=line_mask,noise_mask=noise_annulus_line)
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#MWC_758_SBLB_no_ave_selfcal_time_ave_SiS.contsub_image_natural_4.0sigma.image
#Beam 0.285 arcsec x 0.218 arcsec (-33.60 deg)
#Flux inside disk mask: -195.70 mJy
#Peak intensity of source: 2.02 mJy/beam
#rms: 3.91e-01 mJy/beam
#Peak SNR: 5.16

chanstart = '-17.30km/s'
chanwidth = ' 22.00km/s'
nchan     = 3

imagename = vis_H2CO_303_202[:-3]+'.contsub_image_natural_4.0sigma'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.image.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis=vis_H2CO_303_202+'.contsub',imagename=imagename,
    specmode='cube',restoringbeam='common',nterms=1,
    deconvolver='multiscale',scales=[0,2,4],
    cell='0.016arcsec',imsize=1200,gain=0.1,niter=50000,
    weighting='natural',threshold='1.5440mJy',interactive=False, #~4sigma 
    start=chanstart,width=chanwidth,nchan=nchan,
    outframe='LSRK',veltype='radio',restfreq='{}Hz'.format(rest_freq_H2CO_303_202),
    usemask='user',mask=line_mask,
    #uvtaper='0.1arcsec',#interactive=True, 
)
estimate_SNR(imagename+'.image',disk_mask=line_mask,noise_mask=noise_annulus_line)
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#MWC_758_SBLB_no_ave_selfcal_time_ave_H2CO_303_202.contsub_image_natural_4.0sigma.image
#Beam 0.278 arcsec x 0.217 arcsec (-31.39 deg)
#Flux inside disk mask: 771.37 mJy
#Peak intensity of source: 2.51 mJy/beam
#rms: 3.86e-01 mJy/beam
#Peak SNR: 6.52

chanstart = '-16.50km/s'
chanwidth = ' 22.00km/s'
nchan     = 3

imagename = vis_DCN[:-3]+'.contsub_image_natural_4.0sigma'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.image.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis=vis_DCN+'.contsub',imagename=imagename,
    specmode='cube',restoringbeam='common',nterms=1,
    deconvolver='multiscale',scales=[0,2,4],
    cell='0.013arcsec',imsize=1200,gain=0.1,niter=50000,
    weighting='natural',threshold='1.5520mJy',interactive=False, #~4sigma 
    start=chanstart,width=chanwidth,nchan=nchan,
    outframe='LSRK',veltype='radio',restfreq='{}Hz'.format(rest_freq_DCN_F21),
    usemask='user',mask=line_mask,
    #uvtaper='0.1arcsec',#interactive=True, 
)
estimate_SNR(imagename+'.image',disk_mask=line_mask,noise_mask=noise_annulus_line)
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#MWC_758_SBLB_no_ave_selfcal_time_ave_DCN.contsub_image_natural_4.0sigma.image
#Beam 0.290 arcsec x 0.219 arcsec (-35.68 deg)
#Flux inside disk mask: -413.34 mJy
#Peak intensity of source: 1.90 mJy/beam
#rms: 3.88e-01 mJy/beam
#Peak SNR: 4.90

imagename = vis_SiO[:-3]+'.contsub_image_natural_4.0sigma'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.image.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis=vis_SiO+'.contsub',imagename=imagename,
    specmode='cube',restoringbeam='common',nterms=1,
    deconvolver='multiscale',scales=[0,2,4,8],
    cell='0.009arcsec',imsize=1200,gain=0.1,niter=50000,
    weighting='natural',threshold='0.2404mJy',interactive=False, #~4sigma 
    start=chanstart,width=chanwidth,nchan=nchan,
    outframe='LSRK',veltype='radio',restfreq='{}Hz'.format(rest_freq_SiO),
    usemask='user',mask=line_mask
)
estimate_SNR(imagename+'.image',disk_mask=line_mask,noise_mask=noise_annulus_line)
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#MWC_758_SBLB_no_ave_selfcal_time_ave_SiO.contsub_image_natural_4.0sigma.image
#Beam 0.085 arcsec x 0.072 arcsec (-10.93 deg)
#Flux inside disk mask: 463.20 mJy
#Peak intensity of source: 0.28 mJy/beam
#rms: 6.01e-02 mJy/beam
#Peak SNR: 4.62

imagename = vis_H2S[:-3]+'.contsub_image_natural_4.0sigma'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.image.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis=vis_H2S+'.contsub',imagename=imagename,
    specmode='cube',restoringbeam='common',nterms=1,
    deconvolver='multiscale',scales=[0,2,4,8],
    cell='0.009arcsec',imsize=1200,gain=0.1,niter=50000,
    weighting='natural',threshold='0.2416mJy',interactive=False, #~4sigma 
    start=chanstart,width=chanwidth,nchan=nchan,
    outframe='LSRK',veltype='radio',restfreq='{}Hz'.format(rest_freq_H2S_220_211),
    usemask='user',mask=line_mask
)
estimate_SNR(imagename+'.image',disk_mask=line_mask,noise_mask=noise_annulus_line)
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#MWC_758_SBLB_no_ave_selfcal_time_ave_H2S.contsub_image_natural_4.0sigma.image
#Beam 0.085 arcsec x 0.072 arcsec (-10.68 deg)
#Flux inside disk mask: -351.95 mJy
#Peak intensity of source: 0.26 mJy/beam
#rms: 6.04e-02 mJy/beam
#Peak SNR: 4.33