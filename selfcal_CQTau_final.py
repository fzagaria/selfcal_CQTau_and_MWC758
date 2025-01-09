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

prefix = 'CQ_Tau'

# System properties.
incl  = 35.  # deg, from Ubeira-Gabellini et al. 2019
PA    = 55.  # deg, from Ubeira-Gabellini et al. 2019
v_sys =  6.2 # km/s, from listobs

# Whether to run tclean in parallel or not.
use_parallel = False

data_folderpath = '/data/beegfs/astro-storage/groups/benisty/frzagaria/SO_detections/CQTau/CQTau_Band6_data/'

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

#2013.1.00498.S: ang. res. = 0.260", spec. res. = 0.634 km/s, cont. sens. = 0.0311 mJy/beam (SB0)
path_2013_1_00498_S = './2013.1.00498.S/science_goal.uid___A001_X13a_Xe7/group.uid___A001__X13a_Xe8/member.uid___A001__X13a_Xe9/calibrated/'
data_2013_1_00498_S = ['uid___A002_Xa916fc_X775f.ms']

#2016.A.00026.S: ang. res. = 0.113", spec. res. = 0.635 km/s, cont. sens. = 0.0230 mJy/beam (SB1)
path_2016_A_00026_S = './2016.A.00026.S/science_goal.uid___A001_X1232_X1d1/group.uid___A001_X1232_X1d2/member.uid___A001_X1232_X1d3/calibrated/'
data_2016_A_00026_S = ['uid___A002_Xc3173a_X2204.ms.split.cal']

#2017.1.01404.S: ang. res. = 0.056", spec. res. = 0.635 km/s, cont. sens. = 0.0149 mJy/beam (LB0,1)
path_2017_1_01404_S = './2017.1.01404.S/science_goal.uid___A001_X1284_X141/group.uid___A001_X1284_X142/member.uid___A001_X1284_X143/calibrated/'
data_2017_1_01404_S = ['uid___A002_Xc6ff69_X5d9e.ms.split.cal', 'uid___A002_Xc72427_X1354.ms.split.cal']

PL_calibrated_path = [path_2013_1_00498_S,path_2016_A_00026_S,path_2017_1_01404_S]

PL_calibrated_data = [data_2013_1_00498_S,data_2016_A_00026_S,data_2017_1_01404_S]

PL_calibrated_name = [['SB_EB0'],['SB_EB1'],['LB_EB0','LB_EB1']]

field_name    = ['CQ_Tau', 'CQ_Tau', 'CQ_Tau']

spw_name      = ['0,1,2,3,4,5,6,7', '0,1,2,3,4,5,6,7', '19,23,25,27,29,31,33,35']

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

rest_freq_H2CO_918_919 = 216.5686510e9 #Spw0 (216.63 to 218.62), only 2013.1.boh.S
rest_freq_DCN_F21      = 217.2384000e9 #J=3-2
rest_freq_DCN_F22      = 217.2386307e9 #J=3-2
# rest_freq_SiS_1211    = 217.8176630e9 #J=12-11 (checked later, presumably won't affect the selfcal anyway)
rest_freq_H2CO_303_202 = 218.2221920e9
rest_freq_H2CO_322_221 = 218.4756320e9
rest_freq_H2CO_321_220 = 218.7600660e9 #Spw1 (218.56 to 219.02)
rest_freq_C18O         = 219.5603541e9 #J=2-1, Spw3 (219.49 to 219.96)
rest_freq_SO           = 219.9494420e9 #3Sigma 6(5)-5(4) 
rest_freq_13CO         = 220.3986842e9 #J=2-1 (ignoring splitting), Spw4 (219.96 to 220.43) 
rest_freq_12CO         = 230.5380000e9 #J=2-1, Spw5 (230.49 to 230.96)
# continuum                            #Spw2 (219.03 to 219.49), Spw6 (230.96 to 231.43)

number_of_EBs  = {'SB':2,'LB':2}

#LBs
#0: 231.50, 233.25 continuum*
#1: 216.75, 218.50 rest_freq_DCN_F21,rest_freq_DCN_F22,rest_freq_H2CO_303_202,rest_freq_H2CO_322_221
#2: 231.00, 231.40 continuum
#3: 230.50, 230.96 rest_freq_12CO
#4: 219.05, 219.45 continuum
#5: 218.60, 219.00 rest_freq_H2CO_321_220
#6: 219.96, 220.43 rest_freq_SO,rest_freq_13CO
#7: 219.50, 219.96 rest_freq_C18O,rest_freq_SO

#SB0s
#0: 231.00, 231.40 continuum
#1: 230.50, 230.96 rest_freq_12CO
#2: 231.50, 233.25 continuum*
#3: 219.05, 219.45 continuum
#4: 218.60, 219.00 rest_freq_H2CO_321_220
#5: 219.96, 220.42 rest_freq_SO,rest_freq_13CO
#6: 219.50, 219.96 rest_freq_C18O,rest_freq_SO
#7: 216.20, 217.80 rest_freq_H2CO_918_919,rest_freq_DCN_F21,rest_freq_DCN_F22

#SB1s
#0: 231.50, 233.20 continuum*
#1: 216.75, 218.50 rest_freq_DCN_F21,rest_freq_DCN_F22,rest_freq_H2CO_303_202,rest_freq_H2CO_322_221
#2: 231.00, 231.40 continuum
#3: 230.48, 230.93 rest_freq_12CO
#4: 219.05, 219.45 continuum
#5: 218.58, 218.97 rest_freq_H2CO_321_220
#6: 219.95, 220.42 rest_freq_SO,rest_freq_13CO
#7: 219.48, 219.95 rest_freq_C18O,rest_freq_SO

line_spws_SBs  = [[1, 4, 5,5, 6,6, 7,7,7],[1,1,1,1, 3, 5, 6,6, 7,7]]

line_freqs_SBs = [
    [rest_freq_12CO, rest_freq_H2CO_321_220, rest_freq_SO,rest_freq_13CO, rest_freq_C18O,rest_freq_SO, rest_freq_H2CO_918_919,rest_freq_DCN_F21,rest_freq_DCN_F22],
    [rest_freq_DCN_F21,rest_freq_DCN_F22,rest_freq_H2CO_303_202,rest_freq_H2CO_322_221, rest_freq_12CO, rest_freq_H2CO_321_220, rest_freq_SO,rest_freq_13CO, rest_freq_C18O,rest_freq_SO]
]

spwcont_forplot_SBs = ['2','0']

width_array_SBs = [[60,960,16,30,30,480,480,16],[8,8,30,480,12,12,240,240]]

data_params_LB = {
    f'LB{i}': {
        'vis':          f'{prefix}_LB_EB{i}.ms',
        'name':         f'LB_EB{i}',
        'field':        'CQ_Tau',
        'line_spws':    np.array([1,1,1,1, 3, 5, 6,6, 7,7]), #list of spws containing lines
        'line_freqs':   np.array([
            rest_freq_DCN_F21,rest_freq_DCN_F22,rest_freq_H2CO_303_202,rest_freq_H2CO_322_221,
            rest_freq_12CO,
            rest_freq_H2CO_321_220,
            rest_freq_SO,rest_freq_13CO,
            rest_freq_C18O,rest_freq_SO
        ]), #frequencies (Hz) corresponding to line_spws
        'spwcont_forplot': '0',
        'cont_spws':       '0,1,2,3,4,5,6,7',
        'width_array':     [4,4,15,240,6,6,120,120],
    } for i in range(number_of_EBs['LB']) #phase RMS: 28.985, 22.576 deg
}

data_params_SB = {
    f'SB{i}': {
        'vis':          f'{prefix}_SB_EB{i}.ms',
        'name':         f'SB_EB{i}',
        'field':        'CQ_Tau',
        'line_spws':    line_spws_SBs[i], #list of spws containing lines
        'line_freqs':   line_freqs_SBs[i], #frequencies (Hz) corresponding to line_spws
        'spwcont_forplot': spwcont_forplot_SBs[i],
        'cont_spws':       '0,1,2,3,4,5,6,7',
        'width_array':     width_array_SBs[i],
    } for i in range(number_of_EBs['SB'])
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

#Adjust these plot ranges according to your data, check from weblogs or QA2 report
plotranges = {'SB':[0,3000,0,0.25], 'LB':[0,8550,0,0.25]} #xmin,xmax,ymin,ymax

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
#Then plot amplitude vs frequency to check if flagging worked fine.
for params in data_params.values():
    os.system(f'rm -rf '+prefix+'_flagtest_'+params['name']+'_initcont.ms')

    baseline_key, _ = params['name'].split('_')
    _, idx_key      = params['name'].split('EB')
    contspws        = params['cont_spws']

    flagchannels_string = get_flagchannels(ms_dict=params,output_prefix=prefix,velocity_range=np.array([-15.,15.]) + v_sys)
    # Flagchannels input string for LB_EB0: '1:88~89, 1:88~89, 1:25~26, 1:9~10, 3:112~207, 5:33~35, 6:959~959, 6:53~98, 7:810~855, 7:13~58'
    # Flagchannels input string for LB_EB1: '1:88~89, 1:88~89, 1:25~26, 1:9~10, 3:112~207, 5:33~35, 6:959~959, 6:53~98, 7:810~855, 7:13~58'
    # Flagchannels input string for SB_EB0: '1:128~222, 4:32~35, 5:959~959, 5:45~91, 6:802~847, 6:5~51, 7:90~91, 7:47~48, 7:47~48'
    # Flagchannels input string for SB_EB1: '1:87~89, 1:87~89, 1:24~26, 1:8~9, 3:159~254, 5:31~34, 6:950~959, 6:30~75, 7:787~832, 7:0~36'

    if params['name'] == 'SB_EB0':
        flagchannels_string += ', 0:60~74' #SB0, 0:60~74 to be flagged

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

print(max_chan_width) #72.160, 72.160, 424.431, 166.838
print(f'These datasets require a maximum channel width of {np.amin(max_chan_width)} MHz to avoid bandwidth smearing.')

#LBs [4,4,15,240,6,6,120,120]
#0,  128, 15625.000kHz   4  62.50MHz
#1,  128, 15625.000kHz   4  62.50MHz
#2,  120,  3906.250kHz  15  58.59MHz
#3, 1920,   244.141kHz 240  58.59MHz
#4,   60,  7812.500kHz   6  46.87MHz
#5,   60,  7812.500kHz   6  46.87MHz
#6,  960,   488.281kHz 120  58.59MHz
#7,  960,   488.281kHz 120  58.59MHz

#SB0 [60,960,16,30,30,480,480,16]
#0,  120,  3906.250kHz  60 234.38MHz
#1, 1920,   244.141kHz 960 234.38MHz
#2,  128, 15625.000kHz  16 250.00MHz
#3,   60,  7812.500kHz  30 234.38MHz
#4,   60,  7812.500kHz  30 234.38MHz
#5,  960,   488.281kHz 480 234.38MHz
#6,  960,   488.281kHz 480 234.38MHz
#7,  128, 15625.000kHz  16 250.00MHz

#SB1 [8,8,30,480,12,12,240,240]
#0,  128, 15625.000kHz   8 125.00MHz
#1,  128, 15625.000kHz   8 125.00MHz
#2,  120,  3906.250kHz  30 117.19MHz
#3, 1920,   244.141kHz 480 117.19MHz
#4,   60,  7812.500kHz  12  93.75MHz
#5,   60,  7812.500kHz  12  93.75MHz
#6,  960,   488.281kHz 240 117.18MHz
#7,  960,   488.281kHz 240 117.18MHz

#Do the averaging (using the width_array corresponding to the minimum between the maximum channel width above and 250 MHz, listed above)
for params in data_params.values():
    os.system(f'rm -rf '+prefix+'_'+params['name']+'_initcont.ms')

    flagchannels_string = get_flagchannels(ms_dict=params,output_prefix=prefix,velocity_range=np.array([-15.,15.]) + v_sys)
    #Double-check that the channels idenfied are at the center of the spws, due to a potential issue with the data_desc_id key in the ms table of some programs

    if params['name'] == 'SB_EB0':
        flagchannels_string += ', 0:60~74' #SB0, 0:60~74 to be flagged

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

uv_ranges = {'LB':'125~150m','SB':'125~150m'}

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
mask_pa        = PA #position angle of mask in degrees
mask_semimajor = {'LB':0.75,'SB':0.85} #semimajor axis of mask in arcsec
mask_semiminor = {'LB':0.65,'SB':0.72} #semiminor axis of mask in arcsec #mask_semimajor*np.cos(incl/180.*np.pi)
mask_ra        = {
    'LB':['05h35m58.4723025387s','05h35m58.4723025387s'],'SB':['05h35m58.4655896712s','05h35m58.4678888058s']
} #from imview
mask_dec       = {
    'LB':['24d44m53.5917580665s','24d44m53.5917580665s'],'SB':['24d44m53.7121910554s','24d44m53.6024556007s']
} #from imview

# Cellsize: ~beam/8
cellsize = {'LB':['0.006arcsec','0.006arcsec'],'SB':['0.030arcsec','0.013arcsec']}

# Image size: ~primary beam 1.22*lam/A = 32'' with A=12m (19 arcsec)
imsize = {'LB':[6400,6400],'SB':[1400,3600]}
scales = {'LB':[0,2,4,8,16,24],'SB':[0,2,4,8]} #32 seems too big of a scale for the LBs, same for 16 for the SBs

# Threshold: ~6 sigma for phase cal
thresholds = {'LB':['0.1656mJy','0.1422mJy'],'SB':['1.4640mJy','0.3102mJy']} 

image_png_plot_sizes = [3,10] #sizes in arcsec of the zoomed and overview plots of the pngs

preselfcal_images_png_folder = get_figures_folderpath('3_preselfcal_images')
make_figures_folder(preselfcal_images_png_folder)

for baseline_key,params in zip(('LB','SB'),(data_params_LB,data_params_SB)):
    for EB_key,p in params.items():
        os.system(f'rm -rf '+prefix+'_'+p['name']+'_initcont_image*')
        _, idx_key    = p['name'].split('EB')

        mask = f'ellipse[[{mask_ra[baseline_key][int(idx_key)]},{mask_dec[baseline_key][int(idx_key)]}], [{mask_semimajor[baseline_key]:.3f}arcsec, {mask_semiminor[baseline_key]:.3f}arcsec], {mask_pa:.1f}deg]'
        noise_annulus = f"annulus[[{mask_ra[baseline_key][int(idx_key)]}, {mask_dec[baseline_key][int(idx_key)]}],['4.arcsec', '6.arcsec']]"

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
            cellsize       = cellsize[baseline_key][int(idx_key)],
            imsize         = imsize[baseline_key][int(idx_key)],
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

#CQ_Tau_LB_EB0_initcont_image.image
#Beam 0.073 arcsec x 0.050 arcsec (-25.01 deg)
#Flux inside disk mask: 133.63 mJy
#Peak intensity of source: 1.68 mJy/beam
#rms: 2.77e-02 mJy/beam
#Peak SNR: 60.56

#CQ_Tau_LB_EB1_initcont_image.image
#Beam 0.080 arcsec x 0.048 arcsec (20.61 deg)
#Flux inside disk mask: 141.73 mJy
#Peak intensity of source: 2.04 mJy/beam
#rms: 2.41e-02 mJy/beam
#Peak SNR: 84.67

#CQ_Tau_SB_EB0_initcont_image.image
#Beam 0.273 arcsec x 0.235 arcsec (-14.74 deg)
#Flux inside disk mask: 186.39 mJy
#Peak intensity of source: 27.25 mJy/beam
#rms: 2.51e-01 mJy/beam
#Peak SNR: 108.52

#CQ_Tau_SB_EB1_initcont_image.image
#Beam 0.127 arcsec x 0.099 arcsec (-0.60 deg)
#Flux inside disk mask: 113.73 mJy
#Peak intensity of source: 4.81 mJy/beam
#rms: 5.22e-02 mJy/beam
#Peak SNR: 92.10

for baseline_key,params in zip(('LB','SB'),(data_params_LB,data_params_SB)):
    for EB_key,p in params.items():
        os.system(f'rm -rf '+prefix+'_'+p['name']+'_initcont_statwt.ms')
        shutil.copytree(src=prefix+'_'+p['name']+'_initcont.ms',dst=prefix+'_'+p['name']+'_initcont_statwt.ms')

        statwt(
            vis        = prefix+'_'+p['name']+'_initcont_statwt.ms',
            spw        = p['cont_spws'],
            field      = 'CQ_Tau',
            combine    = 'scan,spw,corr,field',
            timebin    = '10000000s',    #"inf" does not work...
            chanbin    = 'spw',          #"... channel binning occurs within individual spectral windows; bins never span multiple spectral windows..."
            datacolumn = 'DATA',
            intent     = 'OBSERVE_TARGET#ON_SOURCE',
        )

        os.system(f'rm -rf '+prefix+'_'+p['name']+'_initcont_statwt_image*')
        _, idx_key    = p['name'].split('EB')

        mask = f'ellipse[[{mask_ra[baseline_key][int(idx_key)]},{mask_dec[baseline_key][int(idx_key)]}], [{mask_semimajor[baseline_key]:.3f}arcsec, {mask_semiminor[baseline_key]:.3f}arcsec], {mask_pa:.1f}deg]'
        noise_annulus = f"annulus[[{mask_ra[baseline_key][int(idx_key)]}, {mask_dec[baseline_key][int(idx_key)]}],['4.arcsec', '6.arcsec']]"

        imagename = prefix+'_'+p['name']+'_initcont_statwt_image'
        tclean_wrapper(
            vis            = prefix+'_'+p['name']+'_initcont_statwt.ms',
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
            cellsize       = cellsize[baseline_key][int(idx_key)],
            imsize         = imsize[baseline_key][int(idx_key)],
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

#CQ_Tau_LB_EB0_initcont_statwt_image.image
#Beam 0.073 arcsec x 0.050 arcsec (-25.17 deg)
#Flux inside disk mask: 132.49 mJy
#Peak intensity of source: 1.67 mJy/beam
#rms: 2.79e-02 mJy/beam
#Peak SNR: 59.90

#CQ_Tau_LB_EB1_initcont_statwt_image.image
#Beam 0.080 arcsec x 0.048 arcsec (20.69 deg)
#Flux inside disk mask: 141.86 mJy
#Peak intensity of source: 2.00 mJy/beam
#rms: 2.40e-02 mJy/beam
#Peak SNR: 83.30

#CQ_Tau_SB_EB0_initcont_statwt_image.image
#Beam 0.272 arcsec x 0.234 arcsec (-15.01 deg)
#Flux inside disk mask: 186.47 mJy
#Peak intensity of source: 27.18 mJy/beam
#rms: 2.48e-01 mJy/beam
#Peak SNR: 109.37

#CQ_Tau_SB_EB1_initcont_statwt_image.image
#Beam 0.127 arcsec x 0.099 arcsec (-0.57 deg)
#Flux inside disk mask: 113.60 mJy
#Peak intensity of source: 4.78 mJy/beam
#rms: 5.15e-02 mJy/beam
#Peak SNR: 92.71

#Selfcal individual EBs

#Self-calibration parameters
single_EB_contspws = '0~7'
single_EB_spw_mapping = [0,0,0,0,0,0,0,0]

individual_EB_selfcal_shift_folder = get_figures_folderpath('4_individual_EB_selfcal_and_shift_figures')
make_figures_folder(individual_EB_selfcal_shift_folder)

#One round of phase-only self-cal
for params in data_params.values():
    
    #vis = prefix+'_'+params['name']+'_initcont_statwt.ms'
    single_EB_p1 = prefix+'_'+params['name']+'_initcont.p1'
    #os.system(f'rm -rf {single_EB_p1}')
    #gaincal(vis=vis,caltable=single_EB_p1,gaintype='T',spw=single_EB_contspws,combine='scan,spw',calmode='p',solint='inf',minsnr=4,minblperant=3)
       
    #Print calibration png file
    plotfilename = prefix+'_'+params['name']+'_initcont_gain_p1_phase_vs_time.png'
    # plotms(single_EB_p1,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,plotfile=os.path.join(individual_EB_selfcal_shift_folder,plotfilename))
    plotms(
        single_EB_p1,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,
        plotfile=os.path.join(individual_EB_selfcal_shift_folder,plotfilename),
        customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',
        #iteraxis='spw',exprange='all'
    )
    
    if params['name'] == 'SB_EB1':
        flagdata(vis=prefix+'_'+params['name']+'_initcont.p1',mode='manual',antenna='DV12')
    
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
    
#LB_EB0:        1 of 44 solutions flagged due to SNR < 4 in spw=0 at 2017/11/20/07:41:43.3
#LB_EB0_statwt: 3 of 44 solutions flagged due to SNR < 4 in spw=0 at 2017/11/20/07:42:15.2
#LB_EB1:        1 of 47 solutions flagged due to SNR < 4 in spw=0 at 2017/11/23/05:09:07.6
#LB_EB1_statwt: 2 of 47 solutions flagged due to SNR < 4 in spw=0 at 2017/11/23/05:08:25.7

#Image the self-calibrated EBs
for baseline_key,params in zip(('LB','SB'),(data_params_LB,data_params_SB)):
    for EB_key,p in params.items():
        os.system(f'rm -rf '+prefix+'_'+p['name']+'_initcont_selfcal_image*')
        _, idx_key    = p['name'].split('EB')

        mask = f'ellipse[[{mask_ra[baseline_key][int(idx_key)]},{mask_dec[baseline_key][int(idx_key)]}], [{mask_semimajor[baseline_key]:.3f}arcsec, {mask_semiminor[baseline_key]:.3f}arcsec], {mask_pa:.1f}deg]'
        noise_annulus = f"annulus[[{mask_ra[baseline_key][int(idx_key)]}, {mask_dec[baseline_key][int(idx_key)]}],['4.arcsec', '6.arcsec']]"

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
            cellsize       = cellsize[baseline_key][int(idx_key)],
            imsize         = imsize[baseline_key][int(idx_key)],
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

#CQ_Tau_LB_EB0_initcont_selfcal_image.image
#Beam 0.073 arcsec x 0.050 arcsec (-25.17 deg)
#Flux inside disk mask: 132.94 mJy
#Peak intensity of source: 1.69 mJy/beam
#rms: 2.61e-02 mJy/beam
#Peak SNR: 64.83

#CQ_Tau_LB_EB1_initcont_selfcal_image.image
#Beam 0.080 arcsec x 0.048 arcsec (20.69 deg)
#Flux inside disk mask: 142.31 mJy
#Peak intensity of source: 2.03 mJy/beam
#rms: 2.11e-02 mJy/beam
#Peak SNR: 96.25

#CQ_Tau_SB_EB0_initcont_selfcal_image.image
#Beam 0.272 arcsec x 0.234 arcsec (-15.01 deg)
#Flux inside disk mask: 187.39 mJy
#Peak intensity of source: 27.72 mJy/beam
#rms: 2.03e-01 mJy/beam
#Peak SNR: 136.71

#CQ_Tau_SB_EB1_initcont_selfcal_image.image
#Beam 0.127 arcsec x 0.099 arcsec (-0.57 deg)
#Flux inside disk mask: 114.01 mJy
#Peak intensity of source: 4.80 mJy/beam
#rms: 4.86e-02 mJy/beam
#Peak SNR: 98.77

# Cellsize: ~beam/8
cellsize = {'LB':'0.006arcsec','SB':'0.013arcsec'}

# Image size: ~primary beam 1.22*lam/A = 32'' with A=12m (19 arcsec)
imsize   = {'LB':6400,'SB':3600}

#Image the self-calibrated EBs with same cellsize and imsize for ratio
for baseline_key,params in zip(('LB','SB'),(data_params_LB,data_params_SB)):
    for EB_key,p in params.items():
        os.system(f'rm -rf '+prefix+'_'+p['name']+'_initcont_selfcal_aligncomp_image*')
        _, idx_key    = p['name'].split('EB')

        mask = f'ellipse[[{mask_ra[baseline_key][int(idx_key)]},{mask_dec[baseline_key][int(idx_key)]}], [{mask_semimajor[baseline_key]:.3f}arcsec, {mask_semiminor[baseline_key]:.3f}arcsec], {mask_pa:.1f}deg]'
        noise_annulus = f"annulus[[{mask_ra[baseline_key][int(idx_key)]}, {mask_dec[baseline_key][int(idx_key)]}],['4.arcsec', '6.arcsec']]"

        imagename = prefix+'_'+p['name']+'_initcont_selfcal_aligncomp_image'
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

#CQ_Tau_LB_EB0_initcont_selfcal_aligncomp_image.image
#Beam 0.073 arcsec x 0.050 arcsec (-25.17 deg)
#Flux inside disk mask: 132.94 mJy
#Peak intensity of source: 1.69 mJy/beam
#rms: 2.61e-02 mJy/beam
#Peak SNR: 64.83

#CQ_Tau_LB_EB1_initcont_selfcal_aligncomp_image.image
#Beam 0.080 arcsec x 0.048 arcsec (20.69 deg)
#Flux inside disk mask: 142.31 mJy
#Peak intensity of source: 2.03 mJy/beam
#rms: 2.11e-02 mJy/beam
#Peak SNR: 96.25

#CQ_Tau_SB_EB0_initcont_selfcal_aligncomp_image.image
#Beam 0.272 arcsec x 0.234 arcsec (-14.95 deg)
#Flux inside disk mask: 186.77 mJy
#Peak intensity of source: 27.68 mJy/beam
#rms: 2.03e-01 mJy/beam
#Peak SNR: 136.50

#CQ_Tau_SB_EB1_initcont_selfcal_aligncomp_image.image
#Beam 0.127 arcsec x 0.099 arcsec (-0.57 deg)
#Flux inside disk mask: 114.01 mJy
#Peak intensity of source: 4.80 mJy/beam
#rms: 4.86e-02 mJy/beam
#Peak SNR: 98.77

#For the images with same pixel scales, compute the image ratios
for baseline_key,params in zip(('LB','SB'),(data_params_LB,data_params_SB)):
    for p in params.values():
        imagename = prefix+'_'+p['name']+'_initcont_selfcal_aligncomp_image'
        #ratio image, to be compared to the ratio image after alignment
        ref_image = f'{prefix}_{baseline_key}_EB0_initcont_selfcal_aligncomp_image.image'
        ref_rms = params[f'{baseline_key}0']['rms']
        ratio_image = imagename+'.ratio'
        os.system(f'rm -rf {ratio_image}')
        immath(imagename=[ref_image,imagename+'.image'],mode='evalexpr',outfile=ratio_image,expr=f'iif(IM0 > 3*{ref_rms}, IM1/IM0, 0)')
        generate_image_png(
            ratio_image,plot_sizes=[2*mask_semimajor[baseline_key],2*mask_semimajor[baseline_key]],
            color_scale_limits=[0.5,1.5],image_units='ratio',
            save_folder=individual_EB_selfcal_shift_folder
        )

#Align data (go from *initcont_selfcal.ms to *initcont_shift.ms)

#Select the LB EB to act as the reference (usually the best SNR one)
reference_for_LB_alignment = f'{prefix}_LB_EB1_initcont_selfcal.ms'
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
cell_size = 0.01

offset_ms    = prefix+'_LB_EB0_initcont_selfcal.ms'

plotfilename = 'uv_overlap_LB_EB0.png'
offset = alignment.find_offset(
    reference_ms=reference_for_LB_alignment,
    offset_ms=offset_ms,npix=npix,cell_size=cell_size,
    spwid=continuum_spw_id,plot_uv_grid=True,
    uv_grid_plot_filename=os.path.join(individual_EB_selfcal_shift_folder,plotfilename)
)
print(f'#Offset for {offset_ms}: ',offset)
#Offset for CQ_Tau_LB_EB0_initcont_selfcal.ms:  [ 0.00589421 -0.00815621]

for _npix in [256,512,1024,2048]:
    plotfilename = f'uv_overlap_LB_EB0_{_npix}pixels.png'
    offset = alignment.find_offset(
        reference_ms=reference_for_LB_alignment,
        offset_ms=offset_ms,npix=_npix,cell_size=cell_size,
        spwid=continuum_spw_id,plot_uv_grid=True,
        uv_grid_plot_filename=os.path.join(individual_EB_selfcal_shift_folder,plotfilename)
    )
    print(f'#Offset for {offset_ms}: ',offset)
    #Offset for CQ_Tau_LB_EB0_initcont_selfcal.ms:  [ 0.00616088 -0.00698912]
    #Offset for CQ_Tau_LB_EB0_initcont_selfcal.ms:  [ 0.00680972 -0.00868555]
    #Offset for CQ_Tau_LB_EB0_initcont_selfcal.ms:  [ 0.00589421 -0.00815621]
    #Offset for CQ_Tau_LB_EB0_initcont_selfcal.ms:  [ 0.00447982 -0.00548735]

alignment.align_measurement_sets(
    reference_ms=reference_for_LB_alignment,align_ms=offset_LB_EBs,
    npix=npix,cell_size=cell_size,spwid=continuum_spw_id
)
#New coordinates for CQ_Tau_LB_EB0_initcont_selfcal.ms requires a shift of [0.0058942,-0.0081562]
#New coordinates for CQ_Tau_LB_EB1_initcont_selfcal.ms no shift, reference MS.

#Insert offsets from the alignment output
alignment_offsets['LB_EB0'] = [0.0058942,-0.0081562]
alignment_offsets['LB_EB1'] = [0.0,       0.0      ]

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
    #Offset for CQ_Tau_LB_EB0_initcont_selfcal_shift.ms:  [-4.21684156e-05  4.98583148e-07]

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
        npix=npix,cell_size=cell_size,spwid=int(params['spwcont_forplot']),plot_uv_grid=True,
        uv_grid_plot_filename=os.path.join(individual_EB_selfcal_shift_folder,plotfilename)
    )
    print(f'#Offset for {offset_ms}: ',offset)
    #Offset for CQ_Tau_SB_EB0_initcont_selfcal.ms:  [ 0.01896135 -0.01778067]
    #Offset for CQ_Tau_SB_EB1_initcont_selfcal.ms:  [-0.03430092 -0.00735719]
    
    for _npix in [256,512,1024,2048]:
        plotfilename = 'uv_overlap_'+params['name']+f'_{_npix}pixels.png'
        offset = alignment.find_offset(
            reference_ms=reference_for_SB_alignment,
            offset_ms=offset_ms,
            npix=_npix,cell_size=cell_size,spwid=int(params['spwcont_forplot']),plot_uv_grid=True,
            uv_grid_plot_filename=os.path.join(individual_EB_selfcal_shift_folder,plotfilename)
        )
        print(f'#Offset for {offset_ms}: ',offset)
        #Offset for CQ_Tau_SB_EB0_initcont_selfcal.ms:  [ 0.01916959 -0.01579964]
        #Offset for CQ_Tau_SB_EB0_initcont_selfcal.ms:  [ 0.02054213 -0.0149247 ]
        #Offset for CQ_Tau_SB_EB0_initcont_selfcal.ms:  [ 0.01896135 -0.01778067]
        #Offset for CQ_Tau_SB_EB0_initcont_selfcal.ms:  [ 0.01763304 -0.01721871]
        #Offset for CQ_Tau_SB_EB1_initcont_selfcal.ms:  [-0.03291183 -0.00654644]
        #Offset for CQ_Tau_SB_EB1_initcont_selfcal.ms:  [-0.03382735 -0.00793079]
        #Offset for CQ_Tau_SB_EB1_initcont_selfcal.ms:  [-0.03430092 -0.00735719]
        #Offset for CQ_Tau_SB_EB1_initcont_selfcal.ms:  [-0.03444217 -0.00502781]

    alignment.align_measurement_sets(
        reference_ms=reference_for_SB_alignment,align_ms=offset_ms,
        npix=npix,cell_size=cell_size,spwid=int(params['spwcont_forplot'])
    )
    #New coordinates for CQ_Tau_SB_EB0_initcont_selfcal.ms requires a shift of [0.018961,-0.017781]
    #New coordinates for CQ_Tau_SB_EB1_initcont_selfcal.ms requires a shift of [-0.034301,-0.0073572]

alignment_offsets['SB_EB0'] = [ 0.018961,-0.017781 ]
alignment_offsets['SB_EB1'] = [-0.034301,-0.0073572]

shifted_SB_EBs = [EB.replace('.ms','_shift.ms') for EB in offset_SB_EBs]

#Check by calculating offset again
for params in data_params_SB.values():
    shifted_ms = prefix+'_'+params['name']+'_initcont_selfcal_shift.ms'
    offset = alignment.find_offset(
        reference_ms=reference_for_SB_alignment,
        offset_ms=shifted_ms,npix=npix,plot_uv_grid=False,
        cell_size=cell_size,spwid=int(params['spwcont_forplot'])
    )
    print(f'#Offset for {shifted_ms}: ',offset)
#Offset for CQ_Tau_SB_EB0_initcont_selfcal_shift.ms:  [-1.50685424e-05  6.99559272e-05]
#Offset for CQ_Tau_SB_EB1_initcont_selfcal_shift.ms:  [-1.90399831e-04 -3.21126209e-05]

#Remove the '_selfcal' part of the names to match the naming convention below.
for shifted_EB in shifted_LB_EBs+shifted_SB_EBs:
    os.system('mv {} {}'.format(shifted_EB, shifted_EB.replace('_selfcal', '')))

#Check that the images are indeed aligned after the shift
for baseline_key,params in zip(('LB','SB'),(data_params_LB,data_params_SB)):
    for EB_key,p in params.items():
        os.system(f'rm -rf '+prefix+'_'+p['name']+'_initcont_shift_aligncomp_image*')
        _, idx_key    = p['name'].split('EB')

        mask = f"ellipse[[{mask_ra['LB'][1]},{mask_dec['LB'][1]}], [{mask_semimajor['LB']:.3f}arcsec, {mask_semiminor['LB']:.3f}arcsec], {mask_pa:.1f}deg]"
        noise_annulus = f"annulus[[{mask_ra['LB'][1]}, {mask_dec['LB'][1]}],['4.arcsec', '6.arcsec']]"

        imagename = prefix+'_'+p['name']+'_initcont_shift_aligncomp_image'
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

#CQ_Tau_LB_EB0_initcont_shift_aligncomp_image.image
#Beam 0.073 arcsec x 0.050 arcsec (-25.17 deg)
#Flux inside disk mask: 132.66 mJy
#Peak intensity of source: 1.70 mJy/beam
#rms: 2.62e-02 mJy/beam
#Peak SNR: 64.71

#CQ_Tau_LB_EB1_initcont_shift_aligncomp_image.image
#Beam 0.080 arcsec x 0.048 arcsec (20.69 deg)
#Flux inside disk mask: 142.31 mJy
#Peak intensity of source: 2.01 mJy/beam
#rms: 2.11e-02 mJy/beam
#Peak SNR: 95.49

#CQ_Tau_SB_EB0_initcont_shift_aligncomp_image.image
#Beam 0.272 arcsec x 0.234 arcsec (-14.96 deg)
#Flux inside disk mask: 187.87 mJy
#Peak intensity of source: 27.65 mJy/beam
#rms: 2.03e-01 mJy/beam
#Peak SNR: 136.13

#CQ_Tau_SB_EB1_initcont_shift_aligncomp_image.image
#Beam 0.127 arcsec x 0.099 arcsec (-0.57 deg)
#Flux inside disk mask: 116.28 mJy
#Peak intensity of source: 4.82 mJy/beam
#rms: 4.87e-02 mJy/beam
#Peak SNR: 99.06

for baseline_key,params in zip(('LB','SB'),(data_params_LB,data_params_SB)):
    for p in params.values():
        imagename = prefix+'_'+p['name']+'_initcont_shift_aligncomp_image'
        ref_image = f'{prefix}_{baseline_key}_EB0_initcont_shift_aligncomp_image.image'
        ref_rms = params[f'{baseline_key}0']['rms']
        ratio_image = imagename+'.ratio'
        os.system(f'rm -rf {ratio_image}')
        immath(imagename=[ref_image,imagename+'.image'],mode='evalexpr',outfile=ratio_image,expr=f'iif(IM0 > 3*{ref_rms}, IM1/IM0, 0)')
        generate_image_png(
            ratio_image,plot_sizes=[2*mask_semimajor[baseline_key],2*mask_semimajor[baseline_key]],
            color_scale_limits=[0.5,1.5],image_units='ratio',
            save_folder=individual_EB_selfcal_shift_folder
        )
#Also by checking the different images in casaviewer, the different EBs are well aligned

for i,params in enumerate(data_params_LB.values()):
    print('#'+params['name'])
    mask = f"ellipse[[{mask_ra['LB'][i]},{mask_dec['LB'][i]}], [{mask_semimajor['LB']:.3f}arcsec, {mask_semiminor['LB']:.3f}arcsec], {mask_pa:.1f}deg]"
    fit_gaussian(prefix+'_'+params['name']+'_initcont_selfcal_aligncomp_image.image',region=mask)

    #LB_EB0, 05h35m58.472389s +24d44m53.58122s
    #Pixel coordinates of peak: x = 3196.577 y = 3191.452
    #LB_EB1, 05h35m58.471897s +24d44m53.59290s
    #Pixel coordinates of peak: x = 3197.630 y = 3193.436

    mask = f"ellipse[[{mask_ra['LB'][1]},{mask_dec['LB'][1]}], [{mask_semimajor['LB']:.3f}arcsec, {mask_semiminor['LB']:.3f}arcsec], {mask_pa:.1f}deg]"
    fit_gaussian(prefix+'_'+params['name']+'_initcont_shift_aligncomp_image.image',region=mask)

    #LB_EB0, 05h35m58.472414s +24d44m53.58068s
    #Pixel coordinates of peak: x = 3200.699 y = 3194.540
    #LB_EB1, 05h35m58.472402s +24d44m53.58384s
    #Pixel coordinates of peak: x = 3200.725 y = 3195.065

    #Differences are sub-pixels!

for i,params in enumerate(data_params_SB.values()):
    print('#'+params['name'])
    mask = f"ellipse[[{mask_ra['SB'][i]},{mask_dec['SB'][i]}], [{mask_semimajor['SB']:.3f}arcsec, {mask_semiminor['SB']:.3f}arcsec], {mask_pa:.1f}deg]"
    fit_gaussian(prefix+'_'+params['name']+'_initcont_selfcal_aligncomp_image.image',region=mask)

    #SB_EB0, 05h35m58.470699s +24d44m53.65246s
    #Pixel coordinates of peak: x = 1798.674 y = 1796.675
    #SB_EB1, 05h35m58.469115s +24d44m53.58912s
    #Pixel coordinates of peak: x = 1802.021 y = 1796.108
    
    mask = f"ellipse[[{mask_ra['LB'][1]},{mask_dec['LB'][1]}], [{mask_semimajor['SB']:.3f}arcsec, {mask_semiminor['SB']:.3f}arcsec], {mask_pa:.1f}deg]"
    fit_gaussian(prefix+'_'+params['name']+'_initcont_shift_aligncomp_image.image',region=mask)

    #SB_EB0, 05h35m58.472540s +24d44m53.58913s
    #Pixel coordinates of peak: x = 1800.191 y = 1798.129
    #SB_EB1, 05h35m58.471915s +24d44m53.58031s
    #Pixel coordinates of peak: x = 1800.846 y = 1797.451

    #Differences are sub-pixels!

#Now that everything is aligned, we inspect the flux calibration
for params in data_params.values():
    msfile = prefix+'_'+params['name']+'_initcont_shift.ms'
    #Export MS contents into numpy save files
    export_MS(msfile)
    #Measurement set exported to CQ_Tau_LB_EB0_initcont_shift.vis.npz
    #Measurement set exported to CQ_Tau_LB_EB1_initcont_shift.vis.npz
    #Measurement set exported to CQ_Tau_SB_EB0_initcont_shift.vis.npz
    #Measurement set exported to CQ_Tau_SB_EB1_initcont_shift.vis.npz

#Plot deprojected visibility profiles for all data together
list_npz_files = []
for baseline_key,n_EB in number_of_EBs.items():
    list_npz_files += [f'{prefix}_{baseline_key}_EB{i}_initcont_shift.vis.npz' for i in range(n_EB)]

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

#There are traces of waterfall features, but they're rather mild. The SB_EB with the
#highest SNR is the oldest (2013) one, so better not to use it as a guide.
#LB_EB1 and SB_EB1 have rougly the same SNR, go for LB_EB1 that shows less signs
#of decoherence, a higher SNR and lower phase RMS (22.576 vs 28.985 deg) than LB_EB0.

#Decoherence in SB_EB1 is clear from plot of the flux ratios.
#Flux offset of SB_EB0 is clear from the deprojected visibilities.
flux_ref_EB = 'LB_EB1' 

for params in data_params.values():
    plot_label = os.path.join(flux_comparison_folder,'flux_comparison_'+params['name']+f'_to_{flux_ref_EB}.png')
    estimate_flux_scale(
        reference=f'{prefix}_{flux_ref_EB}_initcont_shift.vis.npz',
        comparison=prefix+'_'+params['name']+'_initcont_shift.vis.npz',
        incl=incl,PA=PA,plot_label=plot_label
    )

#The ratio of the fluxes of CQ_Tau_LB_EB0_initcont_shift.vis.npz to
#CQ_Tau_LB_EB1_initcont_shift.vis.npz is 0.95160
#The scaling factor for gencal is 0.976 for your comparison measurement
#The error on the weighted mean ratio is 7.605e-04, although it's likely that
#the weights in the measurement sets are too off by some constant factor

#The ratio of the fluxes of CQ_Tau_LB_EB1_initcont_shift.vis.npz to
#CQ_Tau_LB_EB1_initcont_shift.vis.npz is 1.00000
#The scaling factor for gencal is 1.000 for your comparison measurement
#The error on the weighted mean ratio is 6.636e-04, although it's likely that
#the weights in the measurement sets are too off by some constant factor

#The ratio of the fluxes of CQ_Tau_SB_EB0_initcont_shift.vis.npz to
#CQ_Tau_LB_EB1_initcont_shift.vis.npz is 1.36589
#The scaling factor for gencal is 1.169 for your comparison measurement
#The error on the weighted mean ratio is 1.066e-03, although it's likely that
#the weights in the measurement sets are too off by some constant factor

#The ratio of the fluxes of CQ_Tau_SB_EB1_initcont_shift.vis.npz to
#CQ_Tau_LB_EB1_initcont_shift.vis.npz is 0.91086
#The scaling factor for gencal is 0.954 for your comparison measurement
#The error on the weighted mean ratio is 8.782e-04, although it's likely that
#the weights in the measurement sets are too off by some constant factor

#The flux offsets are all >4%, for SB_EB1 decoherence is a likely culprit, LB_EB0 also shows strange trends with uv-distance
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
# For CQ Tau, I see 'waterfalls' in the SB EBs and descending trends in the flux ratios, indicating decoherence. 
# Therefore, selfcal (phase-only) should fix it. I'll proceed by self-calibrating the SBs without rescaling
# Then I'll check the flux offsets at the end of the process...

#Begin of SB self-cal - iteration 1
#For phase self-cal, clean down to ~6sigma
SB_selfcal_folder = get_figures_folderpath('7_selfcal_SB_figures')
make_figures_folder(SB_selfcal_folder)

SB_cont_p0 = prefix+'_SB_contp0'
os.system('rm -rf %s.ms*' %SB_cont_p0)

concat(
    vis=[f'{prefix}_SB_EB{i}_initcont_shift.ms' for i in range(number_of_EBs['SB'])],
    concatvis=SB_cont_p0+'.ms',dirtol='0.1arcsec',copypointing=False
)
listobs(vis=SB_cont_p0+'.ms',listfile=SB_cont_p0+'.ms.listobs.txt',overwrite=True)
#2024-12-10 21:26:44     WARN    MSConcat::copySysCal    
#    /data/beegfs/astro-storage/groups/benisty/frzagaria/SO_detections/CQTau/selfcal_products/CQ_Tau_SB_contp0.ms does not have a valid syscal table,
#    the MS to be appended, however, has one. Result won't have one.
#2024-12-10 21:26:44     WARN    MSConcat::concatenate
#    (file /source/casa6/casatools/casacore/ms/MSOper/MSConcat.cc, line 1000). Could not merge SysCal subtables.

#Define new SB mask using new center read from listobs
mask_pa        = PA
mask_semimajor = 0.85
mask_semiminor = 0.72
mask_ra  = '05h35m58.472722s'
mask_dec = '24d44m53.613450s'

SB_mask= f"ellipse[[{mask_ra},{mask_dec}], [{mask_semimajor:.3f}arcsec, {mask_semiminor:.3f}arcsec], {mask_pa:.1f}deg]"

noise_annulus_SB = f"annulus[[{mask_ra}, {mask_dec}],['4.arcsec', '6.arcsec']]"

SB_tclean_wrapper_kwargs = {
    'deconvolver':'multiscale','scales':[0,2,4,8,12],
    'smallscalebias':0.6,'gain':0.3,'cycleniter':300,
    'niter':1000000,'mask':SB_mask,
    'cellsize':'0.013arcsec','imsize':3600,
    'parallel':use_parallel,'savemodel':'modelcolumn',
    'robust':0.5,'interactive':False,
    'gridder':'standard',
}

tclean_wrapper(
    vis       = SB_cont_p0+'.ms',
    imagename = SB_cont_p0,
    threshold = '0.3834mJy', 
    **SB_tclean_wrapper_kwargs
)
estimate_SNR(SB_cont_p0+'.image',disk_mask=SB_mask,noise_mask=noise_annulus_SB)
rms_SB = imstat(imagename=SB_cont_p0+'.image',region=noise_annulus_SB)['rms'][0]
generate_image_png(
    SB_cont_p0+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_SB,10*rms_SB],
    save_folder=SB_selfcal_folder
)
#CQ_Tau_SB_contp0.image
#Beam 0.138 arcsec x 0.108 arcsec (-2.02 deg)
#Flux inside disk mask: 169.43 mJy
#Peak intensity of source: 7.17 mJy/beam
#rms: 6.36e-02 mJy/beam
#Peak SNR: 112.77

#Skip: no weblogs for the SBs (too old datasets), don't know how to get the reference antennas
#
#Look for references antennas from the weblog, and pick the first ones that are listed, overlapping with all EBs
#Both EBs: DA63, DA52, DA57, DV08
#
#Get station numbers
#for ref_ant in ('DA63','DA52', 'DA57', 'DV08'):
#    get_station_numbers(SB_cont_p0+'.ms',ref_ant)
#    Observation ID 0: DA63@A035
#    Observation ID 1: DA63@A035
#    Observation ID 0: DA52@A002
#    Observation ID 1: DA52@A002
#    Observation ID 1: DA57@A001
#    Observation ID 0: DV08@A036
#    Observation ID 1: DV08@A036

SB_refant = ''

SB_contspws    = '0~15'
SB_spw_mapping = [
    0,0,0,0,0,0,0,0,
    8,8,8,8,8,8,8,8,
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

#Flag problematic antennas
flagdata(vis=SB_p1,mode='manual',spw='8',antenna='DV12')

#Inspect gain tables to check if flagging worked
plotms(
    SB_p1,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(SB_selfcal_folder,f'{prefix}_SB_p1_gain_phase_vs_time_flagged.png'),
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
    threshold = '0.3810mJy',
    **SB_tclean_wrapper_kwargs
)
estimate_SNR(SB_cont_p1+'.image',disk_mask=SB_mask,noise_mask=noise_annulus_SB)
generate_image_png(
    SB_cont_p1+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_SB,10*rms_SB],
    save_folder=SB_selfcal_folder
)
#CQ_Tau_SB_contp1.image
#Beam 0.138 arcsec x 0.108 arcsec (-2.02 deg)
#Flux inside disk mask: 169.56 mJy
#Peak intensity of source: 7.19 mJy/beam
#rms: 6.35e-02 mJy/beam
#Peak SNR: 113.32

#Second round of phase-only self-cal
SB_p2 = SB_p1.replace('p1','p2')
os.system('rm -rf '+SB_p2)
gaincal(
    vis=SB_cont_p1+'.ms',caltable=SB_p2,
    gaintype='T',combine='scan,spw',calmode='p',solint='360s',
    spw=SB_contspws,refant=SB_refant,
    minsnr=3.,minblperant=4
)

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
flagdata(vis=SB_p2,mode='manual',spw='0',timerange='2015/08/30/11:41:00~2015/08/30/11:43:00',antenna='DV03,DV07')
flagdata(vis=SB_p2,mode='manual',spw='0',timerange='2015/08/30/11:54:00~2015/08/30/11:56:00',antenna='DV07')
flagdata(vis=SB_p2,mode='manual',spw='8',timerange='2017/08/07/12:39:00~2017/08/07/12:41:00',antenna='DA53')
flagdata(vis=SB_p2,mode='manual',spw='8',timerange='2017/08/07/12:42:00~2017/08/07/12:44:00',antenna='DV12')
flagdata(vis=SB_p2,mode='manual',spw='8',timerange='2017/08/07/12:46:00~2017/08/07/12:48:00',antenna='DA46,DV07,DV10,DV12')

#Inspect gain tables to check if flagging worked
plotms(
    SB_p2,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(SB_selfcal_folder,f'{prefix}_SB_p2_gain_phase_vs_time_flagged.png'),
)

#Apply the calibration gains and split-off the corrected ms
applycal(
    vis=SB_cont_p1+'.ms',spw=SB_contspws,spwmap=SB_spw_mapping,
    gaintable=[SB_p2],interp='linearPD',calwt=True,applymode='calonly'
)
SB_cont_p2 = SB_cont_p1.replace('p1','p2')
os.system('rm -rf %s.ms*' %SB_cont_p2)
split(vis=SB_cont_p1+'.ms',outputvis=SB_cont_p2+'.ms',datacolumn='corrected')

#CLEAN, again with a ~6sigma threshold
tclean_wrapper(
    vis       = SB_cont_p2+'.ms',
    imagename = SB_cont_p2,
    threshold = '0.3582mJy',
    **SB_tclean_wrapper_kwargs
)
estimate_SNR(SB_cont_p2+'.image',disk_mask=SB_mask,noise_mask=noise_annulus_SB)
generate_image_png(
    SB_cont_p2+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_SB,10*rms_SB],
    save_folder=SB_selfcal_folder
)
#CQ_Tau_SB_contp2.image
#Beam 0.138 arcsec x 0.108 arcsec (-2.02 deg)
#Flux inside disk mask: 171.03 mJy
#Peak intensity of source: 7.39 mJy/beam
#rms: 5.92e-02 mJy/beam
#Peak SNR: 124.91

#Third round of phase-only self-cal
SB_p3 = SB_p2.replace('p2','p3')
os.system('rm -rf '+SB_p3)
gaincal(
    vis=SB_cont_p2+'.ms',caltable=SB_p3,
    gaintype='T',combine='scan,spw',calmode='p',solint='120s', #diff from exoALMA, kept combine='scan' because scans for SB_EB1~60sec
    spw=SB_contspws,refant=SB_refant,
    minsnr=3.,minblperant=4
)
#11 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:39:52.5
# 1 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:41:18.3
# 1 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:42:27.7
# 1 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:43:56.9
# 2 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:48:13.6
# 2 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:58:25.4
# 1 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:59:26.1
# 2 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/13:03:48.8
# 4 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/13:06:29.7

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
flagdata(vis=SB_p3,mode='manual',spw='0',antenna='DV03,DV07')

flagdata(vis=SB_p3,mode='manual',spw='8',timerange='2017/08/07/12:39:20~2017/08/07/12:39:40',antenna='DA53,DV12')
flagdata(vis=SB_p3,mode='manual',spw='8',timerange='2017/08/07/12:42:20~2017/08/07/12:42:40',antenna='DV12')
flagdata(vis=SB_p3,mode='manual',spw='8',timerange='2017/08/07/12:43:50~2017/08/07/12:44:10',antenna='DV11,DV12')
flagdata(vis=SB_p3,mode='manual',spw='8',timerange='2017/08/07/12:46:50~2017/08/07/12:47:10',antenna='DA53,DV07,DV10,DV12')
flagdata(vis=SB_p3,mode='manual',spw='8',timerange='2017/08/07/13:08:00~2017/08/07/13:08:30',antenna='PM03')
flagdata(vis=SB_p3,mode='manual',spw='8',timerange='2017/08/07/13:09:00~2017/08/07/13:09:30',antenna='DA42')
flagdata(vis=SB_p3,mode='manual',spw='8',timerange='2017/08/07/13:10:30~2017/08/07/13:11:00',antenna='DA46')

#spw='0',timerange='2015/08/30/11:38:00~2015/08/30/11:40:00',antenna='DV03'
#spw='0',timerange='2015/08/30/11:41:00~2015/08/30/11:43:00',antenna='DV03,DV07'
#spw='0',timerange='2015/08/30/11:47:00~2015/08/30/11:49:00',antenna='DV03'
#spw='0',timerange='2015/08/30/11:49:00~2015/08/30/11:51:00',antenna='DV07'
#spw='0',timerange='2015/08/30/11:54:00~2015/08/30/11:56:00',antenna='DV07'
#spw='0',timerange='2015/08/30/11:55:00~2015/08/30/11:57:00',antenna='DV07'

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
    threshold = '0.3444mJy',
    **SB_tclean_wrapper_kwargs
)
estimate_SNR(SB_cont_p3+'.image',disk_mask=SB_mask,noise_mask=noise_annulus_SB)
generate_image_png(
    SB_cont_p3+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_SB,10*rms_SB],
    save_folder=SB_selfcal_folder
)
#CQ_Tau_SB_contp3.image
#Beam 0.138 arcsec x 0.108 arcsec (-2.02 deg)
#Flux inside disk mask: 171.63 mJy
#Peak intensity of source: 7.65 mJy/beam
#rms: 5.72e-02 mJy/beam
#Peak SNR: 133.70

#Fourth round of phase-only self-cal
SB_p4 = SB_p3.replace('p3','p4')
os.system('rm -rf '+SB_p4)
gaincal(
    vis=SB_cont_p3+'.ms',caltable=SB_p4,
    gaintype='T',combine='spw',calmode='p',solint='60s',
    spw=SB_contspws,refant=SB_refant,
    minsnr=3.,minblperant=4
)
#10 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:39:52.5
# 2 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:42:07.5
# 1 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/13:07:12.1
# 1 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/13:11:34.6

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
flagdata(vis=SB_p4,mode='manual',spw='0',antenna='DV03,DV07')
flagdata(vis=SB_p4,mode='manual',spw='8',antenna='DV12')
flagdata(vis=SB_p4,mode='manual',spw='8',timerange='2017/08/07/12:39:30~2017/08/07/12:39:40',antenna='DA53')
flagdata(vis=SB_p4,mode='manual',spw='8',timerange='2017/08/07/12:43:20~2017/08/07/12:43:40',antenna='DV11,DA46')
flagdata(vis=SB_p4,mode='manual',spw='8',timerange='2017/08/07/12:46:20~2017/08/07/12:46:40',antenna='DV10,DV07,DA46')
flagdata(vis=SB_p4,mode='manual',spw='8',timerange='2017/08/07/12:47:50~2017/08/07/12:48:00',antenna='DA53')
flagdata(vis=SB_p4,mode='manual',spw='8',timerange='2017/08/07/12:51:10~2017/08/07/12:51:20',antenna='DV11')
flagdata(vis=SB_p4,mode='manual',spw='8',timerange='2017/08/07/12:58:00~2017/08/07/12:58:20',antenna='PM04')
flagdata(vis=SB_p4,mode='manual',spw='8',timerange='2017/08/07/13:02:00~2017/08/07/13:02:20',antenna='DV23')
flagdata(vis=SB_p4,mode='manual',spw='8',timerange='2017/08/07/13:08:50~2017/08/07/13:09:00',antenna='PM03')

##spw='0',timerange='2015/08/30/11:38:15~2015/08/30/11:38:25',antenna=''
# spw='0',timerange='2015/08/30/11:39:15~2015/08/30/11:39:25',antenna='DV07'
# spw='0',timerange='2015/08/30/11:40:15~2015/08/30/11:40:25',antenna='DV07'
# spw='0',timerange='2015/08/30/11:41:15~2015/08/30/11:41:25',antenna='DV03'
# spw='0',timerange='2015/08/30/11:41:50~2015/08/30/11:42:00',antenna='DV03,DV07'
# spw='0',timerange='2015/08/30/11:43:40~2015/08/30/11:43:50',antenna='DV03'
##spw='0',timerange='2015/08/30/11:46:40~2015/08/30/11:46:50',antenna=''
# spw='0',timerange='2015/08/30/11:47:40~2015/08/30/11:47:50',antenna='DV03'
##spw='0',timerange='2015/08/30/11:48:40~2015/08/30/11:48:50',antenna=''
# spw='0',timerange='2015/08/30/11:49:30~2015/08/30/11:49:50',antenna='DV03,DV07'
# spw='0',timerange='2015/08/30/11:54:10~2015/08/30/11:54:30',antenna='DV07'
# spw='0',timerange='2015/08/30/11:55:15~2015/08/30/11:55:30',antenna='DV07'
# spw='0',timerange='2015/08/30/11:56:00~2015/08/30/11:56:20',antenna='DV07'

# spw='8',timerange='2017/08/07/12:39:30~2017/08/07/12:39:40',antenna='DV12,DA53'
# spw='8',timerange='2017/08/07/12:39:50~2017/08/07/12:40:00',antenna='DV12,DV01'
# spw='8',timerange='2017/08/07/12:42:00~2017/08/07/12:42:20',antenna='DV12'
# spw='8',timerange='2017/08/07/12:43:20~2017/08/07/12:43:40',antenna='DV12,DV11,DA46'
# spw='8',timerange='2017/08/07/12:46:20~2017/08/07/12:46:40',antenna='DV12,DV10,DV07,DA46'
# spw='8',timerange='2017/08/07/12:47:50~2017/08/07/12:48:00',antenna='DV12,DA53'
# spw='8',timerange='2017/08/07/13:02:00~2017/08/07/13:02:20',antenna='DV23'
# spw='8',timerange='2017/08/07/13:08:50~2017/08/07/13:09:00',antenna='PM03'

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
    threshold = '0.3300mJy',
    **SB_tclean_wrapper_kwargs
)
estimate_SNR(SB_cont_p4+'.image',disk_mask=SB_mask,noise_mask=noise_annulus_SB)
generate_image_png(
    SB_cont_p4+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_SB,10*rms_SB],
    save_folder=SB_selfcal_folder
)
#CQ_Tau_SB_contp4.image
#Beam 0.138 arcsec x 0.108 arcsec (-2.02 deg)
#Flux inside disk mask: 172.23 mJy
#Peak intensity of source: 7.82 mJy/beam
#rms: 5.49e-02 mJy/beam
#Peak SNR: 142.64

#Fifth round of phase-only selfcal
SB_p5 = SB_p4.replace('p4','p5')
os.system('rm -rf '+SB_p5)
gaincal(
    vis=SB_cont_p4+'.ms',caltable=SB_p5,
    gaintype='T',combine='spw',calmode='p',solint='20s',
    spw=SB_contspws,refant=SB_refant,
    minsnr=3.,minblperant=4
)
#2 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:41:06.4
#1 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:41:50.4
#2 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:42:27.7
#2 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:43:48.3
#2 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:44:49.8
#1 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:46:15.8
#2 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:46:36.0
#2 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:46:53.1
#3 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:48:13.6
#2 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:49:17.0
#4 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:49:34.2
#1 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:51:17.6
#1 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:51:34.8
#1 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:52:18.1
#1 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:52:36.1
#2 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:54:22.8
#2 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:55:07.4
#4 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:55:44.3
#2 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:57:04.9
#1 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:58:25.4
#2 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:59:08.9
#1 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:59:28.8
#2 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/13:01:51.1
#2 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/13:02:11.2
#1 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/13:02:28.4
#2 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/13:03:48.8
#1 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/13:04:32.1
#4 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/13:05:09.3
#1 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/13:05:52.5
#1 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/13:06:12.6
#5 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/13:06:29.7
#1 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/13:07:12.1
#1 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/13:08:36.8
#2 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/13:09:14.0
#1 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/13:10:17.2
#3 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/13:10:34.3
#1 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/13:11:37.5
#2 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/13:13:34.7

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
flagdata(vis=SB_p5,mode='manual',spw='0',antenna='DV03,DV07')
flagdata(vis=SB_p5,mode='manual',spw='8',antenna='DV12')
flagdata(vis=SB_p5,mode='manual',spw='8',timerange='2017/08/07/12:39:20~2017/08/07/12:39:30',antenna='DA53')
flagdata(vis=SB_p5,mode='manual',spw='8',timerange='2017/08/07/12:39:40~2017/08/07/12:39:50',antenna='DA53')
flagdata(vis=SB_p5,mode='manual',spw='8',timerange='2017/08/07/12:40:25~2017/08/07/12:40:35',antenna='DA53,DA46')
flagdata(vis=SB_p5,mode='manual',spw='8',timerange='2017/08/07/12:41:45~2017/08/07/12:41:55',antenna='PM03,DV23')
flagdata(vis=SB_p5,mode='manual',spw='8',timerange='2017/08/07/12:42:05~2017/08/07/12:42:15',antenna='PM03,DV23')
flagdata(vis=SB_p5,mode='manual',spw='8',timerange='2017/08/07/12:43:05~2017/08/07/12:43:15',antenna='DV11,DA46')
flagdata(vis=SB_p5,mode='manual',spw='8',timerange='2017/08/07/12:43:25~2017/08/07/12:43:35',antenna='DV11,DA46')
flagdata(vis=SB_p5,mode='manual',spw='8',timerange='2017/08/07/12:44:45~2017/08/07/12:44:55',antenna='DV10')
flagdata(vis=SB_p5,mode='manual',spw='8',timerange='2017/08/07/12:46:10~2017/08/07/12:46:20',antenna='DV16,DV10,DV07,DA50,DA46')
flagdata(vis=SB_p5,mode='manual',spw='8',timerange='2017/08/07/12:46:30~2017/08/07/12:46:40',antenna='DV10,DV07,DA46')
flagdata(vis=SB_p5,mode='manual',spw='8',timerange='2017/08/07/12:46:50~2017/08/07/12:47:00',antenna='DV10,DV07,DA46')
flagdata(vis=SB_p5,mode='manual',spw='8',timerange='2017/08/07/12:47:30~2017/08/07/12:47:40',antenna='DV07,DA53')
flagdata(vis=SB_p5,mode='manual',spw='8',timerange='2017/08/07/12:47:50~2017/08/07/12:48:00',antenna='DA53')
flagdata(vis=SB_p5,mode='manual',spw='8',timerange='2017/08/07/12:48:10~2017/08/07/12:48:20',antenna='DV07')
flagdata(vis=SB_p5,mode='manual',spw='8',timerange='2017/08/07/12:56:25~2017/08/07/12:56:35',antenna='PM04,DA53')
flagdata(vis=SB_p5,mode='manual',spw='8',timerange='2017/08/07/13:05:50~2017/08/07/13:05:55',antenna='DA44')
flagdata(vis=SB_p5,mode='manual',spw='8',timerange='2017/08/07/13:08:35~2017/08/07/13:08:40',antenna='PM03')
flagdata(vis=SB_p5,mode='manual',spw='8',timerange='2017/08/07/13:08:50~2017/08/07/13:09:00',antenna='PM03')
flagdata(vis=SB_p5,mode='manual',spw='8',timerange='2017/08/07/13:09:55~2017/08/07/13:10:00',antenna='PM03')
flagdata(vis=SB_p5,mode='manual',spw='8',timerange='2017/08/07/13:11:15~2017/08/07/13:11:20',antenna='DA46')
flagdata(vis=SB_p5,mode='manual',spw='8',timerange='2017/08/07/13:11:35~2017/08/07/13:11:40',antenna='DA46')

# spw='8',timerange='2017/08/07/12:39:40~2017/08/07/12:39:45',antenna='DV12'
# spw='8',timerange='2017/08/07/12:39:50~2017/08/07/12:39:55',antenna='DV12'
# spw='8',timerange='2017/08/07/12:40:25~2017/08/07/12:40:35',antenna='DA53,DA46'
# spw='8',timerange='2017/08/07/12:41:45~2017/08/07/12:41:55',antenna='PM03,DV23'
# spw='8',timerange='2017/08/07/12:42:05~2017/08/07/12:42:15',antenna='PM03,DV23,DV12'
# spw='8',timerange='2017/08/07/12:42:25~2017/08/07/12:42:35',antenna='DV12'
# spw='8',timerange='2017/08/07/12:43:05~2017/08/07/12:43:15',antenna='DV12,DV11,DA46'
# spw='8',timerange='2017/08/07/12:43:25~2017/08/07/12:43:35',antenna='DV12,DV11,DA46'
# spw='8',timerange='2017/08/07/12:44:45~2017/08/07/12:44:55',antenna='DV12,DV10'
# spw='8',timerange='2017/08/07/12:46:10~2017/08/07/12:46:20',antenna='DV16,DV12,DV10,DV07,DA50,DA46'
# spw='8',timerange='2017/08/07/12:46:30~2017/08/07/12:46:40',antenna='DV12,DV10,DV07,DA46'
# spw='8',timerange='2017/08/07/12:46:50~2017/08/07/12:47:00',antenna='DV12,DV10,DV07,DA46'
# spw='8',timerange='2017/08/07/12:47:30~2017/08/07/12:47:40',antenna='DV12,DV07,DA53'
# spw='8',timerange='2017/08/07/12:47:50~2017/08/07/12:48:00',antenna='DA53'
# spw='8',timerange='2017/08/07/12:48:10~2017/08/07/12:48:20',antenna='DV07'
# spw='8',timerange='2017/08/07/13:05:50~2017/08/07/13:05:55',antenna='DA44'
# spw='8',timerange='2017/08/07/13:07:10~2017/08/07/13:07:15',antenna='DV12'
# spw='8',timerange='2017/08/07/13:08:35~2017/08/07/13:08:40',antenna='PM03'
# spw='8',timerange='2017/08/07/13:08:50~2017/08/07/13:09:00',antenna='PM03'
# spw='8',timerange='2017/08/07/13:09:55~2017/08/07/13:10:00',antenna='PM03'
# spw='8',timerange='2017/08/07/13:11:15~2017/08/07/13:11:20',antenna='DA46'
# spw='8',timerange='2017/08/07/13:11:35~2017/08/07/13:11:40',antenna='DA46'

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
    threshold = '0.3198mJy',
    **SB_tclean_wrapper_kwargs
)
estimate_SNR(SB_cont_p5+'.image',disk_mask=SB_mask,noise_mask=noise_annulus_SB)
generate_image_png(
    SB_cont_p5+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_SB,10*rms_SB],
    save_folder=SB_selfcal_folder
)
#CQ_Tau_SB_contp5.image
#Beam 0.138 arcsec x 0.108 arcsec (-2.02 deg)
#Flux inside disk mask: 172.83 mJy
#Peak intensity of source: 8.05 mJy/beam
#rms: 5.32e-02 mJy/beam
#Peak SNR: 151.32

#For each step of the self-cal, check how it improved things
non_self_caled_SB_vis = SB_cont_p0
self_caled_SB_visibilities = {
    'p1':SB_cont_p1,
    'p2':SB_cont_p2,
    'p3':SB_cont_p3,
    'p4':SB_cont_p4,
    'p5':SB_cont_p5
}

#SB_EBs = ('EB0','EB1')
#SB_EB_spws = ('0,1,2,3,4,5,6,7','8,9,10,11,12,13,14,15') #fill out by referring to listobs output
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

for self_cal_step,vis_name in all_SB_visibilities.items():
    #Split out SB EBs
    vis_ms = vis_name+'.ms'
    nametemplate = vis_ms.replace('.ms','_EB')
    split_all_obs(msfile=vis_ms,nametemplate=nametemplate)
    #Saving observation 0 of CQ_Tau_SB_contp0.ms to CQ_Tau_SB_contp0_EB0.ms
    #Saving observation 1 of CQ_Tau_SB_contp0.ms to CQ_Tau_SB_contp0_EB1.ms
    #...
    #Saving observation 0 of CQ_Tau_SB_contp5.ms to CQ_Tau_SB_contp5_EB0.ms
    #Saving observation 1 of CQ_Tau_SB_contp5.ms to CQ_Tau_SB_contp5_EB1.ms

    exported_ms = []
    for i in range(number_of_EBs['SB']):
        EB_vis = f'{nametemplate}{i}.ms'
        listobs(vis=EB_vis,listfile=EB_vis+'.listobs.txt',overwrite=True)
        #Export MS contents into numpy save files
        export_MS(EB_vis)
        exported_ms.append(EB_vis.replace('.ms','.vis.npz'))
        #Measurement set exported to CQ_Tau_SB_contp0_EB0.vis.npz
        #Measurement set exported to CQ_Tau_SB_contp0_EB1.vis.npz
        #...
        #Measurement set exported to CQ_Tau_SB_contp5_EB0.vis.npz
        #Measurement set exported to CQ_Tau_SB_contp5_EB1.vis.npz

    for i,exp_ms in enumerate(exported_ms):
        png_filename = f'flux_comparison_SB_EB{i}_{self_cal_step}_to_{flux_ref_EB}.png'
        plot_label = os.path.join(SB_selfcal_folder,png_filename)
        estimate_flux_scale(
            reference=f'{prefix}_{flux_ref_EB}_initcont_shift.vis.npz',
            comparison=exp_ms,incl=incl,PA=PA,plot_label=plot_label
        )

    fluxscale = [1.,]*number_of_EBs['SB']
    plot_label = os.path.join(SB_selfcal_folder,f'deprojected_vis_profiles_SB_{self_cal_step}.png')
    plot_deprojected(filelist=exported_ms,fluxscale=fluxscale,PA=PA,incl=incl,show_err=True,plot_label=plot_label)

#ratio          = [1.36589,1.36531,1.37131,1.37702,1.37817,1.37872] #CQ_Tau_SB_contp0...p5_EB0.vis.npz vs CQ_Tau_LB_EB1_initcont_shift.vis.npz
#scaling_factor = [1.169  ,1.168  ,1.171  ,1.173  ,1.174  ,1.174  ]
#ratio          = [0.91086,0.91196,0.94078,0.97047,0.99275,1.02373] #CQ_Tau_SB_contp0...p5_EB1.vis.npz vs CQ_Tau_LB_EB1_initcont_shift.vis.npz
#scaling_factor = [0.954  ,0.955  ,0.970  ,0.985  ,0.996  ,1.012  ]

# SB_EB0 shows flux offsets >4% and no decoherence -> correct for flux-scaling and re-run SB selfcal
# SB_EB1 shows flux offsets <4% and the decoherence improves up to p3 -> no flux-scaling and re-run SB selfcal

#Begin of SB+LB self-cal - iteration 1
#For phase self-cal, clean down to ~6sigma
LB_selfcal_folder = get_figures_folderpath('8_selfcal_SBLB_figures')
make_figures_folder(LB_selfcal_folder)

LB_cont_p0 = prefix+'_SBLB_contp0'
os.system('rm -rf %s.ms*' %LB_cont_p0)
concat(
    vis=[SB_cont_p5+'.ms']+[f'{prefix}_LB_EB{i}_initcont_shift.ms' for i in range(number_of_EBs['LB'])],
    concatvis=LB_cont_p0+'.ms',dirtol='0.1arcsec',copypointing=False
)
listobs(vis=LB_cont_p0+'.ms',listfile=LB_cont_p0+'.ms.listobs.txt',overwrite=True)
#2024-12-10 21:26:44     WARN    MSConcat::copySysCal    
#    /data/beegfs/astro-storage/groups/benisty/frzagaria/SO_detections/CQTau/selfcal_products/CQ_Tau_SB_contp0.ms does not have a valid syscal table,
#    the MS to be appended, however, has one. Result won't have one.
#2024-12-10 21:26:44     WARN    MSConcat::concatenate
#    (file /source/casa6/casatools/casacore/ms/MSOper/MSConcat.cc, line 1000). Could not merge SysCal subtables.
#2024-12-13 21:49:15     SEVERE  getcell::TIME   Exception Reported: TableProxy::getCell: no such row
#2024-12-13 21:49:16     WARN    concat::::casa  Some but not all of the input MSs are lacking a populated POINTING table:
#    0: CQ_Tau_SB_contp5.ms The joint dataset will not have a valid POINTING table.
#2024-12-13 21:50:00     WARN    MSConcat::copySysCal    /data/beegfs/astro-storage/groups/benisty/frzagaria/SO_detections/CQTau/selfcal_products/CQ_Tau_SBLB_contp0.ms does not have a valid syscal table,
#    the MS to be appended, however, has one. Result won't have one.
#2024-12-13 21:50:00     WARN    MSConcat::concatenate (file /source/casa6/casatools/casacore/ms/MSOper/MSConcat.cc, line 1000)        Could not merge SysCal subtables 
#    /data/beegfs/astro-storage/groups/benisty/frzagaria/SO_detections/CQTau/selfcal_products/CQ_Tau_SBLB_contp0.ms does not have a valid syscal table,
#    the MS to be appended, however, has one. Result won't have one.
#2024-12-13 21:53:06     WARN    MSConcat::concatenate (file /source/casa6/casatools/casacore/ms/MSOper/MSConcat.cc, line 1000)        Could not merge SysCal subtables

#Define new SB mask using new center read from listobs
mask_pa        = PA
mask_semimajor = 0.85
mask_semiminor = 0.72
mask_ra  = '05h35m58.472722s'
mask_dec = '24d44m53.613450s'

LB_mask = f'ellipse[[{mask_ra},{mask_dec}], [{mask_semimajor:.3f}arcsec, {mask_semiminor:.3f}arcsec], {mask_pa:.1f}deg]'

noise_annulus_LB = f"annulus[[{mask_ra}, {mask_dec}],['4.arcsec', '6.arcsec']]"

LB_tclean_wrapper_kwargs = {
    'deconvolver':'multiscale','scales':[0,2,4,8,16,24],
    'smallscalebias':0.6,'gain':0.3,'cycleniter':300,
    'niter':1000000,'mask':LB_mask,
    'cellsize':'0.006arcsec','imsize':6400,
    'parallel':use_parallel,'savemodel':'modelcolumn',
    'robust':0.5,'interactive':False,
    'gridder':'standard',
}

tclean_wrapper(
    vis       = LB_cont_p0+'.ms', 
    imagename = LB_cont_p0,
    threshold = '0.1056mJy',
    **LB_tclean_wrapper_kwargs
)
estimate_SNR(LB_cont_p0+'.image',disk_mask=LB_mask,noise_mask=noise_annulus_LB)
rms_LB = imstat(imagename=LB_cont_p0+'.image',region=noise_annulus_LB)['rms'][0]
generate_image_png(
    LB_cont_p0+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_LB,10*rms_LB],
    save_folder=LB_selfcal_folder
)
#CQ_Tau_SBLB_contp0.image
#Beam 0.078 arcsec x 0.058 arcsec (7.28 deg)
#Flux inside disk mask: 160.05 mJy
#Peak intensity of source: 2.43 mJy/beam
#rms: 1.75e-02 mJy/beam
#Peak SNR: 138.88

#Skip: no weblogs for the SBs (too old datasets), don't know how to get the reference antennas
#
#Look for references antennas from the weblog, and pick the first ones that are listed, overlapping with all EBs
#Since we are combined SB EBs to LB EBs, ref antennas for both SB and LB are needed, but SB has no recorded refants...
#LB_EB0: DA43, DA57, DV04, DV25, DA42, DA58
#LB_EB1: DA43, DV06, DA60, DA58, DV04, DV01
#
#Get station numbers
#For SB we used:'DA63@A035, DA52@A002, DA57@A001, DV08@A036'
#Bear in mind that here ID0=LB0,ID1=LB1,ID2=LB2,ID3,LB3,ID4=SB0,ID5=SB1
#for ref_ant in ('DA43','DV06','DA57', 'DV04', 'DA63','DA52'):
#    get_station_numbers(LB_cont_p0+'.ms',ref_ant)
#    Observation ID 0: DA43@A035
#    Observation ID 1: DA43@A035
#    Observation ID 2: DA43@A035
#    Observation ID 3: DA43@A073
#    Observation ID 4: DA43@A073
#    Observation ID 5: DA43@A073
#    Observation ID 0: DV06@A024
#    Observation ID 1: DV06@A024
#    Observation ID 2: DV06@A024
#    Observation ID 3: DV06@A024
#    Observation ID 4: DV06@A024
#    Observation ID 5: DV06@A024
#    Observation ID 0: DA57@A043
#    Observation ID 1: DA57@A043
#    Observation ID 2: DA57@A043
#    Observation ID 3: DA57@A104
#    Observation ID 5: DA57@A001
#    Observation ID 0: DV04@A007
#    Observation ID 1: DV04@A007
#    Observation ID 2: DV04@A007
#    Observation ID 3: DV04@A007
#    Observation ID 0: DA63@A123
#    Observation ID 1: DA63@A123
#    Observation ID 2: DA63@A123
#    Observation ID 3: DA63@A116
#    Observation ID 4: DA63@A035
#    Observation ID 5: DA63@A035
#    Observation ID 1: DA52@A089
#    Observation ID 2: DA52@A089
#    Observation ID 4: DA52@A002
#    Observation ID 5: DA52@A002

LB_refant = ''

LB_contspws    = '0~31'
LB_spw_mapping = [
    0,0,0,0,0,0,0,0,
    8,8,8,8,8,8,8,8,
    16,16,16,16,16,16,16,16,
    24,24,24,24,24,24,24,24,
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
    minsnr=2.,minblperant=4
)
#1 of 88 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:42:17.8
#2 of 94 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:08:25.6

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
#flagdata(vis=LB_p1,mode='manual',spw='8',antenna='DV12')
#
#Inspect gain tables to check if flagging worked
#plotms(
#    LB_p1,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
#    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
#    plotfile=os.path.join(LB_selfcal_folder,f'{prefix}_LB_p1_gain_phase_vs_time_flagged.png'),
#)

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
    threshold = '0.1056mJy',
    **LB_tclean_wrapper_kwargs
)
estimate_SNR(LB_cont_p1+'.image',disk_mask=LB_mask,noise_mask=noise_annulus_LB)
generate_image_png(
    LB_cont_p1+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_LB,10*rms_LB],
    save_folder=LB_selfcal_folder
)
#CQ_Tau_SBLB_contp1.image
#Beam 0.078 arcsec x 0.058 arcsec (7.28 deg)
#Flux inside disk mask: 160.00 mJy
#Peak intensity of source: 2.45 mJy/beam
#rms: 1.75e-02 mJy/beam
#Peak SNR: 140.01

#Second round of phase-only self-cal
LB_p2 = LB_p1.replace('p1','p2')
os.system('rm -rf '+LB_p2)
#Pietro: solint='360s' resulted in "mismatched frequencies" error, so I slighlty changed it
gaincal(
    vis=LB_cont_p1+'.ms',caltable=LB_p2,
    gaintype='T',combine='scan,spw',calmode='p',solint='360s',
    spw=LB_contspws,refant=LB_refant,
    minsnr=2.,minblperant=4
)
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:17:38.3
#3 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:21:20.5
#7 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:24:31.1
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:28:26.1
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:31:12.7
#1 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:34:14.5
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:41:45.8
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:47:57.9
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:51:13.6
#1 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:54:42.4
#7 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:58:06.6
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:03:06.9
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:06:37.1
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:44:08.4
#2 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:47:49.5
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:51:04.5
#1 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:55:00.4
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:57:51.7
#2 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:00:58.7
#3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:08:45.3
#3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:14:48.9
#3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:18:10.4
#2 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:21:43.0
#7 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:25:09.5
#2 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:30:11.7
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:33:44.5

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
flagdata(vis=LB_p2,mode='manual',spw='0', timerange='2015/08/30/11:41:00~2015/08/30/11:43:00',antenna='DV03')
flagdata(vis=LB_p2,mode='manual',spw='0', timerange='2015/08/30/11:48:30~2015/08/30/11:50:30',antenna='DV03,DV07')
flagdata(vis=LB_p2,mode='manual',spw='0', timerange='2015/08/30/11:54:00~2015/08/30/11:56:00',antenna='DV07')
flagdata(vis=LB_p2,mode='manual',spw='8', timerange='2017/08/07/12:39:00~2017/08/07/12:40:30',antenna='DV12,DA53')
flagdata(vis=LB_p2,mode='manual',spw='8', timerange='2017/08/07/12:42:00~2017/08/07/12:43:30',antenna='DV12')
flagdata(vis=LB_p2,mode='manual',spw='8', timerange='2017/08/07/12:45:30~2017/08/07/12:48:00',antenna='DV12,DV10,DV07,DA53,DA46')
flagdata(vis=LB_p2,mode='manual',spw='16',timerange='2017/11/20/07:31:00~2017/11/20/07:31:30',antenna='PM04,DV10,DA44')
flagdata(vis=LB_p2,mode='manual',spw='16',timerange='2017/11/20/07:58:00~2017/11/20/07:58:30',antenna='DV13,DA42')
flagdata(vis=LB_p2,mode='manual',spw='16',timerange='2017/11/20/08:06:00~2017/11/20/08:07:00',antenna='DA65')
flagdata(vis=LB_p2,mode='manual',spw='24',timerange='2017/11/23/04:57:30~2017/11/23/04:58:00',antenna='DV05')
flagdata(vis=LB_p2,mode='manual',spw='24',timerange='2017/11/23/05:08:30~2017/11/23/05:09:00',antenna='DV16,DA65')
flagdata(vis=LB_p2,mode='manual',spw='24',timerange='2017/11/23/05:18:00~2017/11/23/05:18:30',antenna='DV16')
flagdata(vis=LB_p2,mode='manual',spw='24',timerange='2017/11/23/05:30:00~2017/11/23/05:30:30',antenna='DV17')

#Inspect gain tables to check if flagging worked
plotms(
    LB_p2,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(LB_selfcal_folder,f'{prefix}_LB_p2_gain_phase_vs_time_flagged.png'),
)

#Apply the calibration gains and split-off the corrected ms
applycal(
    vis=LB_cont_p1+'.ms',spw=LB_contspws,spwmap=LB_spw_mapping,
    gaintable=[LB_p2],interp='linearPD',calwt=True,applymode='calonly'
)
LB_cont_p2 = LB_cont_p1.replace('p1','p2')
os.system('rm -rf %s.ms*' %LB_cont_p2)
split(vis=LB_cont_p1+'.ms',outputvis=LB_cont_p2+'.ms',datacolumn='corrected')

#CLEAN, again with a ~6sigma threshold
tclean_wrapper(
    vis       = LB_cont_p2+'.ms',
    imagename = LB_cont_p2,
    threshold = '0.1008mJy',
    **LB_tclean_wrapper_kwargs
)
estimate_SNR(LB_cont_p2+'.image',disk_mask=LB_mask,noise_mask=noise_annulus_LB)
generate_image_png(
    LB_cont_p2+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_LB,10*rms_LB],
    save_folder=LB_selfcal_folder
)
#CQ_Tau_SBLB_contp2.image
#Beam 0.078 arcsec x 0.058 arcsec (7.28 deg)
#Flux inside disk mask: 160.14 mJy
#Peak intensity of source: 2.49 mJy/beam
#rms: 1.68e-02 mJy/beam
#Peak SNR: 148.36

#Third round of phase-only self-cal
LB_p3 = LB_p2.replace('p2','p3')
os.system('rm -rf '+LB_p3)
gaincal(
    vis=LB_cont_p2+'.ms',caltable=LB_p3,
    gaintype='T',combine='scan,spw',calmode='p',solint='120s', #diff from exoALMA, kept combine='scan' because scans for SB_EB1~60s
    spw=LB_contspws,refant=LB_refant,
    minsnr=2.,minblperant=4
)
#8 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/12:39:52.5
#1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/12:41:18.3
#1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/12:48:13.6
#1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/12:58:25.4
#1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/13:03:48.8
#1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/13:06:29.7
#1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/13:09:12.0
#1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/13:10:45.5
#3 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:17:38.3
#3 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:19:21.7
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:20:30.4
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:22:02.2
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:23:11.5
#3 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:24:40.6
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:26:03.7
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:29:02.4
#2 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:31:23.1
#7 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:32:32.2
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:34:04.1
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:35:13.3
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:36:18.1
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:39:43.7
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:42:04.3
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:43:13.4
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:44:45.4
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:45:54.5
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:46:59.6
#3 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:50:23.9
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:52:44.6
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:53:53.7
#3 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:55:25.4
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:56:34.8
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:57:39.6
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:01:04.6
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:03:25.0
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:04:34.3
#3 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:05:50.4
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:44:08.4
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:45:59.3
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:47:02.6
#3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:48:42.4
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:49:45.9
#2 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:51:16.4
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:52:41.6
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:55:43.3
#3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:58:05.5
#7 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:59:15.5
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:00:48.6
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:01:58.5
#3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:03:04.6
#3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:06:33.8
#2 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:08:55.6
#7 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:10:05.6
#3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:11:38.4
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:12:48.4
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:13:54.4
#3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:17:22.3
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:19:44.2
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:20:54.0
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:22:27.1
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:23:36.6
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:24:42.3
#2 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:28:08.9
#2 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:30:30.7
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:31:40.4
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:32:57.4

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
#flagdata(vis=LB_p3,mode='clip',spw='8', clipminmax=[-50,60],  clipoutside=True,datacolumn='CPARAM')??
#flagdata(vis=LB_p3,mode='clip',spw='16',clipminmax=[-100,150],clipoutside=True,datacolumn='CPARAM',timerange='2017/11/20/07:48:00~2017/11/20/08:06:00')??

flagdata(vis=LB_p3,mode='manual',spw='0', timerange='2015/08/30/11:41:30~2015/08/30/11:42:30',antenna='DV03')
flagdata(vis=LB_p3,mode='manual',spw='0', timerange='2015/08/30/11:49:30~2015/08/30/11:51:00',antenna='DV03,DV07')
flagdata(vis=LB_p3,mode='manual',spw='0', timerange='2015/08/30/11:54:00~2015/08/30/11:55:30',antenna='DV07')
flagdata(vis=LB_p3,mode='manual',spw='0', timerange='2015/08/30/11:55:30~2015/08/30/11:57:00',antenna='DV07')

flagdata(vis=LB_p3,mode='manual',spw='8', timerange='2017/08/07/12:39:30~2017/08/07/12:39:40',antenna='DV12,DA53')
flagdata(vis=LB_p3,mode='manual',spw='8', timerange='2017/08/07/12:42:20~2017/08/07/12:42:40',antenna='DV11')
flagdata(vis=LB_p3,mode='manual',spw='8', timerange='2017/08/07/13:03:40~2017/08/07/13:04:00',antenna='PM03')
flagdata(vis=LB_p3,mode='manual',spw='8', timerange='2017/08/07/13:10:40~2017/08/07/13:11:00',antenna='PM03')
flagdata(vis=LB_p3,mode='manual',spw='8', timerange='2017/08/07/12:39:50~2017/08/07/12:40:00',antenna='DV12')
flagdata(vis=LB_p3,mode='manual',spw='8', timerange='2017/08/07/12:42:20~2017/08/07/12:42:30',antenna='PM03,DV12,DA46')
flagdata(vis=LB_p3,mode='manual',spw='8', timerange='2017/08/07/12:43:50~2017/08/07/12:44:00',antenna='DV12')
flagdata(vis=LB_p3,mode='manual',spw='8', timerange='2017/08/07/12:47:00~2017/08/07/12:47:10',antenna='DV12,DV10,DV07,DA53,DA46')
flagdata(vis=LB_p3,mode='manual',spw='8', timerange='2017/08/07/12:48:00~2017/08/07/12:48:20',antenna='DV23,DV07')
flagdata(vis=LB_p3,mode='manual',spw='8', timerange='2017/08/07/13:06:00~2017/08/07/13:07:00',antenna='PM04,DA53,DA46')
flagdata(vis=LB_p3,mode='manual',spw='8', timerange='2017/08/07/13:08:00~2017/08/07/13:08:30',antenna='PM03')

flagdata(vis=LB_p3,mode='manual',spw='16',timerange='2017/11/20/07:17:30~2017/11/20/07:18:00',antenna='DV17,DV05')
flagdata(vis=LB_p3,mode='manual',spw='16',timerange='2017/11/20/07:28:40~2017/11/20/07:29:20',antenna='DA64,DA65')
flagdata(vis=LB_p3,mode='manual',spw='16',timerange='2017/11/20/07:39:40~2017/11/20/07:39:50',antenna='DA65')
flagdata(vis=LB_p3,mode='manual',spw='16',timerange='2017/11/20/07:50:20~2017/11/20/07:50:30',antenna='DV22')
flagdata(vis=LB_p3,mode='manual',spw='16',timerange='2017/11/20/07:56:30~2017/11/20/07:57:00',antenna='DA65')
flagdata(vis=LB_p3,mode='manual',spw='16',timerange='2017/11/20/07:57:30~2017/11/20/07:58:00',antenna='DV13,DA42')
flagdata(vis=LB_p3,mode='manual',spw='16',timerange='2017/11/20/08:03:00~2017/11/20/08:04:00',antenna='PM03,DA42')

flagdata(vis=LB_p3,mode='manual',spw='24',timerange='2017/11/23/04:47:00~2017/11/23/04:47:10',antenna='DA64,DA65')
flagdata(vis=LB_p3,mode='manual',spw='24',timerange='2017/11/23/04:48:40~2017/11/23/04:48:50',antenna='DV22')
flagdata(vis=LB_p3,mode='manual',spw='24',timerange='2017/11/23/04:51:10~2017/11/23/04:51:20',antenna='DV16,DA64,DA60,PM04')
flagdata(vis=LB_p3,mode='manual',spw='24',timerange='2017/11/23/04:52:40~2017/11/23/04:52:50',antenna='DA46')
flagdata(vis=LB_p3,mode='manual',spw='24',timerange='2017/11/23/04:58:00~2017/11/23/04:58:10',antenna='DV17,DV22')
flagdata(vis=LB_p3,mode='manual',spw='24',timerange='2017/11/23/05:01:50~2017/11/23/05:02:10',antenna='DV05')
flagdata(vis=LB_p3,mode='manual',spw='24',timerange='2017/11/23/05:12:40~2017/11/23/05:13:00',antenna='DV17')
flagdata(vis=LB_p3,mode='manual',spw='24',timerange='2017/11/23/05:20:50~2017/11/23/05:21:00',antenna='DA65,PM03')
flagdata(vis=LB_p3,mode='manual',spw='24',timerange='2017/11/23/05:28:00~2017/11/23/05:28:20',antenna='DV17,DV05')

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
    threshold = '0.0930mJy',
    **LB_tclean_wrapper_kwargs
)
estimate_SNR(LB_cont_p3+'.image',disk_mask=LB_mask,noise_mask=noise_annulus_LB)
generate_image_png(
    LB_cont_p3+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_LB,10*rms_LB],
    save_folder=LB_selfcal_folder
)
#CQ_Tau_SBLB_contp3.image
#Beam 0.078 arcsec x 0.058 arcsec (7.28 deg)
#Flux inside disk mask: 160.64 mJy
#Peak intensity of source: 2.60 mJy/beam
#rms: 1.55e-02 mJy/beam
#Peak SNR: 167.30

#Fourth round of phase-only self-cal
LB_p4 = LB_p3.replace('p3','p4')
os.system('rm -rf '+LB_p4)
gaincal(
    vis=LB_cont_p3+'.ms',caltable=LB_p4,
    gaintype='T',combine='scan,spw',calmode='p',solint='60s', #diff from exoALMA, kept combine='scan' because scans for LB_EB0~30s
    spw=LB_contspws,refant=LB_refant,
    minsnr=2.,minblperant=4
)
#9 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/12:39:52.5
#1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/12:47:53.9
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:17:38.3
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:18:50.4
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:20:10.6
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:21:31.1
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:22:51.8
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:24:12.0
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:25:32.8
#7 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:26:09.7
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:28:19.5
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:29:31.1
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:30:51.7
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:32:12.2
#7 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:33:32.8
#3 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:34:53.3
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:36:13.9
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:36:50.8
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:39:00.7
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:40:12.3
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:41:32.7
#3 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:42:53.4
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:44:14.1
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:45:34.7
#2 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:46:55.5
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:47:32.1
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:49:41.0
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:50:53.0
#3 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:52:13.3
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:53:33.9
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:54:54.2
#3 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:56:14.8
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:57:35.5
#8 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:58:12.6
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:00:21.6
#8 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:01:33.1
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:02:53.7
#3 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:04:14.2
#3 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:05:34.8
#7 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:06:37.1
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:44:08.4
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:45:21.0
#3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:46:42.4
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:48:04.2
#3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:49:25.8
#3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:50:47.4
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:52:08.1
#7 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:52:46.6
#3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:54:59.9
#3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:56:12.3
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:57:33.5
#3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:58:55.3
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:00:16.9
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:01:38.3
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:02:58.9
#8 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:03:37.3
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:05:50.0
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:07:02.6
#3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:08:24.0
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:09:45.4
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:11:06.8
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:12:28.2
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:13:48.7
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:14:27.1
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:16:39.0
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:17:51.2
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:19:11.6
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:20:33.8
#1 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:21:55.4
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:23:16.5
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:24:43.4
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:25:15.6
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:27:25.7
#3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:28:37.8
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:29:59.1
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:31:20.3
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:32:41.6
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:33:44.5

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
flagdata(vis=LB_p4,mode='manual',spw='0', timerange='2015/08/30/11:41:50~2015/08/30/11:42:00',antenna='DV03')
flagdata(vis=LB_p4,mode='manual',spw='0', timerange='2015/08/30/11:49:00~2015/08/30/11:50:00',antenna='DV03,DV07')
flagdata(vis=LB_p4,mode='manual',spw='0', timerange='2015/08/30/11:55:00~2015/08/30/11:55:30',antenna='DV07')
flagdata(vis=LB_p4,mode='manual',spw='0', timerange='2015/08/30/11:56:00~2015/08/30/11:56:30',antenna='DV07')

flagdata(vis=LB_p4,mode='manual',spw='8', timerange='2017/08/07/12:39:50~2017/08/07/12:40:00',antenna='DV12')
flagdata(vis=LB_p4,mode='manual',spw='8', timerange='2017/08/07/12:42:00~2017/08/07/12:42:10',antenna='DV12,PM03,DV23')
flagdata(vis=LB_p4,mode='manual',spw='8', timerange='2017/08/07/12:43:20~2017/08/07/12:43:40',antenna='DV12,DV11,DA46')
flagdata(vis=LB_p4,mode='manual',spw='8', timerange='2017/08/07/12:46:00~2017/08/07/12:47:00',antenna='DV12,DV10,DV07,DA46')
flagdata(vis=LB_p4,mode='manual',spw='8', timerange='2017/08/07/12:47:30~2017/08/07/12:48:30',antenna='DV07,DA53')
flagdata(vis=LB_p4,mode='manual',spw='8', timerange='2017/08/07/13:07:00~2017/08/07/13:08:00',antenna='PM04,DV12')
flagdata(vis=LB_p4,mode='manual',spw='8', timerange='2017/08/07/13:08:00~2017/08/07/13:09:00',antenna='PM03')
flagdata(vis=LB_p4,mode='manual',spw='8', timerange='2017/08/07/13:10:00~2017/08/07/13:11:00',antenna='PM03,DV12')
flagdata(vis=LB_p4,mode='manual',spw='8', timerange='2017/08/07/13:11:00~2017/08/07/13:12:00',antenna='DV12,DA46')

flagdata(vis=LB_p4,mode='manual',spw='16',timerange='2017/11/20/07:17:30~2017/11/20/07:17:40',antenna='DV05')
flagdata(vis=LB_p4,mode='manual',spw='16',timerange='2017/11/20/07:26:00~2017/11/20/07:26:20',antenna='DV17')
flagdata(vis=LB_p4,mode='manual',spw='16',timerange='2017/11/20/07:36:40~2017/11/20/07:37:00',antenna='DV05')
flagdata(vis=LB_p4,mode='manual',spw='16',timerange='2017/11/20/07:38:50~2017/11/20/07:39:10',antenna='DA65')
flagdata(vis=LB_p4,mode='manual',spw='16',timerange='2017/11/20/07:42:30~2017/11/20/07:43:30',antenna='DA60')
flagdata(vis=LB_p4,mode='manual',spw='16',timerange='2017/11/20/07:42:30~2017/11/20/07:43:00',antenna='DA60')
flagdata(vis=LB_p4,mode='manual',spw='16',timerange='2017/11/20/07:44:00~2017/11/20/07:44:30',antenna='DA65')
flagdata(vis=LB_p4,mode='manual',spw='16',timerange='2017/11/20/07:46:50~2017/11/20/07:47:00',antenna='DV17')
flagdata(vis=LB_p4,mode='manual',spw='16',timerange='2017/11/20/07:52:10~2017/11/20/07:52:20',antenna='DV17,DA65')
flagdata(vis=LB_p4,mode='manual',spw='16',timerange='2017/11/20/07:56:00~2017/11/20/07:56:30',antenna='DV17,DA65')
flagdata(vis=LB_p4,mode='manual',spw='16',timerange='2017/11/20/07:58:00~2017/11/20/07:58:30',antenna='DA64,DV23,DV13,DA42')
flagdata(vis=LB_p4,mode='manual',spw='16',timerange='2017/11/20/08:01:00~2017/11/20/08:02:00',antenna='DV22,PM03')
flagdata(vis=LB_p4,mode='manual',spw='16',timerange='2017/11/20/08:02:30~2017/11/20/08:03:00',antenna='PM03,DA42')
flagdata(vis=LB_p4,mode='manual',spw='16',timerange='2017/11/20/08:04:00~2017/11/20/08:04:30',antenna='DA60')

flagdata(vis=LB_p4,mode='manual',spw='24',timerange='2017/11/23/04:52:40~2017/11/23/04:53:00',antenna='DA46')
flagdata(vis=LB_p4,mode='manual',spw='24',timerange='2017/11/23/05:02:50~2017/11/23/05:03:00',antenna='DV16')
flagdata(vis=LB_p4,mode='manual',spw='24',timerange='2017/11/23/05:14:20~2017/11/23/05:14:30',antenna='DV16')
flagdata(vis=LB_p4,mode='manual',spw='24',timerange='2017/11/23/05:16:30~2017/11/23/05:16:40',antenna='DV22')
flagdata(vis=LB_p4,mode='manual',spw='24',timerange='2017/11/23/05:19:10~2017/11/23/05:19:20',antenna='DV16,DV22')
flagdata(vis=LB_p4,mode='manual',spw='24',timerange='2017/11/23/05:20:30~2017/11/23/05:20:40',antenna='DA65,DA60,PM03')
flagdata(vis=LB_p4,mode='manual',spw='24',timerange='2017/11/23/05:24:00~2017/11/23/05:25:00',antenna='DA60')
flagdata(vis=LB_p4,mode='manual',spw='24',timerange='2017/11/23/05:25:00~2017/11/23/05:26:00',antenna='DV16')
flagdata(vis=LB_p4,mode='manual',spw='24',timerange='2017/11/23/05:27:00~2017/11/23/05:27:30',antenna='DV17,DV22')
flagdata(vis=LB_p4,mode='manual',spw='24',timerange='2017/11/23/05:28:30~2017/11/23/05:28:40',antenna='DV17,DV05')
flagdata(vis=LB_p4,mode='manual',spw='24',timerange='2017/11/23/05:33:40~2017/11/23/05:33:50',antenna='DV16')

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
    threshold = '0.0930mJy',
    **LB_tclean_wrapper_kwargs
)
estimate_SNR(LB_cont_p4+'.image',disk_mask=LB_mask,noise_mask=noise_annulus_LB)
generate_image_png(
    LB_cont_p4+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_LB,10*rms_LB],
    save_folder=LB_selfcal_folder
)
#CQ_Tau_SBLB_contp4.image
#Beam 0.078 arcsec x 0.058 arcsec (7.28 deg)
#Flux inside disk mask: 160.79 mJy
#Peak intensity of source: 2.67 mJy/beam
#rms: 1.51e-02 mJy/beam
#Peak SNR: 176.52

#Fifth round of phase-only self-cal
LB_p5 = LB_p4.replace('p4','p5')
os.system('rm -rf '+LB_p5)
gaincal(
    vis=LB_cont_p4+'.ms',caltable=LB_p5,
    gaintype='T',combine='spw',calmode='p',solint='30s',
    spw=LB_contspws,refant=LB_refant,
    minsnr=2.,minblperant=4
)

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
flagdata(vis=LB_p5,mode='manual',spw='0',antenna='DV07')

flagdata(vis=LB_p5,mode='manual',spw='0', timerange='2015/08/30/11:39:50~2015/08/30/11:40:10',antenna='DV03')
flagdata(vis=LB_p5,mode='manual',spw='0', timerange='2015/08/30/11:41:30~2015/08/30/11:41:40',antenna='DV03')
flagdata(vis=LB_p5,mode='manual',spw='0', timerange='2015/08/30/11:41:50~2015/08/30/11:42:00',antenna='DV03')
flagdata(vis=LB_p5,mode='manual',spw='0', timerange='2015/08/30/11:43:20~2015/08/30/11:43:30',antenna='DV03')
flagdata(vis=LB_p5,mode='manual',spw='0', timerange='2015/08/30/11:48:50~2015/08/30/11:49:10',antenna='DV03')
flagdata(vis=LB_p5,mode='manual',spw='0', timerange='2015/08/30/11:49:20~2015/08/30/11:49:40',antenna='DV03')
flagdata(vis=LB_p5,mode='manual',spw='0', timerange='2015/08/30/11:49:50~2015/08/30/11:50:10',antenna='DV03')

flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/08/07/12:39:45~2017/08/07/12:39:55',antenna='DV12')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/08/07/12:41:50~2017/08/07/12:42:00',antenna='DV12,DV23,PM03')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/08/07/12:42:20~2017/08/07/12:42:30',antenna='DV12,PM03')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/08/07/12:43:10~2017/08/07/12:43:20',antenna='DA46,DV11,DV12')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/08/07/12:46:00~2017/08/07/12:46:30',antenna='DV07,DV12')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/08/07/12:46:40~2017/08/07/12:47:00',antenna='DV12')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/08/07/12:47:30~2017/08/07/12:47:50',antenna='DA53,DV07,DV12')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/08/07/12:48:00~2017/08/07/12:48:20',antenna='DV07,DV23')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/08/07/13:04:50~2017/08/07/13:05:10',antenna='DV12')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/08/07/13:05:50~2017/08/07/13:06:10',antenna='DA44')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/08/07/13:07:00~2017/08/07/13:07:20',antenna='DV12,PM04')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/08/07/13:08:30~2017/08/07/13:08:50',antenna='PM03')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/08/07/13:08:50~2017/08/07/13:09:10',antenna='PM03')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/08/07/13:09:50~2017/08/07/13:10:10',antenna='DV12,PM03')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/08/07/13:11:20~2017/08/07/13:11:40',antenna='DA46,DV12')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/08/07/12:46:00~2017/08/07/12:46:30',antenna='DV10,DA46')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/08/07/12:46:30~2017/08/07/12:47:00',antenna='DV10,DV07,DA53,DA46')
flagdata(vis=LB_p5,mode='manual',spw='8', timerange='2017/08/07/12:44:50~2017/08/07/12:45:00',antenna='DV10')

flagdata(vis=LB_p5,mode='manual',spw='16',timerange='2017/11/20/07:42:20~2017/11/20/07:43:00',antenna='DV15,DA60')
flagdata(vis=LB_p5,mode='manual',spw='16',timerange='2017/11/20/08:01:40~2017/11/20/08:01:50',antenna='DV22,DV15')
flagdata(vis=LB_p5,mode='manual',spw='16',timerange='2017/11/20/07:17:50~2017/11/20/07:18:10',antenna='DV17,DA64')
flagdata(vis=LB_p5,mode='manual',spw='16',timerange='2017/11/20/07:19:50~2017/11/20/07:20:10',antenna='DV05')
flagdata(vis=LB_p5,mode='manual',spw='16',timerange='2017/11/20/07:20:20~2017/11/20/07:20:30',antenna='DV05,DA65')
flagdata(vis=LB_p5,mode='manual',spw='16',timerange='2017/11/20/07:23:00~2017/11/20/07:23:10',antenna='DV17,DV22')
flagdata(vis=LB_p5,mode='manual',spw='16',timerange='2017/11/20/07:28:10~2017/11/20/07:28:20',antenna='DA64')
flagdata(vis=LB_p5,mode='manual',spw='16',timerange='2017/11/20/07:28:30~2017/11/20/07:28:40',antenna='DV17')
flagdata(vis=LB_p5,mode='manual',spw='16',timerange='2017/11/20/07:29:10~2017/11/20/07:29:20',antenna='DV17,DV23')
flagdata(vis=LB_p5,mode='manual',spw='16',timerange='2017/11/20/07:30:30~2017/11/20/07:30:50',antenna='DV05,DV23')
flagdata(vis=LB_p5,mode='manual',spw='16',timerange='2017/11/20/07:31:00~2017/11/20/07:31:10',antenna='DV16,DV05,PM04')
flagdata(vis=LB_p5,mode='manual',spw='16',timerange='2017/11/20/07:31:50~2017/11/20/07:32:10',antenna='DA65')
flagdata(vis=LB_p5,mode='manual',spw='16',timerange='2017/11/20/07:34:30~2017/11/20/07:34:50',antenna='DA64')
flagdata(vis=LB_p5,mode='manual',spw='16',timerange='2017/11/20/07:38:50~2017/11/20/07:39:00',antenna='DV23')
flagdata(vis=LB_p5,mode='manual',spw='16',timerange='2017/11/20/07:41:10~2017/11/20/07:41:30',antenna='DV17,DV16')
flagdata(vis=LB_p5,mode='manual',spw='16',timerange='2017/11/20/07:41:30~2017/11/20/07:41:50',antenna='DV15,DA60')
flagdata(vis=LB_p5,mode='manual',spw='16',timerange='2017/11/20/07:43:50~2017/11/20/07:44:10',antenna='DA65')
flagdata(vis=LB_p5,mode='manual',spw='16',timerange='2017/11/20/07:44:20~2017/11/20/07:44:40',antenna='DV22')
flagdata(vis=LB_p5,mode='manual',spw='16',timerange='2017/11/20/07:45:10~2017/11/20/07:45:30',antenna='DV17')
flagdata(vis=LB_p5,mode='manual',spw='16',timerange='2017/11/20/07:46:40~2017/11/20/07:46:50',antenna='DV17')
flagdata(vis=LB_p5,mode='manual',spw='16',timerange='2017/11/20/07:50:35~2017/11/20/07:50:45',antenna='DV22,DV17,DV16')
flagdata(vis=LB_p5,mode='manual',spw='16',timerange='2017/11/20/07:51:00~2017/11/20/07:51:20',antenna='DV16')
flagdata(vis=LB_p5,mode='manual',spw='16',timerange='2017/11/20/07:51:50~2017/11/20/07:52:10',antenna='DV17,DA64')
flagdata(vis=LB_p5,mode='manual',spw='16',timerange='2017/11/20/07:52:20~2017/11/20/07:52:40',antenna='DV17,DA60')
flagdata(vis=LB_p5,mode='manual',spw='16',timerange='2017/11/20/07:53:40~2017/11/20/07:54:00',antenna='DA60')
flagdata(vis=LB_p5,mode='manual',spw='16',timerange='2017/11/20/07:56:20~2017/11/20/07:56:40',antenna='DA65,DA64')
flagdata(vis=LB_p5,mode='manual',spw='16',timerange='2017/11/20/07:57:30~2017/11/20/07:57:50',antenna='DA49,DV13,DA42')
flagdata(vis=LB_p5,mode='manual',spw='16',timerange='2017/11/20/07:58:00~2017/11/20/07:58:20',antenna='DV13,DA42')
flagdata(vis=LB_p5,mode='manual',spw='16',timerange='2017/11/20/08:00:10~2017/11/20/08:00:20',antenna='DV16')
flagdata(vis=LB_p5,mode='manual',spw='16',timerange='2017/11/20/08:01:20~2017/11/20/08:01:40',antenna='PM03')
flagdata(vis=LB_p5,mode='manual',spw='16',timerange='2017/11/20/08:02:40~2017/11/20/08:02:50',antenna='PM03,DA42')
flagdata(vis=LB_p5,mode='manual',spw='16',timerange='2017/11/20/08:03:00~2017/11/20/08:03:10',antenna='DV22,PM03,DA42')
flagdata(vis=LB_p5,mode='manual',spw='16',timerange='2017/11/20/08:04:00~2017/11/20/08:04:10',antenna='DA60')
flagdata(vis=LB_p5,mode='manual',spw='16',timerange='2017/11/20/08:04:20~2017/11/20/08:04:30',antenna='DV05,DA60')

flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/11/23/04:49:40~2017/11/23/04:49:50',antenna='DV16,DA64')
flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/11/23/04:44:20~2017/11/23/04:44:30',antenna='DA64')
flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/11/23/04:49:10~2017/11/23/04:49:20',antenna='DV22')
flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/11/23/04:52:10~2017/11/23/04:52:20',antenna='DA46')
flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/11/23/04:52:40~2017/11/23/04:52:50',antenna='DA65,DA46')
flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/11/23/04:57:10~2017/11/23/04:57:30',antenna='DA60')
flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/11/23/04:57:40~2017/11/23/04:58:00',antenna='DV05,DA60')
flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/11/23/05:00:30~2017/11/23/05:00:40',antenna='DV22,DV05')
flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/11/23/05:01:50~2017/11/23/05:02:00',antenna='DV05')
flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/11/23/05:02:40~2017/11/23/05:03:00',antenna='DV16')
flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/11/23/05:03:00~2017/11/23/05:03:10',antenna='DV16')
flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/11/23/05:03:30~2017/11/23/05:03:40',antenna='DV22,DV15,DA60,PM04')
flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/11/23/05:06:00~2017/11/23/05:06:10',antenna='DV05')
flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/11/23/05:06:40~2017/11/23/05:07:00',antenna='DV22')
flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/11/23/05:07:10~2017/11/23/05:07:20',antenna='DA60')
flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/11/23/05:12:40~2017/11/23/05:13:00',antenna='DV05')
flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/11/23/05:13:40~2017/11/23/05:14:00',antenna='DV17')
flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/11/23/05:16:30~2017/11/23/05:16:40',antenna='DA60')
flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/11/23/05:18:00~2017/11/23/05:18:10',antenna='DV05')
flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/11/23/05:18:50~2017/11/23/05:19:10',antenna='DV22,DV17,DA64')
flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/11/23/05:19:20~2017/11/23/05:19:30',antenna='DA64')
flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/11/23/05:20:20~2017/11/23/05:20:30',antenna='DA65,DA60,PM03')
flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/11/23/05:20:40~2017/11/23/05:20:50',antenna='PM03')
flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/11/23/05:22:10~2017/11/23/05:22:20',antenna='DV05,DA65')
flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/11/23/05:27:20~2017/11/23/05:27:30',antenna='DV22,DV17')
flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/11/23/05:27:40~2017/11/23/05:27:50',antenna='DV22')
flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/11/23/05:31:00~2017/11/23/05:31:10',antenna='DV16,DA64')
flagdata(vis=LB_p5,mode='manual',spw='24',timerange='2017/11/23/05:31:30~2017/11/23/05:31:40',antenna='DA65')

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
    threshold = '0.0906mJy',
    **LB_tclean_wrapper_kwargs
)
estimate_SNR(LB_cont_p5+'.image',disk_mask=LB_mask,noise_mask=noise_annulus_LB)
generate_image_png(
    LB_cont_p5+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_LB,10*rms_LB],
    save_folder=LB_selfcal_folder
)
#CQ_Tau_SBLB_contp5.image
#Beam 0.078 arcsec x 0.058 arcsec (7.28 deg)
#Flux inside disk mask: 161.08 mJy
#Peak intensity of source: 2.69 mJy/beam
#rms: 1.50e-02 mJy/beam
#Peak SNR: 179.52

#Sixth round of phase-only self-cal
LB_p6 = LB_p5.replace('p5','p6')
os.system('rm -rf '+LB_p6)
gaincal(
    vis=LB_cont_p5+'.ms',caltable=LB_p6,
    gaintype='T',combine='spw',calmode='p',solint='18s',
    spw=LB_contspws,refant=LB_refant,
    minsnr=2.,minblperant=4
)
#1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/12:41:49.4
#1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/12:42:25.7
#2 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/12:43:46.3
#1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/12:48:11.6
#1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/12:49:14.0
#1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/12:52:17.1
#1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/12:54:21.8
#1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/12:58:23.4
#1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/12:59:07.9
#1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/12:59:25.8
#1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/13:06:27.7
#1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/13:08:35.8
#1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/13:08:53.8
#1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/13:09:12.0
#1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/13:11:52.7
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:17:29.2
#7 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:17:47.3
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:18:32.2
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:18:49.7
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:19:07.9
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:19:52.4
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:20:10.3
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:20:28.4
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:21:13.0
#7 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:21:30.8
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:21:49.0
#7 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:22:33.6
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:22:51.4
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:23:09.5
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:23:53.8
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:24:11.9
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:24:30.0
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:25:14.4
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:25:32.4
#7 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:26:03.7
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:28:10.4
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:28:28.6
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:29:13.0
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:29:31.0
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:29:49.1
#7 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:30:33.6
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:30:51.5
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:31:09.6
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:31:54.1
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:32:12.0
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:32:30.2
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:33:14.7
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:33:32.6
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:33:50.7
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:34:35.2
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:34:53.1
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:35:11.3
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:35:55.6
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:36:13.7
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:36:44.8
#7 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:38:51.6
#3 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:39:09.8
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:39:54.2
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:40:12.2
#7 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:40:30.3
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:41:14.7
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:41:32.7
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:41:50.6
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:42:35.3
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:42:53.3
#2 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:43:11.4
#7 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:43:55.9
#3 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:44:13.8
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:44:32.0
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:45:16.5
#7 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:45:34.4
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:45:52.5
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:46:37.0
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:46:54.9
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:47:26.0
#7 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:49:31.9
#7 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:49:50.1
#3 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:50:34.9
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:50:52.5
#3 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:51:10.6
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:51:55.2
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:52:13.0
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:52:31.2
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:53:15.7
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:53:33.6
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:53:51.7
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:54:36.1
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:54:54.1
#7 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:55:12.2
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:55:56.7
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:56:14.6
#3 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:56:32.8
#7 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:57:17.0
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:57:35.2
#7 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:58:06.6
#7 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:00:12.5
#3 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:00:30.7
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:01:15.0
#7 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:01:33.1
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:01:51.2
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:02:35.6
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:02:53.6
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:03:11.8
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:03:56.1
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:04:14.2
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:04:32.3
#3 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:05:16.6
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:05:34.7
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:05:52.8
#7 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:06:37.1
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:43:59.4
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:44:17.5
#7 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:45:02.8
#2 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:45:21.0
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:45:39.1
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:46:24.4
#2 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:46:42.6
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:47:00.6
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:47:46.0
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:48:04.2
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:48:22.3
#3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:49:07.6
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:49:25.8
#3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:49:43.9
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:50:29.2
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:50:47.4
#7 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:51:05.5
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:51:50.8
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:52:09.0
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:52:41.6
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:54:50.8
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:55:09.0
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:55:54.2
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:56:12.3
#8 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:56:30.5
#7 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:57:15.7
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:57:33.8
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:57:51.7
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:58:37.2
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:58:55.3
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:59:13.5
#7 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:59:58.8
#7 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:00:16.8
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:00:35.0
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:01:20.2
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:01:38.4
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:01:56.5
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:02:41.7
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:02:59.9
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:03:32.3
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:05:41.2
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:05:59.1
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:06:44.4
#7 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:07:02.6
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:07:20.7
#3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:08:05.9
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:08:24.0
#7 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:08:42.1
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:09:27.3
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:09:45.4
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:10:03.6
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:10:48.7
#7 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:11:06.8
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:11:25.0
#3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:12:10.1
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:12:28.2
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:12:46.4
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:13:31.5
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:13:49.6
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:14:22.0
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:16:29.9
#2 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:16:48.0
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:17:33.1
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:17:51.2
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:18:09.4
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:18:54.4
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:19:12.5
#3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:19:29.9
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:20:15.7
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:20:33.8
#3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:20:52.0
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:21:37.3
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:21:55.2
#2 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:22:13.3
#8 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:22:58.3
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:23:16.5
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:23:34.6
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:24:19.6
#7 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:24:37.8
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:25:09.5
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:27:16.6
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:27:34.8
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:28:19.7
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:28:37.8
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:28:56.0
#7 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:29:41.0
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:29:59.0
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:30:17.2
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:31:02.1
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:31:20.3
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:31:38.4
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:32:23.4
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:32:41.5
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:32:59.6
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:33:44.5

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
flagdata(vis=LB_p6,mode='manual',spw='0',antenna='DV03,DV07')

flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/08/07/12:39:40~2017/08/07/12:39:50',antenna='DV12')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/08/07/12:41:40~2017/08/07/12:41:50',antenna='DV23,DV12,PM03')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/08/07/12:42:00~2017/08/07/12:42:10',antenna='PM03,DV23,DV12')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/08/07/12:42:20~2017/08/07/12:42:30',antenna='PM03,DV12')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/08/07/12:43:00~2017/08/07/12:43:20',antenna='DV12,DV11,DA46')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/08/07/12:43:20~2017/08/07/12:43:40',antenna='DA46,DV11')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/08/07/12:46:10~2017/08/07/12:46:20',antenna='DV12,DV10,DV07,DA46')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/08/07/12:46:30~2017/08/07/12:46:40',antenna='DV10,DV07,DA46')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/08/07/12:46:50~2017/08/07/12:47:00',antenna='DV12,DV10,DV07,DA46')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/08/07/12:47:30~2017/08/07/12:47:40',antenna='DV12,DV07,DA53')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/08/07/12:48:10~2017/08/07/12:48:20',antenna='DV23,DV07')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/08/07/13:05:00~2017/08/07/13:05:10',antenna='DV12')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/08/07/13:05:50~2017/08/07/13:06:00',antenna='DA44')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/08/07/13:07:10~2017/08/07/13:07:20',antenna='PM04,DV12')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/08/07/13:08:30~2017/08/07/13:08:40',antenna='PM03')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/08/07/13:08:50~2017/08/07/13:09:00',antenna='PM03')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/08/07/13:09:50~2017/08/07/13:10:00',antenna='PM03,DV12,DA53')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/08/07/13:10:10~2017/08/07/13:10:20',antenna='PM03,DV12')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/08/07/13:11:10~2017/08/07/13:11:20',antenna='DV12,DA46')
flagdata(vis=LB_p6,mode='manual',spw='8', timerange='2017/08/07/13:11:30~2017/08/07/13:11:40',antenna='PM03,DV12,DA46')

flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/11/20/07:18:30~2017/11/20/07:18:35',antenna='DV17')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/11/20/07:19:50~2017/11/20/07:19:55',antenna='DV05')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/11/20/07:20:05~2017/11/20/07:20:15',antenna='DV22')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/11/20/07:21:30~2017/11/20/07:21:40',antenna='DV05')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/11/20/07:22:45~2017/11/20/07:22:55',antenna='DA65')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/11/20/07:23:05~2017/11/20/07:23:15',antenna='DV17,DA64')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/11/20/07:23:50~2017/11/20/07:23:55',antenna='DV17')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/11/20/07:28:10~2017/11/20/07:28:15',antenna='DA65,DA64')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/11/20/07:29:10~2017/11/20/07:29:20',antenna='DV23')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/11/20/07:29:30~2017/11/20/07:29:40',antenna='DV23')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/11/20/07:29:45~2017/11/20/07:29:55',antenna='DA65')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/11/20/07:31:05~2017/11/20/07:31:15',antenna='DV17')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/11/20/07:31:50~2017/11/20/07:32:00',antenna='DV05')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/11/20/07:32:10~2017/11/20/07:32:20',antenna='DA65')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/11/20/07:38:50~2017/11/20/07:38:55',antenna='DV23')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/11/20/07:39:05~2017/11/20/07:39:15',antenna='DV22,DV17')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/11/20/07:42:30~2017/11/20/07:42:40',antenna='DV15,DA60')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/11/20/07:42:50~2017/11/20/07:43:00',antenna='DV15')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/11/20/07:43:50~2017/11/20/07:44:00',antenna='DA65')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/11/20/07:50:30~2017/11/20/07:50:40',antenna='DV22,DV16')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/11/20/07:50:50~2017/11/20/07:51:00',antenna='DV16')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/11/20/07:53:30~2017/11/20/07:53:40',antenna='DA60')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/11/20/07:53:50~2017/11/20/07:54:00',antenna='DA60')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/11/20/07:55:55~2017/11/20/07:56:00',antenna='DV16')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/11/20/07:56:30~2017/11/20/07:56:40',antenna='DA64,DA65')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/11/20/07:58:05~2017/11/20/07:58:10',antenna='DA64,DV23,DV13,DA42')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/11/20/08:00:30~2017/11/20/08:00:40',antenna='DV16')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/11/20/08:01:10~2017/11/20/08:01:20',antenna='DA60,PM03')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/11/20/08:01:30~2017/11/20/08:01:40',antenna='DA65,DA49,PM03')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/11/20/08:01:50~2017/11/20/08:02:00',antenna='DV22,DV15,PM03')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/11/20/08:02:30~2017/11/20/08:02:40',antenna='PM03,DA42')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/11/20/08:02:50~2017/11/20/08:03:00',antenna='DV05,PM03,DA42')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/11/20/08:03:10~2017/11/20/08:03:20',antenna='DV22,PM03,DA42')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/11/20/08:04:10~2017/11/20/08:04:20',antenna='DA64,DA60')
flagdata(vis=LB_p6,mode='manual',spw='16',timerange='2017/11/20/08:04:30~2017/11/20/08:04:40',antenna='DV05,DA60')

flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/11/23/04:44:10~2017/11/23/04:44:20',antenna='DA65')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/11/23/04:47:00~2017/11/23/04:47:10',antenna='DA64')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/11/23/04:48:20~2017/11/23/04:48:30',antenna='DV22')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/11/23/04:48:15~2017/11/23/04:48:25',antenna='DV22')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/11/23/04:49:05~2017/11/23/04:49:15',antenna='DV22')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/11/23/04:50:40~2017/11/23/04:50:50',antenna='DV17')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/11/23/04:52:00~2017/11/23/04:52:10',antenna='DV16,DA46')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/11/23/04:52:40~2017/11/23/04:52:50',antenna='DA64,DA46')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/11/23/04:57:30~2017/11/23/04:57:40',antenna='DA65')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/11/23/04:57:50~2017/11/23/04:58:00',antenna='DV17,DV05')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/11/23/05:00:30~2017/11/23/05:00:40',antenna='DV17,DV22,DA60')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/11/23/05:01:30~2017/11/23/05:01:40',antenna='DA64')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/11/23/05:01:50~2017/11/23/05:02:00',antenna='DV16,DV05')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/11/23/05:03:30~2017/11/23/05:03:40',antenna='DV22,DV17,DV15,DA60,PM04')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/11/23/05:09:20~2017/11/23/05:09:30',antenna='DV22,DV05,DA64')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/11/23/05:12:05~2017/11/23/05:12:15',antenna='DV05')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/11/23/05:12:40~2017/11/23/05:12:50',antenna='DV17,DV05,DA60')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/11/23/05:13:45~2017/11/23/05:13:55',antenna='DV16,DV17')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/11/23/05:17:45~2017/11/23/05:17:55',antenna='DV16,DA65')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/11/23/05:19:25~2017/11/23/05:19:35',antenna='DV22,DV17,DA64')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/11/23/05:20:10~2017/11/23/05:20:20',antenna='PM03')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/11/23/05:20:30~2017/11/23/05:20:40',antenna='PM03')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/11/23/05:20:50~2017/11/23/05:21:00',antenna='DV22,DA65,DA60,PM03')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/11/23/05:22:10~2017/11/23/05:22:20',antenna='DA64')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/11/23/05:27:10~2017/11/23/05:27:20',antenna='DV17,DA65')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/11/23/05:27:30~2017/11/23/05:27:40',antenna='DV22')
flagdata(vis=LB_p6,mode='manual',spw='24',timerange='2017/11/23/05:33:40~2017/11/23/05:33:50',antenna='DV16')

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
    threshold = '0.0900mJy',
    **LB_tclean_wrapper_kwargs
)
estimate_SNR(LB_cont_p6+'.image',disk_mask=LB_mask,noise_mask=noise_annulus_LB)
generate_image_png(
    LB_cont_p6+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_LB,10*rms_LB],
    save_folder=LB_selfcal_folder
)
#CQ_Tau_SBLB_contp6.image
#Beam 0.078 arcsec x 0.058 arcsec (7.28 deg)
#Flux inside disk mask: 161.20 mJy
#Peak intensity of source: 2.71 mJy/beam
#rms: 1.49e-02 mJy/beam
#Peak SNR: 181.98

#For each step of the self-cal, check how it improved things
non_self_caled_LB_vis = LB_cont_p0
self_caled_LB_visibilities = {
    'p1':LB_cont_p1,
    'p2':LB_cont_p2,
    'p3':LB_cont_p3,
    'p4':LB_cont_p4,
    'p5':LB_cont_p5,
    'p6':LB_cont_p6,
}

#for vis in self_caled_LB_visibilities.values(): 
#    listobs(vis=vis+'.ms',listfile=vis+'.ms.listobs.txt',overwrite=True)
#
#LB_EBs = ('EB0','EB1','EB2','EB3')
#LB_EB_spws = ('0,1,2,3,4,5,6,7', '8,9,10,11,12,13,14,15', '16,17,18,19,20,21,22,23', '24,25,26,27,28,29,30,31') #fill out by referring to listobs output
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
SBLB_flux_ref_EB = 3 #this is LB_EB1

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

#ratio          = [1.36589,1.36531,1.37131,1.37702,1.37817,1.37872] #CQ_Tau_SB_contp0...p5_EB0.vis.npz vs CQ_Tau_LB_EB1_initcont_shift.vis.npz
#scaling_factor = [1.169  ,1.168  ,1.171  ,1.173  ,1.174  ,1.174  ]
#ratio          = [0.91086,0.91196,0.94078,0.97047,0.99275,1.02373] #CQ_Tau_SB_contp0...p5_EB1.vis.npz vs CQ_Tau_LB_EB1_initcont_shift.vis.npz
#scaling_factor = [0.954  ,0.955  ,0.970  ,0.985  ,0.996  ,1.012  ]

#ratio          = [1.37869,1.37799,1.37359,1.34796,1.33492,1.32590,1.32198] #CQ_Tau_SB_contp0...p6_EB0.vis.npz vs CQ_Tau_LB_contp0...p6_EB1.vis.npz
#scaling_factor = [1.174  ,1.174  ,1.172  ,1.161  ,1.155  ,1.151  ,1.150  ]
#ratio          = [1.02371,1.02283,1.01886,0.99424,0.98171,0.97528,0.97514] #CQ_Tau_SB_contp0...p6_EB1.vis.npz vs CQ_Tau_LB_contp0...p6_EB1.vis.npz
#scaling_factor = [1.012  ,1.011  ,1.009  ,0.997  ,0.991  ,0.988  ,0.987  ]
#ratio          = [0.95141,0.95212,0.96038,0.97280,0.98479,0.99312,0.99917] #CQ_Tau_LB_contp0...p6_EB0.vis.npz vs CQ_Tau_LB_contp0...p6_EB1.vis.npz
#scaling_factor = [0.975  ,0.976  ,0.980  ,0.986  ,0.992  ,0.997  ,1.000  ]

#END of COMBINED SB+LB phase-only self-cal iteration 1

#Rescaling flux of the non-self-caled SB and LB EBs
for baseline_key,n_EB in number_of_EBs.items():
    for i in range(n_EB):
        vis = f'{prefix}_{baseline_key}_EB{i}_initcont_shift.ms'
        listobs(vis=vis,listfile=f'{vis}.listobs.txt',overwrite=True)

# Check that you're using the right EB (i.e., no strange stuff happened on the original shifted EBs, and, in terms of flux scales, you do recover the original ones)
# output = f'flux_comparison_SB_EB0_to_LB_EB1.png'
# plot_label = os.path.join(LB_selfcal_folder,output)
# estimate_flux_scale(
#     reference  = prefix+'_LB_EB1_initcont_shift.vis.npz',
#     comparison = prefix+'_SB_EB0_initcont_shift.vis.npz',
#     incl=incl,PA=PA,plot_label=plot_label
# )
#The ratio of the fluxes of CQ_Tau_SB_EB0_initcont_shift.vis.npz to
#CQ_Tau_LB_EB1_initcont_shift.vis.npz is 1.36589
#The scaling factor for gencal is 1.169 for your comparison measurement
#The error on the weighted mean ratio is 1.066e-03, although it's likely that
#the weights in the measurement sets are too off by some constant factor

os.system('rm -rf '+prefix+'_SB_EB0_initcont_shift_rescaled.ms')
rescale_flux(vis=prefix+'_SB_EB0_initcont_shift.ms', gencalparameter=[1.150]) #gencal parameter from SBLB_contp5
#Splitting out rescaled values into new MS: CQ_Tau_SB_EB0_initcont_shift_rescaled.ms
listobs(vis=prefix+'_SB_EB0_initcont_shift_rescaled.ms',listfile=prefix+'_SB_EB0_initcont_shift_rescaled.ms.listobs.txt',overwrite=True)

export_MS(prefix+'_SB_EB0_initcont_shift_rescaled.ms') 
#Measurement set exported to CQ_Tau_SB_EB0_initcont_shift_rescaled.vis.npz

# Check if there is still an offset, indeed the gencalparameters were differente! Selfcal will correct for it!
# output = f'flux_comparison_SB_EB0_rescaled_to_LB_EB1.png'
# plot_label = os.path.join(LB_selfcal_folder,output)
# estimate_flux_scale(
#     reference  = prefix+'_LB_EB1_initcont_shift.vis.npz',
#     comparison = prefix+'_SB_EB0_initcont_shift_rescaled.vis.npz',
#     incl=incl,PA=PA,plot_label=plot_label
# )
#The ratio of the fluxes of CQ_Tau_SB_EB0_initcont_shift_rescaled.vis.npz to
#CQ_Tau_LB_EB1_initcont_shift.vis.npz is 1.03281
#The scaling factor for gencal is 1.016 for your comparison measurement
#The error on the weighted mean ratio is 8.064e-04, although it's likely that
#the weights in the measurement sets are too off by some constant factor

#Begin of SB self-cal - iteration 2
#For phase self-cal, clean down to ~6sigma
SB_selfcal_iteration2_folder = get_figures_folderpath('7.1_selfcal_SB_iteration2_figures')
make_figures_folder(SB_selfcal_iteration2_folder)

SB_iteration2_cont_p0 = prefix+'_SB_iteration2_contp0'
os.system('rm -rf %s.ms*' %SB_iteration2_cont_p0)

concat(
    vis=[f'{prefix}_SB_EB0_initcont_shift_rescaled.ms',f'{prefix}_SB_EB1_initcont_shift.ms'],
    concatvis=SB_iteration2_cont_p0+'.ms',dirtol='0.1arcsec',copypointing=False
)
listobs(vis=SB_iteration2_cont_p0+'.ms',listfile=SB_iteration2_cont_p0+'.ms.listobs.txt',overwrite=True)
#2024-12-16 21:19:22     WARN    concat::::casa  The setup of the input MSs is not fully consistent. The concatenation may fail
#   and/or the affected columns may contain partially only default data.
#   {'CQ_Tau_SB_EB1_initcont_shift.ms': {'Main': {'present_a': True, 'present_b': True, 'missingcol_a': ['MODEL_DATA'], 'missingcol_b': []}}}
#2024-12-16 21:19:55     WARN    MSConcat::copySysCal    /data/beegfs/astro-storage/groups/benisty/frzagaria/SO_detections/CQTau/selfcal_products/CQ_Tau_SB_iteration2_contp0.ms does not have a valid syscal table,
#   the MS to be appended, however, has one. Result won't have one.
#2024-12-16 21:19:55     WARN    MSConcat::concatenate (file /source/casa6/casatools/casacore/ms/MSOper/MSConcat.cc, line 1000)        Could not merge SysCal subtables

#Define new SB mask using the same centre as before (checked and agrees with the listobs one)
mask_pa        = PA
mask_semimajor = 0.85
mask_semiminor = 0.72
mask_ra  = '05h35m58.472722s'
mask_dec = '24d44m53.613450s'

SB_mask= f"ellipse[[{mask_ra},{mask_dec}], [{mask_semimajor:.3f}arcsec, {mask_semiminor:.3f}arcsec], {mask_pa:.1f}deg]"

noise_annulus_SB = f"annulus[[{mask_ra}, {mask_dec}],['4.arcsec', '6.arcsec']]"

SB_tclean_wrapper_kwargs = {
    'deconvolver':'multiscale','scales':[0,2,4,8,12],
    'smallscalebias':0.6,'gain':0.3,'cycleniter':300,
    'niter':1000000,'mask':SB_mask,
    'cellsize':'0.013arcsec','imsize':3600,
    'parallel':use_parallel,'savemodel':'modelcolumn',
    'robust':0.5,'interactive':False,
    'gridder':'standard',
}

tclean_wrapper(
    vis       = SB_iteration2_cont_p0+'.ms',
    imagename = SB_iteration2_cont_p0,
    threshold = '0.2754mJy', 
    **SB_tclean_wrapper_kwargs
)
estimate_SNR(SB_iteration2_cont_p0+'.image',disk_mask=SB_mask,noise_mask=noise_annulus_SB)
rms_iteration2_SB = imstat(imagename=SB_iteration2_cont_p0+'.image',region=noise_annulus_SB)['rms'][0]
generate_image_png(
    SB_iteration2_cont_p0+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_iteration2_SB,10*rms_iteration2_SB],
    save_folder=SB_selfcal_iteration2_folder
)
#CQ_Tau_SB_contp0.image
#Beam 0.138 arcsec x 0.108 arcsec (-2.02 deg)
#Flux inside disk mask: 169.43 mJy
#Peak intensity of source: 7.17 mJy/beam
#rms: 6.36e-02 mJy/beam
#Peak SNR: 112.77

#CQ_Tau_SB_iteration2_contp0.image
#Beam 0.145 arcsec x 0.113 arcsec (-2.87 deg)
#Flux inside disk mask: 142.22 mJy
#Peak intensity of source: 6.79 mJy/beam
#rms: 4.57e-02 mJy/beam
#Peak SNR: 148.65

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

#Flag problematic antennas
flagdata(vis=SB_iteration2_p1,mode='manual',spw='8',antenna='DV12')

#Inspect gain tables to check if flagging worked
plotms(
    SB_iteration2_p1,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(SB_selfcal_iteration2_folder,f'{prefix}_SB_iteration2_p1_gain_phase_vs_time_flagged.png'),
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
    threshold = '0.2730mJy',
    **SB_tclean_wrapper_kwargs
)
estimate_SNR(SB_iteration2_cont_p1+'.image',disk_mask=SB_mask,noise_mask=noise_annulus_SB)
generate_image_png(
    SB_iteration2_cont_p1+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_iteration2_SB,10*rms_iteration2_SB],
    save_folder=SB_selfcal_iteration2_folder
)
#CQ_Tau_SB_contp1.image
#Beam 0.138 arcsec x 0.108 arcsec (-2.02 deg)
#Flux inside disk mask: 169.56 mJy
#Peak intensity of source: 7.19 mJy/beam
#rms: 6.35e-02 mJy/beam
#Peak SNR: 113.32

#CQ_Tau_SB_iteration2_contp1.image
#Beam 0.145 arcsec x 0.113 arcsec (-2.87 deg)
#Flux inside disk mask: 142.29 mJy
#Peak intensity of source: 6.79 mJy/beam
#rms: 4.55e-02 mJy/beam
#Peak SNR: 149.16

#Second round of phase-only self-cal
SB_iteration2_p2 = SB_iteration2_p1.replace('p1','p2')
os.system('rm -rf '+SB_iteration2_p2)
gaincal(
    vis=SB_iteration2_cont_p1+'.ms',caltable=SB_iteration2_p2,
    gaintype='T',combine='scan,spw',calmode='p',solint='360s',
    spw=SB_contspws,refant=SB_refant,
    minsnr=3.,minblperant=4
)

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
flagdata(vis=SB_iteration2_p2,mode='manual',spw='0',timerange='2015/08/30/11:41:00~2015/08/30/11:43:00',antenna='DV03,DV07')
flagdata(vis=SB_iteration2_p2,mode='manual',spw='0',timerange='2015/08/30/11:54:00~2015/08/30/11:56:00',antenna='DV07')
flagdata(vis=SB_iteration2_p2,mode='manual',spw='8',timerange='2017/08/07/12:39:00~2017/08/07/12:41:00',antenna='DA53')
flagdata(vis=SB_iteration2_p2,mode='manual',spw='8',timerange='2017/08/07/12:42:00~2017/08/07/12:44:00',antenna='DV12')
flagdata(vis=SB_iteration2_p2,mode='manual',spw='8',timerange='2017/08/07/12:46:00~2017/08/07/12:48:00',antenna='DA46,DV07,DV10,DV12')

#Inspect gain tables to check if flagging worked
plotms(
    SB_iteration2_p2,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(SB_selfcal_iteration2_folder,f'{prefix}_SB_iteration2_p2_gain_phase_vs_time_flagged.png'),
)

#Apply the calibration gains and split-off the corrected ms
applycal(
    vis=SB_iteration2_cont_p1+'.ms',spw=SB_contspws,spwmap=SB_spw_mapping,
    gaintable=[SB_iteration2_p2],interp='linearPD',calwt=True,applymode='calonly'
)
SB_iteration2_cont_p2 = SB_iteration2_cont_p1.replace('p1','p2')
os.system('rm -rf %s.ms*' %SB_iteration2_cont_p2)
split(vis=SB_iteration2_cont_p1+'.ms',outputvis=SB_iteration2_cont_p2+'.ms',datacolumn='corrected')

#CLEAN, again with a ~6sigma threshold
tclean_wrapper(
    vis       = SB_iteration2_cont_p2+'.ms',
    imagename = SB_iteration2_cont_p2,
    threshold = '0.2490mJy',
    **SB_tclean_wrapper_kwargs
)
estimate_SNR(SB_iteration2_cont_p2+'.image',disk_mask=SB_mask,noise_mask=noise_annulus_SB)
generate_image_png(
    SB_iteration2_cont_p2+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_iteration2_SB,10*rms_iteration2_SB],
    save_folder=SB_selfcal_iteration2_folder
)
#CQ_Tau_SB_contp2.image
#Beam 0.138 arcsec x 0.108 arcsec (-2.02 deg)
#Flux inside disk mask: 171.03 mJy
#Peak intensity of source: 7.39 mJy/beam
#rms: 5.92e-02 mJy/beam
#Peak SNR: 124.91

#CQ_Tau_SB_iteration2_contp2.image
#Beam 0.145 arcsec x 0.113 arcsec (-2.87 deg)
#Flux inside disk mask: 143.33 mJy
#Peak intensity of source: 6.94 mJy/beam
#rms: 4.12e-02 mJy/beam
#Peak SNR: 168.27

#Third round of phase-only self-cal
SB_iteration2_p3 = SB_iteration2_p2.replace('p2','p3')
os.system('rm -rf '+SB_iteration2_p3)
gaincal(
    vis=SB_iteration2_cont_p2+'.ms',caltable=SB_iteration2_p3,
    gaintype='T',combine='scan,spw',calmode='p',solint='120s',
    spw=SB_contspws,refant=SB_refant,
    minsnr=3.,minblperant=4
)
#11 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:39:52.5
# 1 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:41:18.3
# 1 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:43:56.9
# 2 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:48:13.6
# 2 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:58:25.4
# 1 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:59:26.1
# 2 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/13:03:48.8
# 4 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/13:06:29.7

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
flagdata(vis=SB_iteration2_p3,mode='manual',spw='0',antenna='DV03,DV07')

flagdata(vis=SB_iteration2_p3,mode='manual',spw='8',timerange='2017/08/07/12:39:20~2017/08/07/12:39:40',antenna='DA53,DV12')
flagdata(vis=SB_iteration2_p3,mode='manual',spw='8',timerange='2017/08/07/12:42:20~2017/08/07/12:42:40',antenna='DV12')
flagdata(vis=SB_iteration2_p3,mode='manual',spw='8',timerange='2017/08/07/12:43:50~2017/08/07/12:44:10',antenna='DV11,DV12')
flagdata(vis=SB_iteration2_p3,mode='manual',spw='8',timerange='2017/08/07/12:46:50~2017/08/07/12:47:10',antenna='DA53,DV07,DV10,DV12')
flagdata(vis=SB_iteration2_p3,mode='manual',spw='8',timerange='2017/08/07/13:08:00~2017/08/07/13:08:30',antenna='PM03')
flagdata(vis=SB_iteration2_p3,mode='manual',spw='8',timerange='2017/08/07/13:09:00~2017/08/07/13:09:30',antenna='DA42')
flagdata(vis=SB_iteration2_p3,mode='manual',spw='8',timerange='2017/08/07/13:10:30~2017/08/07/13:11:00',antenna='DA46')

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
    threshold = '0.2376mJy',
    **SB_tclean_wrapper_kwargs
)
estimate_SNR(SB_iteration2_cont_p3+'.image',disk_mask=SB_mask,noise_mask=noise_annulus_SB)
generate_image_png(
    SB_iteration2_cont_p3+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_iteration2_SB,10*rms_iteration2_SB],
    save_folder=SB_selfcal_iteration2_folder
)
#CQ_Tau_SB_contp3.image
#Beam 0.138 arcsec x 0.108 arcsec (-2.02 deg)
#Flux inside disk mask: 171.63 mJy
#Peak intensity of source: 7.65 mJy/beam
#rms: 5.72e-02 mJy/beam
#Peak SNR: 133.70

#CQ_Tau_SB_iteration2_contp3.image
#Beam 0.145 arcsec x 0.113 arcsec (-2.87 deg)
#Flux inside disk mask: 144.09 mJy
#Peak intensity of source: 7.21 mJy/beam
#rms: 3.95e-02 mJy/beam
#Peak SNR: 182.82

#Fourth round of phase-only self-cal
SB_iteration2_p4 = SB_iteration2_p3.replace('p3','p4')
os.system('rm -rf '+SB_iteration2_p4)
gaincal(
    vis=SB_iteration2_cont_p3+'.ms',caltable=SB_iteration2_p4,
    gaintype='T',combine='spw',calmode='p',solint='60s',
    spw=SB_contspws,refant=SB_refant,
    minsnr=3.,minblperant=4
)
#2 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:42:07.5
#1 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/13:07:12.1
#1 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/13:11:34.6

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
flagdata(vis=SB_iteration2_p4,mode='manual',spw='0',antenna='DV03,DV07')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='8',antenna='DV12')

flagdata(vis=SB_iteration2_p4,mode='manual',spw='8',timerange='2017/08/07/12:39:30~2017/08/07/12:39:40',antenna='DA53')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='8',timerange='2017/08/07/12:43:20~2017/08/07/12:43:40',antenna='DV11,DA46')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='8',timerange='2017/08/07/12:46:20~2017/08/07/12:46:40',antenna='DV10,DV07,DA46')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='8',timerange='2017/08/07/12:47:50~2017/08/07/12:48:00',antenna='DA53')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='8',timerange='2017/08/07/12:51:10~2017/08/07/12:51:20',antenna='DV11')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='8',timerange='2017/08/07/12:58:00~2017/08/07/12:58:20',antenna='PM04')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='8',timerange='2017/08/07/13:02:00~2017/08/07/13:02:20',antenna='DV23')
flagdata(vis=SB_iteration2_p4,mode='manual',spw='8',timerange='2017/08/07/13:08:50~2017/08/07/13:09:00',antenna='PM03')

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
    threshold = '0.2274mJy',
    **SB_tclean_wrapper_kwargs
)
estimate_SNR(SB_iteration2_cont_p4+'.image',disk_mask=SB_mask,noise_mask=noise_annulus_SB)
generate_image_png(
    SB_iteration2_cont_p4+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_iteration2_SB,10*rms_iteration2_SB],
    save_folder=SB_selfcal_iteration2_folder
)
#CQ_Tau_SB_contp4.image
#Beam 0.138 arcsec x 0.108 arcsec (-2.02 deg)
#Flux inside disk mask: 172.23 mJy
#Peak intensity of source: 7.82 mJy/beam
#rms: 5.49e-02 mJy/beam
#Peak SNR: 142.64

#CQ_Tau_SB_iteration2_contp4.image
#Beam 0.145 arcsec x 0.113 arcsec (-2.87 deg)
#Flux inside disk mask: 144.44 mJy
#Peak intensity of source: 7.39 mJy/beam
#rms: 3.77e-02 mJy/beam
#Peak SNR: 195.85

#Fifth round of phase-only selfcal
SB_iteration2_p5 = SB_iteration2_p4.replace('p4','p5')
os.system('rm -rf '+SB_iteration2_p5)
gaincal(
    vis=SB_iteration2_cont_p4+'.ms',caltable=SB_iteration2_p5,
    gaintype='T',combine='spw',calmode='p',solint='20s',
    spw=SB_contspws,refant=SB_refant,
    minsnr=3.,minblperant=4
)

#1 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:40:29.8
#2 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:41:06.4
#1 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:41:50.4
#2 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:42:27.7
#2 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:43:48.3
#1 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:44:49.8
#1 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:46:15.8
#2 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:46:36.0
#2 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:46:53.1
#1 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:47:56.5
#3 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:48:13.6
#2 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:49:17.0
#5 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:49:34.2
#1 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:51:17.6
#2 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:51:34.8
#1 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:52:18.1
#1 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:52:36.1
#2 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:54:22.8
#1 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:55:07.4
#3 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:55:44.3
#1 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:56:28.1
#2 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:57:04.9
#1 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:58:25.4
#2 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:59:08.9
#1 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/12:59:28.8
#2 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/13:01:51.1
#2 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/13:02:11.2
#1 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/13:02:28.4
#2 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/13:03:48.8
#2 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/13:04:32.1
#2 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/13:05:09.3
#1 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/13:06:12.6
#5 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/13:06:29.7
#1 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/13:07:12.1
#1 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/13:08:36.8
#3 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/13:09:14.0
#2 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/13:10:17.2
#1 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/13:10:34.3
#1 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/13:11:37.5
#2 of 40 solutions flagged due to SNR < 3 in spw=8 at 2017/08/07/13:13:34.7

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
flagdata(vis=SB_iteration2_p5,mode='manual',spw='0',antenna='DV03,DV07')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8',antenna='DV12')

flagdata(vis=SB_iteration2_p5,mode='manual',spw='8',timerange='2017/08/07/12:39:20~2017/08/07/12:39:30',antenna='DA53')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8',timerange='2017/08/07/12:39:40~2017/08/07/12:39:50',antenna='DA53')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8',timerange='2017/08/07/12:40:25~2017/08/07/12:40:35',antenna='DA53,DA46')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8',timerange='2017/08/07/12:41:45~2017/08/07/12:41:55',antenna='PM03,DV23')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8',timerange='2017/08/07/12:42:05~2017/08/07/12:42:15',antenna='PM03,DV23')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8',timerange='2017/08/07/12:43:05~2017/08/07/12:43:15',antenna='DV11,DA46')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8',timerange='2017/08/07/12:43:25~2017/08/07/12:43:35',antenna='DV11,DA46')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8',timerange='2017/08/07/12:44:45~2017/08/07/12:44:55',antenna='DV10,DA46')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8',timerange='2017/08/07/12:46:10~2017/08/07/12:46:20',antenna='DV16,DV10,DV07,DA50,DA46')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8',timerange='2017/08/07/12:46:30~2017/08/07/12:46:40',antenna='DV10,DV07,DA46')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8',timerange='2017/08/07/12:46:50~2017/08/07/12:47:00',antenna='DV10,DV07,DA46')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8',timerange='2017/08/07/12:47:30~2017/08/07/12:47:40',antenna='DV07,DA53')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8',timerange='2017/08/07/12:47:50~2017/08/07/12:48:00',antenna='DA53')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8',timerange='2017/08/07/12:48:10~2017/08/07/12:48:20',antenna='DV07')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8',timerange='2017/08/07/12:56:25~2017/08/07/12:56:35',antenna='PM04,DA53')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8',timerange='2017/08/07/12:58:20~2017/08/07/12:58:30',antenna='PM03')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8',timerange='2017/08/07/13:05:50~2017/08/07/13:05:55',antenna='DA44')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8',timerange='2017/08/07/13:08:35~2017/08/07/13:08:40',antenna='PM03')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8',timerange='2017/08/07/13:08:50~2017/08/07/13:09:00',antenna='PM03')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8',timerange='2017/08/07/13:09:55~2017/08/07/13:10:00',antenna='PM03')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8',timerange='2017/08/07/13:11:15~2017/08/07/13:11:20',antenna='DA46')
flagdata(vis=SB_iteration2_p5,mode='manual',spw='8',timerange='2017/08/07/13:11:35~2017/08/07/13:11:40',antenna='DA46')

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
    threshold = '0.2262mJy',
    **SB_tclean_wrapper_kwargs
)
estimate_SNR(SB_iteration2_cont_p5+'.image',disk_mask=SB_mask,noise_mask=noise_annulus_SB)
generate_image_png(
    SB_iteration2_cont_p5+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_iteration2_SB,10*rms_iteration2_SB],
    save_folder=SB_selfcal_iteration2_folder
)
#CQ_Tau_SB_contp5.image
#Beam 0.138 arcsec x 0.108 arcsec (-2.02 deg)
#Flux inside disk mask: 172.83 mJy
#Peak intensity of source: 8.05 mJy/beam
#rms: 5.32e-02 mJy/beam
#Peak SNR: 151.32

#CQ_Tau_SB_iteration2_contp5.image
#Beam 0.145 arcsec x 0.113 arcsec (-2.87 deg)
#Flux inside disk mask: 144.84 mJy
#Peak intensity of source: 7.63 mJy/beam
#rms: 3.79e-02 mJy/beam
#Peak SNR: 201.44

#Check again how SB phase-only selfcal improved things at each step
non_self_caled_SB_iteration2_vis = SB_iteration2_cont_p0
self_caled_SB_iteration2_visibilities = {
    'p1':SB_iteration2_cont_p1,
    'p2':SB_iteration2_cont_p2,
    'p3':SB_iteration2_cont_p3,
    'p4':SB_iteration2_cont_p4,
    'p5':SB_iteration2_cont_p5
}

#SB_EBs = ('EB0','EB1')
#SB_EB_spws = ('0,1,2,3,4,5,6,7','8,9,10,11,12,13,14,15') #fill out by referring to listobs output
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
        png_filename = f'iteration2_flux_comparison_SB_EB{i}_{self_cal_step}_to_{flux_ref_EB}.png'
        plot_label = os.path.join(SB_selfcal_iteration2_folder,png_filename)
        estimate_flux_scale(
            reference=f'{prefix}_{flux_ref_EB}_initcont_shift.vis.npz',
            comparison=exp_ms,incl=incl,PA=PA,plot_label=plot_label
        )

    fluxscale = [1.,]*number_of_EBs['SB']
    plot_label = os.path.join(SB_selfcal_iteration2_folder,f'deprojected_vis_profiles_SB_iteration2_{self_cal_step}.png')
    plot_deprojected(filelist=exported_ms,fluxscale=fluxscale,PA=PA,incl=incl,show_err=True,plot_label=plot_label)

#iteration_1                
#ratio          = [1.36589,1.36531,1.37131,1.37702,1.37817,1.37872] #CQ_Tau_SB_contp0...p5_EB0.vis.npz vs CQ_Tau_LB_EB1_initcont_shift.vis.npz
#scaling_factor = [1.169  ,1.168  ,1.171  ,1.173  ,1.174  ,1.174  ]
#ratio          = [0.91086,0.91196,0.94078,0.97047,0.99275,1.02373] #CQ_Tau_SB_contp0...p5_EB1.vis.npz vs CQ_Tau_LB_EB1_initcont_shift.vis.npz
#scaling_factor = [0.954  ,0.955  ,0.970  ,0.985  ,0.996  ,1.012  ]

#iteration_2
#ratio          = [1.03281,1.03233,1.03674,1.04108,1.04203,1.04247] #CQ_Tau_SB_contp0...p5_EB0.vis.npz vs CQ_Tau_LB_EB1_initcont_shift.vis.npz
#scaling_factor = [1.016  ,1.016  ,1.018  ,1.020  ,1.021  ,1.021  ]
#ratio          = [0.91086,0.91202,0.94070,0.97043,0.99265,1.02360] #CQ_Tau_SB_contp0...p5_EB1.vis.npz vs CQ_Tau_LB_EB1_initcont_shift.vis.npz
#scaling_factor = [0.954  ,0.955  ,0.970  ,0.985  ,0.996  ,1.012  ]

#Begin of SB+LB self-cal - iteration 2
#For phase self-cal, clean down to ~6sigma; for amplitude self-cal, clean down to ~1 sigma
LB_selfcal_iteration2_folder = get_figures_folderpath('8.1_selfcal_SBLB_iteration2_figures')
make_figures_folder(LB_selfcal_iteration2_folder)

LB_iteration2_cont_p0 = prefix+'_SBLB_iteration2_contp0'
os.system('rm -rf %s.ms*' %LB_iteration2_cont_p0)
concat(
    vis=[SB_iteration2_cont_p5+'.ms',f'{prefix}_LB_EB0_initcont_shift.ms',f'{prefix}_LB_EB1_initcont_shift.ms'],
    concatvis=LB_iteration2_cont_p0+'.ms',dirtol='0.1arcsec',copypointing=False
)
listobs(vis=LB_iteration2_cont_p0+'.ms',listfile=LB_iteration2_cont_p0+'.ms.listobs.txt',overwrite=True)

#Define new SB mask using the same centre as before (checked and agrees with the listobs one)
mask_pa        = PA
mask_semimajor = 0.85
mask_semiminor = 0.72
mask_ra  = '05h35m58.472722s'
mask_dec = '24d44m53.613450s'

LB_mask = f'ellipse[[{mask_ra},{mask_dec}], [{mask_semimajor:.3f}arcsec, {mask_semiminor:.3f}arcsec], {mask_pa:.1f}deg]'

noise_annulus_LB = f"annulus[[{mask_ra}, {mask_dec}],['4.arcsec', '6.arcsec']]"

LB_tclean_wrapper_kwargs = {
    'deconvolver':'multiscale','scales':[0,2,4,8,16,24],
    'smallscalebias':0.6,'gain':0.3,'cycleniter':300,
    'niter':1000000,'mask':LB_mask,
    'cellsize':'0.006arcsec','imsize':6400,
    'parallel':use_parallel,'savemodel':'modelcolumn',
    'robust':0.5,'interactive':False,
    'gridder':'standard',
}

tclean_wrapper(
    vis       = LB_iteration2_cont_p0+'.ms', 
    imagename = LB_iteration2_cont_p0,
    threshold = '0.0888mJy',
    **LB_tclean_wrapper_kwargs
)
estimate_SNR(LB_iteration2_cont_p0+'.image',disk_mask=LB_mask,noise_mask=noise_annulus_LB)
rms_iteration2_LB = imstat(imagename=LB_iteration2_cont_p0+'.image',region=noise_annulus_LB)['rms'][0]
generate_image_png(
    LB_iteration2_cont_p0+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_iteration2_LB,10*rms_iteration2_LB],
    save_folder=LB_selfcal_iteration2_folder
)
#CQ_Tau_SBLB_contp0.image
#Beam 0.078 arcsec x 0.058 arcsec (7.28 deg)
#Flux inside disk mask: 160.05 mJy
#Peak intensity of source: 2.43 mJy/beam
#rms: 1.75e-02 mJy/beam
#Peak SNR: 138.88

#CQ_Tau_SBLB_iteration2_contp0.image
#Beam 0.079 arcsec x 0.060 arcsec (6.77 deg)
#Flux inside disk mask: 146.20 mJy
#Peak intensity of source: 2.42 mJy/beam
#rms: 1.47e-02 mJy/beam
#Peak SNR: 164.87

#First round of phase-only self-cal
#NOTE: you need .p1 instead of _p1 in the caltable name if you want flagdata to work (i.e., flagging problematic antennas non interactively)...
LB_iteration2_p1 = prefix+'_SBLB_iteration2.p1'

os.system('rm -rf '+LB_iteration2_p1)
#If you get many flagged solutions, change gaintype to 'T'
gaincal(
    vis=LB_iteration2_cont_p0+'.ms',caltable=LB_iteration2_p1,
    gaintype='G',combine='scan,spw',calmode='p',solint='inf',
    spw=LB_contspws,refant=LB_refant,
    minsnr=2.,minblperant=4
)
#1 of 88 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:42:17.2
#1 of 94 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:08:25.8

#Inspect gain tables interactively and decide whether to manually flag something
#plotms(LB_iteration2_p1,xaxis='time',yaxis='GainPhase',iteraxis='spw',showgui=True)
#Print calibration png file
#plotms(
#    LB_iteration2_p1,xaxis='time', yaxis='GainPhase',iteraxis='spw',exprange='all',
#    overwrite=True,showgui=False,
#    plotfile=os.path.join(LB_iteration2_cont_p0,f'{prefix}_LB_iteration2_p1_gain_phase_vs_time.png')
#)

#Inspect gain tables and decide whether to flag something
plotms(
    LB_iteration2_p1,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(LB_selfcal_iteration2_folder,f'{prefix}_LB_iteration2_p1_gain_phase_vs_time.png'),
)

#Flag problematic antennas
#flagdata(vis=LB_iteration2_p1,mode='manual',spw='8',antenna='DV12')
#
#Inspect gain tables to check if flagging worked
#plotms(
#    LB_iteration2_p1,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
#    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
#    plotfile=os.path.join(LB_selfcal_iteration2_folder,f'{prefix}_LB_iteration2_p1_gain_phase_vs_time_flagged.png'),
#)

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
    threshold = '0.0882mJy',
    **LB_tclean_wrapper_kwargs
)
estimate_SNR(LB_iteration2_cont_p1+'.image',disk_mask=LB_mask,noise_mask=noise_annulus_LB)
generate_image_png(
    LB_iteration2_cont_p1+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_iteration2_LB,10*rms_iteration2_LB],
    save_folder=LB_selfcal_iteration2_folder
)
#CQ_Tau_SBLB_contp1.image
#Beam 0.078 arcsec x 0.058 arcsec (7.28 deg)
#Flux inside disk mask: 160.00 mJy
#Peak intensity of source: 2.45 mJy/beam
#rms: 1.75e-02 mJy/beam
#Peak SNR: 140.01

#CQ_Tau_SBLB_iteration2_contp1.image
#Beam 0.079 arcsec x 0.060 arcsec (6.77 deg)
#Flux inside disk mask: 146.24 mJy
#Peak intensity of source: 2.43 mJy/beam
#rms: 1.46e-02 mJy/beam
#Peak SNR: 166.40

#Second round of phase-only self-cal
LB_iteration2_p2 = LB_iteration2_p1.replace('p1','p2')
os.system('rm -rf '+LB_iteration2_p2)
gaincal(
    vis=LB_iteration2_cont_p1+'.ms',caltable=LB_iteration2_p2,
    gaintype='T',combine='scan,spw',calmode='p',solint='360s',
    spw=LB_contspws,refant=LB_refant,
    minsnr=2.,minblperant=4
)
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:17:38.3
#1 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:21:20.6
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:24:31.1
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:28:26.0
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:31:12.7
#1 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:34:14.5
#3 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:41:45.8
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:47:57.9
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:51:13.6
#2 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:54:42.4
#8 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:58:06.6
#2 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:03:06.9
#7 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:06:37.1
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:44:08.4
#1 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:47:49.5
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:51:04.5
#2 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:55:00.5
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:57:51.7
#2 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:00:58.7
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:08:45.3
#2 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:14:48.9
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:18:10.4
#2 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:21:43.0
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:25:09.5
#2 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:30:11.7
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:33:44.5

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
flagdata(vis=LB_iteration2_p2,mode='manual',spw='0', timerange='2015/08/30/11:41:00~2015/08/30/11:43:00',antenna='DV03')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='0', timerange='2015/08/30/11:48:30~2015/08/30/11:50:30',antenna='DV03,DV07')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='0', timerange='2015/08/30/11:54:00~2015/08/30/11:56:00',antenna='DV07')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='8', timerange='2017/08/07/12:39:00~2017/08/07/12:40:30',antenna='DV12,DA53')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='8', timerange='2017/08/07/12:42:00~2017/08/07/12:43:30',antenna='DV12')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='8', timerange='2017/08/07/12:45:30~2017/08/07/12:48:00',antenna='DV12,DV10,DV07,DA53,DA46')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='16',timerange='2017/11/20/07:17:30~2017/11/20/07:18:00',antenna='DV17')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='16',timerange='2017/11/20/07:31:00~2017/11/20/07:31:30',antenna='PM04,DV10,DA44')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='16',timerange='2017/11/20/07:58:00~2017/11/20/07:58:30',antenna='DV13,DA42')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='16',timerange='2017/11/20/08:06:00~2017/11/20/08:07:00',antenna='DA65')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='24',timerange='2017/11/23/04:47:30~2017/11/23/04:48:00',antenna='DA64')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='24',timerange='2017/11/23/04:57:30~2017/11/23/04:58:00',antenna='DV05')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='24',timerange='2017/11/23/05:08:30~2017/11/23/05:09:00',antenna='DV22')#DV16,DA65,
flagdata(vis=LB_iteration2_p2,mode='manual',spw='24',timerange='2017/11/23/05:14:30~2017/11/23/05:15:00',antenna='DV05')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='24',timerange='2017/11/23/05:18:00~2017/11/23/05:18:30',antenna='DV16')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='24',timerange='2017/11/23/05:30:00~2017/11/23/05:30:30',antenna='DV17')
flagdata(vis=LB_iteration2_p2,mode='manual',spw='24',timerange='2017/11/23/05:33:40~2017/11/23/05:34:00',antenna='DV16')

#Inspect gain tables to check if flagging worked
plotms(
    LB_iteration2_p2,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(LB_selfcal_iteration2_folder,f'{prefix}_LB_iteration2_p2_gain_phase_vs_time_flagged.png'),
)

#Apply the calibration gains and split-off the corrected ms
applycal(
    vis=LB_iteration2_cont_p1+'.ms',spw=LB_contspws,spwmap=LB_spw_mapping,
    gaintable=[LB_iteration2_p2],interp='linearPD',calwt=True,applymode='calonly'
)
LB_iteration2_cont_p2 = LB_iteration2_cont_p1.replace('p1','p2')
os.system('rm -rf %s.ms*' %LB_iteration2_cont_p2)
split(vis=LB_iteration2_cont_p1+'.ms',outputvis=LB_iteration2_cont_p2+'.ms',datacolumn='corrected')

#CLEAN, again with a ~6sigma threshold
tclean_wrapper(
    vis       = LB_iteration2_cont_p2+'.ms',
    imagename = LB_iteration2_cont_p2,
    threshold = '0.0834mJy',
    **LB_tclean_wrapper_kwargs
)
estimate_SNR(LB_iteration2_cont_p2+'.image',disk_mask=LB_mask,noise_mask=noise_annulus_LB)
generate_image_png(
    LB_iteration2_cont_p2+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_iteration2_LB,10*rms_iteration2_LB],
    save_folder=LB_selfcal_iteration2_folder
)
#CQ_Tau_SBLB_contp2.image
#Beam 0.078 arcsec x 0.058 arcsec (7.28 deg)
#Flux inside disk mask: 160.14 mJy
#Peak intensity of source: 2.49 mJy/beam
#rms: 1.68e-02 mJy/beam
#Peak SNR: 148.36

#CQ_Tau_SBLB_iteration2_contp2.image
#Beam 0.079 arcsec x 0.060 arcsec (6.77 deg)
#Flux inside disk mask: 146.36 mJy
#Peak intensity of source: 2.47 mJy/beam
#rms: 1.39e-02 mJy/beam
#Peak SNR: 178.20

#Third round of phase-only self-cal
LB_iteration2_p3 = LB_iteration2_p2.replace('p2','p3')
os.system('rm -rf '+LB_iteration2_p3)
gaincal(
    vis=LB_iteration2_cont_p2+'.ms',caltable=LB_iteration2_p3,
    gaintype='T',combine='scan,spw',calmode='p',solint='120s',
    spw=LB_contspws,refant=LB_refant,
    minsnr=2.,minblperant=4
)
#9 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/12:39:52.5
#1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/12:41:18.3
#1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/12:42:27.7
#1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/12:48:13.6
#1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/12:58:25.4
#1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/13:03:48.8
#1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/13:06:29.7
#1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/13:09:12.0
#1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/13:11:54.7
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:17:38.3
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:19:21.7
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:20:30.4
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:22:02.2
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:23:11.5
#2 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:24:40.6
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:26:03.7
#7 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:29:02.4
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:31:23.0
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:32:32.2
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:34:04.1
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:35:13.3
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:36:18.1
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:39:43.6
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:42:04.3
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:43:13.4
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:44:45.4
#7 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:45:54.5
#3 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:46:59.6
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:50:23.9
#3 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:52:44.6
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:53:53.7
#3 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:55:25.4
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:56:34.8
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:57:39.6
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:01:04.6
#2 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:03:25.0
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:04:34.3
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:05:50.4
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:44:08.4
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:45:59.3
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:47:02.6
#3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:48:42.4
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:49:45.9
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:51:16.3
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:52:41.6
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:55:43.3
#3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:58:05.5
#7 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:59:15.5
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:00:48.6
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:01:58.5
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:03:04.6
#3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:06:33.8
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:08:55.6
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:10:05.6
#2 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:11:38.4
#7 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:12:48.4
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:13:54.4
#2 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:17:22.3
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:19:44.2
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:20:54.0
#3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:22:27.1
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:23:36.6
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:24:42.3
#3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:28:09.0
#1 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:30:30.7
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:31:40.4
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:32:57.4

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
flagdata(vis=LB_iteration2_p3,mode='manual',spw='0', timerange='2015/08/30/11:41:30~2015/08/30/11:42:30',antenna='DV03')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='0', timerange='2015/08/30/11:49:30~2015/08/30/11:51:00',antenna='DV03,DV07')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='0', timerange='2015/08/30/11:54:00~2015/08/30/11:55:30',antenna='DV07')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='0', timerange='2015/08/30/11:55:30~2015/08/30/11:57:00',antenna='DV07')

flagdata(vis=LB_iteration2_p3,mode='manual',spw='8', timerange='2017/08/07/12:39:30~2017/08/07/12:39:40',antenna='DV12,DA53')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='8', timerange='2017/08/07/12:42:20~2017/08/07/12:42:40',antenna='DV11')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='8', timerange='2017/08/07/13:03:40~2017/08/07/13:04:00',antenna='PM03')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='8', timerange='2017/08/07/13:10:40~2017/08/07/13:11:00',antenna='PM03')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='8', timerange='2017/08/07/12:39:50~2017/08/07/12:40:00',antenna='DV12')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='8', timerange='2017/08/07/12:42:20~2017/08/07/12:42:30',antenna='PM03,DV12,DA46')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='8', timerange='2017/08/07/12:43:50~2017/08/07/12:44:00',antenna='DV12')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='8', timerange='2017/08/07/12:47:00~2017/08/07/12:47:10',antenna='DV12,DV10,DV07,DA53,DA46')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='8', timerange='2017/08/07/12:48:00~2017/08/07/12:48:20',antenna='DV23,DV07')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='8', timerange='2017/08/07/12:58:20~2017/08/07/12:58:30',antenna='PM03')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='8', timerange='2017/08/07/13:06:00~2017/08/07/13:07:00',antenna='PM04,DA53,DA46')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='8', timerange='2017/08/07/13:08:00~2017/08/07/13:08:30',antenna='PM03')

flagdata(vis=LB_iteration2_p3,mode='manual',spw='16',timerange='2017/11/20/07:31:20~2017/11/20/07:31:30',antenna='DV23')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='16',timerange='2017/11/20/07:32:30~2017/11/20/07:32:40',antenna='DV16')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='16',timerange='2017/11/20/07:42:00~2017/11/20/07:42:10',antenna='DV17,DV23')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='16',timerange='2017/11/20/07:46:50~2017/11/20/07:47:00',antenna='DA65')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='16',timerange='2017/11/20/07:53:40~2017/11/20/07:54:00',antenna='DA60,DA44,DA42')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='16',timerange='2017/11/20/07:57:30~2017/11/20/07:57:50',antenna='DA60,PM04,DV13,DA42')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='16',timerange='2017/11/20/08:01:00~2017/11/20/08:01:10',antenna='DV23,DV22,DV15,PM03,DV03')

flagdata(vis=LB_iteration2_p3,mode='manual',spw='24',timerange='2017/11/23/04:45:30~2017/11/23/04:46:30',antenna='DA64')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='24',timerange='2017/11/23/04:47:00~2017/11/23/04:47:30',antenna='DA65')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='24',timerange='2017/11/23/04:52:40~2017/11/23/04:53:00',antenna='DA46')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='24',timerange='2017/11/23/05:01:55~2017/11/23/05:02:05',antenna='DV05')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='24',timerange='2017/11/23/05:06:00~2017/11/23/05:07:00',antenna='DV22')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='24',timerange='2017/11/23/05:08:40~2017/11/23/05:09:00',antenna='DV22,DV11')
flagdata(vis=LB_iteration2_p3,mode='manual',spw='24',timerange='2017/11/23/05:17:20~2017/11/23/05:17:40',antenna='DV22,DA60')

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
    threshold = '0.0756mJy',
    **LB_tclean_wrapper_kwargs
)
estimate_SNR(LB_iteration2_cont_p3+'.image',disk_mask=LB_mask,noise_mask=noise_annulus_LB)
generate_image_png(
    LB_iteration2_cont_p3+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_iteration2_LB,10*rms_iteration2_LB],
    save_folder=LB_selfcal_iteration2_folder
)
#CQ_Tau_SBLB_contp3.image
#Beam 0.078 arcsec x 0.058 arcsec (7.28 deg)
#Flux inside disk mask: 160.64 mJy
#Peak intensity of source: 2.60 mJy/beam
#rms: 1.55e-02 mJy/beam
#Peak SNR: 167.30

#CQ_Tau_SBLB_iteration2_contp3.image
#Beam 0.079 arcsec x 0.060 arcsec (6.77 deg)
#Flux inside disk mask: 146.81 mJy
#Peak intensity of source: 2.58 mJy/beam
#rms: 1.26e-02 mJy/beam
#Peak SNR: 205.64

#Fourth round of phase-only self-cal
LB_iteration2_p4 = LB_iteration2_p3.replace('p3','p4')
os.system('rm -rf '+LB_iteration2_p4)
gaincal(
    vis=LB_iteration2_cont_p3+'.ms',caltable=LB_iteration2_p4,
    gaintype='T',combine='scan,spw',calmode='p',solint='60s',
    spw=LB_contspws,refant=LB_refant,
    minsnr=2.,minblperant=4
)
#9 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/12:39:52.5
#1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/12:42:07.5
#1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/12:47:53.9
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:17:38.3
#3 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:18:50.4
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:20:10.6
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:21:31.1
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:22:51.8
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:24:12.0
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:25:32.8
#8 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:26:09.7
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:28:19.5
#8 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:29:31.1
#2 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:30:51.7
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:32:12.2
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:33:32.8
#3 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:34:53.3
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:36:13.9
#7 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:36:50.8
#3 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:39:00.7
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:40:12.3
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:41:32.6
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:42:53.4
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:44:14.1
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:45:34.7
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:46:55.5
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:47:32.1
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:49:41.0
#3 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:50:53.0
#2 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:52:13.3
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:53:33.9
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:54:54.2
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:56:14.8
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:57:35.5
#8 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:58:12.6
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:00:21.6
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:01:33.1
#3 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:02:53.7
#2 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:04:14.2
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:05:34.8
#7 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:06:37.1
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:44:08.4
#3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:45:21.0
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:46:42.4
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:48:04.2
#2 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:49:25.8
#3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:50:47.4
#2 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:52:08.1
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:52:46.6
#3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:54:59.9
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:56:12.3
#3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:57:33.5
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:58:55.3
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:00:16.9
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:01:38.3
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:02:58.9
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:03:37.3
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:05:50.0
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:07:02.6
#2 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:08:24.0
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:09:45.4
#7 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:11:06.8
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:12:28.2
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:13:48.7
#8 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:14:27.1
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:16:39.0
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:17:51.2
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:19:11.6
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:20:33.8
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:21:55.4
#3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:23:16.5
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:24:43.4
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:25:15.6
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:27:25.7
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:28:37.8
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:29:59.1
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:31:20.3
#7 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:32:41.6
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:33:44.5

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
flagdata(vis=LB_iteration2_p4,mode='manual',spw='0', timerange='2015/08/30/11:41:50~2015/08/30/11:42:00',antenna='DV03')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='0', timerange='2015/08/30/11:49:00~2015/08/30/11:50:00',antenna='DV03,DV07')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='0', timerange='2015/08/30/11:55:00~2015/08/30/11:55:30',antenna='DV07')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='0', timerange='2015/08/30/11:56:00~2015/08/30/11:56:30',antenna='DV07')

flagdata(vis=LB_iteration2_p4,mode='manual',spw='8', timerange='2017/08/07/12:39:50~2017/08/07/12:40:00',antenna='DV12')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='8', timerange='2017/08/07/12:42:00~2017/08/07/12:42:10',antenna='DV12,PM03,DV23')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='8', timerange='2017/08/07/12:43:20~2017/08/07/12:43:40',antenna='DV12,DV11,DA46')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='8', timerange='2017/08/07/12:46:00~2017/08/07/12:47:00',antenna='DV12,DV10,DV07,DA46')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='8', timerange='2017/08/07/12:47:30~2017/08/07/12:48:30',antenna='DV07,DA53')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='8', timerange='2017/08/07/13:07:00~2017/08/07/13:08:00',antenna='PM04,DV12')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='8', timerange='2017/08/07/13:08:00~2017/08/07/13:09:00',antenna='PM03')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='8', timerange='2017/08/07/13:10:00~2017/08/07/13:11:00',antenna='PM03,DV12')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='8', timerange='2017/08/07/13:11:00~2017/08/07/13:12:00',antenna='DV12,DA46')

flagdata(vis=LB_iteration2_p4,mode='manual',spw='16',timerange='2017/11/20/07:17:30~2017/11/20/07:18:00',antenna='DV05,DA64')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='16',timerange='2017/11/20/07:18:30~2017/11/20/07:19:00',antenna='DV22')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='16',timerange='2017/11/20/07:22:40~2017/11/20/07:23:00',antenna='DV16')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='16',timerange='2017/11/20/07:28:10~2017/11/20/07:28:30',antenna='DA51')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='16',timerange='2017/11/20/07:30:30~2017/11/20/07:31:00',antenna='DV17,PM04,DV23')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='16',timerange='2017/11/20/07:32:00~2017/11/20/07:32:30',antenna='DV17')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='16',timerange='2017/11/20/07:33:30~2017/11/20/07:33:40',antenna='DV05')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='16',timerange='2017/11/20/07:38:50~2017/11/20/07:39:10',antenna='DV22,DA64')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='16',timerange='2017/11/20/07:40:00~2017/11/20/07:40:30',antenna='DA65')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='16',timerange='2017/11/20/07:42:50~2017/11/20/07:43:00',antenna='DV17,DV15,DA60,DV23')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='16',timerange='2017/11/20/07:46:50~2017/11/20/07:47:00',antenna='DV17,DV03,DA65')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='16',timerange='2017/11/20/07:47:20~2017/11/20/07:47:40',antenna='DV22,DV03')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='16',timerange='2017/11/20/07:50:30~2017/11/20/07:51:00',antenna='DV22')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='16',timerange='2017/11/20/07:52:00~2017/11/20/07:52:30',antenna='DV17')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='16',timerange='2017/11/20/07:53:30~2017/11/20/07:53:40',antenna='DA44')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='16',timerange='2017/11/20/07:56:00~2017/11/20/07:56:30',antenna='DV17,DA65')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='16',timerange='2017/11/20/07:57:30~2017/11/20/07:58:00',antenna='DV17,DA64,PM04,DV13,DA42,DA65')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='16',timerange='2017/11/20/07:58:00~2017/11/20/07:58:30',antenna='DV22,DV13,DA42')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='16',timerange='2017/11/20/08:00:00~2017/11/20/08:01:00',antenna='DV16')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='16',timerange='2017/11/20/08:01:00~2017/11/20/08:02:00',antenna='DV22,DV15,PM03')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='16',timerange='2017/11/20/08:02:40~2017/11/20/08:03:00',antenna='DV05,PM03')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='16',timerange='2017/11/20/08:04:00~2017/11/20/08:04:20',antenna='DA64,DA60,PM03')

flagdata(vis=LB_iteration2_p4,mode='manual',spw='24',timerange='2017/11/23/04:45:20~2017/11/23/04:45:30',antenna='DV17')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='24',timerange='2017/11/23/04:50:40~2017/11/23/04:51:00',antenna='DV05,DA60')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='24',timerange='2017/11/23/04:52:00~2017/11/23/04:52:10',antenna='DV16,DA46')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='24',timerange='2017/11/23/04:52:40~2017/11/23/04:52:50',antenna='DA65,DA46')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='24',timerange='2017/11/23/04:57:10~2017/11/23/04:57:50',antenna='DV17,DV05,DA65')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='24',timerange='2017/11/23/04:58:50~2017/11/23/04:59:00',antenna='DV17,DV05')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='24',timerange='2017/11/23/05:02:50~2017/11/23/05:03:10',antenna='DA65')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='24',timerange='2017/11/23/05:16:35~2017/11/23/05:16:45',antenna='DV22,DA44')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='24',timerange='2017/11/23/05:17:50~2017/11/23/05:18:00',antenna='DA65,DA60,DA44')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='24',timerange='2017/11/23/05:19:10~2017/11/23/05:19:15',antenna='DV16')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='24',timerange='2017/11/23/05:20:30~2017/11/23/05:20:40',antenna='DV16,DA65')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='24',timerange='2017/11/23/05:23:10~2017/11/23/05:23:20',antenna='DV22')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='24',timerange='2017/11/23/05:25:10~2017/11/23/05:25:20',antenna='DA64')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='24',timerange='2017/11/23/05:28:30~2017/11/23/05:28:40',antenna='DV22')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='24',timerange='2017/11/23/05:31:10~2017/11/23/05:31:30',antenna='DA65,DA60')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='24',timerange='2017/11/23/05:32:30~2017/11/23/05:33:00',antenna='DV15,PM04,DV19,DA44')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='24',timerange='2017/11/23/05:33:40~2017/11/23/05:34:00',antenna='DV16')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='24',timerange='2017/11/23/05:23:10~2017/11/23/05:23:20',antenna='DA60')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='24',timerange='2017/11/23/05:21:50~2017/11/23/05:22:00',antenna='DV05')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='24',timerange='2017/11/23/05:01:30~2017/11/23/05:01:40',antenna='DA60')
flagdata(vis=LB_iteration2_p4,mode='manual',spw='24',timerange='2017/11/23/05:27:20~2017/11/23/05:27:30',antenna='DV22')

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
    threshold = '0.0732mJy',
    **LB_tclean_wrapper_kwargs
)
estimate_SNR(LB_iteration2_cont_p4+'.image',disk_mask=LB_mask,noise_mask=noise_annulus_LB)
generate_image_png(
    LB_iteration2_cont_p4+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_iteration2_LB,10*rms_iteration2_LB],
    save_folder=LB_selfcal_iteration2_folder
)
#CQ_Tau_SBLB_contp4.image
#Beam 0.078 arcsec x 0.058 arcsec (7.28 deg)
#Flux inside disk mask: 160.79 mJy
#Peak intensity of source: 2.67 mJy/beam
#rms: 1.51e-02 mJy/beam
#Peak SNR: 176.52

#CQ_Tau_SBLB_iteration2_contp4.image
#Beam 0.079 arcsec x 0.060 arcsec (6.77 deg)
#Flux inside disk mask: 146.97 mJy
#Peak intensity of source: 2.63 mJy/beam
#rms: 1.22e-02 mJy/beam
#Peak SNR: 214.82

#Fifth round of phase-only self-cal
LB_iteration2_p5 = LB_iteration2_p4.replace('p4','p5')
os.system('rm -rf '+LB_iteration2_p5)
gaincal(
    vis=LB_iteration2_cont_p4+'.ms',caltable=LB_iteration2_p5,
    gaintype='T',combine='spw',calmode='p',solint='30s',
    spw=LB_contspws,refant=LB_refant,
    minsnr=2.,minblperant=4
)

# 1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/12:39:32.3
# 3 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/12:44:54.9
# 2 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/12:52:41.2
# 3 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/12:54:27.8
# 1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/12:59:14.0
# 1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/13:11:22.5
# 3 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/13:13:39.7
# 6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:17:35.2
# 8 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:17:53.4
# 3 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:18:38.3
# 4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:19:04.8
# 4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:19:58.5
# 6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:20:25.4
# 5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:21:19.0
# 5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:21:45.9
# 6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:22:39.7
# 5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:23:06.5
# 3 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:23:59.9
# 5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:24:27.0
# 6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:25:20.4
# 6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:25:38.5
# 5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:26:03.7
# 5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:28:16.5
#10 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:28:34.6
# 4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:29:19.0
# 5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:29:46.1
# 4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:30:39.6
# 3 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:31:06.6
# 4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:32:00.1
# 5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:32:27.2
# 7 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:33:20.7
# 5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:33:47.7
# 2 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:34:41.2
# 6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:35:08.3
# 6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:36:01.6
# 5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:36:19.7
# 4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:36:44.8
# 6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:38:57.7
#10 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:39:15.8
# 6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:40:00.3
# 4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:40:27.3
# 6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:41:20.8
# 6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:41:47.6
# 4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:42:41.3
# 3 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:43:08.4
# 3 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:44:02.0
# 5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:44:28.9
# 6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:45:22.6
# 5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:45:49.5
# 4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:46:43.1
# 6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:47:01.0
# 5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:47:26.0
# 5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:49:38.0
#10 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:49:56.1
# 3 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:50:40.9
# 6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:51:07.6
# 4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:52:01.2
# 2 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:52:28.1
# 5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:53:21.8
# 5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:53:48.7
# 6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:54:42.1
# 5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:55:09.2
# 5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:56:02.7
# 3 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:56:29.8
# 5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:57:23.1
# 8 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:57:41.2
# 5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:58:06.6
# 7 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:00:18.6
# 9 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:00:36.7
# 6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:01:21.0
# 5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:01:48.2
# 6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:02:41.6
# 4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:03:08.7
# 1 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:04:02.1
# 5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:04:29.3
# 2 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:05:22.7
# 4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:05:49.8
# 7 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:06:37.1
# 6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:44:05.4
# 9 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:44:23.6
# 5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:45:08.9
# 4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:45:36.1
# 4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:46:30.5
# 4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:46:57.6
# 6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:47:52.1
# 7 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:48:19.3
# 4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:49:13.7
# 3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:49:40.9
# 5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:50:35.3
# 6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:51:02.5
# 3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:51:56.9
# 4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:52:15.0
# 5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:52:41.6
# 3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:54:56.9
# 9 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:55:15.0
# 4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:56:00.2
# 7 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:56:27.5
# 6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:57:21.7
# 6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:57:48.6
# 5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:58:43.2
# 4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:59:10.5
# 7 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:00:04.8
# 4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:00:32.0
# 5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:01:26.3
# 5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:01:53.5
# 3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:02:47.8
#10 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:03:05.9
# 5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:03:32.3
# 5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:05:47.2
# 9 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:06:05.2
# 5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:06:50.5
# 4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:07:17.7
# 2 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:08:11.9
# 5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:08:39.1
# 6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:09:33.3
# 6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:10:00.5
# 5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:10:54.7
# 4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:11:21.9
# 4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:12:16.1
# 7 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:12:43.3
# 6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:13:37.5
# 7 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:13:55.7
# 4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:14:22.0
# 6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:16:36.0
# 6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:16:54.1
# 6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:17:39.1
# 3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:18:06.3
# 5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:19:00.4
# 5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:19:26.8
# 6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:20:21.7
# 6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:20:49.0
# 3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:21:43.3
# 6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:22:10.3
# 7 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:23:04.4
# 3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:23:31.6
# 3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:24:25.7
# 8 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:24:43.8
# 7 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:25:09.5
# 5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:27:22.7
# 6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:27:40.8
# 5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:28:25.7
# 4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:28:52.9
# 6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:29:47.0
# 3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:30:14.2
# 2 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:31:08.2
# 5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:31:35.4
# 5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:32:29.5
# 6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:32:56.6
# 6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:33:44.5

#Inspect gain tables interactively and decide whether to manually flag something
#plotms(LB_iteration2_p5,xaxis='time',yaxis='GainPhase',iteraxis='spw',showgui=True)
#Print calibration png file
#plotms(
#    LB_iteration2_p5,xaxis='time', yaxis='GainPhase',iteraxis='spw',exprange='all',
#    overwrite=True,showgui=False,
#    plotfile=os.path.join(LB_selfcal_iteration2_folder,f'{prefix}_LB_iteration2_p5_gain_phase_vs_time.png')
#)

#Inspect gain tables and decide whether to flag something
plotms(
    LB_iteration2_p5,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(LB_selfcal_iteration2_folder,f'{prefix}_LB_iteration2_p5_gain_phase_vs_time.png'),
)

#Flag problematic antennas
flagdata(vis=LB_iteration2_p5,mode='manual',spw='0',antenna='DV07')

flagdata(vis=LB_iteration2_p5,mode='manual',spw='0', timerange='2015/08/30/11:39:50~2015/08/30/11:40:10',antenna='DV03')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='0', timerange='2015/08/30/11:41:30~2015/08/30/11:41:40',antenna='DV03')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='0', timerange='2015/08/30/11:41:50~2015/08/30/11:42:00',antenna='DV03')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='0', timerange='2015/08/30/11:43:20~2015/08/30/11:43:30',antenna='DV03')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='0', timerange='2015/08/30/11:48:50~2015/08/30/11:49:10',antenna='DV03')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='0', timerange='2015/08/30/11:49:20~2015/08/30/11:49:40',antenna='DV03')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='0', timerange='2015/08/30/11:49:50~2015/08/30/11:50:10',antenna='DV03')

flagdata(vis=LB_iteration2_p5,mode='manual',spw='8', timerange='2017/08/07/12:39:45~2017/08/07/12:39:55',antenna='DV12')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='8', timerange='2017/08/07/12:41:50~2017/08/07/12:42:00',antenna='DV12,DV23,PM03')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='8', timerange='2017/08/07/12:42:20~2017/08/07/12:42:30',antenna='DV12,PM03')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='8', timerange='2017/08/07/12:43:10~2017/08/07/12:43:20',antenna='DA46,DV11,DV12')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='8', timerange='2017/08/07/12:46:00~2017/08/07/12:46:30',antenna='DV07,DV12')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='8', timerange='2017/08/07/12:46:40~2017/08/07/12:47:00',antenna='DV12')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='8', timerange='2017/08/07/12:47:30~2017/08/07/12:47:50',antenna='DA53,DV07,DV12')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='8', timerange='2017/08/07/12:48:00~2017/08/07/12:48:20',antenna='DV07,DV23')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='8', timerange='2017/08/07/13:04:50~2017/08/07/13:05:10',antenna='DV12')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='8', timerange='2017/08/07/13:05:50~2017/08/07/13:06:10',antenna='DA44')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='8', timerange='2017/08/07/13:07:00~2017/08/07/13:07:20',antenna='DV12,PM04')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='8', timerange='2017/08/07/13:08:30~2017/08/07/13:08:50',antenna='PM03')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='8', timerange='2017/08/07/13:08:50~2017/08/07/13:09:10',antenna='PM03')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='8', timerange='2017/08/07/13:09:50~2017/08/07/13:10:10',antenna='DV12,PM03')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='8', timerange='2017/08/07/13:11:20~2017/08/07/13:11:40',antenna='DA46,DV12')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='8', timerange='2017/08/07/12:46:00~2017/08/07/12:46:30',antenna='DV10,DA46')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='8', timerange='2017/08/07/12:46:30~2017/08/07/12:47:00',antenna='DV10,DV07,DA53,DA46')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='8', timerange='2017/08/07/12:44:50~2017/08/07/12:45:00',antenna='DV10')

flagdata(vis=LB_iteration2_p5,mode='manual',spw='16',timerange='2017/11/20/07:17:50~2017/11/20/07:17:55',antenna='DA64')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='16',timerange='2017/11/20/07:18:30~2017/11/20/07:18:40',antenna='DV17')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='16',timerange='2017/11/20/07:23:00~2017/11/20/07:23:10',antenna='DV22')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='16',timerange='2017/11/20/07:23:50~2017/11/20/07:24:00',antenna='DV22,DV17')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='16',timerange='2017/11/20/07:25:15~2017/11/20/07:25:25',antenna='DA60')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='16',timerange='2017/11/20/07:26:00~2017/11/20/07:26:05',antenna='DA64,DA65')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='16',timerange='2017/11/20/07:28:10~2017/11/20/07:28:20',antenna='DV22,DV15,DA60,DA51')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='16',timerange='2017/11/20/07:28:30~2017/11/20/07:28:40',antenna='DA51')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='16',timerange='2017/11/20/07:29:15~2017/11/20/07:29:25',antenna='DV22,DV23,DV13')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='16',timerange='2017/11/20/07:30:30~2017/11/20/07:30:40',antenna='DV17,DV23')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='16',timerange='2017/11/20/07:31:00~2017/11/20/07:31:10',antenna='DV17,DV16,PM04,DV23')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='16',timerange='2017/11/20/07:32:20~2017/11/20/07:32:30',antenna='DV17')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='16',timerange='2017/11/20/07:33:45~2017/11/20/07:33:50',antenna='DV17,DV15,DA64,DV10')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='16',timerange='2017/11/20/07:34:40~2017/11/20/07:34:50',antenna='DV22')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='16',timerange='2017/11/20/07:35:00~2017/11/20/07:35:10',antenna='DV17')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='16',timerange='2017/11/20/07:36:40~2017/11/20/07:36:50',antenna='DV22')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='16',timerange='2017/11/20/07:38:55~2017/11/20/07:39:00',antenna='DV23')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='16',timerange='2017/11/20/07:41:40~2017/11/20/07:41:50',antenna='DV15,DA60,DV23')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='16',timerange='2017/11/20/07:42:40~2017/11/20/07:42:50',antenna='DV17,DV15,DA60,DV23')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='16',timerange='2017/11/20/07:43:05~2017/11/20/07:43:10',antenna='DV15,DV05,DA65,DA60,DV23')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='16',timerange='2017/11/20/07:45:45~2017/11/20/07:45:50',antenna='DA65,DA60')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='16',timerange='2017/11/20/07:46:40~2017/11/20/07:46:50',antenna='DV17,DV03,DA65')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='16',timerange='2017/11/20/07:46:50~2017/11/20/07:47:10',antenna='DV03')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='16',timerange='2017/11/20/07:47:20~2017/11/20/07:47:30',antenna='DV15,DV03,PM04,DV23')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='16',timerange='2017/11/20/07:49:35~2017/11/20/07:49:40',antenna='DV22')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='16',timerange='2017/11/20/07:50:40~2017/11/20/07:50:45',antenna='DV22,DV17,DV16')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='16',timerange='2017/11/20/07:52:00~2017/11/20/07:52:10',antenna='DV17,DV05')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='16',timerange='2017/11/20/07:52:20~2017/11/20/07:52:30',antenna='DV17,DA65,DA64')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='16',timerange='2017/11/20/07:53:20~2017/11/20/07:53:30',antenna='DV17,DA60,PM04,DV13,DA44')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='16',timerange='2017/11/20/07:53:45~2017/11/20/07:53:50',antenna='DA60,DA44')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='16',timerange='2017/11/20/07:56:00~2017/11/20/07:56:10',antenna='DV17')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='16',timerange='2017/11/20/07:56:25~2017/11/20/07:56:35',antenna='DV05,DA65')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='16',timerange='2017/11/20/07:57:20~2017/11/20/07:57:30',antenna='DV22,DA60,PM04,DV13,DA42')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='16',timerange='2017/11/20/07:57:40~2017/11/20/07:57:50',antenna='DA49,PM04,DV23,DV13')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='16',timerange='2017/11/20/07:58:00~2017/11/20/07:58:20',antenna='DV05,DA65,DA64,DV13,DA42')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='16',timerange='2017/11/20/08:00:30~2017/11/20/08:00:40',antenna='DV05')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='16',timerange='2017/11/20/08:01:20~2017/11/20/08:01:30',antenna='DV15,DA60,DA49,PM03')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='16',timerange='2017/11/20/08:01:40~2017/11/20/08:01:50',antenna='DV22,DV15,DA49,PM03')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='16',timerange='2017/11/20/08:02:40~2017/11/20/08:02:50',antenna='PM03')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='16',timerange='2017/11/20/08:03:00~2017/11/20/08:03:10',antenna='PM03')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='16',timerange='2017/11/20/08:04:00~2017/11/20/08:04:10',antenna='DV16,DA64,DA60,PM03')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='16',timerange='2017/11/20/08:04:20~2017/11/20/08:04:30',antenna='DV22,DA60')

flagdata(vis=LB_iteration2_p5,mode='manual',spw='24',timerange='2017/11/23/04:44:00~2017/11/23/04:44:10',antenna='DA64')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='24',timerange='2017/11/23/04:45:00~2017/11/23/04:45:20',antenna='DV16,DA60')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='24',timerange='2017/11/23/04:46:50~2017/11/23/04:47:00',antenna='DV17')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='24',timerange='2017/11/23/04:49:10~2017/11/23/04:49:20',antenna='DV22')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='24',timerange='2017/11/23/04:49:40~2017/11/23/04:49:50',antenna='DV16,DA64')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='24',timerange='2017/11/23/04:50:30~2017/11/23/04:50:40',antenna='DA60')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='24',timerange='2017/11/23/04:51:50~2017/11/23/04:52:00',antenna='DV16,DA46')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='24',timerange='2017/11/23/04:52:10~2017/11/23/04:52:20',antenna='DV17,DA64,DA46')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='24',timerange='2017/11/23/04:52:40~2017/11/23/04:52:50',antenna='DA46')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='24',timerange='2017/11/23/05:00:30~2017/11/23/05:00:40',antenna='DV22,DA64')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='24',timerange='2017/11/23/05:01:20~2017/11/23/05:01:30',antenna='DA60')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='24',timerange='2017/11/23/05:02:40~2017/11/23/05:03:00',antenna='DV16,DA60')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='24',timerange='2017/11/23/05:03:20~2017/11/23/05:03:40',antenna='DV15,DA60')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='24',timerange='2017/11/23/05:05:45~2017/11/23/05:05:50',antenna='DV22,DV17')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='24',timerange='2017/11/23/05:06:05~2017/11/23/05:06:10',antenna='DV17')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='24',timerange='2017/11/23/05:06:50~2017/11/23/05:07:00',antenna='DV22')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='24',timerange='2017/11/23/05:07:10~2017/11/23/05:07:20',antenna='DV16,DV15,PM03')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='24',timerange='2017/11/23/05:09:30~2017/11/23/05:09:40',antenna='DV22')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='24',timerange='2017/11/23/05:11:20~2017/11/23/05:11:25',antenna='DV22')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='24',timerange='2017/11/23/05:13:30~2017/11/23/05:13:40',antenna='DA60')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='24',timerange='2017/11/23/05:13:50~2017/11/23/05:14:00',antenna='DV17')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='24',timerange='2017/11/23/05:16:35~2017/11/23/05:16:40',antenna='DA44')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='24',timerange='2017/11/23/05:16:35~2017/11/23/05:16:40',antenna='DA44')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='24',timerange='2017/11/23/05:17:30~2017/11/23/05:18:00',antenna='DA60,DA44')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='24',timerange='2017/11/23/05:18:00~2017/11/23/05:18:30',antenna='DV22,DV05,DA60')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='24',timerange='2017/11/23/05:18:50~2017/11/23/05:19:10',antenna='DV05')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='24',timerange='2017/11/23/05:19:20~2017/11/23/05:19:40',antenna='DV17,DV15,DV10')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='24',timerange='2017/11/23/05:20:20~2017/11/23/05:20:40',antenna='DA65')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='24',timerange='2017/11/23/05:21:30~2017/11/23/05:22:00',antenna='DA65')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='24',timerange='2017/11/23/05:22:00~2017/11/23/05:22:30',antenna='DA64')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='24',timerange='2017/11/23/05:23:30~2017/11/23/05:23:40',antenna='DV22,DA60')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='24',timerange='2017/11/23/05:27:20~2017/11/23/05:27:30',antenna='DV17,DV05')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='24',timerange='2017/11/23/05:27:40~2017/11/23/05:27:50',antenna='DV22,DV17,DV05')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='24',timerange='2017/11/23/05:28:20~2017/11/23/05:28:30',antenna='DA60')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='24',timerange='2017/11/23/05:28:50~2017/11/23/05:29:00',antenna='DV22,DA60')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='24',timerange='2017/11/23/05:31:05~2017/11/23/05:31:10',antenna='DV16,DA64')
flagdata(vis=LB_iteration2_p5,mode='manual',spw='24',timerange='2017/11/23/05:32:25~2017/11/23/05:32:35',antenna='DV05,DA65')

#Inspect gain tables to check if flagging worked
plotms(
    LB_iteration2_p5,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(LB_selfcal_iteration2_folder,f'{prefix}_LB_iteration2_p5_gain_phase_vs_time_flagged.png'),
)

#Apply the calibration gains and split-off the corrected ms
applycal(
    vis=LB_iteration2_cont_p4+'.ms',spw=LB_contspws,spwmap=LB_spw_mapping,
    gaintable=[LB_iteration2_p5],interp='linearPD',calwt=True,applymode='calonly'
)
LB_iteration2_cont_p5 = LB_iteration2_cont_p4.replace('p4','p5')
os.system('rm -rf %s.ms*' %LB_iteration2_cont_p5)
split(vis=LB_iteration2_cont_p4+'.ms',outputvis=LB_iteration2_cont_p5+'.ms',datacolumn='corrected')

#CLEAN, again with a ~6sigma threshold
tclean_wrapper(
    vis       = LB_iteration2_cont_p5+'.ms',
    imagename = LB_iteration2_cont_p5,
    threshold = '0.0732mJy',
    **LB_tclean_wrapper_kwargs
)
estimate_SNR(LB_iteration2_cont_p5+'.image',disk_mask=LB_mask,noise_mask=noise_annulus_LB)
generate_image_png(
    LB_iteration2_cont_p5+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_iteration2_LB,10*rms_iteration2_LB],
    save_folder=LB_selfcal_iteration2_folder
)
#CQ_Tau_SBLB_contp5.image
#Beam 0.078 arcsec x 0.058 arcsec (7.28 deg)
#Flux inside disk mask: 161.08 mJy
#Peak intensity of source: 2.69 mJy/beam
#rms: 1.50e-02 mJy/beam
#Peak SNR: 179.52

#CQ_Tau_SBLB_iteration2_contp5.image
#Beam 0.079 arcsec x 0.060 arcsec (6.77 deg)
#Flux inside disk mask: 147.16 mJy
#Peak intensity of source: 2.66 mJy/beam
#rms: 1.22e-02 mJy/beam
#Peak SNR: 218.24

#Sixth round of phase-only self-cal
LB_iteration2_p6 = LB_iteration2_p5.replace('p5','p6')
os.system('rm -rf '+LB_iteration2_p6)
gaincal(
    vis=LB_iteration2_cont_p5+'.ms',caltable=LB_iteration2_p6,
    gaintype='T',combine='spw',calmode='p',solint='18s',
    spw=LB_contspws,refant=LB_refant,
    minsnr=2.,minblperant=4
)
#1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/12:41:04.4
#1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/12:41:49.4
#1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/12:42:25.7
#2 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/12:43:46.3
#1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/12:48:11.6
#1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/12:49:14.0
#1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/12:52:17.1
#1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/12:54:21.8
#1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/12:58:23.4
#1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/12:59:07.9
#1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/12:59:25.8
#1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/13:08:53.8
#1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/13:09:12.0
#1 of 40 solutions flagged due to SNR < 2 in spw=8  at 2017/08/07/13:11:52.7
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:17:29.2
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:17:47.3
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:18:32.2
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:18:49.7
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:19:07.9
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:19:52.4
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:20:10.3
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:20:28.4
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:21:13.0
#7 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:21:30.8
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:21:49.0
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:22:33.6
#3 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:22:51.4
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:23:09.5
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:23:53.8
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:24:11.9
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:24:30.0
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:25:14.4
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:25:32.4
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:26:03.7
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:28:10.4
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:28:28.6
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:29:13.0
#8 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:29:31.0
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:29:49.1
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:30:33.6
#2 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:30:51.5
#2 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:31:09.6
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:31:54.1
#3 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:32:12.0
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:32:30.2
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:33:14.7
#3 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:33:32.6
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:33:50.7
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:34:35.2
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:34:53.1
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:35:11.3
#7 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:35:55.6
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:36:13.7
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:36:44.8
#7 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:38:51.6
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:39:09.8
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:39:54.2
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:40:12.2
#7 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:40:30.3
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:41:14.7
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:41:32.7
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:41:50.6
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:42:35.3
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:42:53.3
#3 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:43:11.4
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:43:55.9
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:44:13.8
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:44:32.0
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:45:16.5
#7 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:45:34.4
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:45:52.5
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:46:37.0
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:46:54.9
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:47:26.0
#7 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:49:31.9
#8 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:49:50.1
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:50:34.9
#7 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:50:52.5
#2 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:51:10.6
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:51:55.2
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:52:13.0
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:52:31.2
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:53:15.7
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:53:33.6
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:53:51.7
#7 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:54:36.1
#7 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:54:54.1
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:55:12.2
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:55:56.7
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:56:14.6
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:56:32.8
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:57:17.0
#7 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:57:35.2
#7 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/07:58:06.6
#7 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:00:12.5
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:00:30.7
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:01:15.0
#7 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:01:33.1
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:01:51.2
#6 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:02:35.6
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:02:53.6
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:03:11.8
#3 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:03:56.1
#4 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:04:14.2
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:04:32.3
#2 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:05:16.6
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:05:34.7
#5 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:05:52.8
#7 of 44 solutions flagged due to SNR < 2 in spw=16 at 2017/11/20/08:06:37.1
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:43:59.4
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:44:17.5
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:45:02.8
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:45:21.0
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:45:39.1
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:46:24.4
#3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:46:42.6
#3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:47:00.6
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:47:46.0
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:48:04.2
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:48:22.3
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:49:07.6
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:49:25.8
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:49:43.9
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:50:29.2
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:50:47.4
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:51:05.5
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:51:50.8
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:52:09.0
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:52:41.6
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:54:50.8
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:55:09.0
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:55:54.2
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:56:12.3
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:56:30.5
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:57:15.7
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:57:33.8
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:57:51.7
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:58:37.2
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:58:55.3
#3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:59:13.5
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/04:59:58.8
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:00:16.8
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:00:35.0
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:01:20.2
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:01:38.4
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:01:56.5
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:02:41.7
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:02:59.9
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:03:32.3
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:05:41.2
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:05:59.1
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:06:44.4
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:07:02.6
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:07:20.7
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:08:05.9
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:08:24.0
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:08:42.1
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:09:27.3
#7 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:09:45.4
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:10:03.6
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:10:48.7
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:11:06.8
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:11:25.0
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:12:10.1
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:12:28.2
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:12:46.4
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:13:31.5
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:13:49.6
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:14:22.0
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:16:29.9
#2 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:16:48.0
#7 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:17:33.1
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:17:51.2
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:18:09.4
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:18:54.4
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:19:12.5
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:19:29.9
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:20:15.7
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:20:33.8
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:20:52.0
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:21:37.2
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:21:55.2
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:22:13.3
#8 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:22:58.3
#3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:23:16.5
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:23:34.6
#4 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:24:19.6
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:24:37.8
#7 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:25:09.5
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:27:16.6
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:27:34.8
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:28:19.7
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:28:37.8
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:28:56.0
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:29:41.0
#2 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:29:59.0
#3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:30:17.2
#3 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:31:02.1
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:31:20.3
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:31:38.4
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:32:23.4
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:32:41.5
#5 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:32:59.6
#6 of 47 solutions flagged due to SNR < 2 in spw=24 at 2017/11/23/05:33:44.5

#Inspect gain tables interactively and decide whether to manually flag something
#plotms(LB_iteration2_p6,xaxis='time',yaxis='GainPhase',iteraxis='spw',showgui=True)
#Print calibration png file
#plotms(
#    LB_iteration2_p6,xaxis='time', yaxis='GainPhase',iteraxis='spw',exprange='all',
#    overwrite=True,showgui=False,
#    plotfile=os.path.join(LB_selfcal_iteration2_folder,f'{prefix}_LB_iteration2_p6_gain_phase_vs_time.png')
#)

#Inspect gain tables and decide whether to flag something
plotms(
    LB_iteration2_p6,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(LB_selfcal_iteration2_folder,f'{prefix}_LB_iteration2_p6_gain_phase_vs_time.png'),
)

#Flag problematic antennas
flagdata(vis=LB_iteration2_p6,mode='manual',spw='0',antenna='DV03,DV07')

flagdata(vis=LB_iteration2_p6,mode='manual',spw='8', timerange='2017/08/07/12:39:40~2017/08/07/12:39:50',antenna='DV12')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='8', timerange='2017/08/07/12:41:40~2017/08/07/12:41:50',antenna='DV23,DV12,PM03')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='8', timerange='2017/08/07/12:42:00~2017/08/07/12:42:10',antenna='PM03,DV23,DV12')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='8', timerange='2017/08/07/12:42:20~2017/08/07/12:42:30',antenna='PM03,DV12')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='8', timerange='2017/08/07/12:43:00~2017/08/07/12:43:20',antenna='DV12,DV11,DA46')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='8', timerange='2017/08/07/12:43:20~2017/08/07/12:43:40',antenna='DA46,DV11')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='8', timerange='2017/08/07/12:46:10~2017/08/07/12:46:20',antenna='DV12,DV10,DV07,DA46')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='8', timerange='2017/08/07/12:46:30~2017/08/07/12:46:40',antenna='DV10,DV07,DA46')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='8', timerange='2017/08/07/12:46:50~2017/08/07/12:47:00',antenna='DV12,DV10,DV07,DA46')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='8', timerange='2017/08/07/12:47:30~2017/08/07/12:47:40',antenna='DV12,DV07,DA53')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='8', timerange='2017/08/07/12:48:10~2017/08/07/12:48:20',antenna='DV23,DV07')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='8', timerange='2017/08/07/13:05:00~2017/08/07/13:05:10',antenna='DV12')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='8', timerange='2017/08/07/13:05:50~2017/08/07/13:06:00',antenna='DA44')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='8', timerange='2017/08/07/13:07:10~2017/08/07/13:07:20',antenna='PM04,DV12')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='8', timerange='2017/08/07/13:08:30~2017/08/07/13:08:40',antenna='PM03')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='8', timerange='2017/08/07/13:08:50~2017/08/07/13:09:00',antenna='PM03')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='8', timerange='2017/08/07/13:09:50~2017/08/07/13:10:00',antenna='PM03,DV12,DA53')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='8', timerange='2017/08/07/13:10:10~2017/08/07/13:10:20',antenna='PM03,DV12')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='8', timerange='2017/08/07/13:11:10~2017/08/07/13:11:20',antenna='DV12,DA46')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='8', timerange='2017/08/07/13:11:30~2017/08/07/13:11:40',antenna='PM03,DV12,DA46')

flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/07:17:40~2017/11/20/07:17:50',antenna='DV17')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/07:18:25~2017/11/20/07:18:35',antenna='DV17,DV05,PM04')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/07:20:05~2017/11/20/07:20:15',antenna='DV22')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/07:22:30~2017/11/20/07:22:40',antenna='DV16')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/07:22:45~2017/11/20/07:22:55',antenna='DV16,DA65,DA64')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/07:23:05~2017/11/20/07:23:15',antenna='DV17')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/07:25:10~2017/11/20/07:25:20',antenna='DV16')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/07:25:30~2017/11/20/07:25:40',antenna='DV17,DV15,DA65')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/07:28:05~2017/11/20/07:28:15',antenna='DA65,DA64,DA51')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/07:28:25~2017/11/20/07:28:35',antenna='DV17,DA51')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/07:29:05~2017/11/20/07:29:15',antenna='DV17,DV23')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/07:29:25~2017/11/20/07:29:35',antenna='DV23')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/07:29:45~2017/11/20/07:29:55',antenna='DV22,DA65')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/07:30:25~2017/11/20/07:30:35',antenna='DV05,DA60,DV23')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/07:30:45~2017/11/20/07:30:55',antenna='DV05,PM04,DV23')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/07:31:05~2017/11/20/07:31:15',antenna='DV17,DV05,PM04,DV23')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/07:31:50~2017/11/20/07:32:00',antenna='DA64')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/07:32:05~2017/11/20/07:32:15',antenna='DV22,DA65')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/07:32:25~2017/11/20/07:32:35',antenna='DV17,DA65')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/07:33:45~2017/11/20/07:33:55',antenna='DV17,DV15,PM04,DV10')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/07:36:10~2017/11/20/07:36:20',antenna='DA65')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/07:36:40~2017/11/20/07:36:50',antenna='DV22,DV17')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/07:38:45~2017/11/20/07:38:55',antenna='DV23')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/07:39:05~2017/11/20/07:39:15',antenna='DV22')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/07:41:10~2017/11/20/07:41:20',antenna='DV17')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/07:41:45~2017/11/20/07:41:55',antenna='DV15,DA60')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/07:42:30~2017/11/20/07:42:40',antenna='DV17,DV15,DA60,DV23')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/07:42:50~2017/11/20/07:43:00',antenna='DV17,DV15,DV23')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/07:43:10~2017/11/20/07:43:20',antenna='DV22,DA64,DV23')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/07:43:50~2017/11/20/07:44:00',antenna='DV15,DV23')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/07:44:10~2017/11/20/07:44:20',antenna='DV05')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/07:45:50~2017/11/20/07:46:00',antenna='DV05,DA60')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/07:46:30~2017/11/20/07:46:40',antenna='DV05,DA65,DV03')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/07:46:50~2017/11/20/07:47:00',antenna='DV17,DV03,PM03')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/07:47:20~2017/11/20/07:47:30',antenna='DV15,DV03,DA64,PM04,DV23')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/07:51:05~2017/11/20/07:51:15',antenna='DV16')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/07:51:50~2017/11/20/07:52:00',antenna='DA64')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/07:52:25~2017/11/20/07:52:35',antenna='DA65')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/07:53:30~2017/11/20/07:53:35',antenna='DA44')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/07:53:50~2017/11/20/07:53:55',antenna='DA60,DA44')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/07:56:30~2017/11/20/07:56:40',antenna='DA65')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/07:57:10~2017/11/20/07:57:20',antenna='DV22,DA42')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/07:57:30~2017/11/20/07:57:40',antenna='DA49,PM04,DV13,DA42')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/07:58:00~2017/11/20/07:58:30',antenna='DV16,DV23,DV13,DA42')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/08:00:10~2017/11/20/08:00:20',antenna='DA51,DV13,DA42')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/08:00:25~2017/11/20/08:00:35',antenna='DV16,DV15,DV23')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/08:01:10~2017/11/20/08:01:20',antenna='DV15,DA64,DA60,PM03,DV23')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/08:01:30~2017/11/20/08:01:40',antenna='DV15,DA49,PM03')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/08:01:50~2017/11/20/08:02:00',antenna='DV22,DV15,DA49,PM03')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/08:02:30~2017/11/20/08:02:40',antenna='PM03')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/08:02:50~2017/11/20/08:03:00',antenna='PM03')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/08:03:10~2017/11/20/08:03:20',antenna='DV05,PM03')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/08:03:50~2017/11/20/08:04:00',antenna='DA64,PM03')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/08:04:10~2017/11/20/08:04:20',antenna='DA64,DA60,PM03')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/08:04:30~2017/11/20/08:04:40',antenna='DV22,DA60')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/08:05:30~2017/11/20/08:05:40',antenna='DV16')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='16',timerange='2017/11/20/08:05:50~2017/11/20/08:06:00',antenna='DA64')

flagdata(vis=LB_iteration2_p6,mode='manual',spw='24',timerange='2017/11/23/04:44:10~2017/11/23/04:44:20',antenna='DA65,DA64')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='24',timerange='2017/11/23/04:46:40~2017/11/23/04:46:50',antenna='DV22,DV05')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='24',timerange='2017/11/23/04:47:00~2017/11/23/04:47:10',antenna='DA64')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='24',timerange='2017/11/23/04:48:00~2017/11/23/04:48:10',antenna='DA65')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='24',timerange='2017/11/23/04:50:45~2017/11/23/04:50:50',antenna='DV22,DV05,DA60')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='24',timerange='2017/11/23/04:51:50~2017/11/23/04:51:55',antenna='DV05,DA46')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='24',timerange='2017/11/23/04:52:05~2017/11/23/04:52:15',antenna='DV17,DV16,DA46')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='24',timerange='2017/11/23/04:52:40~2017/11/23/04:52:45',antenna='DA46')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='24',timerange='2017/11/23/04:56:30~2017/11/23/04:56:35',antenna='DA65')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='24',timerange='2017/11/23/04:57:30~2017/11/23/04:57:35',antenna='DV22,DA65')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='24',timerange='2017/11/23/04:58:35~2017/11/23/04:58:40',antenna='DV05')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='24',timerange='2017/11/23/05:00:15~2017/11/23/05:00:20',antenna='DA60')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='24',timerange='2017/11/23/05:00:30~2017/11/23/05:00:40',antenna='DV22,DA65')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='24',timerange='2017/11/23/05:01:15~2017/11/23/05:01:25',antenna='DA60')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='24',timerange='2017/11/23/05:01:50~2017/11/23/05:02:00',antenna='DV16,DA65')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='24',timerange='2017/11/23/05:02:35~2017/11/23/05:02:45',antenna='DV16,DA60')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='24',timerange='2017/11/23/05:02:55~2017/11/23/05:03:00',antenna='DA60')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='24',timerange='2017/11/23/05:03:30~2017/11/23/05:03:35',antenna='DV15,DA65,DA60')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='24',timerange='2017/11/23/05:05:40~2017/11/23/05:05:45',antenna='DV05')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='24',timerange='2017/11/23/05:05:55~2017/11/23/05:06:00',antenna='DV16')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='24',timerange='2017/11/23/05:06:40~2017/11/23/05:06:45',antenna='DV22,PM03')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='24',timerange='2017/11/23/05:12:40~2017/11/23/05:12:50',antenna='DA60')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='24',timerange='2017/11/23/05:13:20~2017/11/23/05:13:40',antenna='DV05')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='24',timerange='2017/11/23/05:13:45~2017/11/23/05:13:55',antenna='DV17')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='24',timerange='2017/11/23/05:16:25~2017/11/23/05:16:35',antenna='DA44')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='24',timerange='2017/11/23/05:18:05~2017/11/23/05:18:15',antenna='DV22,DV16')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='24',timerange='2017/11/23/05:18:50~2017/11/23/05:19:00',antenna='DV17,DA65')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='24',timerange='2017/11/23/05:20:10~2017/11/23/05:20:20',antenna='DA65')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='24',timerange='2017/11/23/05:20:30~2017/11/23/05:20:40',antenna='DV22')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='24',timerange='2017/11/23/05:20:50~2017/11/23/05:21:00',antenna='DV05')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='24',timerange='2017/11/23/05:21:50~2017/11/23/05:22:00',antenna='DV05')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='24',timerange='2017/11/23/05:22:10~2017/11/23/05:22:15',antenna='DA65,DA60')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='24',timerange='2017/11/23/05:23:10~2017/11/23/05:23:20',antenna='DV22,DA64')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='24',timerange='2017/11/23/05:23:30~2017/11/23/05:23:40',antenna='DA60')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='24',timerange='2017/11/23/05:24:15~2017/11/23/05:24:25',antenna='DV17,DA65,DA64')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='24',timerange='2017/11/23/05:24:35~2017/11/23/05:24:40',antenna='DV17')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='24',timerange='2017/11/23/05:27:10~2017/11/23/05:27:20',antenna='DV22,DV17')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='24',timerange='2017/11/23/05:27:30~2017/11/23/05:27:40',antenna='DV22,DV17')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='24',timerange='2017/11/23/05:28:15~2017/11/23/05:28:25',antenna='DV05')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='24',timerange='2017/11/23/05:28:50~2017/11/23/05:29:00',antenna='DV22')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='24',timerange='2017/11/23/05:29:35~2017/11/23/05:29:45',antenna='DV17,DA60')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='24',timerange='2017/11/23/05:29:55~2017/11/23/05:30:05',antenna='DV22,DV05,DA60')
flagdata(vis=LB_iteration2_p6,mode='manual',spw='24',timerange='2017/11/23/05:33:40~2017/11/23/05:33:50',antenna='DV16')

#Inspect gain tables to check if flagging worked
plotms(
    LB_iteration2_p6,xaxis='time',yaxis='GainPhase',overwrite=True,showgui=False,exprange='all',
    customsymbol=True,customflaggedsymbol=True,flaggedsymbolshape='autoscaling',iteraxis='spw',
    plotfile=os.path.join(LB_selfcal_iteration2_folder,f'{prefix}_LB_iteration2_p6_gain_phase_vs_time_flagged.png'),
)

#Apply the calibration gains and split-off the corrected ms
applycal(
    vis=LB_iteration2_cont_p5+'.ms',spw=LB_contspws,spwmap=LB_spw_mapping,
    gaintable=[LB_iteration2_p6],interp='linearPD',calwt=True,applymode='calonly'
)
LB_iteration2_cont_p6 = LB_iteration2_cont_p5.replace('p5','p6')
os.system('rm -rf %s.ms*' %LB_iteration2_cont_p6)
split(vis=LB_iteration2_cont_p5+'.ms',outputvis=LB_iteration2_cont_p6+'.ms',datacolumn='corrected')

#CLEAN, again with a ~6sigma threshold
tclean_wrapper(
    vis       = LB_iteration2_cont_p6+'.ms',
    imagename = LB_iteration2_cont_p6,
    threshold = '0.0732mJy',
    **LB_tclean_wrapper_kwargs
)
estimate_SNR(LB_iteration2_cont_p6+'.image',disk_mask=LB_mask,noise_mask=noise_annulus_LB)
generate_image_png(
    LB_iteration2_cont_p6+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_iteration2_LB,10*rms_iteration2_LB],
    save_folder=LB_selfcal_iteration2_folder
)
#CQ_Tau_SBLB_contp6.image
#Beam 0.078 arcsec x 0.058 arcsec (7.28 deg)
#Flux inside disk mask: 161.20 mJy
#Peak intensity of source: 2.71 mJy/beam
#rms: 1.49e-02 mJy/beam
#Peak SNR: 181.98

#CQ_Tau_SBLB_iteration2_contp6.image
#Beam 0.079 arcsec x 0.060 arcsec (6.77 deg)
#Flux inside disk mask: 147.34 mJy
#Peak intensity of source: 2.71 mJy/beam
#rms: 1.22e-02 mJy/beam
#Peak SNR: 222.13

#CLEAN with a ~1sigma threshold before amplitude self-cal
tclean_wrapper(
    vis       = LB_iteration2_cont_p6+'.ms',
    imagename = LB_iteration2_cont_p6,
    threshold = '0.0122mJy',
    **LB_tclean_wrapper_kwargs
)
estimate_SNR(LB_iteration2_cont_p6+'.image',disk_mask=LB_mask,noise_mask=noise_annulus_LB)
generate_image_png(
    LB_iteration2_cont_p6+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_iteration2_LB,10*rms_iteration2_LB],
    save_folder=LB_selfcal_iteration2_folder
)
#CQ_Tau_SBLB_iteration2_contp6.image
#Beam 0.079 arcsec x 0.060 arcsec (6.77 deg)
#Flux inside disk mask: 148.27 mJy
#Peak intensity of source: 2.66 mJy/beam
#rms: 1.20e-02 mJy/beam
#Peak SNR: 222.35

#First round of amplitude self-cal
LB_iteration2_ap0 = prefix+'_SBLB_iteration2.ap0'
os.system('rm -rf '+LB_iteration2_ap0)
gaincal(
    vis=LB_iteration2_cont_p6+'.ms',caltable=LB_iteration2_ap0,
    gaintype='T',combine='scan,spw',calmode='ap',solint='inf',
    spw=LB_contspws,refant=LB_refant,
    minsnr=5.0,minblperant=4,solnorm=False
)

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
flagdata(vis=LB_iteration2_ap0,mode='manual',spw='0', antenna='DV03,DV07')
flagdata(vis=LB_iteration2_ap0,mode='manual',spw='8', antenna='DA53,DV12,DV23,PM03,PM04')
flagdata(vis=LB_iteration2_ap0,mode='manual',spw='16',antenna='DV22,DV17,DV16')
flagdata(vis=LB_iteration2_ap0,mode='manual',spw='24',antenna='DV16,DA60')

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
    vis=LB_iteration2_cont_p6+'.ms',spw=LB_contspws,spwmap=LB_spw_mapping,
    gaintable=[LB_iteration2_ap0],interp='linearPD',calwt=True,applymode='calonly'
)
LB_iteration2_cont_ap0 = prefix+'_SBLB_iteration2_contap0'
os.system('rm -rf %s.ms*' %LB_iteration2_cont_ap0)
split(vis=LB_iteration2_cont_p6+'.ms',outputvis=LB_iteration2_cont_ap0+'.ms',datacolumn='corrected')

#CLEAN with a ~1sigma threshold before amplitude self-cal
tclean_wrapper(
    vis=LB_iteration2_cont_ap0+'.ms',
    imagename = LB_iteration2_cont_ap0,
    threshold = '0.0120mJy',
    **LB_tclean_wrapper_kwargs
)
estimate_SNR(LB_iteration2_cont_ap0+'.image',disk_mask=LB_mask,noise_mask=noise_annulus_LB)
generate_image_png(
    LB_iteration2_cont_ap0+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_iteration2_LB,10*rms_iteration2_LB],
    save_folder=LB_selfcal_iteration2_folder
)
#CQ_Tau_SBLB_iteration2_contap0.image
#Beam 0.081 arcsec x 0.059 arcsec (12.15 deg)
#Flux inside disk mask: 147.67 mJy
#Peak intensity of source: 2.67 mJy/beam
#rms: 1.19e-02 mJy/beam
#Peak SNR: 223.66

#Try ampl self-cal on scan length intervals
LB_iteration2_ap1 = prefix+'_SBLB_iteration2.ap1'
os.system('rm -rf '+LB_iteration2_ap1)
gaincal(
    vis=LB_iteration2_cont_ap0+'.ms',caltable=LB_iteration2_ap1,
    gaintype='T',combine='spw',calmode='ap',solint='inf',
    spw=LB_contspws,refant=LB_refant,
    minsnr=5.0,minblperant=4,solnorm=False
)
#1 of 40 solutions flagged due to SNR < 5 in spw=8  at 2017/08/07/12:39:35.3
#2 of 40 solutions flagged due to SNR < 5 in spw=8  at 2017/08/07/12:42:07.5
#1 of 40 solutions flagged due to SNR < 5 in spw=8  at 2017/08/07/12:43:28.3
#1 of 40 solutions flagged due to SNR < 5 in spw=8  at 2017/08/07/12:44:39.9
#3 of 40 solutions flagged due to SNR < 5 in spw=8  at 2017/08/07/12:47:54.0
#1 of 40 solutions flagged due to SNR < 5 in spw=8  at 2017/08/07/12:52:26.2
#1 of 40 solutions flagged due to SNR < 5 in spw=8  at 2017/08/07/12:55:24.5
#1 of 40 solutions flagged due to SNR < 5 in spw=8  at 2017/08/07/12:58:05.7
#1 of 40 solutions flagged due to SNR < 5 in spw=8  at 2017/08/07/12:59:26.1
#1 of 40 solutions flagged due to SNR < 5 in spw=8  at 2017/08/07/13:02:08.2
#1 of 40 solutions flagged due to SNR < 5 in spw=8  at 2017/08/07/13:07:12.0
#1 of 40 solutions flagged due to SNR < 5 in spw=8  at 2017/08/07/13:10:14.3
#7 of 44 solutions flagged due to SNR < 5 in spw=16 at 2017/11/20/07:17:38.3
#7 of 44 solutions flagged due to SNR < 5 in spw=16 at 2017/11/20/07:18:50.4
#7 of 44 solutions flagged due to SNR < 5 in spw=16 at 2017/11/20/07:20:10.6
#7 of 44 solutions flagged due to SNR < 5 in spw=16 at 2017/11/20/07:21:31.1
#7 of 44 solutions flagged due to SNR < 5 in spw=16 at 2017/11/20/07:22:51.8
#7 of 44 solutions flagged due to SNR < 5 in spw=16 at 2017/11/20/07:24:12.0
#7 of 44 solutions flagged due to SNR < 5 in spw=16 at 2017/11/20/07:25:23.4
#9 of 44 solutions flagged due to SNR < 5 in spw=16 at 2017/11/20/07:26:03.7
#7 of 44 solutions flagged due to SNR < 5 in spw=16 at 2017/11/20/07:28:19.5
#8 of 44 solutions flagged due to SNR < 5 in spw=16 at 2017/11/20/07:29:31.1
#9 of 44 solutions flagged due to SNR < 5 in spw=16 at 2017/11/20/07:30:51.8
#7 of 44 solutions flagged due to SNR < 5 in spw=16 at 2017/11/20/07:32:12.2
#8 of 44 solutions flagged due to SNR < 5 in spw=16 at 2017/11/20/07:33:32.8
#7 of 44 solutions flagged due to SNR < 5 in spw=16 at 2017/11/20/07:34:53.3
#7 of 44 solutions flagged due to SNR < 5 in spw=16 at 2017/11/20/07:36:04.6
#8 of 44 solutions flagged due to SNR < 5 in spw=16 at 2017/11/20/07:36:44.8
#8 of 44 solutions flagged due to SNR < 5 in spw=16 at 2017/11/20/07:39:00.7
#7 of 44 solutions flagged due to SNR < 5 in spw=16 at 2017/11/20/07:40:12.3
#7 of 44 solutions flagged due to SNR < 5 in spw=16 at 2017/11/20/07:41:32.7
#9 of 44 solutions flagged due to SNR < 5 in spw=16 at 2017/11/20/07:42:53.4
#7 of 44 solutions flagged due to SNR < 5 in spw=16 at 2017/11/20/07:44:14.2
#7 of 44 solutions flagged due to SNR < 5 in spw=16 at 2017/11/20/07:45:34.7
#7 of 44 solutions flagged due to SNR < 5 in spw=16 at 2017/11/20/07:46:46.2
#8 of 44 solutions flagged due to SNR < 5 in spw=16 at 2017/11/20/07:47:26.0
#7 of 44 solutions flagged due to SNR < 5 in spw=16 at 2017/11/20/07:49:41.0
#7 of 44 solutions flagged due to SNR < 5 in spw=16 at 2017/11/20/07:50:53.1
#7 of 44 solutions flagged due to SNR < 5 in spw=16 at 2017/11/20/07:52:13.3
#7 of 44 solutions flagged due to SNR < 5 in spw=16 at 2017/11/20/07:53:33.9
#7 of 44 solutions flagged due to SNR < 5 in spw=16 at 2017/11/20/07:54:54.2
#6 of 44 solutions flagged due to SNR < 5 in spw=16 at 2017/11/20/07:56:14.8
#8 of 44 solutions flagged due to SNR < 5 in spw=16 at 2017/11/20/07:57:26.1
#7 of 44 solutions flagged due to SNR < 5 in spw=16 at 2017/11/20/07:58:06.6
#7 of 44 solutions flagged due to SNR < 5 in spw=16 at 2017/11/20/08:00:21.6
#9 of 44 solutions flagged due to SNR < 5 in spw=16 at 2017/11/20/08:01:33.1
#7 of 44 solutions flagged due to SNR < 5 in spw=16 at 2017/11/20/08:02:53.7
#7 of 44 solutions flagged due to SNR < 5 in spw=16 at 2017/11/20/08:04:14.2
#6 of 44 solutions flagged due to SNR < 5 in spw=16 at 2017/11/20/08:05:34.8
#7 of 44 solutions flagged due to SNR < 5 in spw=16 at 2017/11/20/08:06:37.1
#7 of 47 solutions flagged due to SNR < 5 in spw=24 at 2017/11/23/04:44:08.4
#5 of 47 solutions flagged due to SNR < 5 in spw=24 at 2017/11/23/04:45:21.0
#5 of 47 solutions flagged due to SNR < 5 in spw=24 at 2017/11/23/04:46:42.4
#7 of 47 solutions flagged due to SNR < 5 in spw=24 at 2017/11/23/04:48:04.2
#5 of 47 solutions flagged due to SNR < 5 in spw=24 at 2017/11/23/04:49:25.8
#7 of 47 solutions flagged due to SNR < 5 in spw=24 at 2017/11/23/04:50:47.4
#7 of 47 solutions flagged due to SNR < 5 in spw=24 at 2017/11/23/04:51:59.9
#7 of 47 solutions flagged due to SNR < 5 in spw=24 at 2017/11/23/04:52:41.6
#5 of 47 solutions flagged due to SNR < 5 in spw=24 at 2017/11/23/04:54:59.9
#7 of 47 solutions flagged due to SNR < 5 in spw=24 at 2017/11/23/04:56:12.3
#7 of 47 solutions flagged due to SNR < 5 in spw=24 at 2017/11/23/04:57:33.4
#5 of 47 solutions flagged due to SNR < 5 in spw=24 at 2017/11/23/04:58:55.3
#7 of 47 solutions flagged due to SNR < 5 in spw=24 at 2017/11/23/05:00:16.9
#5 of 47 solutions flagged due to SNR < 5 in spw=24 at 2017/11/23/05:01:38.3
#7 of 47 solutions flagged due to SNR < 5 in spw=24 at 2017/11/23/05:02:50.8
#9 of 47 solutions flagged due to SNR < 5 in spw=24 at 2017/11/23/05:03:32.3
#7 of 47 solutions flagged due to SNR < 5 in spw=24 at 2017/11/23/05:05:49.9
#7 of 47 solutions flagged due to SNR < 5 in spw=24 at 2017/11/23/05:07:02.6
#6 of 47 solutions flagged due to SNR < 5 in spw=24 at 2017/11/23/05:08:24.0
#7 of 47 solutions flagged due to SNR < 5 in spw=24 at 2017/11/23/05:09:45.4
#7 of 47 solutions flagged due to SNR < 5 in spw=24 at 2017/11/23/05:11:06.8
#7 of 47 solutions flagged due to SNR < 5 in spw=24 at 2017/11/23/05:12:28.2
#7 of 47 solutions flagged due to SNR < 5 in spw=24 at 2017/11/23/05:13:40.6
#9 of 47 solutions flagged due to SNR < 5 in spw=24 at 2017/11/23/05:14:22.0
#7 of 47 solutions flagged due to SNR < 5 in spw=24 at 2017/11/23/05:16:39.0
#7 of 47 solutions flagged due to SNR < 5 in spw=24 at 2017/11/23/05:17:51.2
#7 of 47 solutions flagged due to SNR < 5 in spw=24 at 2017/11/23/05:19:11.6
#7 of 47 solutions flagged due to SNR < 5 in spw=24 at 2017/11/23/05:20:33.8
#7 of 47 solutions flagged due to SNR < 5 in spw=24 at 2017/11/23/05:21:55.4
#6 of 47 solutions flagged due to SNR < 5 in spw=24 at 2017/11/23/05:23:16.5
#7 of 47 solutions flagged due to SNR < 5 in spw=24 at 2017/11/23/05:24:28.7
#8 of 47 solutions flagged due to SNR < 5 in spw=24 at 2017/11/23/05:25:09.5
#7 of 47 solutions flagged due to SNR < 5 in spw=24 at 2017/11/23/05:27:25.7
#7 of 47 solutions flagged due to SNR < 5 in spw=24 at 2017/11/23/05:28:37.8
#7 of 47 solutions flagged due to SNR < 5 in spw=24 at 2017/11/23/05:29:59.2
#7 of 47 solutions flagged due to SNR < 5 in spw=24 at 2017/11/23/05:31:20.3
#7 of 47 solutions flagged due to SNR < 5 in spw=24 at 2017/11/23/05:32:41.6
#8 of 47 solutions flagged due to SNR < 5 in spw=24 at 2017/11/23/05:33:44.5

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
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='0', timerange='2015/08/30/11:45:00~2015/08/30/11:48:00',antenna='DV09,DV07')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='0', timerange='2015/08/30/11:54:00~2015/08/30/11:57:00',antenna='DV07')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='8', timerange='2017/08/07/12:43:25~2017/08/07/12:43:35',antenna='DA46,DV11')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='8', timerange='2017/08/07/12:46:30~2017/08/07/12:46:40',antenna='DA46,DV07,DV10')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='8', timerange='2017/08/07/13:07:10~2017/08/07/13:07:20',antenna='DV12,PM04')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='8', timerange='2017/08/07/13:10:10~2017/08/07/13:10:20',antenna='DV12')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='16',timerange='2017/11/20/07:28:10~2017/11/20/07:28:30',antenna='DV15,DA51')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='16',timerange='2017/11/20/07:30:40~2017/11/20/07:31:00',antenna='DV23')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='16',timerange='2017/11/20/07:47:20~2017/11/20/07:47:40',antenna='DV15')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='16',timerange='2017/11/20/07:53:20~2017/11/20/07:54:00',antenna='DA44')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='16',timerange='2017/11/20/07:57:00~2017/11/20/07:58:00',antenna='DA42')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='16',timerange='2017/11/20/07:58:00~2017/11/20/07:58:20',antenna='PM04,DV13,DA42')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='16',timerange='2017/11/20/08:02:40~2017/11/20/08:03:00',antenna='PM03')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='24',timerange='2017/11/23/04:45:10~2017/11/23/04:45:30',antenna='DV16')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='24',timerange='2017/11/23/04:49:20~2017/11/23/04:49:30',antenna='DV16')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='24',timerange='2017/11/23/04:51:40~2017/11/23/04:52:20',antenna='DA46')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='24',timerange='2017/11/23/04:52:20~2017/11/23/04:53:00',antenna='DA46')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='24',timerange='2017/11/23/04:54:50~2017/11/23/04:55:10',antenna='DA65')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='24',timerange='2017/11/23/05:01:30~2017/11/23/05:01:40',antenna='DV16')
flagdata(vis=LB_iteration2_ap1,mode='manual',spw='24',timerange='2017/11/23/05:16:20~2017/11/23/05:17:00',antenna='DA44')

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
    threshold = '0.0118mJy',
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

#Check again how LB phase-only selfcal improved things at each step
non_self_caled_LB_iteration2_vis = LB_iteration2_cont_p0
self_caled_LB_iteration2_visibilities = {
    'p1' :LB_iteration2_cont_p1 ,
    'p2' :LB_iteration2_cont_p2 ,
    'p3' :LB_iteration2_cont_p3 ,
    'p4' :LB_iteration2_cont_p4 ,
    'p5' :LB_iteration2_cont_p5 ,
    'p6' :LB_iteration2_cont_p6 ,
    'ap0':LB_iteration2_cont_ap0,
    'ap1':LB_iteration2_cont_ap1,
}

#for vis in self_caled_LB_iteration2_visibilities.values(): 
#    listobs(vis=vis+'.ms',listfile=vis+'.ms.listobs.txt',overwrite=True)
#
#LB_EBs = ('EB0','EB1','EB2','EB3')
#LB_EB_spws = ('0,1,2,3,4,5,6,7', '8,9,10,11,12,13,14,15', '16,17,18,19,20,21,22,23', '24,25,26,27,28,29,30,31') #fill out by referring to listobs output
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
SBLB_flux_ref_EB = 3 #this is LB_EB1

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

#iteration_1
#ratio          = [1.36589,1.36531,1.37131,1.37702,1.37817,1.37872] #CQ_Tau_SB_contp0...p5_EB0.vis.npz vs CQ_Tau_LB_EB1_initcont_shift.vis.npz
#scaling_factor = [1.169  ,1.168  ,1.171  ,1.173  ,1.174  ,1.174  ]
#ratio          = [0.91086,0.91196,0.94078,0.97047,0.99275,1.02373] #CQ_Tau_SB_contp0...p5_EB1.vis.npz vs CQ_Tau_LB_EB1_initcont_shift.vis.npz
#scaling_factor = [0.954  ,0.955  ,0.970  ,0.985  ,0.996  ,1.012  ]

#ratio          = [1.37869,1.37799,1.37359,1.34796,1.33492,1.32590,1.32198] #CQ_Tau_SB_contp0...p6_EB0.vis.npz vs CQ_Tau_LB_contp0...p6_EB1.vis.npz
#scaling_factor = [1.174  ,1.174  ,1.172  ,1.161  ,1.155  ,1.151  ,1.150  ]
#ratio          = [1.02371,1.02283,1.01886,0.99424,0.98171,0.97528,0.97514] #CQ_Tau_SB_contp0...p6_EB1.vis.npz vs CQ_Tau_LB_contp0...p6_EB1.vis.npz
#scaling_factor = [1.012  ,1.011  ,1.009  ,0.997  ,0.991  ,0.988  ,0.987  ]
#ratio          = [0.95141,0.95212,0.96038,0.97280,0.98479,0.99312,0.99917] #CQ_Tau_LB_contp0...p6_EB0.vis.npz vs CQ_Tau_LB_contp0...p6_EB1.vis.npz
#scaling_factor = [0.975  ,0.976  ,0.980  ,0.986  ,0.992  ,0.997  ,1.000  ]

#iteration_2
#ratio          = [1.03281,1.03233,1.03674,1.04108,1.04203,1.04247] #CQ_Tau_SB_contp0...p5_EB0.vis.npz vs CQ_Tau_LB_EB1_initcont_shift.vis.npz
#scaling_factor = [1.016  ,1.016  ,1.018  ,1.020  ,1.021  ,1.021  ]
#ratio          = [0.91086,0.91202,0.94070,0.97043,0.99265,1.02360] #CQ_Tau_SB_contp0...p5_EB1.vis.npz vs CQ_Tau_LB_EB1_initcont_shift.vis.npz
#scaling_factor = [0.954  ,0.955  ,0.970  ,0.985  ,0.996  ,1.012  ]

#ratio          = [1.04245,1.04198,1.03864,1.01923,1.00943,1.00258,0.99957,0.99209,0.99104] #CQ_Tau_SB_contp0...ap1_EB0.vis.npz vs CQ_Tau_LB_contp0...p6_EB1.vis.npz
#scaling_factor = [1.021  ,1.021  ,1.019  ,1.010  ,1.005  ,1.001  ,1.000  ,0.996  ,0.996  ]
#ratio          = [1.02359,1.02280,1.01886,0.99413,0.98169,0.97538,0.97518,0.99748,0.99367] #CQ_Tau_SB_contp0...ap1_EB1.vis.npz vs CQ_Tau_LB_contp0...p6_EB1.vis.npz
#scaling_factor = [1.012  ,1.011  ,1.009  ,0.997  ,0.991  ,0.988  ,0.988  ,0.999  ,0.997  ]
#ratio          = [0.95141,0.95208,0.96012,0.97206,0.98394,0.99174,0.99767,1.00281,0.99899] #CQ_Tau_LB_contp0...ap1_EB0.vis.npz vs CQ_Tau_LB_contp0...p6_EB1.vis.npz
#scaling_factor = [0.975  ,0.976  ,0.980  ,0.986  ,0.992  ,0.996  ,0.999  ,1.001  ,0.999  ]

#All the EBs show flux offsets <0.9% after ap1 and <0.8 after ap0. ap0 and ap1 have the same SNR and noise structure. I chose ap1.

#END of COMBINED SB+LB phase-only self-cal iteration 2

#Split out final continuum ms table, with a 30s timebin
LB_iteration2_cont_averaged = f'{prefix}_time_ave_continuum'
os.system(f'rm -rf {LB_iteration2_cont_averaged}.ms*')
split(vis=LB_iteration2_cont_ap1+'.ms',outputvis=LB_iteration2_cont_averaged+'.ms',datacolumn='data',keepflags=False,timebin='30s')

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

#Apply statwt to the individual EBs
for params in data_params.values():
    shutil.copytree(src=prefix+'_'+params['name']+'.ms',dst=prefix+'_'+params['name']+'_statwt.ms')
    statwt(
        vis        = prefix+'_'+params['name']+'_statwt.ms',
        spw        = params['cont_spws'],
        field      = 'CQ_Tau',
        combine    = 'scan,spw,corr,field',
        timebin    = '10000000s',    #"inf" does not work...
        chanbin    = 'spw',          #"... channel binning occurs within individual spectral windows; bins never span multiple spectral windows..."
        datacolumn = 'DATA',
        intent     = 'OBSERVE_TARGET#ON_SOURCE',
    )

#Apply the gaintables of individual EBs
for params in data_params.values():
    single_EB_p1 = prefix+'_'+params['name']+'_initcont.p1'
    vis          = prefix+'_'+params['name']+'_statwt.ms'
    applycal(vis=vis,spw=single_EB_contspws,spwmap=single_EB_spw_mapping,gaintable=[single_EB_p1],interp='linearPD',applymode='calonly',calwt=True)
    split(vis=vis,outputvis=prefix+'_'+params['name']+'_no_ave_selfcal.ms',datacolumn='corrected')

#Align the data
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

#If you have re-scaled fluxes, you need to re-scale the shifted *no_ave* EBs as well
rescale_flux(vis=prefix+'_SB_EB0_no_ave_selfcal_shift.ms', gencalparameter=[1.150])
listobs(vis=prefix+'_SB_EB0_no_ave_selfcal_shift_rescaled.ms',listfile=prefix+'_SB_EB0_no_ave_selfcal_shift_rescaled.ms.listobs.txt',overwrite=True)

#Concat the non-averaged SB data
SB_combined = f'{prefix}_SB_no_ave_concat'
os.system('rm -rf %s.ms*' %SB_combined)
concat(
    vis=[f'{prefix}_SB_EB0_no_ave_selfcal_shift_rescaled.ms',f'{prefix}_SB_EB1_no_ave_selfcal_shift.ms'],
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
    gaintable  = [SB_iteration2_p1,SB_iteration2_p2,SB_iteration2_p3,SB_iteration2_p4,SB_iteration2_p5],
    spw        = SB_contspws,
    spwmap     = [SB_spw_mapping]*5,
    interp     = ['linearPD']*5, 
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
    vis=[SB_no_ave_selfcal]+[f'{prefix}_LB_EB0_no_ave_selfcal_shift.ms',f'{prefix}_LB_EB1_no_ave_selfcal_shift.ms'],
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
    gaintable  = [LB_iteration2_p1,LB_iteration2_p2,LB_iteration2_p3,LB_iteration2_p4,LB_iteration2_p5,LB_iteration2_p6,LB_iteration2_ap0,LB_iteration2_ap1],
    spw        = LB_contspws,
    spwmap     = [LB_spw_mapping]*8,
    interp     = ['linearPD']*8,
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
    'field':     'CQ_Tau',
    'line_spws': np.array([
        1, 4, 5,5, 6,6, 7,7,7,            #1, 4, 5,5, 6,6, 7,7,7,
        9,9,9,9, 11, 13, 14,14, 15,15,    #1,1,1,1, 3, 5, 6,6, 7,7,
        17,17,17,17, 19, 21, 22,22, 23,23,#1,1,1,1, 3, 5, 6,6, 7,7,
        25,25,25,25, 27, 29, 30,30, 31,31,#1,1,1,1, 3, 5, 6,6, 7,7,
    ]), # list of spws containing lines
    'line_freqs':  np.array([
        rest_freq_12CO, 
        rest_freq_H2CO_321_220, 
        rest_freq_SO,rest_freq_13CO, 
        rest_freq_C18O,rest_freq_SO,
        rest_freq_H2CO_918_919,rest_freq_DCN_F21,rest_freq_DCN_F22,

        rest_freq_DCN_F21,rest_freq_DCN_F22,rest_freq_H2CO_303_202,rest_freq_H2CO_322_221, 
        rest_freq_12CO, 
        rest_freq_H2CO_321_220, 
        rest_freq_SO,rest_freq_13CO, 
        rest_freq_C18O,rest_freq_SO,
        
        rest_freq_DCN_F21,rest_freq_DCN_F22,rest_freq_H2CO_303_202,rest_freq_H2CO_322_221,
        rest_freq_12CO,
        rest_freq_H2CO_321_220,
        rest_freq_SO,rest_freq_13CO,
        rest_freq_C18O,rest_freq_SO,
        
        rest_freq_DCN_F21,rest_freq_DCN_F22,rest_freq_H2CO_303_202,rest_freq_H2CO_322_221,
        rest_freq_12CO,
        rest_freq_H2CO_321_220,
        rest_freq_SO,rest_freq_13CO,
        rest_freq_C18O,rest_freq_SO,
    ]), #frequencies (Hz) corresponding to line_spws
    'spwcont_forplot': ['2','0','0','0'],
    'cont_spws':       '0~31',
    'width_array':     [
        60,960,16,30,30,480,480,16,
        8,8,30,480,12,12,240,240,
        4,4,15,240,6,6,120,120,
        4,4,15,240,6,6,120,120,
    ],
}

#Use the output of get_flagchannels at the beginning of the script to define fitspw
# Flagchannels input string for LB_EB0: '1:88~89, 1:88~89, 1:25~26, 1:9~10, 3:112~207, 5:33~35, 6:959~959, 6:53~98, 7:810~855, 7:13~58'
# Flagchannels input string for LB_EB1: '1:88~89, 1:88~89, 1:25~26, 1:9~10, 3:112~207, 5:33~35, 6:959~959, 6:53~98, 7:810~855, 7:13~58'
# Flagchannels input string for SB_EB0: '1:128~222, 4:32~35, 5:959~959, 5:45~91, 6:802~847, 6:5~51, 7:90~91, 7:47~48, 7:47~48'
# Flagchannels input string for SB_EB1: '1:87~89, 1:87~89, 1:24~26, 1:8~9, 3:159~254, 5:31~34, 6:950~959, 6:30~75, 7:787~832, 7:0~36'

#Remember how the spws are organized in the SBLB files (seen from listobs)
# [0,1,2,3,4,5,6,7]         -> SB EB0
# [8,9,10,11,12,13,14,15]   -> SB EB1
# [16,17,18,19,20,21,22,23] -> LB EB0
# [24,25,26,27,28,29,30,31] -> LB EB1
fitspw =  '0:60~74, 1:128~222, 2:0, 3:0, 4:32~35, 5:959~959, 5:45~91, 6:802~847, 6:5~51, 7:90~91, 7:47~48, 7:47~48, '\
         +'8:0, 9:87~89, 9:87~89, 9:24~26, 9:8~9, 10:0, 11:159~254, 12:0, 13:31~34, 14:950~959, 14:30~75, 15:787~832, 15:0~36, '\
         +'16:0, 17:88~89, 17:88~89, 17:25~26, 17:9~10, 18:0, 19:112~207, 20:0, 21:33~35, 22:959~959, 22:53~98, 23:810~855, 23:13~58, '\
         +'24:0, 25:88~89, 25:88~89, 25:25~26, 25:9~10, 26:0, 27:112~207, 28:0, 29:33~35, 30:959~959, 30:53~98, 31:810~855, 31:13~58'

avg_cont(
    ms_dict=complete_dataset_dict,output_prefix=prefix,flagchannels=fitspw,
    contspws=complete_dataset_dict['cont_spws'],width_array=complete_dataset_dict['width_array']
)

#Image the avg then cal with same parameters as in last step of self-cal
tclean_wrapper(
    vis       = LB_iteration2_cont_averaged+'.ms',
    imagename = LB_iteration2_cont_averaged+'_image', 
    threshold = '0.0118mJy',
    **LB_tclean_wrapper_kwargs
)
estimate_SNR(LB_iteration2_cont_averaged+'_image'+'.image',disk_mask=LB_mask,noise_mask=noise_annulus_LB)
generate_image_png(
    LB_iteration2_cont_averaged+'_image'+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_iteration2_LB,10*rms_iteration2_LB],save_folder=calibrate_linedata_folder
)
#CQ_Tau_SBLB_iteration2_contap1.image
#Beam 0.081 arcsec x 0.059 arcsec (11.53 deg)
#Flux inside disk mask: 147.09 mJy
#Peak intensity of source: 2.65 mJy/beam
#rms: 1.18e-02 mJy/beam
#Peak SNR: 224.17

#CQ_Tau_time_ave_continuum_image.image
#Beam 0.081 arcsec x 0.059 arcsec (11.53 deg)
#Flux inside disk mask: 147.11 mJy
#Peak intensity of source: 2.65 mJy/beam
#rms: 1.18e-02 mJy/beam
#Peak SNR: 224.58

#Image the cal then avg with same parameters as in last step of self-cal
complete_dataset_image = prefix+'_'+complete_dataset_dict['name']+'_initcont_image'
tclean_wrapper(
    vis       = prefix+'_'+complete_dataset_dict['name']+'_initcont.ms',
    imagename = complete_dataset_image, 
    threshold = '0.0118mJy',
    **LB_tclean_wrapper_kwargs
)
estimate_SNR(complete_dataset_image+'.image',disk_mask=LB_mask,noise_mask=noise_annulus_LB)
generate_image_png(
    complete_dataset_image+'.image',plot_sizes=image_png_plot_sizes,
    color_scale_limits=[-3*rms_iteration2_LB,10*rms_iteration2_LB],save_folder=calibrate_linedata_folder
)
#CQ_Tau_time_ave_continuum_image.image
#Beam 0.081 arcsec x 0.059 arcsec (11.53 deg)
#Flux inside disk mask: 147.11 mJy
#Peak intensity of source: 2.65 mJy/beam
#rms: 1.18e-02 mJy/beam
#Peak SNR: 224.58

#CQ_Tau_SBLB_concat_initcont_image.image
#Beam 0.081 arcsec x 0.060 arcsec (11.65 deg)
#Flux inside disk mask: 147.34 mJy
#Peak intensity of source: 2.70 mJy/beam
#rms: 1.22e-02 mJy/beam
#Peak SNR: 221.84

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
#Check OK: the ratio is equal to ~1.01-1.02
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
#Complaints:
#2024-12-27 16:35:36     SEVERE  uvcontsub::::casa       Task uvcontsub raised an exception of class ValueError with the following message: combine must include 'spw' when the fit is being applied to spws outside fitspw.
#2024-12-27 16:35:36     SEVERE  uvcontsub::::casa       Exception Reported: Error in uvcontsub: combine must include 'spw' when the fit is being applied to spws outside fitspw.
#Solution: added null channels for the other spws in fitspw

#2024-12-27 16:44:24     WARN    calibrater::setvi(bool,bool)    Forcing use of OLD VisibilityIterator.
#2024-12-27 17:14:05     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [231.380, 233.364] (GHz).
#     The frequency range used for the continuum fit was [231.396, 233.364] (GHz).
#2024-12-27 17:14:06     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [219.038, 219.498] (GHz).
#     The frequency range used for the continuum fit was [219.038, 219.491] (GHz).
#2024-12-27 17:14:08     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [219.971, 220.440] (GHz).
#     The frequency range used for the continuum fit was [219.972, 220.44] (GHz).
#2024-12-27 17:14:24     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [231.380, 233.365] (GHz).
#     The frequency range used for the continuum fit was [231.396, 233.365] (GHz).
#2024-12-27 17:14:25     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [230.960, 231.425] (GHz).
#     The frequency range used for the continuum fit was [230.964, 231.425] (GHz).
#2024-12-27 17:14:25     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [219.026, 219.487] (GHz).
#     The frequency range used for the continuum fit was [219.026, 219.479] (GHz).
#2024-12-27 17:14:26     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [219.960, 220.428] (GHz).
#     The frequency range used for the continuum fit was [219.965, 220.428] (GHz).
#2024-12-27 17:14:26     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [219.491, 219.960] (GHz).
#     The frequency range used for the continuum fit was [219.491, 219.942] (GHz).
#2024-12-27 17:15:30     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [231.384, 233.368] (GHz).
#     The frequency range used for the continuum fit was [231.4, 233.368] (GHz).
#2024-12-27 17:15:31     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [230.964, 231.429] (GHz).
#     The frequency range used for the continuum fit was [230.968, 231.429] (GHz).
#2024-12-27 17:15:32     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [219.030, 219.491] (GHz).
#     The frequency range used for the continuum fit was [219.03, 219.483] (GHz).
#2024-12-27 17:15:32     WARN    VBContinuumSubtractor::apply    Extrapolating to cover [219.964, 220.432] (GHz).
#     The frequency range used for the continuum fit was [219.964, 220.432] (GHz).
"""
#Split final ms table into separate spws for each targeted line (added SiS, double velocity range coverage, just in case, but check that no overlap occurs)
rest_freq_12CO         = 230.5380000e9 #J=2-1
    #1:81~270,11:112~301,19:65~254,27:65~254
rest_freq_13CO         = 220.3986842e9 #J=2-1 (ignoring splitting)
    #5:23~113,14:8~98,22:30~121,30:30~121
rest_freq_C18O         = 219.5603541e9 #J=2-1
    #6:780~870,15:765~855,23:787~877,31:787~877
rest_freq_H2CO_303_202 = 218.2221920e9
    #9:23~26,17:24~27,25:24~27
rest_freq_H2CO_321_220 = 218.7600660e9
    #4:31~36,13:30~35,21:31~37,29:31~37
rest_freq_H2CO_322_221 = 218.4756320e9
    #9:7~10,17:8~11,25:8~11
rest_freq_H2CO_918_919 = 216.5686510e9
    #7:89~92,
rest_freq_DCN_F21      = 217.2384000e9 #J=3-2
rest_freq_DCN_F22      = 217.2386307e9 #J=3-2
    #7:46~49,9:86~89,17:87~90,25:87~90
rest_freq_SiS_1211     = 217.8176630e9 #J=12-11
    #7:9~12,9:49~52,17:50~53,25:50~53
rest_freq_SO           = 219.9494420e9 #3Sigma 6(5)-5(4)
    #5:943~959,6:0~73,14:928~959,15:0~58,22:950~959,23:0~81,30:950~959,31:0~81

line_spws_SBs  = [[1, 4, 5,5, 6,6, 7,7,7,7],[1,1,1,1,1, 3, 5, 6,6, 7,7]]

line_freqs_SBs = [
    [rest_freq_12CO, rest_freq_H2CO_321_220, rest_freq_SO,rest_freq_13CO, rest_freq_C18O,rest_freq_SO, rest_freq_H2CO_918_919,rest_freq_DCN_F21,rest_freq_DCN_F22,rest_freq_SiS_1211],
    [rest_freq_DCN_F21,rest_freq_DCN_F22,rest_freq_SiS_1211,rest_freq_H2CO_303_202,rest_freq_H2CO_322_221, rest_freq_12CO, rest_freq_H2CO_321_220, rest_freq_SO,rest_freq_13CO, rest_freq_C18O,rest_freq_SO]
]

spwcont_forplot_SBs = ['2','0']

width_array_SBs = [[60,960,16,30,30,480,480,16],[8,8,30,480,12,12,240,240]]

data_params_LB = {
    f'LB{i}': {
        'vis':          f'{prefix}_LB_EB{i}.ms',
        'name':         f'LB_EB{i}',
        'field':        'CQ_Tau',
        'line_spws':    np.array([1,1,1,1,1, 3, 5, 6,6, 7,7]), #list of spws containing lines
        'line_freqs':   np.array([
            rest_freq_DCN_F21,rest_freq_DCN_F22,rest_freq_SiS_1211,rest_freq_H2CO_303_202,rest_freq_H2CO_322_221,
            rest_freq_12CO,
            rest_freq_H2CO_321_220,
            rest_freq_SO,rest_freq_13CO,
            rest_freq_C18O,rest_freq_SO
        ]), #frequencies (Hz) corresponding to line_spws
        'spwcont_forplot': '0',
        'cont_spws':       '0,1,2,3,4,5,6,7',
        'width_array':     [4,4,15,240,6,6,120,120],
    } for i in range(number_of_EBs['LB'])#Phase RMS: 28.985, 22.576 deg
}

data_params_SB = {
    f'SB{i}': {
        'vis':          f'{prefix}_SB_EB{i}.ms',
        'name':         f'SB_EB{i}',
        'field':        'CQ_Tau',
        'line_spws':    line_spws_SBs[i], #list of spws containing lines
        'line_freqs':   line_freqs_SBs[i], #frequencies (Hz) corresponding to line_spws
        'spwcont_forplot': spwcont_forplot_SBs[i],
        'cont_spws':       '0,1,2,3,4,5,6,7',
        'width_array':     width_array_SBs[i],
    } for i in range(number_of_EBs['SB'])
}

data_params = data_params_LB.copy()
data_params.update(data_params_SB)

for params in data_params.values():
    flagchannels_string = get_flagchannels(ms_dict=params,output_prefix=prefix,velocity_range=np.array([-30.,30.]) + v_sys)

# Flagchannels input string for LB_EB0: '1:87~90, 1:87~90, 1:50~53, 1:24~27, 1:8~11, 3:65~254, 5:31~37, 6:950~959, 6:30~121, 7:787~877, 7:0~81'
# Flagchannels input string for LB_EB1: '1:87~90, 1:87~90, 1:50~53, 1:24~27, 1:8~11, 3:65~254, 5:31~37, 6:950~959, 6:30~121, 7:787~877, 7:0~81'
# Flagchannels input string for SB_EB0: '1:81~270, 4:31~36, 5:943~959, 5:23~113, 6:780~870, 6:0~73, 7:89~92, 7:46~49, 7:46~49, 7:9~12'
# Flagchannels input string for SB_EB1: '1:86~89, 1:86~89, 1:49~52, 1:23~26, 1:7~10, 3:112~301, 5:30~35, 6:928~959, 6:8~98, 7:765~855, 7:0~58'

# 1 : 81~270, rest_freq_12CO
# 4 : 31~ 36, rest_freq_H2CO_321_220
# 5 :943~959, rest_freq_SO
# 5 : 23~113, rest_freq_13CO
# 6 :780~870, rest_freq_C18O
# 6 :  0~ 73, rest_freq_SO
# 7 : 89~ 92, rest_freq_H2CO_918_919
# 7 : 46~ 49, rest_freq_DCN_F21
# 7 : 46~ 49, rest_freq_DCN_F22
# 7 :  9~ 12, rest_freq_SiS_1211
# 9 : 86~ 89, rest_freq_DCN_F21
# 9 : 86~ 89, rest_freq_DCN_F22
# 9 : 49~ 52, rest_freq_SiS_1211
# 9 : 23~ 26, rest_freq_H2CO_303_202
# 9 :  7~ 10, rest_freq_H2CO_322_221
# 11:112~301, rest_freq_12CO
# 13: 30~ 35, rest_freq_H2CO_321_220
# 14:928~959, rest_freq_SO
# 14:  8~ 98, rest_freq_13CO
# 15:765~855, rest_freq_C18O
# 15:  0~ 58, rest_freq_SO
# 17: 87~ 90, rest_freq_DCN_F21
# 17: 87~ 90, rest_freq_DCN_F22
# 17: 50~ 53, rest_freq_SiS_1211
# 17: 24~ 27, rest_freq_H2CO_303_202
# 17:  8~ 11, rest_freq_H2CO_322_221
# 19: 65~254, rest_freq_12CO
# 21: 31~ 37, rest_freq_H2CO_321_220
# 22:950~959, rest_freq_SO
# 22: 30~121, rest_freq_13CO
# 23:787~877, rest_freq_C18O
# 23:  0~ 81, rest_freq_SO
# 25: 87~ 90, rest_freq_DCN_F21
# 25: 87~ 90, rest_freq_DCN_F22
# 25: 50~ 53, rest_freq_SiS_1211
# 25: 24~ 27, rest_freq_H2CO_303_202
# 25:  8~ 11, rest_freq_H2CO_322_221
# 27: 65~254, rest_freq_12CO
# 29: 31~ 37, rest_freq_H2CO_321_220
# 30:950~959, rest_freq_SO
# 30: 30~121, rest_freq_13CO
# 31:787~877, rest_freq_C18O
# 31:  0~ 81, rest_freq_SO

SBLB_no_ave_selfcal = f'{prefix}_SBLB_no_ave_selfcal_time_ave.ms'
contsub_vis = f'{SBLB_no_ave_selfcal}.contsub'

#12CO
vis_12CO = SBLB_no_ave_selfcal[:-3]+'_12CO.ms'
os.system(f'rm -rf {vis_12CO}*')
spw_12CO = '1:81~270,11:112~301,19:65~254,27:65~254'

split(vis=SBLB_no_ave_selfcal,outputvis=vis_12CO,spw=spw_12CO,datacolumn='data',keepflags=False)
listobs(vis=vis_12CO,listfile=vis_12CO+'.listobs.txt',overwrite=True)

split(vis=contsub_vis,outputvis=f'{vis_12CO}.contsub',spw=spw_12CO,datacolumn='data',keepflags=False)
listobs(vis=vis_12CO+'.contsub',listfile=vis_12CO+'.contsub.listobs.txt',overwrite=True)

#13CO
vis_13CO = SBLB_no_ave_selfcal[:-3]+'_13CO.ms'
os.system(f'rm -rf {vis_13CO}*')
spw_13CO = '5:23~113,14:8~98,22:30~121,30:30~121'

split(vis=SBLB_no_ave_selfcal,outputvis=vis_13CO,spw=spw_13CO,datacolumn='data',keepflags=False)
listobs(vis=vis_13CO,listfile=vis_13CO+'.listobs.txt',overwrite=True)

split(vis=contsub_vis,outputvis=f'{vis_13CO}.contsub',spw=spw_13CO,datacolumn='data',keepflags=False)
listobs(vis=vis_13CO+'.contsub',listfile=vis_13CO+'.contsub.listobs.txt',overwrite=True)

#C18O
vis_C18O = SBLB_no_ave_selfcal[:-3]+'_C18O.ms'
os.system(f'rm -rf {vis_C18O}*')
spw_C18O = '6:780~870,15:765~855,23:787~877,31:787~877'

split(vis=SBLB_no_ave_selfcal,outputvis=vis_C18O,spw=spw_C18O,datacolumn='data',keepflags=False)
listobs(vis=vis_C18O,listfile=vis_C18O+'.listobs.txt',overwrite=True)

split(vis=contsub_vis,outputvis=f'{vis_C18O}.contsub',spw=spw_C18O,datacolumn='data',keepflags=False)
listobs(vis=vis_C18O+'.contsub',listfile=vis_C18O+'.contsub.listobs.txt',overwrite=True)

#H2CO_303_202
vis_H2CO_303_202 = SBLB_no_ave_selfcal[:-3]+'_H2CO_303_202.ms'
os.system(f'rm -rf {vis_H2CO_303_202}*')
spw_H2CO_303_202 = '9:23~26,17:24~27,25:24~27'

split(vis=SBLB_no_ave_selfcal,outputvis=vis_H2CO_303_202,spw=spw_H2CO_303_202,datacolumn='data',keepflags=False)
listobs(vis=vis_H2CO_303_202,listfile=vis_H2CO_303_202+'.listobs.txt',overwrite=True)

split(vis=contsub_vis,outputvis=f'{vis_H2CO_303_202}.contsub',spw=spw_H2CO_303_202,datacolumn='data',keepflags=False)
listobs(vis=vis_H2CO_303_202+'.contsub',listfile=vis_H2CO_303_202+'.contsub.listobs.txt',overwrite=True)

#H2CO_321_220
vis_H2CO_321_220 = SBLB_no_ave_selfcal[:-3]+'_H2CO_321_220.ms'
os.system(f'rm -rf {vis_H2CO_321_220}*')
spw_H2CO_321_220 = '4:31~36,13:30~35,21:31~37,29:31~37'

split(vis=SBLB_no_ave_selfcal,outputvis=vis_H2CO_321_220,spw=spw_H2CO_321_220,datacolumn='data',keepflags=False)
listobs(vis=vis_H2CO_321_220,listfile=vis_H2CO_321_220+'.listobs.txt',overwrite=True)

split(vis=contsub_vis,outputvis=f'{vis_H2CO_321_220}.contsub',spw=spw_H2CO_321_220,datacolumn='data',keepflags=False)
listobs(vis=vis_H2CO_321_220+'.contsub',listfile=vis_H2CO_321_220+'.contsub.listobs.txt',overwrite=True)

#H2CO_322_221
vis_H2CO_322_221 = SBLB_no_ave_selfcal[:-3]+'_H2CO_322_221.ms'
os.system(f'rm -rf {vis_H2CO_322_221}*')
spw_H2CO_322_221 = '9:7~10,17:8~11,25:8~11'

split(vis=SBLB_no_ave_selfcal,outputvis=vis_H2CO_322_221,spw=spw_H2CO_322_221,datacolumn='data',keepflags=False)
listobs(vis=vis_H2CO_322_221,listfile=vis_H2CO_322_221+'.listobs.txt',overwrite=True)

split(vis=contsub_vis,outputvis=f'{vis_H2CO_322_221}.contsub',spw=spw_H2CO_322_221,datacolumn='data',keepflags=False)
listobs(vis=vis_H2CO_322_221+'.contsub',listfile=vis_H2CO_322_221+'.contsub.listobs.txt',overwrite=True)

#H2CO_918_919
vis_H2CO_918_919 = SBLB_no_ave_selfcal[:-3]+'_H2CO_918_919.ms'
os.system(f'rm -rf {vis_H2CO_918_919}*')
spw_H2CO_918_919 = '7:89~92'

split(vis=SBLB_no_ave_selfcal,outputvis=vis_H2CO_918_919,spw=spw_H2CO_918_919,datacolumn='data',keepflags=False)
listobs(vis=vis_H2CO_918_919,listfile=vis_H2CO_918_919+'.listobs.txt',overwrite=True)

split(vis=contsub_vis,outputvis=f'{vis_H2CO_918_919}.contsub',spw=spw_H2CO_918_919,datacolumn='data',keepflags=False)
listobs(vis=vis_H2CO_918_919+'.contsub',listfile=vis_H2CO_918_919+'.contsub.listobs.txt',overwrite=True)

#DCN
vis_DCN = SBLB_no_ave_selfcal[:-3]+'_DCN.ms'
os.system(f'rm -rf {vis_DCN}*')
spw_DCN = '7:46~49,9:86~89,17:87~90,25:87~90'

split(vis=SBLB_no_ave_selfcal,outputvis=vis_DCN,spw=spw_DCN,datacolumn='data',keepflags=False)
listobs(vis=vis_DCN,listfile=vis_DCN+'.listobs.txt',overwrite=True)

split(vis=contsub_vis,outputvis=f'{vis_DCN}.contsub',spw=spw_DCN,datacolumn='data',keepflags=False)
listobs(vis=vis_DCN+'.contsub',listfile=vis_DCN+'.contsub.listobs.txt',overwrite=True)

#SiS
vis_SiS = SBLB_no_ave_selfcal[:-3]+'_SiS.ms'
os.system(f'rm -rf {vis_SiS}*')
spw_SiS = '7:9~12,9:49~52,17:50~53,25:50~53'

split(vis=SBLB_no_ave_selfcal,outputvis=vis_SiS,spw=spw_SiS,datacolumn='data',keepflags=False)
listobs(vis=vis_SiS,listfile=vis_SiS+'.listobs.txt',overwrite=True)

split(vis=contsub_vis,outputvis=f'{vis_SiS}.contsub',spw=spw_SiS,datacolumn='data',keepflags=False)
listobs(vis=vis_SiS+'.contsub',listfile=vis_SiS+'.contsub.listobs.txt',overwrite=True)

#SO
vis_SO = SBLB_no_ave_selfcal[:-3]+'_SO.ms'
os.system(f'rm -rf {vis_SO}*')
spw_SO = '5:943~959,6:0~73,14:928~959,15:0~58,22:950~959,23:0~81,30:950~959,31:0~81'

split(vis=SBLB_no_ave_selfcal,outputvis=vis_SO,spw=spw_SO,datacolumn='data',keepflags=False)
listobs(vis=vis_SO,listfile=vis_SO+'.listobs.txt',overwrite=True)

split(vis=contsub_vis,outputvis=f'{vis_SO}.contsub',spw=spw_SO,datacolumn='data',keepflags=False)
listobs(vis=vis_SO+'.contsub',listfile=vis_SO+'.contsub.listobs.txt',overwrite=True)
"""
for _vis,_freq in zip(
    [vis_12CO,vis_13CO,vis_C18O,vis_H2CO_303_202,vis_H2CO_321_220,vis_H2CO_322_221,vis_H2CO_918_919,vis_DCN,vis_SiS,vis_SO],
    [rest_freq_12CO,rest_freq_13CO,rest_freq_C18O,rest_freq_H2CO_303_202,rest_freq_H2CO_321_220,rest_freq_H2CO_322_221,rest_freq_H2CO_918_919,rest_freq_DCN_F21,rest_freq_SiS_1211,rest_freq_SO]
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
mask_semimajor = 1.7
mask_semiminor = 1.5
mask_ra  = '05h35m58.472722s'
mask_dec = '24d44m53.613450s'

line_mask = f'ellipse[[{mask_ra},{mask_dec}], [{mask_semimajor:.3f}arcsec, {mask_semiminor:.3f}arcsec], {mask_pa:.1f}deg]'

noise_annulus_line = f"annulus[[{mask_ra}, {mask_dec}],['4.arcsec', '6.arcsec']]"

chanstart = '-11.30km/s'
chanwidth = '  0.35km/s'
nchan     = 100

imagename = vis_12CO[:-3]+'.contsub_image_1.0robust_4.0sigma'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis=vis_12CO+'.contsub',imagename=imagename,
    specmode='cube',restoringbeam='common',nterms=1,
    deconvolver='multiscale',scales=[0,2,4,8],
    cell='0.011arcsec',imsize=1200,gain=0.1,niter=50000,
    weighting='briggs',robust=1.0,threshold='4.0800mJy',interactive=False, #~4sigma
    start=chanstart,width=chanwidth,nchan=nchan,
    outframe='LSRK',veltype='radio',restfreq='{}Hz'.format(rest_freq_12CO),
    usemask='user',mask=line_mask, 
)
estimate_SNR(imagename+'.image',disk_mask=line_mask,noise_mask=noise_annulus_line)
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#CQ_Tau_SBLB_no_ave_selfcal_time_ave_12CO.contsub_image_1.0robust_4.0sigma.image
#Beam 0.105 arcsec x 0.087 arcsec (0.06 deg)
#Flux inside disk mask: 2601.95 mJy
#Peak intensity of source: 33.59 mJy/beam
#rms: 1.02e+00 mJy/beam
#Peak SNR: 32.93

chanstart = '-11.3km/s'
chanwidth = '  0.7km/s'
nchan     = 50

imagename = vis_13CO[:-3]+'.contsub_image_1.0robust_4.0sigma'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis=vis_13CO+'.contsub',imagename=imagename,
    specmode='cube',restoringbeam='common',nterms=1,
    deconvolver='multiscale',scales=[0,2,4,8],
    cell='0.012arcsec',imsize=1200,gain=0.1,niter=50000,
    weighting='briggs',robust=1.0,threshold='3.0440mJy',interactive=False, #~4sigma
    start=chanstart,width=chanwidth,nchan=nchan,
    outframe='LSRK',veltype='radio',restfreq='{}Hz'.format(rest_freq_13CO),
    usemask='user',mask=line_mask,
    #uvtaper='0.1arcsec',#interactive=True, 
)
estimate_SNR(imagename+'.image',disk_mask=line_mask,noise_mask=noise_annulus_line)
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#CQ_Tau_SBLB_no_ave_selfcal_time_ave_13CO.contsub_image_1.0robust_4.0sigma.image
#Beam 0.109 arcsec x 0.092 arcsec (0.08 deg)
#Flux inside disk mask: 678.83 mJy
#Peak intensity of source: 15.21 mJy/beam
#rms: 7.61e-01 mJy/beam
#Peak SNR: 20.00

chanstart = '-11.3km/s'
chanwidth = '  0.7km/s'
nchan     = 50

imagename = vis_C18O[:-3]+'.contsub_image_1.0robust_4.0sigma'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis=vis_C18O+'.contsub',imagename=imagename,
    specmode='cube',restoringbeam='common',nterms=1,
    deconvolver='multiscale',scales=[0,2,4,8],
    cell='0.012arcsec',imsize=1200,gain=0.1,niter=50000,
    weighting='briggs',robust=1.0,threshold='2.3320mJy',interactive=False, #~4sigma
    start=chanstart,width=chanwidth,nchan=nchan,
    outframe='LSRK',veltype='radio',restfreq='{}Hz'.format(rest_freq_C18O),
    usemask='user',mask=line_mask,
    #uvtaper='0.1arcsec',#interactive=True, 
)
estimate_SNR(imagename+'.image',disk_mask=line_mask,noise_mask=noise_annulus_line)
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#CQ_Tau_SBLB_no_ave_selfcal_time_ave_C18O.contsub_image_1.0robust_4.0sigma.image
#Beam 0.110 arcsec x 0.093 arcsec (0.13 deg)
#Flux inside disk mask: 364.37 mJy
#Peak intensity of source: 10.06 mJy/beam
#rms: 5.83e-01 mJy/beam
#Peak SNR: 17.24

chanstart = '-11.3km/s'
chanwidth = '  0.7km/s'
nchan     = 50

imagename = vis_SO[:-3]+'.contsub_image_natural_4.0sigma'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis=vis_SO+'.contsub',imagename=imagename,
    specmode='cube',restoringbeam='common',nterms=1,
    deconvolver='multiscale',scales=[0,2,4],
    cell='0.013arcsec',imsize=1200,gain=0.1,niter=50000,
    weighting='natural',threshold='2.8080mJy',interactive=False, #~4sigma 
    start=chanstart,width=chanwidth,nchan=nchan,
    outframe='LSRK',veltype='radio',restfreq='{}Hz'.format(rest_freq_SO),
    usemask='user',mask=line_mask,
    #uvtaper='0.1arcsec',#interactive=True, 
)
estimate_SNR(imagename+'.image',disk_mask=line_mask,noise_mask=noise_annulus_line)
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#CQ_Tau_SBLB_no_ave_selfcal_time_ave_SO.contsub_image_natural_4.0sigma.image
#Beam 0.126 arcsec x 0.105 arcsec (-1.10 deg)
#Flux inside disk mask: -100.71 mJy
#Peak intensity of source: 6.81 mJy/beam
#rms: 7.02e-01 mJy/beam
#Peak SNR: 9.71

imagename = vis_SO[:-3]+'.contsub_image_0.5robust_4.0sigma'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis=vis_SO+'.contsub',imagename=imagename,
    specmode='cube',restoringbeam='common',nterms=1,
    deconvolver='multiscale',scales=[0,2,4],
    cell='0.008arcsec',imsize=1200,gain=0.1,niter=50000,
    weighting='briggs',robust=0.5,threshold='3.2120mJy',interactive=False, #~4sigma 
    start=chanstart,width=chanwidth,nchan=nchan,
    outframe='LSRK',veltype='radio',restfreq='{}Hz'.format(rest_freq_SO),
    usemask='user',mask=line_mask,
    #uvtaper='0.1arcsec',#interactive=True, 
)
estimate_SNR(imagename+'.image',disk_mask=line_mask,noise_mask=noise_annulus_line)
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#CQ_Tau_SBLB_no_ave_selfcal_time_ave_SO.contsub_image_0.5robust_4.0sigma.image
#Beam 0.083 arcsec x 0.063 arcsec (9.97 deg)
#Flux inside disk mask: -32.62 mJy
#Peak intensity of source: 4.62 mJy/beam
#rms: 8.03e-01 mJy/beam
#Peak SNR: 5.75

imagename = vis_SO[:-3]+'.contsub_image_0.75robust_4.0sigma'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis=vis_SO+'.contsub',imagename=imagename,
    specmode='cube',restoringbeam='common',nterms=1,
    deconvolver='multiscale',scales=[0,2,4],
    cell='0.010arcsec',imsize=1200,gain=0.1,niter=50000,
    weighting='briggs',robust=0.75,threshold='2.9600mJy',interactive=False, #~4sigma 
    start=chanstart,width=chanwidth,nchan=nchan,
    outframe='LSRK',veltype='radio',restfreq='{}Hz'.format(rest_freq_SO),
    usemask='user',mask=line_mask,
    #uvtaper='0.1arcsec',#interactive=True, 
)
estimate_SNR(imagename+'.image',disk_mask=line_mask,noise_mask=noise_annulus_line)
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#CQ_Tau_SBLB_no_ave_selfcal_time_ave_SO.contsub_image_0.75robust_4.0sigma.image
#Beam 0.098 arcsec x 0.078 arcsec (0.21 deg)
#Flux inside disk mask: -114.59 mJy
#Peak intensity of source: 5.27 mJy/beam
#rms: 7.40e-01 mJy/beam
#Peak SNR: 7.12

imagename = vis_SO[:-3]+'.contsub_image_1.0robust_4.0sigma'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis=vis_SO+'.contsub',imagename=imagename,
    specmode='cube',restoringbeam='common',nterms=1,
    deconvolver='multiscale',scales=[0,2,4],
    cell='0.011arcsec',imsize=1200,gain=0.1,niter=50000,
    weighting='briggs',robust=1.0,threshold='2.8480mJy',interactive=False, #~4sigma 
    start=chanstart,width=chanwidth,nchan=nchan,
    outframe='LSRK',veltype='radio',restfreq='{}Hz'.format(rest_freq_SO),
    usemask='user',mask=line_mask,
    #uvtaper='0.1arcsec',#interactive=True, 
)
estimate_SNR(imagename+'.image',disk_mask=line_mask,noise_mask=noise_annulus_line)
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#CQ_Tau_SBLB_no_ave_selfcal_time_ave_SO.contsub_image_1.0robust_4.0sigma.image
#Beam 0.110 arcsec x 0.093 arcsec (0.10 deg)
#Flux inside disk mask: -82.41 mJy
#Peak intensity of source: 6.11 mJy/beam
#rms: 7.12e-01 mJy/beam
#Peak SNR: 8.58

chanstart = '-30.00km/s'
chanwidth = ' 22.00km/s'
nchan     = 4

imagename = vis_SiS[:-3]+'.contsub_image_natural_4.0sigma'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis=vis_SiS+'.contsub',imagename=imagename,
    specmode='cube',restoringbeam='common',nterms=1,
    deconvolver='multiscale',scales=[0,2,4],
    cell='0.016arcsec',imsize=1200,gain=0.1,niter=50000,
    weighting='natural',threshold='0.5800mJy',interactive=False, #~4sigma 
    start=chanstart,width=chanwidth,nchan=nchan,
    outframe='LSRK',veltype='radio',restfreq='{}Hz'.format(rest_freq_SiS_1211),
    usemask='user',mask=line_mask,
    #uvtaper='0.1arcsec',#interactive=True, 
)
estimate_SNR(imagename+'.image',disk_mask=line_mask,noise_mask=noise_annulus_line)
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#CQ_Tau_SBLB_no_ave_selfcal_time_ave_SiS.contsub_image_natural_4.0sigma.image
#Beam 0.162 arcsec x 0.127 arcsec (-1.23 deg)
#Flux inside disk mask: -77.79 mJy
#Peak intensity of source: 1.04 mJy/beam
#rms: 1.45e-01 mJy/beam
#Peak SNR: 7.17

chanstart = '-26.00km/s'
chanwidth = ' 11.00km/s'
nchan     = 6

imagename = vis_H2CO_321_220[:-3]+'.contsub_image_natural_4.0sigma'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis=vis_H2CO_321_220+'.contsub',imagename=imagename,
    specmode='cube',restoringbeam='common',nterms=1,
    deconvolver='multiscale',scales=[0,2,4,8],
    cell='0.014arcsec',imsize=1200,gain=0.1,niter=50000,
    weighting='natural',threshold='1.624mJy',interactive=False, #~4sigma 
    start=chanstart,width=chanwidth,nchan=nchan,
    outframe='LSRK',veltype='radio',restfreq='{}Hz'.format(rest_freq_H2CO_321_220),
    usemask='user',mask=line_mask,
    #uvtaper='0.1arcsec',#interactive=True, 
)
estimate_SNR(imagename+'.image',disk_mask=line_mask,noise_mask=noise_annulus_line)
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#CQ_Tau_SBLB_no_ave_selfcal_time_ave_H2CO_321_220.contsub_image_natural_4.0sigma.image
#Beam 0.128 arcsec x 0.110 arcsec (-0.17 deg)
#Flux inside disk mask: 282.63 mJy
#Peak intensity of source: 4.46 mJy/beam
#rms: 4.09e-01 mJy/beam
#Peak SNR: 10.91

chanstart = '-30.00km/s'
chanwidth = ' 22.00km/s'
nchan     = 4

imagename = vis_H2CO_303_202[:-3]+'.contsub_image_natural_4.0sigma'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis=vis_H2CO_303_202+'.contsub',imagename=imagename,
    specmode='cube',restoringbeam='common',nterms=1,
    deconvolver='multiscale',scales=[0,2,4,8],
    cell='0.016arcsec',imsize=1200,gain=0.1,niter=50000,
    weighting='natural',threshold='0.5880mJy',interactive=False, #~4sigma 
    start=chanstart,width=chanwidth,nchan=nchan,
    outframe='LSRK',veltype='radio',restfreq='{}Hz'.format(rest_freq_H2CO_303_202),
    usemask='user',mask=line_mask,
    #uvtaper='0.1arcsec',#interactive=True, 
)
estimate_SNR(imagename+'.image',disk_mask=line_mask,noise_mask=noise_annulus_line)
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#CQ_Tau_SBLB_no_ave_selfcal_time_ave_H2CO_303_202.contsub_image_natural_4.0sigma.image
#Beam 0.162 arcsec x 0.127 arcsec (-1.23 deg)
#Flux inside disk mask: 33.35 mJy
#Peak intensity of source: 0.77 mJy/beam
#rms: 1.47e-01 mJy/beam
#Peak SNR: 5.27

chanstart = '-1.45km/s' #off-centre (because it's at the spw edge and only ~half of the line is covered, and channel number doesn't go to 0!? maybe due to pipeline flagging...)
chanwidth = '22.00km/s'
nchan     = 2

imagename = vis_H2CO_322_221[:-3]+'.contsub_image_natural_4.0sigma'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis=vis_H2CO_322_221+'.contsub',imagename=imagename,
    specmode='cube',restoringbeam='common',nterms=1,
    deconvolver='multiscale',scales=[0,2,4,8],
    cell='0.012arcsec',imsize=1200,gain=0.1,niter=50000,
    weighting='natural',threshold='0.6640mJy',interactive=False, #~4sigma 
    start=chanstart,width=chanwidth,nchan=nchan,
    outframe='LSRK',veltype='radio',restfreq='{}Hz'.format(rest_freq_H2CO_322_221),
    usemask='user',mask=line_mask,
    #uvtaper='0.1arcsec',#interactive=True, 
)
estimate_SNR(imagename+'.image',disk_mask=line_mask,noise_mask=noise_annulus_line)
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#CQ_Tau_SBLB_no_ave_selfcal_time_ave_H2CO_322_221.contsub_image_natural_4.0sigma.image
#Beam 0.117 arcsec x 0.098 arcsec (0.08 deg)
#Flux inside disk mask: 54.28 mJy
#Peak intensity of source: 0.81 mJy/beam
#rms: 1.66e-01 mJy/beam
#Peak SNR: 4.86

chanstart = '-15.80km/s'
chanwidth = ' 22.00km/s'
nchan     = 3

imagename = vis_H2CO_918_919[:-3]+'.contsub_image_natural_4.0sigma'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis=vis_H2CO_918_919+'.contsub',imagename=imagename,
    specmode='cube',restoringbeam='common',nterms=1,
    deconvolver='multiscale',scales=[0,2,4],
    cell='0.036arcsec',imsize=1200,gain=0.1,niter=50000,
    weighting='natural',threshold='1.3320mJy',interactive=False, #~4sigma 
    start=chanstart,width=chanwidth,nchan=nchan,
    outframe='LSRK',veltype='radio',restfreq='{}Hz'.format(rest_freq_H2CO_918_919),
    usemask='user',mask=line_mask,
    #uvtaper='0.1arcsec',#interactive=True, 
)
estimate_SNR(imagename+'.image',disk_mask=line_mask,noise_mask=noise_annulus_line)
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#2025-01-06 17:23:21     WARN    SIImageStore::restore (file /source/casa6/casatools/src/code/synthesis/ImagerObjects/SIImageStore.cc, line 2284)      
#    Restoring with an empty model image. Only residuals will be processed to form the output restored image.
#CQ_Tau_SBLB_no_ave_selfcal_time_ave_H2CO_918_919.contsub_image_natural_4.0sigma.image
#Beam 0.333 arcsec x 0.287 arcsec (-15.32 deg)
#Flux inside disk mask: -210.62 mJy
#Peak intensity of source: 0.90 mJy/beam
#rms: 3.34e-01 mJy/beam
#Peak SNR: 2.69

chanstart = '-15.80km/s'
chanwidth = ' 22.00km/s'
nchan     = 3

imagename = vis_DCN[:-3]+'.contsub_image_natural_4.0sigma'
for ext in ['.image','.mask','.model','.pb','.psf','.residual','.sumwt','.fits']:
    os.system('rm -rf '+ imagename + ext)
tclean(
    vis=vis_DCN+'.contsub',imagename=imagename,
    specmode='cube',restoringbeam='common',nterms=1,
    deconvolver='multiscale',scales=[0,2,4],
    cell='0.013arcsec',imsize=1200,gain=0.1,niter=50000,
    weighting='natural',threshold='0.3804mJy',interactive=False, #~4sigma 
    start=chanstart,width=chanwidth,nchan=nchan,
    outframe='LSRK',veltype='radio',restfreq='{}Hz'.format(rest_freq_DCN_F21),
    usemask='user',mask=line_mask,
    #uvtaper='0.1arcsec',#interactive=True, 
)
estimate_SNR(imagename+'.image',disk_mask=line_mask,noise_mask=noise_annulus_line)
exportfits(imagename=imagename+'.image',fitsimage=imagename+'.image.fits',overwrite=True)
#CQ_Tau_SBLB_no_ave_selfcal_time_ave_DCN.contsub_image_natural_4.0sigma.image
#Beam 0.127 arcsec x 0.106 arcsec (0.70 deg)
#Flux inside disk mask: 135.44 mJy
#Peak intensity of source: 0.38 mJy/beam
#rms: 9.51e-02 mJy/beam
#Peak SNR: 3.96