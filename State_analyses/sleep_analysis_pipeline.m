function sleep_analysis_pipeline(sData)

% Written by Christoffer Berge || Vervaeke lab

% Pipeline for sleep-swr project. Input has to be a sData struct (see sData 
% template sleep-swr project)

% This pipeline (1) filters various ephys channels in different frequency bands, 
%               (2) opens GUI for manual brain state scoring, merges NREM
%               bouts separated by microarousals;
%               (3) opens manual sorting of putative SWRs
%               (4) opens manual sorting of putative sleep spindles, and
%               (5) classifies SWRs as awake, spindle-coupled, and spindle-uncoupled. 

%% Locate file name in folder (for saving file)
folder_name = dir;
folder_cell = strfind({folder_name.name}, '.mat');
for i = 1:length(folder_cell)
    if folder_cell{1,i} > 0
    idx = i;
    end
end

if length(idx) > 1
    error('Idx is > 1, could not locate sData file name')
end

fileName = folder_name(idx).name;

%% Filtering the ephys signals

% Filter the raw electrophysiology signals into different frequency bands
sData = filter_ephys_sig(sData);
% 
clearvars -except sData fileName
save(fileName)

%% Sleep scoring
% Uses Begonia software from the Letten center


% ECOG used for brain state scoring (RSC channel)
LFP2          = sData.ephysdata2.lfp;
EMG           = sData.ephysdata3.lfp;
run_speed     = sData.daqdata.runSpeed;
signal_length = length(sData.ephysdata.lfp);
srate         = 2500;
delta_t       = 1/srate; % Sampling interval
tt            = delta_t:delta_t:signal_length*delta_t;

% Start Begonia manual sleep scorer
sData.episodes = begonia.processing.mark_sleep.mark_sleep(LFP2, tt, EMG, tt,run_speed,tt,[]);

clearvars -except sData fileName
save(fileName)

%% Merge NREM episodes separated by micro awakenings/microarousals

sData = merge_nrem(sData);
%% SWR analysis
% Run Anna's code for marking ripples. (Have adjusted SD for treshholding
% and window size for the moving mean threshold to better capture SWRs
% occurring in clusters.

sData = markRipples(sData,1,[]);

% Merge overlapping SWRs
sData = merge_ripples(sData);
% Remove duplicate SWRs
sData = delete_duplicate_swr(sData);

save(fileName)

%% Spindle analysis

% Inspect and mark spindles
sData = mark_spindle(sData,1);

% Sleep spindles are sometimes duplicated... need to fix. For now this code
% removes them. 
sData = deleteduplicate(sData);

% Find all sleep spindles occurring in NREM sleep episodes
sData = detect_NREM_spindles(sData);

save(fileName)
%% SWR classification
sData = swrspindle(sData);

clearvars -except sData fileName
save(fileName)
