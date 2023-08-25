function plot_swrs_spindles(sData)

% Written by Christoffer Berge || Vervaeke lab

% Plots a figure consisting of multiple subplots:(1) raw and (2) 100-250 Hz 
% filtered CA1 LFP and all detected SWRs; and (3) raw and (4) 10-16 Hz 
% (sigma) filtered RSC ECoG, and the center of detected spindles. 

NREMep = sData.episodes.state == 'NREM' & sData.episodes.state_duration > 30;
% check that there is REM episodes in recording
if ~isempty(NREMep)
    NREMarr       = table2array(sData.episodes(NREMep,2:3));
    NREM_start_end = NREMarr;
end
sessionID      = sData.sessionInfo.sessionID;
REM_episodes   = rem_sleep(sData);
REM_start_end  = REM_episodes./2500;


%% Select spindle freq
prompt = sprintf('8-18Hz (1) or 10-16hz (2)? ');
spindle_band_select = input(prompt);

if spindle_band_select == 1
    ephys_select   = [];
    spindle_select = [];
elseif spindle_band_select == 2
    ephys_select   = num2str(2);
    spindle_select = '1016';
end

unclassified_swr_str           = strcat('unclassified_swr', spindle_select);
NREM_spindle_uncoupled_swr_str = strcat('NREM_spindle_uncoupled_swr', spindle_select);
spindle_coupled_swr_str        = strcat('spindle_coupled_swr', spindle_select);
spin_idx_str                   = strcat('NREMAbsSpindleIdx', spindle_select);
nr_of_spindles_in_analysis     = length(sData.ephysdata2.(spin_idx_str) );

%% Select SWRs for analysis
prompt = sprintf('All SWRs (1) | Rem. locom. SWRs (2) | Rem. clust. SWRs (3) | Rem- locom. & clust. SWRs (4) ');
swr_for_analysis = input(prompt);

frames  = sData.daqdata.frame_onset_reference_frame;
swr_idx = sData.ephysdata.absRipIdx;

if swr_for_analysis == 1
    unclassified_SWR_idx = sData.ephysdata.(unclassified_swr_str);
    NREM_spinU_SWR_idx   = sData.ephysdata.(NREM_spindle_uncoupled_swr_str);
    NREM_spinC_SWR_idx   = sData.ephysdata.(spindle_coupled_swr_str);
elseif swr_for_analysis == 2
    [RippleIdx,~] = riprun2(sData, swr_idx);
elseif swr_for_analysis == 3
    [RippleIdx] = RemoveRip(swr_idx);
elseif swr_for_analysis == 4
    [RippleIdx,~] = riprun2(sData, swr_idx);
    [RippleIdx]   = RemoveRip(RippleIdx);
end

if swr_for_analysis ~= 1 
    temp_varx            = ismember( RippleIdx,sData.ephysdata.(unclassified_swr_str));
    unclassified_SWR_idx = RippleIdx(temp_varx);

    temp_vary          = ismember(RippleIdx,sData.ephysdata.(NREM_spindle_uncoupled_swr_str));
    NREM_spinU_SWR_idx = RippleIdx(temp_vary);

    temp_varz          = ismember(RippleIdx,sData.ephysdata.(spindle_coupled_swr_str));
    NREM_spinC_SWR_idx = RippleIdx(temp_varz);
end

%% Find SWR onset/offset
[swr_start_stop,~] = mark_ripple_onset_offset(sData);

unclassified_swr_idx = ismember(sData.ephysdata.absRipIdx, unclassified_SWR_idx);
spindleC_swr_idx     = ismember(sData.ephysdata.absRipIdx, NREM_spinC_SWR_idx);
spindleUC_swr_idx    = ismember(sData.ephysdata.absRipIdx, NREM_spinU_SWR_idx);

temp_var1 = repmat(unclassified_swr_idx',1,2);
temp_var2 = repmat(spindleC_swr_idx',1,2);
temp_var3 = repmat(spindleUC_swr_idx',1,2);

swr_times1 = swr_start_stop(temp_var1);
swr_times2 = swr_start_stop(temp_var2);
swr_times3 = swr_start_stop(temp_var3);

resize_idx1 = length(swr_times1)/2;
resize_idx2 = length(swr_times2)/2;
resize_idx3 = length(swr_times3)/2;

unclassified_swr_onset_offset = reshape(swr_times1, resize_idx1,2 );
spindleC_swr_onset_offset     = reshape(swr_times2, resize_idx2,2 );
spindleUC_swr_onset_offset    = reshape(swr_times3, resize_idx3,2 );


%% Get ephys data
temp_str          = strcat('sigmaband', num2str(ephys_select) );
raw_lfp           = sData.ephysdata.lfp;
filt_lfp          = sData.ephysdata.ripplefreq;
raw_ECoG          = sData.ephysdata2.lfp;
ECoG_filt         = sData.ephysdata2.lfpFilt;
sigma_ECoG        = sData.ephysdata2.(temp_str);
delta_ECoG        = sData.ephysdata2.deltaband;
spindle_start_end = strcat('NREMspindleStartEnd', spindle_select);
ECoG_spindle      = sData.ephysdata2.(spindle_start_end)/2500;
time              = (0:length(sData.ephysdata.lfp)-1)/2500;

%% Plot figures
figure,
hAx(1) = subplot(211);
plot(time, raw_lfp*0.2), hold on,
plot(time, filt_lfp-0.3)
title(sessionID)
set(gca, 'xlim', [0 max(time)]);
ylabel('mV')
title('CA1 LFP')

% Loop over unclassified SWRs and plot them as individual patches
for i = 1:length(unclassified_swr_onset_offset)
    x = [unclassified_swr_onset_offset(i,1) unclassified_swr_onset_offset(i,1) unclassified_swr_onset_offset(i,2) unclassified_swr_onset_offset(i,2)]/2500;
    y = [-1 1 1 -1];
    h1 = patch(x, y, [0.8500 0.3250 0.0980], 'edgecolor', 'none', 'FaceAlpha', .4,'LineWidth',2);
end
% Loop over spindle-coupled SWRs and plot them as individual patches
for i = 1:length(spindleC_swr_onset_offset)
    x = [spindleC_swr_onset_offset(i,1) spindleC_swr_onset_offset(i,1) spindleC_swr_onset_offset(i,2) spindleC_swr_onset_offset(i,2)]/2500;
    y = [-1 1 1 -1];
    h2 = patch(x, y, [0.4940 0.1840 0.5560], 'edgecolor', 'none', 'FaceAlpha', .4,'LineWidth',2);
end
% Loop over spindle-uncoupled SWRs and plot them as individual patches
for i = 1:length(spindleUC_swr_onset_offset)
    x = [spindleUC_swr_onset_offset(i,1) spindleUC_swr_onset_offset(i,1) spindleUC_swr_onset_offset(i,2) spindleUC_swr_onset_offset(i,2)]/2500;
    y = [-1 1 1 -1];
    h3 = patch(x, y, [0.9290 0.6940 0.1250], 'edgecolor', 'none', 'FaceAlpha', .5,'LineWidth',2);
end
set(gca, 'ylim', [-.6 .2])


% NREM bouts
for i = 1:length(NREM_start_end)
    x = [NREM_start_end(i,1) NREM_start_end(i,1) NREM_start_end(i,2) NREM_start_end(i,2)];
    y = [-1 1 1 -1];
    h4 = patch(x, y, 'black', 'edgecolor', 'none', 'FaceAlpha', .1);
end
set(gca, 'ylim', [-.6 .4])

 % REM bouts
if length(REM_start_end) == 2
    a = [REM_start_end(1,1) REM_start_end(1,1) REM_start_end(1,2) REM_start_end(1,2)];
    b = [-1 1 1 -1];
    h5 = patch(a, b, 'red', 'edgecolor', 'none', 'FaceAlpha', .2);
else
    for i = 1:length(REM_start_end)
        a = [REM_start_end(i,1) REM_start_end(i,1) REM_start_end(i,2) REM_start_end(i,2)];
        b = [-1 1 1 -1];
        h5 = patch(a, b, 'red', 'edgecolor', 'none', 'FaceAlpha', .2);
    end
end

hAx(2) = subplot(212);
plot(time, raw_ECoG*0.3), hold on,
plot(time, sigma_ECoG-0.3)
plot(time, delta_ECoG+0.2)
title(sessionID)
set(gca, 'xlim', [0 max(time)]);
ylabel('mV')
xlabel('Time (s)')
title('RSC ECoG')
%% Plot sleep spindles 
for k = 1:length(sData.ephysdata2.(spindle_start_end))
    z = [ECoG_spindle(k,1) ECoG_spindle(k,1) ECoG_spindle(k,2) ECoG_spindle(k,2) ];
    v = [-1 1 1 -1];
    h6 = patch(z,v, 'magenta', 'edgecolor', 'none', 'FaceAlpha', .2);
end

% NREM bouts
for i = 1:length(NREM_start_end)
    x = [NREM_start_end(i,1) NREM_start_end(i,1) NREM_start_end(i,2) NREM_start_end(i,2)];
    y = [-1 1 1 -1];
    h4 = patch(x, y, 'black', 'edgecolor', 'none', 'FaceAlpha', .1);
end
 set(gca, 'ylim', [-.6 .4])

 % REM bouts
if length(REM_start_end) == 2
    a = [REM_start_end(1,1) REM_start_end(1,1) REM_start_end(1,2) REM_start_end(1,2)];
    b = [-1 1 1 -1];
    h5 = patch(a, b, 'red', 'edgecolor', 'none', 'FaceAlpha', .2);
else
    for i = 1:length(REM_start_end)
        a = [REM_start_end(i,1) REM_start_end(i,1) REM_start_end(i,2) REM_start_end(i,2)];
        b = [-1 1 1 -1];
        h5 = patch(a, b, 'red', 'edgecolor', 'none', 'FaceAlpha', .2);
    end
end

legend( [h1, h2, h3, h4, h5, h6], 'Unclassified','Coupled', 'Uncoupled', ...
    'NREM', 'REM', 'Spindle');

 
linkaxes(hAx, 'x')