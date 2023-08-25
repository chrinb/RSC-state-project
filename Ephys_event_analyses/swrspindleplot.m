function plot_swrs_spindles(sData)

% Written by Christoffer Berge || Vervaeke lab

% Plots a figure consisting of multiple subplots:(1) raw and (2) 100-250 Hz 
% filtered CA1 LFP and all detected SWRs; and (3) raw and (4) 8-18 Hz 
% (sigma) filtered RSC ECoG, and the center of detected spindles. 

% Optional prompt to only select a specific part of the e-phys traces. If
% 'y' is chosen as input, specify first the start of the interval in
% seconds and then the end of the interval. 

% prompt = sprintf('Select time interval? ');
% signalSelect = input(prompt,'s');
% 
% if strcmp(signalSelect,'y')
%     prompt = sprintf('Start (sec): ');
%     start = input(prompt);
%     prompt = sprintf('End (sec): ');
%     stop = input(prompt);
% else
%     start = -1;
% end
ECoG_spindle = sData.ephysdata2.spindleStartEnd/2500;

[~, NREM_start_end] = swrspindle(sData);
NREM_start_end = NREM_start_end./2500;
sessionID = sData.sessionInfo.sessionID;
%% Select SWRs for analysis
prompt = sprintf('All ripples? (y = yes | everything else = no) ');
allrip = input(prompt,'s');

if strcmp(allrip,'y') %keep all ripples
    awakeSWRidx = sData.ephysdata.awake_swr;
    NREMspindleUncoupledSWRidx = sData.ephysdata.NREM_spindle_uncoupled_swr;
    NREMspindleCoupledSWRidx = sData.ephysdata.spindle_coupled_swr;
    
else
    prompt = sprintf('Remove locomotion SWR? (y = yes | everything else = no) ');
    riprun = input(prompt, 's');
    
    prompt = sprintf('Remove temporally close SWR? (y = yes | everything else = no) ');
    removerip = input(prompt, 's');
    
    % if remove locomotion SWR but not temporally close SWR
    if strcmp(riprun, 'y') && ~strcmp(removerip, 'y')
        [awakeSWRidx, NREMspindleUncoupledSWRidx, NREMspindleCoupledSWRidx ] = ripRunAn(sData); 
    % if remove temporally close SWR but not locomotion SWR
    elseif strcmp(removerip, 'y') && ~strcmp(riprun, 'y')
        [awakeSWRidx, NREMspindleUncoupledSWRidx, NREMspindleCoupledSWRidx] = removeCloseRip(sData);

    % if remove both temporally close and locotion SWR
    elseif strcmp(removerip, 'y') && strcmp(riprun, 'y')
        [awakeSWRidx, NREMspindleUncoupledSWRidx, NREMspindleCoupledSWRidx ] = ripRunAn(sData,1);
    end
end

%% Find SWR onset/offset
[swr_start_stop,~] = mark_ripple_onset_offset(sData);

awake_swr_idx = ismember(sData.ephysdata.absRipIdx, awakeSWRidx);
spindleC_swr_idx = ismember(sData.ephysdata.absRipIdx, NREMspindleCoupledSWRidx);
spindleUC_swr_idx = ismember(sData.ephysdata.absRipIdx, NREMspindleUncoupledSWRidx);

temp_var1 = repmat(awake_swr_idx',1,2);
temp_var2 = repmat(spindleC_swr_idx',1,2);
temp_var3 = repmat(spindleUC_swr_idx',1,2);

swr_times1 = swr_start_stop(temp_var1);
swr_times2 = swr_start_stop(temp_var2);
swr_times3 = swr_start_stop(temp_var3);

resize_idx1 = length(swr_times1)/2;
resize_idx2 = length(swr_times2)/2;
resize_idx3 = length(swr_times3)/2;

awake_swr_onset_offset = reshape(swr_times1, resize_idx1,2 );
spindleC_swr_onset_offset = reshape(swr_times2, resize_idx2,2 );
spindleUC_swr_onset_offset = reshape(swr_times3, resize_idx3,2 );


%% Get ephys data

raw_lfp = sData.ephysdata.lfp;
filt_lfp = sData.ephysdata.ripplefreq;
raw_ECoG = sData.ephysdata2.lfp;
filt_ECoG = sData.ephysdata2.sigmaband;

SWR_freq_power = abs(sData.ephysdata.ripplefreq).^2;
spindle_freq_power = abs(sData.ephysdata2.sigmaband).^2;

swrIdx = zeros(1, length(sData.ephysdata.lfp));
swrIdx(sData.ephysdata.absRipIdx) = 0.02;
spindleIdx = zeros(1, length(sData.ephysdata.lfp));
spindleIdx(round(sData.ephysdata2.absSpindleIdx)) = max(spindle_freq_power);

time = (0:length(sData.ephysdata.lfp)-1)/2500;

figure,
hAx(1) = subplot(211);
plot(time, raw_lfp*0.2), hold on,
plot(time, filt_lfp-0.3)
title(sessionID)
set(gca, 'xlim', [0 max(time)]);
ylabel('mV')
title('CA1 LFP')

for i = 1:length(awake_swr_onset_offset)
    x = [awake_swr_onset_offset(i,1) awake_swr_onset_offset(i,1) awake_swr_onset_offset(i,2) awake_swr_onset_offset(i,2)]/2500;
    y = [-1 1 1 -1];
    patch(x, y, [0.8500 0.3250 0.0980], 'edgecolor', 'none', 'FaceAlpha', .4,'LineWidth',2);
end
for i = 1:length(spindleC_swr_onset_offset)
    x = [spindleC_swr_onset_offset(i,1) spindleC_swr_onset_offset(i,1) spindleC_swr_onset_offset(i,2) spindleC_swr_onset_offset(i,2)]/2500;
    y = [-1 1 1 -1];
    patch(x, y, [0.4940 0.1840 0.5560], 'edgecolor', 'none', 'FaceAlpha', .4,'LineWidth',2);
end

for i = 1:length(spindleUC_swr_onset_offset)
    x = [spindleUC_swr_onset_offset(i,1) spindleUC_swr_onset_offset(i,1) spindleUC_swr_onset_offset(i,2) spindleUC_swr_onset_offset(i,2)]/2500;
    y = [-1 1 1 -1];
    patch(x, y, [0.9290 0.6940 0.1250], 'edgecolor', 'none', 'FaceAlpha', .5,'LineWidth',2);
end
 set(gca, 'ylim', [-.6 .2])



hAx(2) = subplot(212);
plot(time, raw_ECoG*0.3), hold on,
plot(time, filt_ECoG-0.3)
title(sessionID)
set(gca, 'xlim', [0 max(time)]);
ylabel('mV')
title('RSC ECoG')
%% Plot sleep spindles 
for k = 1:length(sData.ephysdata2.spindleStartEnd)
    z = [ECoG_spindle(k,1) ECoG_spindle(k,1) ECoG_spindle(k,2) ECoG_spindle(k,2) ];
    v = [-1 1 1 -1];
    patch(z,v, 'magenta', 'edgecolor', 'none', 'FaceAlpha', .2);
end
 set(gca, 'ylim', [-.6 .2])

linkaxes(hAx, 'x')