function plot_deconv_events_awake(sData, params)

% Written by Christoffer Berge | Vervaeke Lab

% Function that plots deconvolved dF/F data and marks SWRs. 

try 
    ecog          = sData.ephysdata2.lfp; 
catch
    ecog          = sData.ephysdata.lfp; 
end

hpc           = sData.ephysdata.lfp;
if nargin > 1
    if strcmp( params.exp_type, 'default')
        signal        = sData.imdata.roiSignals(2).ciaDeconvolved;
    elseif strcmp( params.exp_type, 'axon')
        signal        = sData.imdata.roiSignals(2).mergedAxonsDec;
    end
else
    signal        = sData.imdata.roiSignals(2).ciaDeconvolved;
end
nr_rois       = size(signal,1);
frames        = sData.daqdata.frame_onset_reference_frame;
signal_length = length(signal);
sessionID     = sData.sessionInfo.sessionID;

% find frames per second for imaging time vector
nr_of_frames_in_rec = length(signal);
recording_length_sec = length(sData.ephysdata.lfp)/2500;
imag_srate = nr_of_frames_in_rec/recording_length_sec;

time_imaging = linspace(0, length(signal), length(signal))/imag_srate;
time_ephys   = linspace(0, length(ecog), length(ecog))/2500;

% time_imaging = (0:length(signal)-1)/frames_per_sec;
time_vec = (1:(imag_srate*20):signal_length);

SWR_idx = sData.ephysdata.absRipIdx;
 
%% Find SWR onset/offset
[swr_start_stop,~]     = mark_ripple_onset_offset(sData);
% temp_var1              = repmat(SWR_idx',1,2);
% swr_times1             = swr_start_stop(temp_var1);
% resize_idx1            = length(swr_times1)/2;
% awake_swr_onset_offset = reshape(swr_times1, resize_idx1,2 );
awake_swr_onset_offset = swr_start_stop;
%% Plot ECoG signal
figure,
sgtitle(sessionID,'Interpreter', 'none') 

hAx(1) = subplot(412);
plot(time_ephys, hpc);
ylabel('Vm', 'FontSize',14)
set(gca, 'xlim',[time_ephys(1), time_ephys(end)])
title('HPC CA1 LFP', 'FontSize',10)

hAx(2) = subplot(411);
plot(time_ephys, ecog);
ylabel('Vm', 'FontSize',14)
set(gca, 'xlim',[time_ephys(1), time_ephys(end)])
title('ECoG', 'FontSize',10)
% Plot roi x deconvolved events

hAx(3) = subplot(413);
y1     = [1 size(signal,1)];
imagesc(time_imaging, y1, signal), colormap(flipud(gray)), caxis([0 .05]),
ylabel('# ROIs', 'FontSize',14)

hold on 

% Plot SWRs
for i = 1:length(awake_swr_onset_offset)
    x = [awake_swr_onset_offset(i,1) awake_swr_onset_offset(i,1) awake_swr_onset_offset(i,2) awake_swr_onset_offset(i,2)]/2500;
    y = [1 nr_rois nr_rois 1];
    patch(x, y, [0.8500 0.3250 0.0980], 'edgecolor', 'none', 'FaceAlpha', .4,'LineWidth',2);
end
%% Plot smoothed mean deconvolved events

hAx(4) = subplot(414);
y_max = max( smooth(mean(signal, 'omitnan')));
plot(time_imaging, smooth(mean(signal, 'omitnan')))
xlabel('Time (s)', 'FontSize',14);
ylabel('Spike rate (Hz)', 'FontSize',14);

set(gca, 'ylim',[0, y_max])
set(gca, 'xlim',[time_imaging(1), time_imaging(end)])
title('Mean deconv. DF/F (smoothed)', 'FontSize',10)
hold on,


% Plot SWRs
for i = 1:length(awake_swr_onset_offset)
    x = [awake_swr_onset_offset(i,1) awake_swr_onset_offset(i,1) awake_swr_onset_offset(i,2) awake_swr_onset_offset(i,2)]/2500;
    y = [0 y_max y_max 0];
    patch(x, y, [0.8500 0.3250 0.0980], 'edgecolor', 'none', 'FaceAlpha', .4,'LineWidth',2);
end

linkaxes(hAx, 'x');
