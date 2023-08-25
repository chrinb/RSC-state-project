function plot_slow_wave_events(sData)

% Written by Christoffer Berge | Vervaeke Lab

% Function that plots deconvolved dF/F data and marks sleep episodes, slow
% oscillations, and delta waves. 

signal = sData.imdata.roiSignals(2).ciaDeconvolved;
nr_rois = size(signal,1);
%% Find sleep episodes
sessionID = sData.sessionInfo.sessionID;
% ECoG_spindle = sData.ephysdata2.spindleStartEnd/2500;

REM_episodes = rem_sleep(sData);
NREM_episodes = nrem_sleep(sData);

NREM_start_end = NREM_episodes./2500;
REM_start_end = REM_episodes./2500;

% %% Select spindle freq
% prompt = sprintf('8-18Hz (1) or 10-16hz (2)? ');
% spindle_band_select = input(prompt);
% 
% if spindle_band_select == 1
%     spindle_select = [];
% elseif spindle_band_select == 2
%     spindle_select = '1016';
% end
% 
% spin_start_end_str = strcat('spindleStartEnd', spindle_select);
% ECoG_spindle = sData.ephysdata2.(spin_start_end_str)/2500;
% 
% awake_swr_str                  = strcat('awake_swr', spindle_select);
% NREM_spindle_uncoupled_swr_str = strcat('NREM_spindle_uncoupled_swr', spindle_select);
% spindle_coupled_swr_str        = strcat('spindle_coupled_swr', spindle_select);

% %% Select SWRs for analysis
% prompt = sprintf('All ripples? (y = yes | everything else = no) ');
% allrip = input(prompt,'s');
% 
frames = sData.daqdata.frame_onset_reference_frame;
% 
signal_length = length(signal);

% find frames per second for imaging time vector
nr_of_frames_in_rec = length(signal);
recording_length_sec = length(sData.ephysdata.lfp)/2500;
frames_per_sec = nr_of_frames_in_rec/recording_length_sec;

time_imaging = (0:length(signal)-1)/frames_per_sec;
time_vec = (1:(frames_per_sec*20):signal_length);
tick_labels = time_vec/frames_per_sec;
% 
% if strcmp(allrip,'y') %keep all ripples
%     awakeSWRidx = sData.ephysdata.(awake_swr_str);
%     NREMspindleUncoupledSWRidx = sData.ephysdata.(NREM_spindle_uncoupled_swr_str);
%     NREMspindleCoupledSWRidx = sData.ephysdata.(spindle_coupled_swr_str);
%     
% else
%     prompt = sprintf('Remove locomotion SWR? (y = yes | everything else = no) ');
%     riprun = input(prompt, 's');
%     
%     prompt = sprintf('Remove temporally close SWR? (y = yes | everything else = no) ');
%     removerip = input(prompt, 's');
%     
%     % if remove locomotion SWR but not temporally close SWR
%     if strcmp(riprun, 'y') && ~strcmp(removerip, 'y')
%         [awakeSWRidx, NREMspindleUncoupledSWRidx, NREMspindleCoupledSWRidx ] = ripRunAn(sData); 
%     % if remove temporally close SWR but not locomotion SWR
%     elseif strcmp(removerip, 'y') && ~strcmp(riprun, 'y')
%         [awakeSWRidx, NREMspindleUncoupledSWRidx, NREMspindleCoupledSWRidx] = removeCloseRip(sData);
% 
%     % if remove both temporally close and locotion SWR
%     elseif strcmp(removerip, 'y') && strcmp(riprun, 'y')
%         [awakeSWRidx, NREMspindleUncoupledSWRidx, NREMspindleCoupledSWRidx ] = ripRunAn(sData,1);
%     end
% end

% %% Find SWR onset/offset
% [swr_start_stop,~] = mark_ripple_onset_offset(sData);
% 
% awake_swr_idx = ismember(sData.ephysdata.absRipIdx, awakeSWRidx);
% spindleC_swr_idx = ismember(sData.ephysdata.absRipIdx, NREMspindleCoupledSWRidx);
% spindleUC_swr_idx = ismember(sData.ephysdata.absRipIdx, NREMspindleUncoupledSWRidx);
% 
% temp_var1 = repmat(awake_swr_idx',1,2);
% temp_var2 = repmat(spindleC_swr_idx',1,2);
% temp_var3 = repmat(spindleUC_swr_idx',1,2);
% 
% swr_times1 = swr_start_stop(temp_var1);
% swr_times2 = swr_start_stop(temp_var2);
% swr_times3 = swr_start_stop(temp_var3);
% 
% resize_idx1 = length(swr_times1)/2;
% resize_idx2 = length(swr_times2)/2;
% resize_idx3 = length(swr_times3)/2;
% 
% awake_swr_onset_offset = reshape(swr_times1, resize_idx1,2 );
% spindleC_swr_onset_offset = reshape(swr_times2, resize_idx2,2 );
% spindleUC_swr_onset_offset = reshape(swr_times3, resize_idx3,2 );
[nrem_SOs, nrem_delta_waves] = mark_slow_wave(sData);


%% Plotting
figure,
sgtitle(sessionID),
% time_ephys = (1:length(sData.ephysdata.lfp))/2500;
% hAx(3) = subplot(313);
% plot(sData.ephysdata.lfp*0.2), hold on,
% plot(sData.ephysdata.ripplefreq-0.4);

hAx(1) = subplot(211);
x1 = [time_imaging(1) time_imaging(end)];
y1 = [1 size(signal,1)];
imagesc(x1, y1, signal), colormap(flipud(gray)), caxis([0 .05]),
hold on 

%% Plot sleep 
% NREM bouts
for i = 1:length(NREM_start_end)
    x = [NREM_start_end(i,1) NREM_start_end(i,1) NREM_start_end(i,2) NREM_start_end(i,2)];
    y = [1 nr_rois nr_rois 1];
    patch(x, y, 'green', 'edgecolor', 'none', 'FaceAlpha', .2);
end
% REM bouts
if length(REM_start_end) == 2
    a = [REM_start_end(1,1) REM_start_end(1,1) REM_start_end(1,2) REM_start_end(1,2)];
    b = [1 nr_rois nr_rois 1];
    patch(a, b, 'red', 'edgecolor', 'none', 'FaceAlpha', .2);
else
    for i = 1:length(REM_start_end)
        a = [REM_start_end(i,1) REM_start_end(i,1) REM_start_end(i,2) REM_start_end(i,2)];
        b = [1 nr_rois nr_rois 1];
        patch(a, b, 'red', 'edgecolor', 'none', 'FaceAlpha', .2);
    end
end



%% Plot slow oscillations and delta waves
for i = 1:length(nrem_SOs)
    x = [nrem_SOs(i,1) nrem_SOs(i,1) nrem_SOs(i,2) nrem_SOs(i,2)]/2500;
    y = [1 nr_rois nr_rois 1];
    patch(x, y, [0.8500 0.3250 0.0980], 'edgecolor', 'none', 'FaceAlpha', .4,'LineWidth',2);
end

for i = 1:length(nrem_delta_waves)
    x = [nrem_delta_waves(i,1) nrem_delta_waves(i,1) nrem_delta_waves(i,2) nrem_delta_waves(i,2)]/2500;
    y = [1 nr_rois nr_rois 1];
    patch(x, y, [0.4940 0.1840 0.5560], 'edgecolor', 'none', 'FaceAlpha', .4,'LineWidth',2);
end



%% Second plot
hAx(2) = subplot(212);
plot(time_imaging, smooth(nanmean(signal)))
set(gca, 'ylim',[0, 0.005])
set(gca, 'xlim',[time_imaging(1), time_imaging(end)])
hold on,

%% Plot sleep
% NREM bouts
for i = 1:length(NREM_start_end)
    x = [NREM_start_end(i,1) NREM_start_end(i,1) NREM_start_end(i,2) NREM_start_end(i,2)];
    y = [0 0.005 0.005 0];
    patch(x, y, 'green', 'edgecolor', 'none', 'FaceAlpha', .2);
end
% REM bouts
if length(REM_start_end) == 2
    a = [REM_start_end(1,1) REM_start_end(1,1) REM_start_end(1,2) REM_start_end(1,2)];
    b = [0 0.005 0.005 0];
    patch(a, b, 'red', 'edgecolor', 'none', 'FaceAlpha', .2);
else
    for i = 1:length(REM_start_end)
        a = [REM_start_end(i,1) REM_start_end(i,1) REM_start_end(i,2) REM_start_end(i,2)];
        b = [0 0.005 0.005 0];
        patch(a, b, 'red', 'edgecolor', 'none', 'FaceAlpha', .2);
    end
end


%% Plot slow oscillations and delta waves
for i = 1:length(nrem_SOs)
    x = [nrem_SOs(i,1) nrem_SOs(i,1) nrem_SOs(i,2) nrem_SOs(i,2)]/2500;
    y = [0 0.005 0.005 0];
    patch(x, y, [0.8500 0.3250 0.0980], 'edgecolor', 'none', 'FaceAlpha', .4,'LineWidth',2);
end

for i = 1:length(nrem_delta_waves)
    x = [nrem_delta_waves(i,1) nrem_delta_waves(i,1) nrem_delta_waves(i,2) nrem_delta_waves(i,2)]/2500;
    y = [0 0.005 0.005 0];
    patch(x, y, [0.4940 0.1840 0.5560], 'edgecolor', 'none', 'FaceAlpha', .4,'LineWidth',2);
end

linkaxes(hAx, 'x');