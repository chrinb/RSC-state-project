function plot_missing_frames(sData)

cam_samples_sec = sData.daqdata.cam_frame_times;

%% Plot camera frame times and diff to check for missing frames
[val, val_idx] = unique(diff(cam_samples_sec));

frame_count_vector = 1:length( diff(cam_samples_sec));


med_val = median( diff(cam_samples_sec));

% Check for frame intervals above or below median sample rate, + - a small
% difference

above_median_idx = diff(cam_samples_sec) > med_val + 0.005;
below_median_idx = diff(cam_samples_sec) < med_val - 0.005;

if sum(above_median_idx) > 0
    idx1 = frame_count_vector(above_median_idx);
else
    idx1 = [];
end

if sum(below_median_idx) > 0
    idx2 = frame_count_vector(below_median_idx);
else
    idx2 = [];
end

missing_frame_idx = [idx1, idx2];

srate = find_imaging_framerate(sData);

time_vec = (0:length(cam_samples_sec)-1)/srate;

%% Plot
figure(1), clf
sgtitle(sData.sessionInfo.sessionID, 'interpreter', 'none')
h(1) = subplot(211);
plot( cam_samples_sec)

% missing_frame_idx = missing_frame_idx./srate;
if ~isempty(missing_frame_idx)

    for i = 1:numel(missing_frame_idx)
        xline(missing_frame_idx(i), 'r--')
    end
end
% try xline(missing_frame_idx, 'r--'), catch, end
title('Camera frame times')
legend({'frame times', 'Identified time stamps of missing frames'})

h(2) = subplot(212);
plot(diff(cam_samples_sec))

title('Diff of camera frame times')
xlabel('samples')
linkaxes(h, 'x' )

