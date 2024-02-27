function sData = calc_cam_motion(sessionObject)

%% Get path to facemap output folder
root_path     = sessionObject.DataLocation(2).RootPath;
folder_path   = sessionObject.DataLocation(2).Subfolders;
mouse_folder  = [root_path,'\', folder_path];
facemap_folder_path     = [mouse_folder, '\FaceCam'];
files = dir(fullfile(facemap_folder_path, '*proc.mat'));
 % Check if any files were found
if isempty(files)
    disp('No proc.mat files found in the specified folder.');
else
    % Get the first file (assuming there's only one .tdms file in the folder)
    face_file_name = files(1).name;
end

sData = sessionObject.Data.sData;

%% Load Facemap output
facemap_data_path = [facemap_folder_path,'\' face_file_name];

data = load(facemap_data_path);

start_frame = sData.daqdata.camera_start_frame  ;


threshold  = 5*std(abs( data.motSVD_4(start_frame:start_frame+100, 1)));
    
%% Plot eye ROI 
figure(1), clf 
hold on
plot( abs(data.motSVD_4(start_frame:end, 1)),'LineWidth',2)
yline(threshold, 'r--', 'LineWidth',2)
legend({'1st PC, eye ROI', 'Threshold'})
set(gca, 'xlim', [0 100])
xlabel('Frames', 'FontSize',14)

[~,locs] = findpeaks( abs( data.motSVD_4(start_frame:start_frame+100, 1)), "MinPeakHeight",threshold);

% If more than one peak, go with the first one
if numel(locs) > 1
    locs = locs(1);
end

% Correct for 1 frame difference
locs_corr = locs-1;

true_start_frame = start_frame+locs_corr;

sData.daqdata.true_camera_start_frame = true_start_frame;

% Store the mean of the absolute values (bc components can have opposite
% signs) of the top 3 components
sData.analysis.paw_data     = double( mean( abs( data.motSVD_1(true_start_frame:end, 1:3) ),2) );
sData.analysis.face_data    = double( mean( abs( data.motSVD_2(true_start_frame:end, 1:3) ),2) );
sData.analysis.whisker_data = double( mean( abs( data.motSVD_3(true_start_frame:end, 1:3) ),2) );

%% Find closest frame betweeen face cam and two-photon 

% Use the frame when the 2P light shines through the eye ("locs") as the
% first frame in the "cam_frame_times" variable. Then find closest
% corresponding 2P frame for rest of recording.

cam_samples_sec            = sData.daqdata.cam_frame_times;
cam_frame_time_at_2P_onset = cam_samples_sec(true_start_frame);

corrected_t = cam_samples_sec-cam_frame_time_at_2P_onset;
corrected_t_trim = corrected_t(true_start_frame:end);



corr_cam_frame_times = corrected_t_trim(2:end);
sData.daqdata.cam_frame_times_corr = corr_cam_frame_times;

closest_matches = zeros(size(corr_cam_frame_times));

% Loop over each timestamp in vector1
for i = 1:length(corr_cam_frame_times)

    % Compute the absolute differences between the current timestamp in vector1 and all timestamps in vector2
    time_diffs = abs(corr_cam_frame_times(i) - sData.daqdata.two_photon_frame_times);
    
    % Find the index of the minimum difference
    [~, min_index] = min(time_diffs);
    
    % Store the closest match
    closest_matches(i) = sData.daqdata.two_photon_frame_times(min_index);
end

sData.daqdata.camera_2p_closest_match  = closest_matches;

%% Plot 1st component motion SVDs
cam_srate = sData.daqdata.faceCam_fps;

time_vec = (0:length(data.motSVD_2(true_start_frame:end, 1))-1)/cam_srate;

figure(5), clf 
hold on
sgtitle('Mean absolute values top 3 components (SVD)')
h(1) = subplot(313);
% plot(time_vec, data.motSVD_1(true_start_frame:end, 1))
plot(time_vec, mean(abs(data.motSVD_1(true_start_frame:end, 1:3)),2))

title('Paw')
xlabel('Time (s)', FontSize=14)

h(2) = subplot(312);
% plot(time_vec,  data.motSVD_2(true_start_frame:end, 1))
plot(time_vec, mean(abs(data.motSVD_2(true_start_frame:end, 1:3)),2))

title('Face')
h(3) = subplot(311);
% plot(time_vec, data.motSVD_3(true_start_frame:end, 1))
plot(time_vec, mean(abs(data.motSVD_3(true_start_frame:end, 1:3)),2))

title('Whisker')
linkaxes(h, 'x')
set(gca, 'xlim', [time_vec(1) time_vec(end)])


%% Calculate movement thresholds
srate = find_imaging_framerate(sData);


face_threshold    = movstd(sData.analysis.face_data, 31) > 1;
paw_threshold     = movstd(sData.analysis.paw_data, 31) > 1;
whisker_threshold = movstd(sData.analysis.whisker_data, 31) > 1;

all_threshold = sum( [face_threshold' ;paw_threshold' ;whisker_threshold'], 1);
all_threshold = all_threshold > 0;


% Add 0.5 seconds to all movement bouts as safety margin
[mov_start, mov_stop] = findTransitions(all_threshold);
adjusted_threshold = false(1, length(all_threshold));

for movement_nr = 1:length(mov_start)

    tmp_mov_start = mov_start(movement_nr) - round(srate/2);
    tmp_mov_stop  = mov_stop(movement_nr) + round(srate/2);

    if tmp_mov_start < 0 
        tmp_mov_start = 1;
    end

    if  tmp_mov_stop > length(paw_threshold)
        tmp_mov_stop = length(paw_threshold);
    end
    
    adjusted_threshold(tmp_mov_start:tmp_mov_stop) = true;
end


sData.analysis.movement_threshold = adjusted_threshold;

[~, idx] = min(abs(sData.daqdata.camera_2p_closest_match(end)-sData.daqdata.two_photon_frame_times  ));

% Correction for those (few) sessions where face camera sample rate was
% ~60Hz OR dealing with missing face camera frames
if sData.daqdata.faceCam_fps < 50 && idx < numel(sData.daqdata.camera_2p_closest_match)
    sData.analysis.cutoff_frame = idx;
else
    [~, idx] = min(abs(sData.daqdata.two_photon_frame_times(idx) - sData.daqdata.camera_2p_closest_match ));
        sData.analysis.cutoff_frame = idx;
end

% Calculate QW duration in rec.
qw_frames = sData.analysis.movement_threshold(1:sData.analysis.cutoff_frame);
qw_frames(qw_frames > 0) = [];

sData.analysis.qw_duration = length(qw_frames)/sData.daqdata.faceCam_fps;

%% Plot more data
% 
% root_path     = sessionObjects.DataLocation(2).RootPath;
% folder_path   = sessionObjects.DataLocation(2).Subfolders;
% mouse_folder  = [root_path,'\', folder_path];
% facemap_folder_path     = [mouse_folder, '\motion_corrected'];
% files = dir(fullfile(facemap_folder_path, '*stats.mat'));
%  % Check if any files were found
% if isempty(files)
%     disp('No proc.mat files found in the specified folder.');
% else
%     % Get the first file (assuming there's only one .tdms file in the folder)
%     motion_corr_file_name = files(1).name;
% end
% 
% motion_corr_data_path = [facemap_folder_path,'\' motion_corr_file_name];
% 
% motion_corr_stats = load(motion_corr_data_path);
% 
% figure(3), clf
% 
% time_imaging  = linspace(1, size(dff,2), size(dff,2) )/31;
% y1 = [1 size(dff,1)];
% 
% h(1) = subplot(311); hold on
% imagesc(time_imaging, y1, dffD)
% caxis([0 .4])
% title('DF/F')
% 
% y_data = mean(dffD);
% y_data_scaled = (y_data - min(y_data))/(max(y_data)-min(y_data));
% nTrials = 700;
% new_min = nTrials;
% new_max = nTrials-200;
% y_data_scaled = nTrials + y_data_scaled*(new_max-new_min);
% plot(time_imaging, y_data_scaled, 'w', 'LineWidth',1)
% 
% h(2) = subplot(312);
% plot(time_vec, sData.analysis.paw_data)
% hold on
% plot(time_vec, sData.analysis.face_data)
% legend({'Paw', 'Face'})
% title('Face camera data')
% h(3) = subplot(313);
% 
% plot(time_imaging, motion_corr_stats.MotionCorrectionStats{1, 1}.offsetX  )
% hold on
% plot(time_imaging, motion_corr_stats.MotionCorrectionStats{1, 1}.offsetY  )
% linkaxes(h, 'x')
% set(gca, 'xlim', [time_vec(1) time_vec(end)])
% legend({'X offset', 'Y offset'})
% title('Motion correction')