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

%% Find closest frame betweeen face cam and two-photon 

% Use the frame when the 2P light shines through the eye ("locs") as the
% first frame in the "cam_frame_times" variable. Then find closest
% corresponding 2P frame for rest of recording.

cam_samples_sec            = sData.daqdata.cam_frame_times;
cam_frame_time_at_2P_onset = cam_samples_sec(true_start_frame);

adjusted_missing_frame_idx = sData.daqdata.missing_cam_frames > true_start_frame;
adjusted_missing_frame_idx = sData.daqdata.missing_cam_frames(adjusted_missing_frame_idx)-true_start_frame;

corrected_t      = cam_samples_sec-cam_frame_time_at_2P_onset;
corrected_t_trim = corrected_t(true_start_frame:end);

corr_cam_frame_times = corrected_t_trim;
sData.daqdata.cam_frame_times_corr = corr_cam_frame_times;

n_2p_samples = size( sData.imdata.roiSignals(2).newdff, 2);

% closest_matches = zeros(size(corr_cam_frame_times));
[idx_two_photon, facecam_match_with_2p] = deal( zeros(n_2p_samples,1  ));

% Loop over face cam frame times and find the closest 2P frame
for i = 1:length(corr_cam_frame_times)
    
    if isnan( corr_cam_frame_times(i))
        facecam_match_with_2p(i) = NaN;
        idx_two_photon(i) = NaN;
    else

    % Find the index in the 2P vector that is closest to the current face cam frame vector    
    time_diffs         = abs(corr_cam_frame_times(i) - sData.daqdata.two_photon_frame_times);    
    [minval, min_index]     = min(time_diffs);
    % facecam_match_with_2p(i) = sData.daqdata.two_photon_frame_times(min_index);
    % facecam_match_with_2p(i) =corr_cam_frame_times(min_index);
    
    % IMPORTANT: this variable contains indices pointing to the two-photon
    % frame at the same time as the face cam. For example, if face cam
    % frame 10000 = 323.0760, the index is 9992, corresponding to 2P data
    % frame 9992 = 323.0638
    idx_two_photon(i) = min_index;
    end
end
idx_two_photon(idx_two_photon==0) = NaN;

%% Test 

idx_motion = zeros(1, n_2p_samples);

for i = 1:n_2p_samples
    
    % if sData.daqdata.two_photon_frame_times(i) < min_sig(end)
    %     i
        % Find the index in the 2P vector that is closest to the current face cam frame vector    
        time_diffs         = abs(sData.daqdata.two_photon_frame_times(i) - corr_cam_frame_times);    
        [min_val, min_index]     = min(time_diffs);
        
        % While loop index is smaller than nr of frames in face cam, check
        % if corresponding face cam frame is NaN
        if i < numel(corr_cam_frame_times)
            if isnan( corr_cam_frame_times(i))
                idx_motion(i) = NaN;
            else
                idx_motion(i) = min_index;
            end
        end

    % end

end
idx_motion(idx_motion==0)= NaN;


%% Align motion vectors from face camp with 2p data

% First insert NaNs into motion vector based on face cam frame time vector
nan_indices = isnan(corr_cam_frame_times);

paw_data_adjusted     = double( mean( abs( data.motSVD_1(true_start_frame:end, 1:3) ),2) );
face_data_adjusted    = double( mean( abs( data.motSVD_2(true_start_frame:end, 1:3) ),2) );
whisker_data_adjusted = double( mean( abs( data.motSVD_3(true_start_frame:end, 1:3) ),2) );

% paw_data_adjusted(nan_indices) = NaN;

nan_stretches = find( [false; diff(nan_indices) == 1]);

% Update corresponding indices in Y
for i = 1:length(nan_stretches)
    start_idx = nan_stretches(i);
    end_idx   = find(~nan_indices(start_idx:end), 1, 'first') + start_idx - 1;

    tmp_nan_vec = NaN( length( start_idx:end_idx-1), 1);

    paw_data_adjusted     = [ paw_data_adjusted(1:start_idx-1); tmp_nan_vec; paw_data_adjusted(start_idx+1:end)];
    face_data_adjusted    = [ face_data_adjusted(1:start_idx-1); tmp_nan_vec; face_data_adjusted(start_idx+1:end)];
    whisker_data_adjusted = [ whisker_data_adjusted(1:start_idx-1); tmp_nan_vec; whisker_data_adjusted(start_idx+1:end)];

end

%% Match motion with 2p vector
[paw_motion_2p_match, face_motion_2p_match, whisker_motion_2p_match] = deal( zeros(n_2p_samples, 1));

for i = 1:n_2p_samples 
    
    
    tmp_idx = idx_motion(i);
    if isnan(tmp_idx)       
        paw_motion_2p_match(i)     = NaN;
        face_motion_2p_match(i)    = NaN;
        whisker_motion_2p_match(i) = NaN;

    elseif tmp_idx > length(paw_data_adjusted)
        paw_motion_2p_match(i)     = NaN;
        face_motion_2p_match(i)    = NaN;
        whisker_motion_2p_match(i) = NaN;
    else
        paw_motion_2p_match(i)     = paw_data_adjusted( tmp_idx);
        face_motion_2p_match(i)    = face_data_adjusted( tmp_idx);
        whisker_motion_2pt_match(i) = whisker_data_adjusted( tmp_idx);
    end

end

% Store the mean of the absolute values (bc components can have opposite
% signs) of the top 3 components
sData.analysis.paw_data     = paw_motion_2p_match;
sData.analysis.face_data    = face_motion_2p_match;
sData.analysis.whisker_data = whisker_motion_2p_match;

%% Plot mean of the absolute values of top 3 SVD components
cam_srate = sData.daqdata.faceCam_fps;

time_vec = (0:length(sData.analysis.paw_data)-1)/cam_srate;

figure(5), clf 
hold on
sgtitle('Mean absolute values top 3 components (SVD)')
h(1) = subplot(313);
plot(time_vec, sData.analysis.paw_data)
title('Paw')

h(2) = subplot(312);
plot(time_vec, sData.analysis.face_data)
title('Face')

h(3) = subplot(311);
plot(time_vec, sData.analysis.whisker_data)
title('Whisker')

linkaxes(h, 'x')
set(gca, 'xlim', [time_vec(1) time_vec(end)])
xlabel('Time (s)', FontSize=14)

%% Calculate movement thresholds
srate = find_imaging_framerate(sData);

window = sData.daqdata.faceCam_fps;

if window > 40

    face_threshold    = ( movstd(sData.analysis.face_data, window) - mean(movstd(sData.analysis.face_data, window)) ) > 1;
    paw_threshold     = movstd(sData.analysis.paw_data, window) - mean(movstd(sData.analysis.paw_data, window)) > 1;
    whisker_threshold = movstd(sData.analysis.whisker_data, window) - mean(movstd(sData.analysis.whisker_data, window)) > 1;
    msgbox([sData.sessionInfo.sessionID, ' window is larger than 40 samples, investigate'])
else
    face_threshold    = movstd(sData.analysis.face_data, window)   > 1;
    paw_threshold     = movstd(sData.analysis.paw_data, window) > 1;
    whisker_threshold = movstd(sData.analysis.whisker_data, window)  > 1;
end
all_threshold = sum( [face_threshold' ;paw_threshold' ;whisker_threshold'], 1);
all_threshold = all_threshold > 0;

% Add 0.5 seconds to all movement bouts as safety margin
[mov_start, mov_stop] = findTransitions(all_threshold);
adjusted_threshold = false(1, length(all_threshold));

for movement_nr = 1:length(mov_start)

    tmp_mov_start = mov_start(movement_nr) - round(window/2);
    tmp_mov_stop  = mov_stop(movement_nr) + round(window/2);

    if tmp_mov_start < 0 
        tmp_mov_start = 1;
    end

    if  tmp_mov_stop > length(paw_threshold)
        tmp_mov_stop = length(paw_threshold);
    end
    
    adjusted_threshold(tmp_mov_start:tmp_mov_stop) = true;
end


sData.analysis.movement_threshold = adjusted_threshold;
sData.analysis.paw_threshold      = paw_threshold;

%% Identify quiet wakefulness duration
log_vec                    = ~isnan( sData.analysis.paw_data);
threshold_vec_trimmed      = adjusted_threshold(log_vec);
qw_length_vec              = 1:sum(threshold_vec_trimmed < 1);

length(qw_length_vec)/sData.daqdata.faceCam_fps

sData.analysis.qw_duration = length(qw_length_vec)/sData.daqdata.faceCam_fps;

