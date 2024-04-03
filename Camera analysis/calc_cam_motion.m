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
    
%% Find 2P laser onset and set that as first frame in face cam data
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

% Use the frame when the 2P light shines through the eye ("locs") as the
% first frame in the "cam_frame_times" variable. Then find closest
% corresponding 2P frame for rest of recording.

cam_samples_sec            = sData.daqdata.cam_frame_times;
cam_frame_time_at_2P_onset = cam_samples_sec(true_start_frame);

% adjusted_missing_frame_idx = sData.daqdata.missing_cam_frames > true_start_frame;
% adjusted_missing_frame_idx = sData.daqdata.missing_cam_frames(adjusted_missing_frame_idx)-true_start_frame;

corrected_t          = cam_samples_sec-cam_frame_time_at_2P_onset;
corr_cam_frame_times = corrected_t(true_start_frame:end);
% corr_cam_frame_times = corr_cam_frame_times(2:end);

% For one session, cam frame rate is ~60Hz. So for that session delete every second sample.
if sData.daqdata.faceCam_fps > 40
    corr_cam_frame_times1 = corr_cam_frame_times(1:2:end);
else 
    corr_cam_frame_times1 = corr_cam_frame_times;
end


%% Align the face camera frames with two-photon frames
n_2p_samples = size( sData.imdata.roiSignals(2).newdff, 2);

% [idx_motion, facecam_match_with_2p] = deal( zeros(n_2p_samples,1  ));
% 
% % Loop over face cam frame times and find the closest 2P frame
% for i = 1:length(corr_cam_frame_times)
% 
%     if isnan( corr_cam_frame_times(i))
%         % facecam_match_with_2p(i) = NaN;
%         idx_motion(i) = NaN;
%     else
% 
%     % Find the index in the 2P vector that is closest to the current face cam frame vector    
%     time_diffs         = abs(corr_cam_frame_times(i) - sData.daqdata.two_photon_frame_times);    
%     [minval, min_index]     = min(time_diffs);
%     % facecam_match_with_2p(i) = sData.daqdata.two_photon_frame_times(min_index);
%     % facecam_match_with_2p(i) =corr_cam_frame_times(min_index);
% 
%     % IMPORTANT: this variable contains indices pointing to the two-photon
%     % frame at the same time as the face cam. For example, if face cam
%     % frame 10000 = 323.0760, the index is 9992, corresponding to 2P data
%     % frame 9992 = 323.0638
%     idx_motion(i) = min_index;
%     end
% end
% idx_motion(idx_motion==0) = NaN;

%% Test 
n_2p_samples  = size( sData.imdata.roiSignals(2).newdff, 2);

if ~sData.daqdata.two_photon_frame_times(1) == 0
    sData.daqdata.two_photon_frame_times = [0; sData.daqdata.two_photon_frame_times];
end

idx_motion = zeros(1, n_2p_samples);
%% Loop
for i = 1:n_2p_samples

    if i <= numel( sData.daqdata.two_photon_frame_times)
         
        % sample_time_diff = abs(sData.daqdata.two_photon_frame_times(i) - corr_cam_frame_times1(i));    
        % % Find the index in the 2P vector that is closest to the current face cam frame vector    
        % all_time_diffs           = abs(sData.daqdata.two_photon_frame_times(i) - corr_cam_frame_times1);    
        % [min_val, min_index] = min(all_time_diffs);
        % 
        % While loop index is smaller than nr of frames in face cam, check
        % if corresponding face cam frame is NaN
        if i < numel(corr_cam_frame_times1)
            
            % Compare 
            % sample_time_diff     = abs( sData.daqdata.two_photon_frame_times(i) - corr_cam_frame_times1(i) );    
            all_time_diffs       = abs( sData.daqdata.two_photon_frame_times(i) - corr_cam_frame_times1 );    
            [min_val, min_index] = min(all_time_diffs);


            if isnan( corr_cam_frame_times1(i))
                idx_motion(i) = NaN;
                
            else

                % if abs( sample_time_diff - min_val ) < 0.01
                %     idx_motion(i) = i;
                % else
                    idx_motion(i) = min_index;
                % end
            end
        else
            idx_motion(i) = NaN;
        end
    end
end

% IMPORTANT: this variable contains indices pointing to the face cam sample 
% closest in time to the corresponding 2P frame. For example, if 2P sample
% frame 10000 = 323.3196, inputting 10000 into idx_motion = 10009, showing that 
% the nearest face cam frame occurs at index 10009 (corr_cam_frame_times(
% idx_motion(10000)) = 323.3320
idx_motion(idx_motion==0)= NaN;

% % Quick fix
% inter_sample_times = diff(idx_motion);
% inter_sample_times = [1 inter_sample_times];
% for i = 1:n_2p_samples
% 
%     % When the index goes to the next frame, check whether there is
%     % meaningful difference
%     if ~( inter_sample_times(i) == 1)
% 
%         if ~(idx_motion(i) == i)
%             pre_frame_idx_val = abs( sData.daqdata.two_photon_frame_times(i) - corr_cam_frame_times1( idx_motion(i) - 1) );
%             cur_frame_idx_val  = abs( sData.daqdata.two_photon_frame_times(i) - corr_cam_frame_times1( idx_motion(i) ) );
%             post_frame_idx_val = abs( sData.daqdata.two_photon_frame_times(i) - corr_cam_frame_times1( idx_motion(i) + 1) );
%             % potential_val = abs( sData.daqdata.two_photon_frame_times(i)-corr_cam_frame_times1(i) );
%             % true_val      = abs( sData.daqdata.two_photon_frame_times(i)-corr_cam_frame_times1( idx_motion(i)) );
% 
%             current_to_pre   = abs( cur_frame_idx_val - pre_frame_idx_val);
%             current_to_post  = abs( cur_frame_idx_val - post_frame_idx_val);
% 
%             % Find which frame idx is closer, pre or post
%             combined_idx = [current_to_pre, current_to_post];
%             [~, min_idx] = min( combined_idx);
% 
%             if min_idx == 1
%                 idx_motion(i) = idx_motion(i)-1 ;
%             elseif min_idx == 2
%                 idx_motion(i) = idx_motion(i)+1 ;
%             end
%         end
%     end
% end

sData.daqdata.aligned_2p_face_cam_idx = idx_motion;
sData.daqdata.cam_frame_times_corr    = corr_cam_frame_times1;
% rmfield(sData.daqdata,"facecam_2p_idx");

%% Align motion vectors from face camp with 2p data

% First insert NaNs into motion vector based on face cam frame time vector
nan_indices = isnan(corr_cam_frame_times1);

paw_data_adjusted     = double( mean( abs( data.motSVD_1(true_start_frame:end, 1:3) ),2) );
face_data_adjusted    = double( mean( abs( data.motSVD_2(true_start_frame:end, 1:3) ),2) );
whisker_data_adjusted = double( mean( abs( data.motSVD_3(true_start_frame:end, 1:3) ),2) );

if sData.daqdata.faceCam_fps > 40
    paw_data_adjusted = paw_data_adjusted(1:2:end);
    face_data_adjusted = face_data_adjusted(1:2:end);
    whisker_data_adjusted = whisker_data_adjusted(1:2:end);
end
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
        whisker_motion_2p_match(i) = whisker_data_adjusted( tmp_idx);
    end

end

% Store the mean of the absolute values (bc components can have opposite
% signs) of the top 3 components
sData.analysis.paw_data     = paw_motion_2p_match;
sData.analysis.face_data    = face_motion_2p_match;
sData.analysis.whisker_data = whisker_motion_2p_match;

%% Plot mean of the absolute values of top 3 SVD components

time_vec = (0:length(sData.analysis.paw_data)-1)/31;

figure(2), clf 
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

if sData.daqdata.faceCam_fps > 40
    window = 31;
else
    window = sData.daqdata.faceCam_fps;
end

if sData.daqdata.faceCam_fps > 40

    face_threshold    = ( movstd(sData.analysis.face_data, window) - mean(movstd(sData.analysis.face_data, window),'omitnan') ) > 1;
    paw_threshold     = movstd(sData.analysis.paw_data, window) - mean(movstd(sData.analysis.paw_data, window),'omitnan') > 1;
    whisker_threshold = movstd(sData.analysis.whisker_data, window) - mean(movstd(sData.analysis.whisker_data, window),'omitnan') > 1;
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

