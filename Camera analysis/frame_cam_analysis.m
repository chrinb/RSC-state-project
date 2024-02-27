function sData = frame_cam_analysis(sessionObject)


sData = sessionObject.Data.sData;

% Get path to TDMS file
root_path     = sessionObject.DataLocation(2).RootPath;
folder_path   = sessionObject.DataLocation(2).Subfolders;
mouse_folder  = [root_path,'\', folder_path];
facemap_folder_path     = [mouse_folder, '\FaceCam'];

% files = dir(fullfile(facemap_folder_path, '*.tdms'));
%  % Check if any files were found
% if isempty(files)
%     disp('No .tdms files found in the specified folder.');
% else
%     % Get the first file (assuming there's only one .tdms file in the folder)
%     tdms_file_name = files(1).name;
% end


files = dir(fullfile(facemap_folder_path));

 % Check if any files were found
for i = 1:size(files,1)

    if contains(files(i).name, '.ini' )
        config_file = files(i).name;
    end

    if contains(files(i).name, '.tdms' )
        tdms_file_name = files(i).name;
    end
end

% Define parameters
config_data_path = [facemap_folder_path,'\' config_file];
ini_contents     = textread(config_data_path, '%s', 'delimiter', '\n', 'whitespace', '');
cam_fps_char     = ini_contents{3};
numbersAsNumeric = str2double(cam_fps_char);
fps_num          = regexp(cam_fps_char, '\d+', 'match');

fps                       = str2double( [fps_num{1} '.' fps_num{2}]);
sData.daqdata.faceCam_fps = fps;


% Load behavior-rig TDMS file to get frame times
sData = SA_loadTDMSdata(sData, facemap_folder_path, tdms_file_name);  

if ~sum(sData.daqdataB.faceCam) == 0
        msgbox(['There is data in face cam field for session: ', sData.sessionInfo.sessionID]);
        sData.daqdata.faceCam = sData.daqdataB.faceCam;
end

sData.daqdata.frameSignal_behavior_rig = sData.daqdataB.frameSignal;
sData.daqdata.frameIndex_behavior_rig  = sData.daqdataB.frameIndex;
% sData                                     = rmfield(sData, 'daqdataB');

%% Load camera frame times 

files = dir(fullfile(facemap_folder_path, '*.txt'));

% % Initialize variables to store file name and path
% cam_frames_path = '';

% Check if any files were found
if isempty(files)
    disp('No .txt files found in the specified folder.');
else
    % Loop through each file in the folder
    for i = 1:length(files)
        % Check if the file name contains "frametimes"
        if contains(files(i).name, 'frametimes')
            % Get the file name and path
            cam_file_name = files(i).name;
            cam_frames_path = fullfile(facemap_folder_path, cam_file_name);
            break; % Stop searching once the first matching file is found
        end
    end
end

cam_time_stamps = vr.importPupilFrameTimes(cam_frames_path);
%% Get imaging frame rate
srate = find_imaging_framerate(sData);

%% Find closest match between camera and two-photon frame times

% Convert from datetime object to 
t_pos = posixtime(cam_time_stamps);

% Could try to "make" the camera frame time begin at the first frame when
% 2P light comes on, but it actually makes no difference for the zero index
% calculated below. Here, the 2P light came on in frame 53.
% t = t_pos-t_pos(53);

% Remove the posixtime 
cam_samples_sec = t_pos-t_pos(1);

% Find time stamps of 2P frames in seconds
frame_idx        = sData.daqdata.frameIndex_behavior_rig;
frame_idx_length = 1:length(frame_idx);

two_photon_frame_times = frame_idx_length*(1/srate);
two_photon_frame_times = two_photon_frame_times';

% msgbox(sprintf('Number of camera frames and 2P frames are: %d vs %d',numel(cam_samples_sec), numel(two_photon_frame_times) ));

%Find difference in 2P and face cam recording length
rec_diff = cam_samples_sec(end)-two_photon_frame_times(end);

corrected_t = cam_samples_sec - rec_diff;

%Find the camera frame index where the recording should start
[~, zero_idx] = min(abs(corrected_t));
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

% missing_frame_idx = zeros(1, numel(val));
% mf = [];
% for i = 1:numel(val)
%     % if val(i) > med_val 
% 
%         missing_frame_idx(i) = val_idx(i);
% 
%         k = l == val(i);
%         mf = [mf, l2(k)];
% 
%     % Also look for frame interval smaller than frame rate (if interval is
%     % large, the next interval can be shorter than the normal ~32ms interval) 
%     % elseif val(i) < med_val
%         missing_frame_idx(i) = val_idx(i);
%          k = l == val(i);
%         mf = [mf, l2(k)];
% 
%     % end
%     clear k
% end
% missing_frame_idx(missing_frame_idx == 0) = [];

% Check of camera FPS from ini file differs a lot from empirical fps. If
% so, use empirical
% mod_val = mode(val);

if abs(fps - med_val*1000) > 1.5
    % fps_corrected = med_val*1000;
    fps_corrected = 30.9;
    sData.daqdata.faceCam_fps = fps_corrected;
    % msgbox( [num2str(abs(fps - med_val*1000)),' difference in Hz between ini file data and empirical data. Corrected to ', num2str(fps_corrected)])
    msgbox( ['ini file data is ', num2str(fps) 'Hz, now hardcoded to 30.9 Hz'])

end
%%
figure(1), clf
sgtitle(sData.sessionInfo.sessionID, 'interpreter', 'none')
h(1) = subplot(211);
plot(cam_samples_sec)

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

sData.daqdata.missing_cam_frames = missing_frame_idx;

% 
% msgbox(sprintf(['Difference in camera vs 2P frames in session ',  sData.sessionInfo.sessionID, ' is: ', num2str(numel(corrected_t(zero_idx:end))-numel(two_photon_frame_times)) ]));
% 
% figure, 
% plot(corrected_t(zero_idx:end))
% hold on
% plot(two_photon_frame_times)
% xlabel('Samples')
% ylabel('Time (s)')
% 
% camera_frame_times = corrected_t(zero_idx:end);

%% Find closest match between the two vectors

% Initialize an array to store the closest matches
% closest_matches = zeros(size(camera_frame_times));
% 
% % Loop over each timestamp in vector1
% for i = 1:length(camera_frame_times)
%     % Compute the absolute differences between the current timestamp in vector1 and all timestamps in vector2
%     time_diffs = abs(camera_frame_times(i) - two_photon_frame_times);
% 
%     % Find the index of the minimum difference
%     [~, min_index] = min(time_diffs);
% 
%     % Store the closest match
%     closest_matches(i) = two_photon_frame_times(min_index);
% end

%% Remove unnecessary fields from struct to save data
% if isfield(sData, 'behavior')
%     sData = rmfield(sData, 'behavior');
% end

if isfield(sData, 'ephysdata3')
    sData = rmfield(sData, 'ephysdata3');
end

sData.daqdata.cam_frame_times        = cam_samples_sec;
sData.daqdata.two_photon_frame_times = two_photon_frame_times;
sData.daqdata.camera_start_frame     = zero_idx;
try
    msgbox(sprintf('Difference in two-photon frames in ephys and behavior rig: %d ', sum(diff(sData.daqdata.frameSignal  )==32) - sum(diff(sData.daqdata.frameSignal_behavior_rig  )==1) ));
catch
    msgbox(sprintf('Difference in two-photon frames in ephys and behavior rig: %d ', sum(diff(sData.daqdataB.frameSignal  )==1) - sum(diff(sData.daqdata.frameSignal_behavior_rig  )==1) ));
end