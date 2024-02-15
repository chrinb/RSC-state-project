function sData = frame_cam_analysis(sData, tdms_path)


files = dir(fullfile(tdms_path, '*.tdms'));
 % Check if any files were found
if isempty(files)
    disp('No .tdms files found in the specified folder.');
else
    % Get the first file (assuming there's only one .tdms file in the folder)
    tdms_file_name = files(1).name;
end

% Load behavior-rig TDMS file to get frame times
sData = SA_loadTDMSdata(sData, tdms_path, tdms_file_name);  

if ~sum(sData.daqdataB.faceCam) == 0
        msgbox(['There is data in face cam field for session: ', sData.sessionInfo.sessionID]);
        sData.daqdata.faceCam = sData.daqdataB.faceCam;
end

sData.daqdata.frameSignal_behavior_rig = sData.daqdataB.frameSignal;
sData.daqdata.frameIndex_behavior_rig  = sData.daqdataB.frameIndex;
% sData                                     = rmfield(sData, 'daqdataB');

%% Load camera frame times 

files = dir(fullfile(tdms_path, '*.txt'));

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
            cam_frames_path = fullfile(tdms_path, cam_file_name);
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

% Find difference in 2P and face cam recording length
rec_diff = cam_samples_sec(end)-two_photon_frame_times(end);

corrected_t = cam_samples_sec - rec_diff;

% Find the camera frame index where the recording should start
[~, zero_idx] = min(abs(corrected_t));

msgbox(sprintf(['Difference in camera vs 2P frames in session ',  sData.sessionInfo.sessionID, ' is: ', num2str(numel(corrected_t(zero_idx:end))-numel(two_photon_frame_times)) ]));

figure, 
plot(corrected_t(zero_idx:end))
hold on
plot(two_photon_frame_times)
xlabel('Samples')
ylabel('Time (s)')

camera_frame_times = corrected_t(zero_idx:end);

%% Find closest match between the two vectors

% Initialize an array to store the closest matches
closest_matches = zeros(size(camera_frame_times));

% Loop over each timestamp in vector1
for i = 1:length(camera_frame_times)
    % Compute the absolute differences between the current timestamp in vector1 and all timestamps in vector2
    time_diffs = abs(camera_frame_times(i) - two_photon_frame_times);
    
    % Find the index of the minimum difference
    [~, min_index] = min(time_diffs);
    
    % Store the closest match
    closest_matches(i) = two_photon_frame_times(min_index);
end

%% Remove unnecessary fields from struct to save data
% if isfield(sData, 'behavior')
%     sData = rmfield(sData, 'behavior');
% end

if isfield(sData, 'ephysdata3')
    sData = rmfield(sData, 'ephysdata3');
end


sData.daqdata.cam_frame_times          = camera_frame_times;
sData.daqdata.camera_2p_closest_match  = closest_matches;
sData.daqdata.camera_start_frame       = zero_idx;

try
    msgbox(sprintf('Difference in two-photon frames in ephys and behavior rig: %d ', sum(diff(sData.daqdata.frameSignal  )==32) - sum(diff(sData.daqdata.frameSignal_behavior_rig  )==1) ));
catch
    msgbox(sprintf('Difference in two-photon frames in ephys and behavior rig: %d ', sum(diff(sData.daqdataB.frameSignal  )==1) - sum(diff(sData.daqdata.frameSignal_behavior_rig  )==1) ));
end