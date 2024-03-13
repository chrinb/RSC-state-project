function convert_raw_to_avi(sessionObjects)

%% Load raw file
root_path     = sessionObjects.DataLocation(2).RootPath;
folder_path   = sessionObjects.DataLocation(2).Subfolders;
mouse_folder  = [root_path,'\', folder_path];
facemap_folder_path     = [mouse_folder, '\FaceCam'];
% files = dir(fullfile(facemap_folder_path, '*.raw'));
files = dir(fullfile(facemap_folder_path));

 % Check if any files were found
for i = 1:size(files,1)

    if contains(files(i).name, '.ini' )
        config_file = files(i).name;
    end

    if contains(files(i).name, '.raw' )
        raw_file_name = files(i).name;
    end
end

% Load raw file
raw_data_path = [facemap_folder_path,'\' raw_file_name];
data          = fopen(raw_data_path, 'r');
raw_data      = fread(data, inf, 'uint8');
fclose(data);

% Split the file name and extension
[~, recording_name, ~] = fileparts(raw_file_name);

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
t_pos           = posixtime(cam_time_stamps);
cam_samples_sec = t_pos-t_pos(1);
cam_samples_sec = cam_samples_sec(2:end);
med_val = median( diff(cam_samples_sec));

% Define parameters
config_data_path = [facemap_folder_path,'\' config_file];
ini_contents     = textread(config_data_path, '%s', 'delimiter', '\n', 'whitespace', '');
cam_fps_char     = ini_contents{3};
fps_num          = regexp(cam_fps_char, '\d+', 'match');

fps = str2double( [fps_num{1} '.' fps_num{2}]);

if abs(fps - med_val*1000) > 1.5
    fps_corrected = 30.9;
    msgbox( ['ini file data is ', num2str(fps) 'Hz, now hardcoded to 30.9 Hz'])
else 
    fps_corrected = fps;
end

output_AVI_Name = [recording_name, '.avi'];  % Output .avi file name
width = 250;  % Width of the image
height = 250; % Height of the image
% fps = 31;     % Frames per second for the AVI file


% Reshape the raw data to image dimensions
num_frames = numel(raw_data) / (width * height);
video_data = reshape(raw_data, [width, height, num_frames]);

% Construct the full path to the output AVI file
outputAVIPath = fullfile(facemap_folder_path, output_AVI_Name);

% Create VideoWriter object with specified output path
% v = VideoWriter(outputAVIPath, 'Grayscale AVI');
v = VideoWriter(outputAVIPath, 'MPEG-4');
v.FrameRate = fps_corrected;

% Open the VideoWriter object
open(v);

% Write each frame to the video file
for frameIdx = 1:num_frames
    frame = video_data(:, :, frameIdx);
    frame_normalized = frame / 255;
    writeVideo(v, frame_normalized);
end

% Close the VideoWriter object
close(v);

disp('Video conversion completed.');
