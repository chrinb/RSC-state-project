function align_imaging_ephys_data(sData)

% Written by Christoffer Berge | Vervaeke lab

% Align 2P imaging data and ephys data. In the behavioral experiment data,
% the onset of ephys data is already aligned to imaging onset, but the
% ephys data outlasts the imaging data by a few seconds. Therefore the last
% imaging frame has to be found and used as an index to trim away excess
% ephys data. 

% Find last frame nr 
nr_of_frames = max(sData.daqdata.frame_onset_reference_frame);

log_idx = (sData.daqdata.frame_onset_reference_frame == nr_of_frames );

% Find index of last frame in ephys time
last_frame_ephys_idx = find(log_idx,1);

