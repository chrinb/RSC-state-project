function output = rand_nrem_time_comp(sData)

% Written by Christoffer Berge || Vervaeke lab

% Function that extracts random NREM windows to compare with sleep spindle
% data.

% Create different variables
signal        = mean(sData.imdata.roiSignals(2).newdff);
nr_of_seconds = 3; % nr of seconds before/after spindle onset
nr_of_frames  = (nr_of_seconds*31*2)+1; % nr of recording frames in spindle window
time          = linspace(-nr_of_seconds,nr_of_seconds,nr_of_frames); % time vector in seconds
frames        = sData.daqdata.frame_onset_reference_frame;
sessionID     = sData.sessionInfo.sessionID;
rand_times    = rand_nrem_times(sData);
rand_times_2ptime = frames(rand_times);

% find frames per second for imaging time vector
nrofFramesInRec    = length(sData.imdata.roiSignals(2).newdff);
RecordingLengthSec = length(sData.ephysdata.lfp)/2500;
framesPerSec       = nrofFramesInRec/RecordingLengthSec;

% loop over spindle onset times
for rand_nr = 1:length(rand_times_2ptime) % creates vector of nr of ripples during recording
    
    rand_window_start = rand_times_2ptime(rand_nr) - (nr_of_seconds*31); 
    rand_window_end   = rand_times_2ptime(rand_nr) + (nr_of_seconds*31);
    
    rand_activity(rand_nr, :) = signal(rand_window_start:rand_window_end); 

    baselineDFF = mean(rand_activity(rand_nr,1:31), 'omitnan');
        rand_activity(rand_nr,:) = rand_activity(rand_nr,:)-baselineDFF;

end

rand_activity_zscore = zscore(rand_activity,0,2);

output = {rand_activity, rand_activity_zscore};