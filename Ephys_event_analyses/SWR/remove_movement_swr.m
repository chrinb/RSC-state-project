function [RippleIdx,RunSpeedDuringRipples, locomotion_vector] = remove_movement_swr(sData, swr_idx)

% Written by Christoffer Berge | Vervaeke Lab

% Removes SWR occuring during or within 2s of locomotion bouts, defined as 
% either periods when instantaneous running speed is > 1cm/s. Periods of
% immobility lasting less than 2s also counts as locomotion. SWRs occurring
% within 3s of an locomotion bout is also removed. 

% THIS FUNCTION CAN BE (AND IS) USED FOR ALL RUNNING ANALYSES (awake, sleep etc.)FIX
% TITLE!! (14.02.22)

RunSpeedDuringRipples = zeros(length(swr_idx),10001); 

RippleIdx         = swr_idx; 
runSpeed          = sData.daqdata.runSpeed;
threshold         = 3; % cm/s
locomotion_vector = zeros(1,length(runSpeed));
srate             = 2500;
[~, locs] = findpeaks(runSpeed, 'MinPeakHeight', threshold);

% Check if recording contains running speed equal to or above 1 cm/sec. 
if max(sData.daqdata.runSpeed) >= 1
    
    % loop over nr of peaks in recording where running speed exceeds
    % threshold (= 1cm/sec). 
    for i = 1:length(locs)
        
        % calculate a window 2s before and after each peak exceeding
        % threshold and store each locomotion bout in a matrix
        locomotion_bout_start = locs(i) - (2500*2);
        locomotion_bout_end   = locs(i) + (2500*2);
        locomotion_bout(i,:)  = locomotion_bout_start:locomotion_bout_end;
    
        % check if a locomotion bout window begins before recording start.
        % If so delete indicies corresponding to negative values. Set each 
        % locomotion bout to have a value of 1. 
        if locomotion_bout_start <= 0
            tempVar                    = locomotion_bout_start:locomotion_bout_end;
            tempVar(tempVar <= 0)      = [];
            locomotion_vector(tempVar) = 1;
        else
            locomotion_vector(locomotion_bout(i,:)) = 1;
        end
    end


% calculate duration of immobility bouts
inverse_locomotion_vector = ~locomotion_vector;
a = strfind([0 inverse_locomotion_vector], [0 1])-1;  % find indicies of immobility bout onsets
b = strfind([inverse_locomotion_vector 0], [1 0]);    % find indicies of immobility bout onsets
if a(1) == 0
    a(1) = 1;
end

% loop over immobility bout onsets
for idx = 1:length(a)
    
    % Check of each immobility bout is smaller than 2 sec. If true, set
    % that bout = 1 (to exclude periods of immobility shorter than 2 sec)
    immobility_duration = b(idx)-a(idx);
    if immobility_duration < 2500*2
        locomotion_vector(a(idx):b(idx)) = 1;
    end
end

% calculate onset/offset of locomotion bouts
locomotion_bout_start = strfind([0 locomotion_vector], [0 1])-1;  % gives indices of beginning of groups
locomotion_bout_end   = strfind([locomotion_vector 0], [1 0]);    % gives indices of end of groups


%% Remove SWRs occurring within a locomotion bout

% loop over nr of locomotion bouts
for loc_bouts_n = 1:length(locomotion_bout_start)
    
    % loop over nr of SWRs in recording
    for swr_nr = 1:length(RippleIdx)
        % check if each SWR occurs within a locomotion bout + 3s before and
        % after. If true set that SWR index to 0
        if RippleIdx(swr_nr) > ( locomotion_bout_start(loc_bouts_n) - (srate*3) ) && ...
        RippleIdx(swr_nr) < ( locomotion_bout_end(loc_bouts_n) + (srate*3) )
        RippleIdx(swr_nr) = 0;
        end
    end
    
    % Remove any SWR indicies = 0
    RippleIdx(RippleIdx == 0) = [];

end
    
    
% if running speed does not exceed threshold, do nothing.     
else
    RippleIdx = swr_idx;
end
