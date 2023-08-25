function [awakeSWRidx, NREMspindleUncoupledSWRidx, NREMspindleCoupledSWRidx ] ...
    = ripRunAn(varargin)

% Written by Christoffer Berge || Vervaeke lab

% Removes SWR occuring during or within 2s of locomotion bouts, defined as 
% either periods when instantaneous running speed is > 1cm/s. Periods of
% immobility lasting less than 2s also counts as locomotion. SWRs occurring
% within 3s of an locomotion bout is also removed. 

sData = varargin{1, 1};

unclassified_swr_str                  = (varargin{1,3});
NREM_spindle_uncoupled_swr_str = (varargin{1,4});
spindle_coupled_swr_str        = (varargin{1,5});

removeTemporallyCloseRip  = varargin{1, 2};

frames = sData.daqdata.frame_onset_reference_frame;
runSpeed = sData.daqdata.runSpeed;
threshold = 1; % cm/s
locomotion_vector = zeros(1,length(runSpeed));

[~, locs] = findpeaks(runSpeed, 'MinPeakHeight', threshold);

% check of running speed at any time point is equal to or exceeds 1 cm/sec.
% If not, skip analysis. 
if max(sData.daqdata.runSpeed) >= 1

    % Loop over peaks in running speed exceeding threshold (1cm/sec) and
    % extract segment 2 sec before and after that peak. 
    for i = 1:length(locs)
    locomotion_bout_start = locs(i) - (2500*2);
    locomotion_bout_end = locs(i) + (2500*2);
    locomotion_bout(i,:) = locomotion_bout_start:locomotion_bout_end;
    
    % check if extracted locomotion bout has a negative index (i.e. starts before
    % recording start) delete negative indicies. Store all locomotion bouts
    % in a new matrix (locomotion vector). 
    if locomotion_bout_start <= 0
    tempVar = locomotion_bout_start:locomotion_bout_end;
    tempVar(tempVar <= 0) = [];
    locomotion_vector(tempVar) = 1;
    else
    locomotion_vector(locomotion_bout(i,:)) = 1;
    end
    end
    
% calculate duration of immobility bouts
inverse_locomotion_vector = ~locomotion_vector;
a = strfind([0 inverse_locomotion_vector], [0 1])-1;  %gives indices of beginning of groups
b = strfind([inverse_locomotion_vector 0], [1 0]);    %gives indices of end of groups

for idx = 1:length(a)
    immobility_duration = b(idx)-a(idx);
    if length(immobility_duration) > 2500*2
        locomotion_vector(a(idx):b(idx)) = 1;
    end
end

% calculate onset/offset of locomotion bouts
locomotion_bout_start = strfind([0 locomotion_vector], [0 1])-1;  %gives indices of beginning of groups
locomotion_bout_end = strfind([locomotion_vector 0], [1 0]);    %gives indices of end of groups

if removeTemporallyCloseRip == 0
awakeSWRidx = sData.ephysdata.(unclassified_swr_str);
NREMspindleUncoupledSWRidx = sData.ephysdata.(NREM_spindle_uncoupled_swr_str);
NREMspindleCoupledSWRidx = sData.ephysdata.(spindle_coupled_swr_str);

else
    [awakeSWRidx, NREMspindleUncoupledSWRidx, NREMspindleCoupledSWRidx] ...
        = removeCloseRip(sData);
end


%% check if any SWR-type occurs during a locomotion bout (and within 3s of one)

% awake SWR
if removeTemporallyCloseRip == 0
    for y = 1:length(locomotion_bout_start)
        for l = 1:length(awakeSWRidx)
            if awakeSWRidx(l) > ( locomotion_bout_start(y) - (2500*3) ) && ...
            awakeSWRidx(l) < ( locomotion_bout_end(y) + (2500*3) )
            awakeSWRidx(l) = 0;
            end
        end
        awakeSWRidx(awakeSWRidx == 0) = [];
        % NREM spindle-uncoupled SWR
        for u = 1:length(NREMspindleUncoupledSWRidx)
            if NREMspindleUncoupledSWRidx(u) > ( locomotion_bout_start(y) - (2500*3) ) && ...
               NREMspindleUncoupledSWRidx(u) < ( locomotion_bout_end(y) + (2500*3) )
               NREMspindleUncoupledSWRidx(u) = 0;
            end
        end
        NREMspindleUncoupledSWRidx(NREMspindleUncoupledSWRidx == 0) = [];
        
        % NREM spindle-coupled SWR
        for p = 1:length(NREMspindleCoupledSWRidx)
            if NREMspindleCoupledSWRidx(p) > ( locomotion_bout_end(y) - (2500*3) ) && ...
               NREMspindleCoupledSWRidx(p) < ( locomotion_bout_end(y) + (2500*3) )
               NREMspindleCoupledSWRidx(p) = 0;
            end
        end
        NREMspindleCoupledSWRidx(NREMspindleCoupledSWRidx == 0) = [];
        
    end
    

else
    for y = 1:length(locomotion_bout_start)
        for l = 1:length(awakeSWRidx)
            if awakeSWRidx(l) > ( locomotion_bout_start(y) - (2500*3) ) && ...
            awakeSWRidx(l) < ( locomotion_bout_end(y) + (2500*3) )
            awakeSWRidx(l) = 0;
            end
        end
        awakeSWRidx(awakeSWRidx == 0) = [];
        
        % NREM spindle-uncoupled SWR
        for u = 1:length(NREMspindleUncoupledSWRidx)
            if NREMspindleUncoupledSWRidx(u) > ( locomotion_bout_start(y) - (2500*3) ) && ...
               NREMspindleUncoupledSWRidx(u) < ( locomotion_bout_end(y) + (2500*3) )
               NREMspindleUncoupledSWRidx(u) = 0;
            end
        end
        NREMspindleUncoupledSWRidx(NREMspindleUncoupledSWRidx == 0) = [];
    % NREM spindle-coupled SWR
        for p = 1:length(NREMspindleCoupledSWRidx)
            if NREMspindleCoupledSWRidx(p) > ( locomotion_bout_end(y) - (2500*3) ) && ...
               NREMspindleCoupledSWRidx(p) < ( locomotion_bout_end(y) + (2500*3) )
               NREMspindleCoupledSWRidx(p) = 0;
            end
        end
        NREMspindleCoupledSWRidx(NREMspindleCoupledSWRidx == 0) = [];
    end
end

elseif removeTemporallyCloseRip == 0
    awakeSWRidx = sData.ephysdata.awake_swr;
    NREMspindleUncoupledSWRidx = sData.ephysdata.NREM_spindle_uncoupled_swr;
    NREMspindleCoupledSWRidx = sData.ephysdata.spindle_coupled_swr;
else
    [awakeSWRidx, NREMspindleUncoupledSWRidx, NREMspindleCoupledSWRidx] = removeCloseRip(sData);
end
        
