function output = get_state_logicals(varargin)

% Written by Christoffer Berge | Vervaeke lab

% Function that converts behavioral state vectors (1 for state on, 0 for
% not) for different states (active, quiet wakefulness, NREM, REM) from
% ephys time to two-photon time.

sData = varargin{1,1};

% User can input a 4x1 cell array as second input argument. The function
% will then use those data as state vectors. These data are state vectors
% where states happening close to the beginning or end have been removed. 
if nargin > 1
    all_state_vectors_adjusted = varargin{1,2};
    
    sData.behavior.active_wakefulness = all_state_vectors_adjusted{1, 1};  
    sData.behavior.quiet_wakefulness  = all_state_vectors_adjusted{2, 1};  
    sData.behavior.NREM_vector        = all_state_vectors_adjusted{3, 1};  
    sData.behavior.REM_vector         = all_state_vectors_adjusted{4, 1};  
end

frames_in_rec = size(sData.imdata.roiSignals(2).newdff, 2);
frames        = sData.daqdata.frame_onset_reference_frame;

% Check if multi-plane (piezo) recording. If so, adjustments have to be
% made to account for the fact that 4 frames are taken in one piezo
% cycle, and event starts/stops have to be assigned to a particular cycle. 
% A typical piezo session with 4 planes at 62hz for 10 min > 30000
% frames. OR, if recording is 512 lines, then FPS is ~7Hz and total is less
% than 30000. Then check if there are ROI indices at other planes than the
% first 1.
if frames(end) > 30000 || sum(unique(sData.imdata.plane_indices) > 1) > 0

    % First create a vector going from 1 to max nr of frames in rec
    frame_counter = unique(sData.daqdata.frame_onset_reference_frame);

    % Then create a new vector selecting only every 4th frame. The
    % intervening (missing) frames between two frames, e.g. 13 and 17,  "belongs"
    % to the first frame (13), because frame 14-16 are taken before piezo jumps back to
    % original plane.
    first_frame_per_piezo_cycle = (1:4:frame_counter(end));
    
    % Divide by 4 (nr of planes) to get frames per plane
%     frame_counter_per_plane = ceil( first_frame_per_piezo_cycle./4);
end

%% Find epochs 

output = cell(1,4);
% Loop over states
for state_nr = 1:4

    switch state_nr
        case 1
            state_vec_name = 'NREM_vector';
        case 2
            state_vec_name = 'REM_vector';
        case 3
            state_vec_name = 'quiet_wakefulness';
        case 4
            state_vec_name = 'active_wakefulness';
    end


    if isfield(sData.behavior, state_vec_name) && ~isempty(sData.behavior.(state_vec_name)) && ~(sum(sData.behavior.(state_vec_name)) == 0)
        % Load NREM logical
        state_vector                    = sData.behavior.(state_vec_name);
        [eventStartIdx, eventStopIdx ]  = findTransitions( state_vector );
        
        % convert to 2P time
        eventStartIdx = frames(eventStartIdx);
        eventStopIdx  = frames(eventStopIdx);
        
        % Correction for multi-plane recordings
        if frames(end) > 30000 || sum(unique(sData.imdata.plane_indices) > 1) > 0
           
            % For each event start/stop frame, find the closest match to the first 
            % frame of each piezo cycle
            index_start = dsearchn( first_frame_per_piezo_cycle' , eventStartIdx' );
            index_stop  = dsearchn( first_frame_per_piezo_cycle' , eventStopIdx' );
            
            % If start/stop time is smaller than the closest matching frame, that
            % means that the start/stop time belongs to the previous frame/index (because
            % if the closest match is larger, then that frame is the start of a new
            % piezo cycle)
            start_check = bsxfun(@lt, eventStartIdx, first_frame_per_piezo_cycle(index_start));
            stop_check  = bsxfun(@lt, eventStopIdx, first_frame_per_piezo_cycle(index_stop));
    
    %        corr_frame_index_start = bsxfun(@minus, index_start(start_check), 1)
    %        corr_frame_index_stop  
    
           for i = 1:length(index_start)
               if start_check(i)
                   index_start(i) = index_start(i)-1;
               end
    
               if stop_check(i)
                   index_stop(i) = index_stop(i)-1;
               end
           end
        
            eventStartIdx = index_start;
            eventStopIdx  = index_stop;
        end
    
        % Sometimes ephys is a little bit longer than imaging data. Check if
        % that's the case and correct in necessary.
        if eventStopIdx(end) > frames_in_rec
            eventStopIdx(end) = frames_in_rec;
        end

        % Create new logical (now in 2P time)
        state_vec_2P = false(1, frames_in_rec);
        for i = 1:size(eventStartIdx,2)
            state_vec_2P( eventStartIdx(i):eventStopIdx(i)) = true;
        end
    
        output{1,state_nr} = state_vec_2P;
    end
end
