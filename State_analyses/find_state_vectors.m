function sData = find_state_vectors(sData)

% Written Christoffer Berge | Vervaeke lab

% Create logical vectors marking various behavioral states in recording
% (active/quiet wakefulness, NREM sleep, REM sleep). Active wakefulness is
% defined as epochs in the rec. where running speed > 3 cm/sec. In awake
% recs., quiet wakefulness is everything else. In sleep recordings,
% quiet wakefulness is everything else which *doesn't* overlap with NREM or
% REM bouts.

% TO DO: This code might need some polishing to be certain about the quiet
% wakefulness state (e.g., what about whisker movements? Shrugs and jerks?)
srate             = 2500;
runSpeed          = sData.daqdata.runSpeed;
runSpeed          = smoothdata(sData.daqdata.runSpeed, 'gaussian', srate);

threshold         = 1; % cm/s
time              = linspace(1, size(runSpeed,1), size(runSpeed,1) )/srate;

active_wakefulness = [];
% Check if recording contains running speed equal to or above 1 cm/sec. 
if max(runSpeed) >= threshold
    
        % Find peaks above running speed threshold
        logical_pos_peaks       = runSpeed > threshold;
       
        % Correct for running at imaging start
        if logical_pos_peaks(1) == true
            logical_pos_peaks(1) = false;
        end
        
        [eventStartIdx, eventStopIdx ]  = findTransitions( logical_pos_peaks );
        
        % Loop over run events
        for i = 2:size(eventStartIdx, 1)

            % If inter-run event time is > 2.5 sec, merge run events
            if ( eventStartIdx(i) - eventStopIdx(i-1) ) < srate*2.5
                logical_pos_peaks(eventStopIdx(i-1):eventStartIdx(i)) = ones;
            end
        end

        % Remove locomotion bouts < 0.5 sec
         [eventStartIdx, eventStopIdx ]  = findTransitions( logical_pos_peaks );
         event_start_stop_diff           = eventStopIdx-eventStartIdx;
         log_idx = event_start_stop_diff > srate/2;

         for i = 1:size(log_idx,1)
             if ~any( log_idx(i))
                 logical_pos_peaks(eventStartIdx(i):eventStopIdx(i)) = zeros;
             end
         end


 
        active_wakefulness = logical_pos_peaks';
else
    active_wakefulness = false( 1, length(time));
end

% If sleep recording, distinguish quiet wakefulness, NREM sleep, & REM
% sleep states (if present in rec.). 
if isfield(sData, 'episodes')
    

    NREM_log_vec = zeros(1, size(runSpeed,1));
    REM_log_vec  = zeros(1, size(runSpeed,1));

    % Create logical vector of NREM epochs
    if sum( sData.episodes.state == 'NREM') > 0
        NREM_episodes = nrem_sleep(sData);

        % Loop over NREM episodes
        for i = 1:size(NREM_episodes,1)
            NREM_log_vec( NREM_episodes(i,1):NREM_episodes(i,2)) = ones;
        end
    end
    NREM_log_vec = logical(NREM_log_vec);
    
    % Create logical vector of REM epochs
    if sum( sData.episodes.state == 'REM') > 0
        REM_episodes = rem_sleep(sData);
    
        % Loop over REM episodes
        for i = 1:size(REM_episodes,1)
            REM_log_vec( REM_episodes(i,1):REM_episodes(i,2)) = ones;
        end
    end
    REM_log_vec = logical(REM_log_vec);

    % Concatenate active wakefulness, NREM, and REM vectors to create one
    % logical vector indicating all three states
    cat_states = logical( sum( [active_wakefulness; REM_log_vec; NREM_log_vec], 1));
    
    % Define quiet wakefulness as all remaining time epochs in rec. not active
    % wakefulness, NREM, or REM sleep. 
    quiet_wakefulness = not(cat_states);

    % Correct for ones at imaging start
    if quiet_wakefulness(1) == true
        quiet_wakefulness(1) = false;
    end
    
    sData.behavior.active_wakefulness = active_wakefulness;
    sData.behavior.quiet_wakefulness  = quiet_wakefulness;
    sData.behavior.NREM_vector        = NREM_log_vec;
    sData.behavior.REM_vector         = REM_log_vec;
else
    % If not sleep rec, define quiet wakefulness as everything else besides
    % active wakefulness. 
    quiet_wakefulness = not(active_wakefulness);

    % Correct for ones at imaging start
    if quiet_wakefulness(1) == true
        quiet_wakefulness(1) = false;
    end
    NREM_log_vec = zeros(1, length(runSpeed));
    REM_log_vec = zeros(1, length(runSpeed));

    sData.behavior.active_wakefulness = active_wakefulness;
    sData.behavior.quiet_wakefulness  = quiet_wakefulness;
    sData.behavior.NREM_vector        = NREM_log_vec;
    sData.behavior.REM_vector         = REM_log_vec;

end
