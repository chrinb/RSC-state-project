function state_times = get_state_start_stop_times(sData)

% Written by Christoffer Berge | Vervake lab

%{
Compute start stop times for AW, QW, NREM, and REM episodes in recording in
imaging time
%}

% frames = sData.daqdata.frame_onset_reference_frame;

state_vectors_2p = get_state_logicals(sData);

for state_nr = 1:4

    temp_state_vector = state_vectors_2p{1, state_nr };

    if ~isempty(temp_state_vector)
        [state_start, state_stop] = findTransitions(temp_state_vector);
    
        state_times{state_nr, 1} = [state_start', state_stop'];
    else
        state_times{state_nr, 1} = [];
    end
end

% % Get 2P frame rate
% imaging_sampling_rate = find_imaging_framerate(sData);
% 
% % Get start/stop times for different states
% [active_start, active_stop] = findTransitions(sData.behavior.active_wakefulness); 
% active_times                = [active_start', active_stop']./imaging_sampling_rate;
% 
% [quiet_start, quiet_stop]   = findTransitions(sData.behavior.quiet_wakefulness);
% quiet_times                 = [quiet_start', quiet_stop']./imaging_sampling_rate;
% 
% [NREM_start, NREM_stop]     = findTransitions(sData.behavior.NREM_vector);
% NREM_times                  = round([NREM_start', NREM_stop']./imaging_sampling_rate);
% NREM_times2 =  [frames(NREM_start)', frames(NREM_stop)']
% [REM_start, REM_stop]       = findTransitions(sData.behavior.REM_vector);
% REM_times                   = [REM_start', REM_stop']./imaging_sampling_rate;
% 
% all_state_times = {active_times, quiet_times, NREM_times, REM_times};