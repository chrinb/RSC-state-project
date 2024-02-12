function [event_idx, xlabel_text, window_size] = select_event_type(sData, params, frames)

%{
Select event type to analyze and get associated indices and variables
%}
sleep_win_size = 5;

imaging_sampling_rate = find_imaging_framerate(sData);

signal_end = size(sData.imdata.roiSignals(2).newdff,2);

% Set threshold from recording start/end
threshold       = 3; % seconds
threshold_start = imaging_sampling_rate*threshold;
threshold_stop  = signal_end - (imaging_sampling_rate*threshold);


switch params.event_type

    case 'SWR'

        if strcmp(params.beh_state, 'awake')
            select_swr_idx = sData.ephysdata.absRipIdx; 
        elseif strcmp(params.beh_state, 'sleep')
            select_swr_idx = sort([sData.ephysdata.NREM_spindle_uncoupled_swr, sData.ephysdata.spindle_coupled_swr]);
        end

        % Get the indicies of user specified SWR types
        event_idx = get_swr_idx(params.swr_for_analysis ,sData,select_swr_idx, params);
        
        % Convert SWR time stamps from e-phys to 2P time
        event_idx = frames(round(event_idx));
        
        window_size = 1;
        xlabel_text = ' SWR peak ';
    case 'Spindle'
        event_idx = frames( sData.ephysdata2.NREMspindleStartEnd(:,1));  
        window_size = 1;
        xlabel_text = ' Spindle onset ';
    case 'SWA'
        [SOs, delta_waves] = mark_slow_wave(sData);
        % event_idx     = frames( SOs(:,3)); 
        event_idx  = frames( delta_waves(:,3));
        window_size = 0.33;
        xlabel_text = ' SWA trough ';
    case 'NREM_start'
        nrem_start_stop = nrem_sleep(sData);
        event_idx       = frames(nrem_start_stop(:,1));
        event_idx       = event_idx(event_idx > threshold_start);
        xlabel_text     = ' NREM start ';
        window_size     = sleep_win_size;
    case 'NREM_stop'
        nrem_start_stop = nrem_sleep(sData);
        event_idx       = frames(nrem_start_stop(:,2));
        event_idx       = event_idx(event_idx < threshold_stop);
        xlabel_text     = ' NREM stop ';
        window_size     = sleep_win_size;

     case 'REM_start'
        rem_start_stop = rem_sleep(sData);
        event_idx      = frames(rem_start_stop(:,1));
        event_idx       = event_idx(event_idx > threshold_start);
        xlabel_text    = ' REM start ';
        window_size    = sleep_win_size;

     case 'REM_stop'
        rem_start_stop = rem_sleep(sData);
        event_idx      = frames(rem_start_stop(:,2));   
        event_idx       = event_idx(event_idx < threshold_stop);
        xlabel_text    = ' REM stop ';
        window_size     = sleep_win_size;

end