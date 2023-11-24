function [event_idx, xlabel_text] = select_event_type(sData, params, frames)

%{
Select event type to analyze and get associated indices and variables
%}


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
             
        xlabel_text = ' SWR peak ';
    case 'Spindle'
%         output = sleep_spindle_an(sData, params);
        xlabel_text = ' Spindle onset ';
    case 'SWA'
%         output = so_sig_an(sData, params);
    case 'NREM_start'
        nrem_start_stop = nrem_sleep(sData);
        event_idx       = frames(nrem_start_stop(:,1));
        xlabel_text     = ' NREM start ';
   
    case 'NREM_stop'
        nrem_start_stop = nrem_sleep(sData);
        event_idx       = frames(nrem_start_stop(:,2));
        xlabel_text     = ' NREM stop ';
     case 'REM_start'
        rem_start_stop = rem_sleep(sData);
        event_idx      = frames(rem_start_stop(:,1));
        xlabel_text    = ' REM start ';

     case 'REM_stop'
        rem_start_stop = rem_sleep(sData);
        event_idx      = frames(rem_start_stop(:,2));   
        xlabel_text    = ' REM stop ';
end