function output = peri_event_activity(sData, params)

%{
Plot peri-event activity for used-specified events such as SWRs, sleep
spindles, and slow wave activity.

%}

% Get ROI signals
[signal, ~, ~, ~] = get_roi_signals_from_sData(sData, params );

switch params.cell_type
    case 'pc'
        signal = signal{1,:};
    case 'in'
        signal = signal{2,:};
    case 'axon'
        signal = signal{1,:};
end

% Define variables
nr_of_seconds  = 3;
nr_of_frames   = (nr_of_seconds*31*2)+1;
frames         = sData.daqdata.frame_onset_reference_frame;
time           = (-(31*nr_of_seconds):(31*nr_of_seconds))./31;
sessionID      = sData.sessionInfo.sessionID;

% switch params.event_type
[event_idx, xlabel_text] = select_event_type(sData, params, frames);

mean_event_activity = zeros(size(signal,1), nr_of_frames);
peri_event_activity = zeros(length(event_idx), nr_of_frames) ;
all_data            = zeros(size(signal,1), size(event_idx,2), size(time,2));     
%     case 'SWR'
% 
%         if strcmp(params.beh_state, 'awake')
%             select_swr_idx = sData.ephysdata.absRipIdx; 
%         elseif strcmp(params.beh_state, 'sleep')
%             select_swr_idx = sort([sData.ephysdata.NREM_spindle_uncoupled_swr, sData.ephysdata.spindle_coupled_swr]);
%         end
% 
%         % Get the indicies of user specified SWR types
%         event_idx = get_swr_idx(params.swr_for_analysis ,sData,select_swr_idx, params);
% 
%         % Convert SWR time stamps from e-phys to 2P time
%         event_idx = frames(round(event_idx));
% 
%         % Preallocate 
%         mean_event_activity = zeros(size(signal,1), nr_of_frames);
%         peri_event_activity = zeros(length(event_idx), nr_of_frames) ;
%         all_data            = zeros(size(signal,1), size(event_idx,2), size(time,2));                
% 
%     case 'Spindle'
% %         output = sleep_spindle_an(sData, params);
%     case 'SWA'
% %         output = so_sig_an(sData, params);
%     case 'NREM_start'
%         nrem_start_stop     = nrem_sleep(sData);
%         event_idx           = frames(nrem_start_stop(:,1));
%         mean_event_activity = zeros(size(signal,1), nr_of_frames);
%         peri_event_activity = zeros(length(event_idx), nr_of_frames) ;
%         all_data            = zeros(size(signal,1), size(event_idx,2), size(time,2));      
% 
%     case 'NREM_stop'
%         nrem_start_stop = nrem_sleep(sData);
%         event_idx       = frames(nrem_start_stop(:,2));
% end

%% Run main analysis
t = 1;

% Loop over nr of ROIs
for roi_nr = 1:size(signal,1) 
    
    % Extract signal from currrent ROI
    roi_signal = signal(roi_nr,:);

    % Get all peri-SWR activity windows from current ROI
    peri_event_activity = extract_avg_activity(event_idx, roi_signal, peri_event_activity, nr_of_seconds, params);

    % store all roi x time x SWR data
    all_data(t,:,:)   = peri_event_activity;
    
    % store mean SWR-activity
    mean_event_activity(roi_nr,:)  = mean(peri_event_activity, 'omitnan');

    t = t+1;
end


%% Save output as struct
% output = struct;

% if strcmp(params.exp_type, 'bulk')
% output.mean_event_activity   = mean_event_activity;
% output.params                = params;
% output.time                  = time;
% output.sessionID             = sessionID;
% output.event_idx             = event_idx;

output{1,1} = mean_event_activity;
output{2,1} = params;
output{3,1} = time;
output{4,1} = sessionID;
output{5,1} = event_idx;
output{6,1} = all_data;
output{7,1} = xlabel_text;
