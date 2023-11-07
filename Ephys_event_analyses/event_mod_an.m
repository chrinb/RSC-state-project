function sData = event_mod_an(sData, params)

% Written by Christoffer Berge | Vervaeke Lab

% Function that looks for event-modulated ROIs by creating a Peri-event 
% time histogram (PETH) for each ROI and comparing the mean squared error 
% (MSE) of that PETH to a distribution of MSE's created from random,
% cirular shuffeling of the original PETH. 

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
        
        % Preallocate 
        mean_event_activity = zeros(size(signal,1), nr_of_frames);
        peri_event_activity = zeros(length(event_idx), nr_of_frames) ;
        all_data            = zeros(size(signal,1), size(event_idx,2), size(time,2));                

        time_win = (62:108); % -1s to +500ms after SWR peak
    case 'Spindle'
        spindle_select     = [];
        spin_start_end_str = strcat('NREMspindleStartEnd', spindle_select);
        
        % convert sleep spindle time stamps from e-phys to 2P time
        event_idx = sData.ephysdata2.(spin_start_end_str)(:,1);
    
        % Get time points for sleep spindle onsets
        event_idx = get_spindle_idx(params , event_idx);
        
        time_win             = (93:124); % 0s to +1s after spindle onset
    case 'SWA'
%         output = so_sig_an(sData, params);
end

%% Define variables for analysis
nr_of_seconds       = 3; % Before/after SWR peak
nr_of_shuffles      = 1000; % nr of shuffle iterations
min_nr_events       = 5; 
ephys_to_imag_time  = sData.daqdata.frame_onset_reference_frame;
win_length          = (nr_of_seconds*31)*2+1;

%% Run main analysis

% Preallocate
signal_event_activity             = zeros( length(event_idx), win_length);
% signal_event_activity_zscore      = zeros( length(event_idx), win_length);
signal_event_activity_mean        = zeros( size(signal,1), win_length);
% signal_event_activity_mean_zscore = zeros( size(signal,1), win_length);
event_mod_ROI_idx                 = zeros( size(signal,1) ,1);

t = 1;
% roi_idx = zeros(1, size(signal,1));
% Loop over ROIs
for roi_nr = 1:size(signal,1) 

    % Extract signal from currrent ROI
    roi_signal        = signal(roi_nr, :); 

    % Get all peri-SWR activity windows from current ROI
    peri_event_activity = extract_avg_activity(event_idx, roi_signal, peri_event_activity, nr_of_seconds, params);

    % store mean SWR-activity
    mean_event_activity(roi_nr,:)  = mean(peri_event_activity, 'omitnan');

    % store all roi x time x SWR data
    all_data(t,:,:) = peri_event_activity;
    t = t+1;
    % Shuffle analysis
    event_mod_ROI_idx = shuffle_analysis(win_length, nr_of_shuffles, peri_event_activity, ...
        roi_nr, event_mod_ROI_idx, time_win, min_nr_events);
end

event_mod_ROI_idx(event_mod_ROI_idx ==0) = [];

% Determine positive or negative event modulation
[ROIs_activated_idx, ROIs_suppressed_idx] = split_modulated_rois(mean_event_activity, win_length, params, event_mod_ROI_idx);

% Store results
sData.imdata.([params.event_type, '_', params.cell_type, '_pos_cells']) = ROIs_activated_idx;
sData.imdata.([params.event_type, '_', params.cell_type, '_neg_cells']) = ROIs_suppressed_idx;
sData.imdata.([params.event_type, '_', params.cell_type, '_params'])     = params;