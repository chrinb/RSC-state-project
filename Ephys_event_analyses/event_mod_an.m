function sData = event_mod_an(varargin)

% Written by Christoffer Berge | Vervaeke Lab

% Function that looks for event-modulated ROIs by creating a Peri-event 
% time histogram (PETH) for each ROI and comparing the mean squared error 
% (MSE) of that PETH to a distribution of MSE's created from random,
% cirular shuffeling of the original PETH. 


sData  = varargin{1,1};
params = varargin{1,2};

% Get ROI signals
[signal, params, ~, rois_for_an, ~, ~] = get_roi_signals_from_sData(sData, params );

% Z-score data. 
zscore_signal = bsxfun(@minus, signal, mean(signal,2) ); 
zscore_signal = bsxfun(@rdivide, zscore_signal, std(signal,0,2));


%% Define variables for analysis
nr_of_seconds       = 3; % Before/after SWR peak
nr_of_shuffles      = 1000; % nr of shuffle iterations
min_nr_events       = 5; 
ephys_to_imag_time  = sData.daqdata.frame_onset_reference_frame;
win_length          = (nr_of_seconds*31)*2+1;

%% Select events for analysis

switch params.event_type
    case 'swr'
        select_swr_idx = select_swr_for_mod_an(sData, params);
        
        % Get the indicies of user specified SWR types
        event_idx = get_swr_idx(params.swr_for_analysis ,sData,select_swr_idx, params);
        
        % Convert SWR time stamps from e-phys to 2P time
        event_idx = ephys_to_imag_time(round(event_idx));
    
        time_win = (62:108); % -1s to +500ms after SWR peak

    case 'spindle'
        spindle_select     = [];
        spin_start_end_str = strcat('NREMspindleStartEnd', spindle_select);
        
        % convert sleep spindle time stamps from e-phys to 2P time
        event_idx = sData.ephysdata2.(spin_start_end_str)(:,1);
    
        % Get time points for sleep spindle onsets
        event_idx = get_spindle_idx(params , event_idx);
        
        time_win             = (93:124); % 0s to +1s after spindle onset

    case 'SWA'
end

%% Run main analysis

% Preallocate
signal_event_activity             = zeros( length(event_idx), win_length);
signal_event_activity_zscore      = zeros( length(event_idx), win_length);
signal_event_activity_mean        = zeros( size(signal,1), win_length);
signal_event_activity_mean_zscore = zeros( size(signal,1), win_length);
event_mod_ROI_idx                 = zeros( size(rois_for_an,1) ,1);

% Loop over ROIs
for roi_nr = rois_for_an

    roi_signal        = signal(roi_nr, :); 
    roi_signal_zscore = zscore_signal(roi_nr,:); 

    % Get all peri-SWR activity 
    [signal_event_activity, signal_event_activity_zscore] = extract_avg_activity(event_idx, ...
        roi_signal, roi_signal_zscore, signal_event_activity, signal_event_activity_zscore, ...
        nr_of_seconds, params, roi_nr);

    % Mean
    signal_event_activity_mean(roi_nr,:)        = mean(signal_event_activity, 'omitnan');
    signal_event_activity_mean_zscore(roi_nr,:) = mean(signal_event_activity_zscore, 'omitnan');

    % Shuffle analysis
    event_mod_ROI_idx = shuffle_analysis(win_length, nr_of_shuffles, signal_event_activity, ...
        roi_nr, event_mod_ROI_idx, signal_event_activity_mean, time_win, min_nr_events);
end

event_mod_ROI_idx(event_mod_ROI_idx ==0) = [];

% Determine positive or negative event modulation
[ROIs_activated_idx, ROIs_suppressed_idx] = split_modulated_rois(signal_event_activity_mean_zscore, win_length, params, event_mod_ROI_idx);

% Store results
sData.imdata.([params.event_type, '_', params.cell_type, '_pos_cells']) = ROIs_activated_idx;
sData.imdata.([params.event_type, '_', params.cell_type, '_neg_cells']) = ROIs_suppressed_idx;
sData.imdata.([params.event_type, '_', params.cell_type, '_params'])     = params;