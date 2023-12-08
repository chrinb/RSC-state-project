function slope = calc_rem_slope_best_fit(sData, params)

%{
Calculate the slope from the best fit line of the mean population activity
during REM
%}

params.cell_type                = 'pc';   
% params.cell_type_               = {'pc','in', 'axon'};
params.zscore                   = 'no';
% params.zscore_                  = {'no', 'yes'};
params.beh_state                = {'NREM'};
% params.beh_state_               = {'NREM', 'REM', 'AW', 'QW'};
params.signal_type              = 'Dff';
% params.signal_type_             = {'Dff'; 'deconv'; 'transients'};
params.use_roi_classification   =  'grid';   
% params.use_roi_classification_  = {'grid', 'ch2_across_session', 'no'};
params.filter                   = 'no';
% params.filter_                  = {'yes', 'no'};


% Get signals
[signal_to_plot, ~] = get_roi_signals_from_sData(sData, params);

if strcmp(params.cell_type, 'pc') || strcmp(params.cell_type, 'axon')
    signal = signal_to_plot{1,1};
elseif strcmp( params.cell_type, 'in')
    signal = signal_to_plot{1,2};
end

if strcmp(params.filter, 'yes')
    signal = okada(signal,2);
end

% Get 2P frame rate
% imaging_sampling_rate = find_imaging_framerate(sData);

frames = sData.daqdata.frame_onset_reference_frame;

if isfield(sData, 'episodes')
    REM_episodes = rem_sleep(sData);
end

REM_start_end = zeros( size(REM_episodes,1),2);
slope         = zeros( size(REM_episodes,1),1);
for rem_ep_nr = 1:size(REM_episodes,1)
    REM_start_end(rem_ep_nr,:) = frames( [REM_episodes(rem_ep_nr,1), REM_episodes(rem_ep_nr,2)] );

    % Get REM episode snippet
    signal_rem_snip      = signal(:, REM_start_end(1):REM_start_end(2) );
    signal_mean_rem_snip = mean(signal_rem_snip);
    
    % Best fit line
    x = 1:length(signal_mean_rem_snip);
    p = polyfit(x, signal_mean_rem_snip, 1);
    
    % 
    % figure, 
    % plot(x, signal_mean_rem_snip)
    % hold on
    % plot(x, polyval(p, x), 'LineWidth',2);
    
    slope(rem_ep_nr) = p(1);
end
