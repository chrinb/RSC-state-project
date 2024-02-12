function output = plot_state_sextiles(sData, adjusted_states, params)

[signal_to_plot, ~] = get_roi_signals_from_sData(sData, params);

% state_vectors_2p = get_state_logicals(sData, adjusted_states);
state_vectors_2p = adjusted_states;

% Get 2P frame rate
imaging_sampling_rate = find_imaging_framerate(sData);

if isfield(sData, 'episodes')
    REM_episodes  = rem_sleep(sData);
    NREM_episodes = nrem_sleep(sData);

    % NREM_start_end = round( (NREM_episodes./2500)*imaging_sampling_rate);
    [nrem_start, nrem_stop] = findTransitions(state_vectors_2p{1,1});
    NREM_start_end = [nrem_start', nrem_stop'];

    % REM_start_end  = round( (REM_episodes./2500)*imaging_sampling_rate);
     [rem_start, rem_stop] = findTransitions(state_vectors_2p{1,2});
    REM_start_end = [rem_start', rem_stop'];
end

%% select signal
if strcmp(params.cell_type, 'pc')
    data = signal_to_plot{1,1};
elseif strcmp(params.cell_type, 'in')
        data = signal_to_plot{2,1};
elseif strcmp(params.cell_type, 'axon')
        data = signal_to_plot{1,1};
end


%% Remove NaNs
data = rmmissing(data);

%% Filter
if strcmp(params.filter,'yes') && ~strcmp(params.cell_type, 'axon')
    data = okada(data, 2);
end

%% Sort ROIs based on activity in quiet wakefulness
switch params.sorting
    case 'QW'
    data_to_sort = data(:, state_vectors_2p{1,3});
    roiStat      = getRoiActivityStats(data_to_sort);
    [~, sortID]  = sort(roiStat.activityLevel, 'descend');

    % sum_sorted_data = sum(data_to_sort,2);
    
    case 'NREM'
    data_to_sort = data(:, state_vectors_2p{1,1});
    roiStat      = getRoiActivityStats(data_to_sort);
    [~, sortID]  = sort(roiStat.activityLevel, 'descend');

    % sum_sorted_data = sum(data_to_sort,2);

    case 'Session'
    % sum_sorted_data = sum(data,2);
     roiStat   = getRoiActivityStats(data);
    [~, sortID] = sort(roiStat.activityLevel, 'descend');

    case 'QW+NREM'
    qw_nrem_combined = sum( [state_vectors_2p{1,1}; state_vectors_2p{1,3}]) > 0;
    data_to_sort    = data(:, qw_nrem_combined);
    roiStat    = getRoiActivityStats(data_to_sort);
    [~, sortID] = sort(roiStat.activityLevel, 'descend');

    case 'AW'
    data_to_sort = data(:, state_vectors_2p{1,4});
    roiStat      = getRoiActivityStats(data_to_sort);
    [~, sortID]  = sort(roiStat.activityLevel, 'descend');

end

activity_sorted = roiStat.activityLevel(sortID);
peak_sorted     = roiStat.peakDff(sortID);
snr_sorted      = roiStat.signalToNoise(sortID);
% [mean_qw_sorted, sortID] = sort(sum_sorted_data, 'descend');

data_sorted = data(sortID,:);

% Compute ROI statistics and alternative sorting indices
% roiStat = getRoiActivityStats(data_sorted);
roiStat_all_session = getRoiActivityStats(data);

[session_activity_sorted, activity_sortID] = sort(roiStat_all_session.activityLevel, 'descend');
[session_snr_sorted, snr_sortID]           = sort(roiStat_all_session.signalToNoise, 'descend');
[session_peak_sorted, peak_sortID]         = sort(roiStat_all_session.peakDff, 'descend');

% data_sorted = data_sorted(activity_sortID,:);
% 
% roiStat = getRoiActivityStats(data_sorted);
% 
% [activity_sorted, activity_sortID] = sort(roiStat.activityLevel, 'descend');
% [snr_sorted, snr_sortID]           = sort(roiStat.signalToNoise, 'descend');
% [peak_sorted, peak_sortID]         = sort(roiStat.peakDff, 'descend');
%% Sort into sextiles
% Define the row indices for sextiles
sextileBoundaries = round(linspace(1, size(data, 1) + 1, 7));
 
% Initialize cell array to store sextile data
[sextiles, activity_level, peak_dff, signal_to_noise] = deal( cell(1, 6));

% Split the matrix into sextiles
for i = 1:6
    rowIndices         = sextileBoundaries(i):sextileBoundaries(i + 1) - 1;
    sextiles{i}        = data_sorted(rowIndices, :);
    activity_level{i}  = activity_sorted(rowIndices);
    peak_dff{i}        = peak_sorted(rowIndices);
    signal_to_noise{i} = snr_sorted(rowIndices);

    sorted_activity{i} = session_activity_sorted(rowIndices);
    sorted_snr{i}      = session_snr_sorted(rowIndices);
    sorted_peak{i}     = session_peak_sorted(rowIndices);
end

%% Calculate activity per sextile

% REM
if ~isempty(REM_start_end)
    rem_top_sextile    = arrayfun(@(x) sextiles{1,1}(:,  REM_start_end(x,1):REM_start_end(x,2)), (1:size(REM_start_end,1)), 'UniformOutput', false);
    rem_bottom_sextile = arrayfun(@(x) sextiles{1,6}(:,  REM_start_end(x,1):REM_start_end(x,2)), (1:size(REM_start_end,1)), 'UniformOutput', false);
    
    n_bins_rem = repmat({10}, 1, numel(rem_top_sextile)); 
    
    [norm_rem_mean_top, norm_rem_SEM_top]       = cellfun(@episode_third_analysis, rem_top_sextile, n_bins_rem, 'uni', 0);
    [norm_rem_mean_bottom, norm_rem_SEM_bottom]  = cellfun(@episode_third_analysis, rem_bottom_sextile, n_bins_rem,'uni', 0);
    
    if size(norm_rem_SEM_top,1) > 1
        rem_binned_mean_top    = mean( vertcat(norm_rem_mean_top{:}));
        rem_binned_mean_bottom = mean( vertcat(norm_rem_mean_bottom{:}));
    
        rem_binned_sem_top    = mean( vertcat(norm_rem_SEM_top{:}));
        rem_binned_sem_bottom = mean( vertcat(norm_rem_SEM_bottom{:}));
    else
        rem_binned_mean_top    = norm_rem_mean_top{:};
        rem_binned_mean_bottom = norm_rem_mean_bottom{:};
    
        rem_binned_sem_top    = norm_rem_SEM_top{:};
        rem_binned_sem_bottom = norm_rem_SEM_bottom{:};
    end
else
        rem_binned_mean_top    = NaN(1,10);
        rem_binned_mean_bottom = NaN(1,10);
        rem_binned_sem_top     = NaN(1,10);
        rem_binned_sem_bottom  = NaN(1,10);
end

% NREM
if ~isempty(NREM_start_end)

    nrem_top_sextile    = arrayfun(@(x) sextiles{1,1}(:,  NREM_start_end(x,1):NREM_start_end(x,2)), (1:size(NREM_start_end,1)), 'UniformOutput', false);
    nrem_bottom_sextile = arrayfun(@(x) sextiles{1,6}(:,  NREM_start_end(x,1):NREM_start_end(x,2)), (1:size(NREM_start_end,1)), 'UniformOutput', false);
    
    n_bins_nrem = repmat({10}, 1, numel(nrem_top_sextile)); 
    
    [norm_nrem_mean_top, norm_nrem_SEM_top]        = cellfun(@episode_third_analysis, nrem_top_sextile, n_bins_nrem, 'uni', 0);
    [norm_nrem_mean_bottom, norm_nrem_SEM_bottom]  = cellfun(@episode_third_analysis, nrem_bottom_sextile, n_bins_nrem,'uni', 0);
    
    if size(norm_nrem_SEM_bottom, 1) > 1
        nrem_binned_mean_top    = mean( vertcat(norm_nrem_mean_top{:}));
        nrem_binned_mean_bottom = mean( vertcat(norm_nrem_mean_bottom{:}));
    
        nrem_binned_sem_top    = mean( vertcat(norm_nrem_SEM_top{:}));
        nrem_binned_sem_bottom = mean( vertcat(norm_nrem_SEM_bottom{:}));
    else
        nrem_binned_mean_top    = norm_nrem_mean_top{:};
        nrem_binned_mean_bottom =norm_nrem_mean_bottom{:};
    
        nrem_binned_sem_top    = norm_nrem_SEM_top{:};
        nrem_binned_sem_bottom = norm_nrem_SEM_bottom{:};
    end
else
        nrem_binned_mean_top    = NaN(1,10);
        nrem_binned_mean_bottom = NaN(1,10);
        nrem_binned_sem_top     = NaN(1,10);
        nrem_binned_sem_bottom  = NaN(1,10);
end

output{1,1} = rem_binned_mean_top;
output{1,2} = rem_binned_mean_bottom;
output{1,3} = rem_binned_sem_top;
output{1,4} = rem_binned_sem_bottom;

output{1,5} = nrem_binned_mean_top;
output{1,6} = nrem_binned_mean_bottom;
output{1,7} = nrem_binned_sem_top;
output{1,8} = nrem_binned_sem_bottom;

output{1,9}  = activity_level{1};
output{1,10} = activity_level{6};
output{1,11} = peak_dff{1};
output{1,12} = peak_dff{6};  
output{1,13} = signal_to_noise{1};
output{1,14} = signal_to_noise{6};

output{1,15} = sorted_activity{1};
output{1,16} = sorted_activity{6};
output{1,17} = sorted_peak{1};
output{1,18} = sorted_peak{6};
output{1,19} = sorted_snr{1};
output{1,20} = sorted_snr{6};
%% Plot session
time_vector = linspace(0, length(sData.imdata.roiSignals(2).newdff), length(sData.imdata.roiSignals(2).newdff))/ imaging_sampling_rate;

top_sextile_n_cells    = size( sextiles{1,1},1);
bottom_sextile_n_cells = size( sextiles{1,6}, 1);

scale_factor = 6;

figure, 
hold on
plot( time_vector, sextiles{1,1} + (1:top_sextile_n_cells)'*scale_factor, 'color', 'r')
plot( time_vector, sextiles{1,6} + ((2:bottom_sextile_n_cells+1)+bottom_sextile_n_cells)'*scale_factor, 'color', 'b')

plot(time_vector, mean( sextiles{1,1})-8, 'color', 'r')
plot(time_vector, mean( sextiles{1,6})-10, 'Color', 'b')

temp_y_max = sum(top_sextile_n_cells + bottom_sextile_n_cells)*scale_factor+10;
% temp_y_max = max( 1:20'*scale_factor)+40 ;
temp_y_min = -20;
 if isfield(sData, 'episodes')

        % NREM bouts
        for ep_nr = 1:length(NREM_start_end)
            x = [NREM_start_end(ep_nr,1) NREM_start_end(ep_nr,1) NREM_start_end(ep_nr,2) NREM_start_end(ep_nr,2)]./imaging_sampling_rate;
            y = [temp_y_min temp_y_max temp_y_max temp_y_min];
            patch(x, y, 'green', 'edgecolor', 'none', 'FaceAlpha', .2);
        end

        % REM bouts
        if size(REM_start_end(:,1),1) == 1

            a = [REM_start_end(1,1) REM_start_end(1,1) REM_start_end(1,2) REM_start_end(1,2)]./imaging_sampling_rate;
            b =[temp_y_min temp_y_max temp_y_max temp_y_min];
            patch(a, b, 'red', 'edgecolor', 'none', 'FaceAlpha', .2);
        else
            for ep_nr = 1:length(REM_start_end)
                a = [REM_start_end(ep_nr,1) REM_start_end(ep_nr,1) REM_start_end(ep_nr,2) REM_start_end(ep_nr,2)]./imaging_sampling_rate;
                b = [temp_y_min temp_y_max temp_y_max temp_y_min];
                patch(a, b, 'red', 'edgecolor', 'none', 'FaceAlpha', .2);
            end
        end
 end

 title(['Sorting: ', params.sorting, ' activity'])