function output = roi_state_preference_idx(sData)

%{
Compare each ROI mean activity during quiet wakefulness against other
states to determine state preference.
%}
params.signal_type   = 'Dff'; 
params.zscore        = 'no';
params.cell_type     = 'pc';   

% mean_activity = sData.analysis.state_activity;
% [signal_to_plot, cmap, pc_rois, in_rois] = get_roi_signals_from_sData(sData, params);

i = 1;
try AW_data_all_rois   = sData.analysis.state_activity.(['mean_active_', params.signal_type, '_', params.cell_type]); catch, end
try QW_data_all_rois   = sData.analysis.state_activity.(['mean_quiet_' params.signal_type, '_', params.cell_type]); catch, end
try NREM_data_all_rois = sData.analysis.state_activity.(['mean_NREM_', params.signal_type, '_', params.cell_type]); catch, end
try REM_data_all_rois  = sData.analysis.state_activity.(['mean_REM_', params.signal_type, '_', params.cell_type]); catch, end

[pc_idx, in_idx]     = remove_cells_longitudinal(sData);

%% 

% QW_data_all_rois = session_state_data{i,2};

% Check if sessions contain 
if isempty(QW_data_all_rois)
    msgbox(['No quiet wakefulness data for session: ', num2str(sData.sessionInfo.sessionID  )])
    return
end

% Check if correct version of mean activity analysis have been ran
if size(QW_data_all_rois,2) == 1
    msgbox('Only one data point in mean activity field. Rerun mean activity analysis!')
    return
end

% QW_data_roi = QW_data_all_rois(roi_nr);

roi_state_averages_cat = horzcat(QW_data_all_rois', NREM_data_all_rois', REM_data_all_rois');

% Determine which states have max values
% [max_val, max_id] = max(roi_state_averages_cat, [], 2);

% Ratio of activity relative to QW (high ratio = higher activity in that
% state than QW)
try AW_QW_ratio = AW_data_all_rois./QW_data_all_rois; catch, end
try NREM_QW_ratio = NREM_data_all_rois./QW_data_all_rois; catch, end
try REM_QW_ratio = REM_data_all_rois./QW_data_all_rois; catch, end


% [~, QW_ratio_sort_idx] = sortrows(QW_ratio)

%% Plot 
jitter_factor = 0.005;
rand_vec_for_plot_NREM = (randn(size(roi_state_averages_cat(:,1),1), 1)*jitter_factor) + 1;
rand_vec_for_plot_REM  = (randn(size(roi_state_averages_cat(:,1),1), 1)*jitter_factor) + 1.2;
rand_vec_for_plot_AW   = randn(size(roi_state_averages_cat(:,1),1), 1)*jitter_factor;


[sorted_REM_ratio, sort_idx_REM ] = sort(REM_QW_ratio);
[sorted_NREM_ratio, sort_idx_NREM ] = sort(NREM_QW_ratio);

below_one_idx_NREM = NREM_QW_ratio < 1;
above_one_idx_NREM = NREM_QW_ratio > 1;

below_one_idx_REM = REM_QW_ratio < 1;
above_one_idx_REM = REM_QW_ratio > 1;

figure,
scatter(rand_vec_for_plot_NREM(below_one_idx_NREM), NREM_QW_ratio((below_one_idx_NREM)),'filled'), hold on
scatter(rand_vec_for_plot_NREM(above_one_idx_NREM), NREM_QW_ratio((above_one_idx_NREM)),'filled')
plot(1, mean(NREM_QW_ratio), 'ko', 'LineWidth', 3)

scatter(rand_vec_for_plot_REM(below_one_idx_REM), REM_QW_ratio((below_one_idx_REM)),'filled'), hold on
scatter(rand_vec_for_plot_REM(above_one_idx_REM), REM_QW_ratio((above_one_idx_REM)),'filled')
plot(1.2, mean(REM_QW_ratio), 'ko', 'LineWidth', 3)

plot_middle = mean( [1,1.2]);
plot_edgeL  = plot_middle-0.3;
plot_edgeR  = plot_middle+0.3;
set(gca, 'xlim', [plot_edgeL plot_edgeR])
