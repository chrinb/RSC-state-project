function ratios_rel_to_QW = roi_state_preference_idx(sData, params)

%{
Compare each ROI mean activity during quiet wakefulness against other
states to determine state preference.
%}

mean_data_all_rois = cell( 1, 4);
% Load mean activity per ROI per state
try mean_data_all_rois{:,1} = sData.analysis.state_activity.(['mean_active_', params.signal_type, '_', params.cell_type]); catch, end
try mean_data_all_rois{:,2} = sData.analysis.state_activity.(['mean_quiet_' params.signal_type, '_', params.cell_type]); catch, end
try mean_data_all_rois{:,3} = sData.analysis.state_activity.(['mean_NREM_', params.signal_type, '_', params.cell_type]); catch, end
try mean_data_all_rois{:,4} = sData.analysis.state_activity.(['mean_REM_', params.signal_type, '_', params.cell_type]); catch, end

% Check if sessions contains quiet wakefulness data
if isempty(mean_data_all_rois{1,2} )
    msgbox(['No quiet wakefulness data for session: ', num2str(sData.sessionInfo.sessionID  )])
    return
end

% Check if quiet wakefulness data contains one value or mean activity per
% ROI
if size(mean_data_all_rois{1,2} ,2) == 1
    msgbox('Only one data point in mean activity field. Rerun mean activity analysis!')
    return
end

%% Compute ratio of mean state activity relative to mean quiet wakefulness
ratios_rel_to_QW = zeros( size(mean_data_all_rois{1,2}, 2), 3);
state_list       = [1 3 4];

for state_nr = 1:3
   
    idx = state_list(state_nr);
    if sum( mean_data_all_rois{1, idx}  > 0) > 0 && ~isempty(mean_data_all_rois{1, idx})
        
        ratios_rel_to_QW(:, state_nr) = abs( mean_data_all_rois{1, idx} ./ mean_data_all_rois{1, 2} );
    end
end
