function output = roi_state_preference_idx(sData)

%{
Compare each ROI mean activity during quiet wakefulness against other
states to determine state preference.
%}
params.signal_type   = 'dff'; 
params.zscore        = 'no';
params.cell_type     = 'axon';   

mean_activity = sData.analysis.state_activity;
[signal_to_plot, cmap, pc_rois, in_rois] = get_roi_signals_from_sData(sData, params);

n_rois = size( signal_to_plot{1,1},1)

for roi_nr = 1:n_rois
    