function sData = state_activity(sData, params)

% Written by Christoffer Berge | Vervaeke lab

%{
Function computes (1) mean DF/F per ROI and on average per state, (2) 
calcium transient & deconvolved event rate (Hz) per state, (3) proportion of 
transient and deconvolved events per state, and (4) percentages of cells
with maximal DF/F per state
%}


%% Get exctitatory and inhibitory indices
[pc_rois, in_rois] = remove_cells_longitudinal(sData);

% Select data
switch params.cell_type
    case 'axon'
     dff    = sData.imdata.roiSignals(2).mergedAxonsDffFilt;
     dec = sData.imdata.roiSignals(2).mergedAxonsDec;

    case 'in'
    dff = sData.imdata.roiSignals(2).newdff(in_rois,:);
    dec = sData.imdata.roiSignals(2).ciaDeconvolved(in_rois,:);

    case 'pc'
    dff = sData.imdata.roiSignals(2).newdff(pc_rois,:);
    dec = sData.imdata.roiSignals(2).ciaDeconvolved(pc_rois,:);
end
% Z-score data
dff_zscore = zscore(dff, [], 2);
dec_zscore = zscore(dec, [], 2);

n_rois = size(dff,1);

% Convert state logical vector from ephys to 2P time
state_vectors = get_state_logicals(sData);

%% Compute mean DF/F
n_states    = 4;
state_label = {'NREM', 'REM', 'quiet', 'active'};

imaging_sampling_rate = find_imaging_framerate(sData);

% Loop over the different states
for state_nr = 1:n_states

    % Get state vector (in 2P time) for current state
    current_state_vec_2P = state_vectors{1, state_nr};
    
    state_name = ['mean_', state_label{state_nr}];

    % Check that sessions does in fact contain state
    if ~isempty(current_state_vec_2P)

        % Extract DF/F during 
        state_dff = dff(:, current_state_vec_2P);
        state_dec = dec(:, current_state_vec_2P);

        total_state_length_sec = size(state_dff,2)/ imaging_sampling_rate;

        state_dff_z = dff_zscore(:, current_state_vec_2P);
        state_dec_z = dec_zscore(:, current_state_vec_2P);
        
        % Compute mean DF/F per ROI
        roi_mean_dff   = mean(state_dff, 2, 'omitnan');
        roi_mean_dff_z = mean(state_dff_z, 2, 'omitnan');

        % Compute mean DF/F across all ROIS during session
        mean_dff        =  mean(roi_mean_dff,'omitnan');
        mean_dff_zscore =  mean(roi_mean_dff_z,'omitnan');

        % Compute mean firing rate (Hz) per ROI

        % Consider binarizing deconvolved signal based on some threshold (as
        % Grosmark et al 2021)
        roi_dec_peaks         = arrayfun(@(x) findpeaks(state_dec(x,:)), (1:n_rois), 'uni', 0);
        roi_dec_num_events    = cell2mat( cellfun(@numel, roi_dec_peaks, 'uni', 0));
        mean_deconvolved_rate = roi_dec_num_events/total_state_length_sec;

        % Compute mean transient event rate, duration, and peak DF/F
        for roi_nr = 1:n_rois   
            roi_transient_idx{roi_nr} = current_state_vec_2P( sData.analysis.transients.([params.cell_type, '_transient_times']){roi_nr, 1} );
        end
        anon                     = @(x) sum(x > 0);
        roi_transient_num        = cell2mat( cellfun(anon, roi_transient_idx, 'uni', 0));
        mean_transient_event_rate = roi_transient_num/total_state_length_sec;
        
        anon2                   = @(x,y) x(y);
        transient_durations     = cellfun(anon2, sData.analysis.transients.([params.cell_type, '_transient_dur_sec']), roi_transient_idx', 'uni', 0); 
        transient_peak_dff      = cellfun(anon2, sData.analysis.transients.([params.cell_type, '_transient_ampl']), roi_transient_idx', 'uni', 0); 
        mean_transient_duration = cell2mat( cellfun(@mean, transient_durations, 'uni', 0));
        
        % Store in sData
        sData.analysis.state_activity.([state_name, '_Dff_', params.cell_type])              = mean_dff;
        sData.analysis.state_activity.([state_name, '_Dff_zscore_', params.cell_type])       = mean_dff_zscore;
        sData.analysis.state_activity.([state_name,'_deconv_rate_', params.cell_type])       = mean_deconvolved_rate;
        sData.analysis.state_activity.([state_name,'_transient_rate_', params.cell_type])     = mean_transient_event_rate;
        sData.analysis.state_activity.([state_name,'_transient_duration_', params.cell_type]) = mean_transient_duration;
        sData.analysis.state_activity.([state_name,'_transient_peakDFF_', params.cell_type])  = transient_peak_dff;
    end
end

