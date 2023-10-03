function sData = state_activity(varargin)

% Written by Christoffer Berge | Vervaeke lab

%{
Function computes (1) mean DF/F per ROI and on average per state, (2) 
calcium transient & deconvolved event rate (Hz) per state, (3) proportion of 
transient and deconvolved events per state, and (4) percentages of cells
with maximal DF/F per state
%}

sData = varargin{1,1};

[pc_rois, in_rois] = remove_cells(sData);

if isfield(sData.imdata, 'roi_classification')
    pc_roi_idx = sData.imdata.roi_classification(pc_rois); % of all ROIs, index out cell type
    in_roi_idx = sData.imdata.roi_classification(in_rois); % of all ROIs, index out cell type
    log_idx_pc = pc_roi_idx == 1;
    log_idx_in = in_roi_idx == 1;
    pc_rois    = pc_rois(log_idx_pc);
    in_rois    = in_rois(log_idx_in);
end

% Load data
if nargin > 1 
    switch varargin{1,2}
        case 'axon'
             dff = sData.imdata.roiSignals(2).mergedAxonsDffFilt;
             dec = sData.imdata.roiSignals(2).mergedAxonsDec;
             cell_type = 'axon';
        case 'in'
            dff = sData.imdata.roiSignals(2).newdff(in_rois,:);
            dec = sData.imdata.roiSignals(2).ciaDeconvolved(in_rois,:);
            cell_type = 'in';
        case 'pc'
            dff = sData.imdata.roiSignals(2).newdff(pc_rois,:);
            dec = sData.imdata.roiSignals(2).ciaDeconvolved(pc_rois,:);
            cell_type = 'pc';
    end
else
    dff = sData.imdata.roiSignals(2).newdff;
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
    
    field_name = ['mean_', state_label{state_nr}];

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
        state_mean_dff        =  mean(roi_mean_dff,'omitnan');
        state_mean_dff_zscore =  mean(roi_mean_dff_z,'omitnan');

        % Compute mean firing rate (Hz) per ROI

        % Consider binarizing deconvolved signal based on some threshold (as
        % Grosmark et al 2021)
%         nrem_dec = dec;
%         nrem_dec(dec > 0) = 1;
        roi_dec_peaks      = arrayfun(@(x) findpeaks(state_dec(x,:)), (1:n_rois), 'uni', 0);
        roi_dec_num_events = cell2mat( cellfun(@numel, roi_dec_peaks, 'uni', 0));
        roi_dec_event_rate = roi_dec_num_events/total_state_length_sec;

        % Compute NREM event rate (Ca2+ transients and deconvolved) and total nr of
        % transients in NREM.
        
        % Preallocate
        state_event_rate     = zeros(n_rois,1 );
        state_deconv_rate    = zeros(n_rois,1 );
        max_ampl_time       = zeros(n_rois,1);
        [all_and_nrem_events, all_and_nrem_events_dec] = deal( zeros(n_rois,2 ));
        
        % Loop over ROIs
        for roi_nr = 1:n_rois
        
            % Create zero vector equal to nr of frames in recording. Get transient
            % times for a ROI and insert them as 1s in zero vector
    
    %         nrem_log_vec      = zeros(1, size(dff,2));
    %         transient_times = sData.analysis.transients.transient_times{roi_nr, 1};
    %         nrem_log_vec(transient_times) = 1; 
    %     
            % Repeat for deconvolved data
            nrem_log_vec_dec     = zeros(1, size(dff,2));
            deconv_times         = 1:size(dff,2);
            
            % quick fix for ROIs with NaN values:
            try
                deconv_times_logical = logical( nrem_dec(roi_nr,:));
                deconv_times  = deconv_times(deconv_times_logical);
                nrem_log_vec_dec(deconv_times) = 1;
            catch
                nrem_log_vec_dec(deconv_times) = NaN;
            end
        
            % Conatenate NREM logical vector (NREM frames = 1) and zeros vector with 
            % transient times = 1s. Sum them and them remove values = 1 (indicating events outside 
            % of NREM epochs or just NREM epochs) and set values = 2 
            % (indicating concurrent transient and NREM frames) equal to 1. The
            % resulting n_nrem_events contains all NREM events.

    %         temp_vec = [nrem_log_vec; NREM_vec_2P]; 
    %         n_nrem_events = sum(temp_vec,1);
    %         n_nrem_events(n_nrem_events == 1) = 0;
    %         n_nrem_events(n_nrem_events == 2) = 1;
        
            % Repeat for deconvolved events
            temp_vec_dec = [nrem_log_vec_dec; current_state_vec_2P]; 
            n_nrem_dec = sum(temp_vec_dec,1);
            n_nrem_dec(n_nrem_dec == 1) = 0;
            n_nrem_dec(n_nrem_dec == 2) = 1;
        
            % Compute fraction of events in NREM
    %         all_and_nrem_events(roi_nr, 1) = sum(n_nrem_events) ;
    %         all_and_nrem_events(roi_nr, 2) = numel(transient_times);
        
            % And for deconvolved
            all_and_nrem_events_dec(roi_nr, 2) = sum(n_nrem_dec) ;
            all_and_nrem_events_dec(roi_nr, 1) = numel(deconv_times);
        
            % Sum NREM logical vector to obtain NREM duration (in frames). Sum n_nrem_events
            % and divide by NREM duration to get average event rate during NREM. 
            nrem_duration            = sum( current_state_vec_2P);
    %         nrem_event_rate(roi_nr)  = sum(n_nrem_events)/nrem_duration;
            state_deconv_rate(roi_nr) = sum(n_nrem_dec)/(nrem_duration/31);
        
            % Find the frame nr of the max DF/F Ca2+ transient for each ROI
    %         roi_ampl              = sData.analysis.transients.all_amplitudes{roi_nr, 1};
    %         [max_ampl, max_id]    = max(roi_ampl); 
%             max_ampl_time(roi_nr) = transient_times(max_id);
        end
        
        % Find which of the ROI's max DF/F transients occurred in NREM. First
        % create zeros vector where the max amplitude tranients are indicates
        % as 1s. Then concatenate the zeros vector with NREM logical vector
        % (NREM = 1s), sum, and remove values = 1 (indicating NREM sleep or max
        % transient but not both) and set values = 2 = 1 (indicating overlap
        % between max transients and NREM sleep

    %     find_max                = zeros(1, size(dff,2));
    %     find_max(max_ampl_time) = 1;
    %     % Correct for transients occurring at same frame
    %     temp     = max_ampl_time;
    %     [~, idx,~] = unique(temp);
    %     temp(idx) = [];
    %     for i = 1:size(temp,1)
    %         find_max(temp(i)) = find_max(temp(i))+1;
    %     end
    % 
    %     if ~(sum(find_max) == roi_nr)
    %         error('Nr of max DF/F transients does not match nr of ROIs')
    %     end
    
        % TO DO: need to fix transients that overlap!!!
    %     cat_vec               = [find_max; NREM_vec_2P]; 
    %     sum_vec               = sum(cat_vec,1);
    %     sum_vec(sum_vec == 1) = 0;
    %     sum_vec(sum_vec == 2) = 1;
    %     
    %     sum_max_in_nrem       = sum(sum_vec);
        % fraction_of_events_in_nrem = sum( all_and_nrem_events(:,1)) / sum( all_and_nrem_events(:,2));
        
    % Store mean DF/F during NREM sleep

    sData.analysis.state_activity.([field_name, '_Dff_', cell_type])                  = state_mean_dff;
    sData.analysis.state_activity.([field_name, '_Dff_zscore_', cell_type])           = state_mean_dff_zscore;

    %     sData.analysis.state_activity.NREM_transient_rate           = nrem_event_rate;
    sData.analysis.state_activity.([field_name,'_deconv_rate_', cell_type])              = state_deconv_rate;
    %     sData.analysis.state_activity.all_and_NREM_transient_events = all_and_nrem_events;
%     sData.analysis.state_activity.(['all_and_NREM_deconv_events_', cell_type])    = all_and_nrem_events_dec;
%         sData.analysis.state_activity.sum_max_in_NREM               = sum_max_in_nrem;
    else
            % Store mean DF/F during NREM sleep
        sData.analysis.state_activity.([field_name, '_Dff_', cell_type])                 = [];
        sData.analysis.state_activity.([field_name, '_Dff_zscore_', cell_type])           = [];

        %     sData.analysis.state_activity.NREM_transient_rate           = nrem_event_rate;
        sData.analysis.state_activity.([field_name,'_deconv_rate_', cell_type])              = [];
        %     sData.analysis.state_activity.all_and_NREM_transient_events = all_and_nrem_events;
%         sData.analysis.state_activity.(['all_and_NREM_deconv_events_', cell_type])    = all_and_nrem_events_dec;
%         sData.analysis.state_activity.sum_max_in_NREM               = sum_max_in_nrem;
    end
end

