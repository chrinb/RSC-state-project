function [sig_corr_mat, true_corr_coefs, n_rois] = significant_state_correlations(varargin)

% Written by Christoffer Berge | Vervaeke lab

% Compute significant correlations between cell activity and state vectors.
% Use a shuffle approach to determine significance 

% (All code until shuffle procedure is copied from "plot_session_activity".


sData  = varargin{1,1};
params = varargin{1,2};

%% Get signal data
[signal_to_plot, ~, pc_rois, in_rois] = get_roi_signals_from_sData(sData, params);

if strcmp(params.cell_type, 'pc')
    signal_to_plot = signal_to_plot{1,:};
    txt = 'Excitatory cells';
elseif strcmp(params.cell_type, 'in')
    signal_to_plot = signal_to_plot{2,:};
    txt = 'Inhibitory cells';
elseif strcmp(params.cell_type, 'axon')
    txt = 'Axons';
end

state_vectors_2p = get_state_logicals(sData);

%% Create hypnogram vector

% Initialize empty vector
hypnogram_vector = zeros(1, size(signal_to_plot,2) );

% Create hypnogram vector
hypnogram_vector(state_vectors_2p{1,1} == 1)  = 2;
hypnogram_vector(state_vectors_2p{1,2} == 1)  = 3;
hypnogram_vector(state_vectors_2p{1,3} == 1)  = 1;
hypnogram_vector(state_vectors_2p{1,4}  == 1) = 0;

%% State correlation and shuffle analysis
nr_of_shuffles     = 1000;
signal             = signal_to_plot;
n_rois             = size(signal,1);
true_corr_coefs    = zeros(n_rois, 4  );
true_corr_coefs_p  = zeros(n_rois, 4);

sig_corr_mat = zeros(size(true_corr_coefs));


for roi_nr = 1:size(signal)
    tic;
    roi_data = signal(roi_nr,:);
%     roi_data = roi_data./max(roi_data);

    for state_nr = 1:4
        current_state_vec = state_vectors_2p{1, state_nr};
        
        if ~isempty(current_state_vec)
            [temp_corr, temp_corr_p]       = corrcoef(current_state_vec, roi_data);
            true_corr_coefs(roi_nr, state_nr )  = temp_corr(2,1);
            true_corr_coefs_p(roi_nr, state_nr) = temp_corr_p(2,1);

        % Shuffle analysis
        
        % Preallocate
        roi_shuffles = repmat(roi_data, nr_of_shuffles, 1);
        % Matrix specifying a randow displacement for each row of the
        % ROI shuffle matrix
        shuffle_idx = randi(size(roi_shuffles,2),1,size(roi_shuffles,1));

        shuffled_event_activity = cell2mat( arrayfun(@(x) circshift(roi_shuffles(x,:), shuffle_idx(x) ), (1:numel(shuffle_idx))', 'UniformOutput', false));
        current_state_vec_rep  = double( repmat(current_state_vec, nr_of_shuffles,1));

        for i = 1:nr_of_shuffles
            shuf_corr{state_nr,i}  = corrcoef(shuffled_event_activity(i,:), double(current_state_vec));
            shuf_corr1{i} =  shuf_corr{state_nr, i}(2,1);
        end

%         shuffled_corr_sorted = sort( cell2mat(shuf_corr1));
%         if true_corr_coefs(roi_nr, :) > perctile(shuffle_corr_sorted, 99)
        
        threshold_higher = prctile( cell2mat(shuf_corr1), 95);
        threshold_lower  = prctile( cell2mat(shuf_corr1), 5);
        
        % Check of ROI correlation with current state exceeds either
        % threshold
        if true_corr_coefs(roi_nr, state_nr) < threshold_lower  || true_corr_coefs(roi_nr, state_nr) > threshold_higher
            sig_corr_mat(roi_nr, state_nr) = true;
        else
            sig_corr_mat(roi_nr, state_nr) = false;
        end
        
%         figure, histogram( cell2mat(shuf_corr1), 20)
%         xline(threshold_lower)
%         xline(threshold_higher)
        end
    end

t = toc;
fprintf('\n Calculating ROI state correlations in # %d in %.1f seconds', roi_nr,t) 
end
%% Sort
state_indicies = {'NREM', 'REM', 'AW', 'QW'};

if strcmp(params.beh_state, 'NREM')
    idx = 1;
    txt2 = 'NREM';
elseif strcmp(params.beh_state, 'REM')
    idx = 2;
    txt2 = 'REM';
elseif strcmp(params.beh_state, 'QW')
    idx = 4;
    txt2 = 'QW';
elseif strcmp(params.beh_state, 'AW')
    idx = 3;
    txt2 = 'AW';
end
[sorted_coefs, sort_coef_indx_all] = sort(true_corr_coefs,1, 'descend');
sorted_coefs_p_all                 = sig_corr_mat(sort_coef_indx_all);

coefs_to_analyze = sorted_coefs(:,idx);
sort_coef_indx   = sort_coef_indx_all(:, idx);
sorted_coefs_p   = sorted_coefs_p_all(:, idx);


%% Histogram 
figure,
histogram(coefs_to_analyze, 20)
set(gca, 'xlim', [-1 1])
xline(0, 'r--', 'LineWidth',1)
ylabel('Counts')
xlabel('Correlation coefficients')
title([txt, ' - ', txt2])
axis square


