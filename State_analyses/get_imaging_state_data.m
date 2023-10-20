function [state_data, pc_rois, in_rois] = get_imaging_state_data(sData, params)

% Written by Christoffer Berge | Vervaeke lab

%{
THIS FUNCTION NEEDS FIX: it SHOULD extract 2P for each state in a
session, and exclude episodes that are (1) ongoing at recording start, 
(2) continues past recording end, or (3) are < certain duration threshold.
%}

% Function that extracts the indicies of different states, excluding states
% that happens in the very beginning or end of the recording (as these most
% likely are ongoing states).

% Get imaging data
[signal, ~, pc_rois, in_rois] = get_roi_signals_from_sData(sData, params);

if strcmp(params.zscore, 'yes')
    signal = zscore(signal{1,1}, [], 2);
else
    signal = signal{1,1};
end

% Get 2P frame rate
imaging_sampling_rate = find_imaging_framerate(sData);

% Get start/stop times for NREM, REM, QW, and AW (in that order!)
all_state_times = get_state_start_stop_times(sData);

signal_end = size(sData.imdata.roiSignals(2).newdff,2);

% Set threshold from recording start/end
threshold       = 3; % seconds
threshold_start = imaging_sampling_rate*threshold;
threshold_stop  = signal_end - (imaging_sampling_rate*threshold);

% Minimum duration in seconds for AW, QW, NREM, and REM episodes
min_ep_duration = [1 1 20 30];

% Loop over the four different state  (AW, QW, NREM, REM)
for state_nr = 1:4
    
    % Get start/stop times for current state
    temp_state_times = all_state_times{state_nr, 1 };

    % Check if current state contains data
    if any( sum( temp_state_times))
        
        % Loop over episodes for current state
        for state_ep_nr = 1:size(temp_state_times,1)
            
            % Get start/stop time for current episode
            temp_ep_times = [temp_state_times(state_ep_nr, 1), temp_state_times(state_ep_nr, 2)];
            % Check if episode fits criteria
            ep_fits_criteria = @(x, y, z, a, b) x(1) > y && x(2) < z && numel(x(1):x(2)) > a(b);

            if ep_fits_criteria(temp_ep_times, threshold_start, threshold_stop, min_ep_duration, state_nr )
                
                temp_state_snippet_times{state_ep_nr,1} = temp_ep_times(1):temp_ep_times(2);
                temp_state_snippet{state_ep_nr, 1}      = signal(:, temp_state_snippet_times{state_ep_nr,1});
            else
                temp_state_snippet_times{state_ep_nr,1} = [];
                temp_state_snippet{state_ep_nr,1}       = [];
            end
        end

        % Put adjusted state vectors into new cell array
        state_times_2use{state_nr,1}    = temp_state_snippet_times;
        state_snippets_2use{state_nr,1} = temp_state_snippet ;
    else
        state_times_2use{state_nr,1}    = [];
        state_snippets_2use{state_nr,1} = [];
    end
    clear temp_state_snippet_times temp_state_snippet
end

% Store output in struct
state_data = struct();
state_data.all_state_times    = state_times_2use;
state_data.all_state_snippets = state_snippets_2use;
% 
% % Convert state logical vector from ephys to 2P time (NOTE: states are
% % changed: 1 = NREM, 2 = REM, 3 = active, 4 = quiet)
% state_vectors_2p = get_state_logicals(sData, all_state_vectors_adjusted);
% 
% %% Select signal and cell type data to be extracted
% 
% % NOTE: this function will remove z-drift affected ROIs!!
% [pc_rois, in_rois] = remove_cells(sData);
% 
% % Check if session is multi-day recording and only select cells present on
% % current day 
% % if isfield(sData.imdata, 'roi_classification')
% %     pc_roi_idx = sData.imdata.roi_classification(pc_rois); % of all ROIs, index out cell type
% %     in_roi_idx = sData.imdata.roi_classification(in_rois); % of all ROIs, index out cell type
% %     log_idx_pc = pc_roi_idx == 1;
% %     log_idx_in = in_roi_idx == 1;
% %     pc_rois    = pc_rois(log_idx_pc);
% %     in_rois    = in_rois(log_idx_in);
% % end
% 
% % For now, extract ALL the data from a cell type (exc/in) and later set ROI
% % data where cell is not present in that session to NaNs.
% 
% %%  Load data
% 
% % Get signal data
% if strcmp(params.signal_type, 'dff')
%     dff = sData.imdata.roiSignals(2).newdff;
%     if strcmp(params.cell_type, 'axon')
%         dff = sData.imdata.roiSignals(2).mergedAxonsDff;
%     end
% elseif strcmp(params.signal_type, 'deconv')
%     dec = sData.imdata.roiSignals(2).ciaDeconvolved;
%     if strcmp(params.cell_type, 'axon')
%         dec = sData.imdata.roiSignals(2).mergedAxonsDec;
%     end
% end
% % if nargin > 1 
% %     switch varargin{1,2}.cell_type
% %         case 'axon'
% %              dff = sData.imdata.roiSignals(2).mergedAxonsDffFilt;
% %              dec = sData.imdata.roiSignals(2).mergedAxonsDec;
% %         case 'in'
% %             dff = sData.imdata.roiSignals(2).newdff(in_rois,:);
% %             dec = sData.imdata.roiSignals(2).ciaDeconvolved(in_rois,:);
% %         case 'pc'
% %             dff = sData.imdata.roiSignals(2).newdff(pc_rois,:);
% %             dec = sData.imdata.roiSignals(2).ciaDeconvolved(pc_rois,:);
% %     end
% % else
% %     dff = sData.imdata.roiSignals(2).newdff;
% % end
% %% Loop over the epochs of different states within a session and extract imaging data 
% 
% % Loop over nr of states
% for n_state = 1:4
% 
%     current_state_2p_vector = state_vectors_2p{1,n_state};
%     
%     % Check that sum of vector is not zero (meaning it's not present in
%     % the session)
%     if any( sum( current_state_2p_vector))
% 
%         % Find the start and stop indicies of the state vector
%         [eventStartIdx, eventStopIdx ]  = findTransitions( current_state_2p_vector);
% 
%         % Loop over nr of bouts and for each bout extract the DF/F activity of
%         % that bout and store in cell array
%         for n_bouts = 1:numel(eventStartIdx)
%             temp = false(1, length(dff));
%             temp(eventStartIdx(n_bouts):eventStopIdx(n_bouts)) = true;
%             state_data_cell{n_bouts,1} = dff(:, temp);
%         end
% 
%         % Store un-padded state activity
%         all_data{n_state} = state_data_cell;
%         clear state_data_cell
% 
% %         % Loop over active bouts and find longest bout
% %         for i = 1:size(test_cell)
% %             state_cell_size(i) = size(test_cell{i,1},2);
% %         end
% %         largest_bout = max(state_cell_size);
% %         
% %         % Loop over active bouts and pad bouts shorter than the longest one with NaNs
% %         state_padded = test_cell;
% %         for i = 1:size(test_cell)
% %             state_padded{i,1}(:, state_cell_size(i)+1:largest_bout) = NaN;
% %         end
% %     % Convert cell array to matrix
% %     padded_mat{n_state} = cell2mat(state_padded);
% % 
% %     clear test_cell state_cell_size state_padded
% %     else
% %         padded_mat{n_state} = NaN;
% 
%     end
% 
% end
% 
% %% Average over individual cells 
% 
% % n_rois = size(dff,1);
% % 
% % % for roi_nr 1:n_rois
% % 
% % roi1 = all_data{1,1}{1,1}(1,:);
% % 
% % roi2 = all_data{1,1}{1,2}(1,:);
% % 
% % roi3 = all_data{1,1}{1,3}(1,:);
% % 
% % 
% % data_max = size(roi1,2);
% % 
% % roi2_i = interp1(1:size(roi2,2), roi2, linspace(1, size(roi2,2), data_max) );
% % roi3_i = interp1(1:size(roi3,2), roi3, linspace(1, size(roi3,2), data_max) );
% % 
% % test3 = interp1(1:4, a, linspace(1, 4, 7));
% 
% %% Loop over bouts in different states, find longest, pad shorter bout with NaNs
% % active_DFF = dff(:, output{1,3});
% % quiet_DFF  = dff(:, output{1,4});
% % NREM_DFF   = dff(:, output{1,1});
% % REM_DFF    = dff(:, output{1,2});
