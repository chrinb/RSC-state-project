function SWR_ephys = swr_aligned_ephys(varargin)

% Written by Christoffer Berge | Vervaeke Lab

% Calculate mean SWR-aligned ephys activity. This could for example be
% either ECoG or EMG signals. Works with awake and sleep sessions. 

sData     = varargin{1,1};
params    = varargin{1,2};
sleep_idx = varargin{1,3};

% Load and assign variables
nr_of_seconds  = 1;
win_length     = (nr_of_seconds*2500*2)+1;
time           = (-(nr_of_seconds*2500):(nr_of_seconds*2500))./2500;
sessionID      = sData.sessionInfo.sessionID;
select_swr_idx = sData.ephysdata.absRipIdx; 

% Select signal
if strcmp(params.signal_type, 'ECoG')
    signal = sData.ephysdata2.lfp';
elseif strcmp(params.signal_type, 'EMG')
    signal = sData.ephysdata3.lfp';
elseif strcmp(params.signal_type,'Running speed')
    signal = sData.daqdata.runSpeed';
elseif strcmp(params.signal_type,'LFP')
    signal = sData.ephysdata.lfp';
end

% Z-score signal (not running speed)
if strcmp(params.zscore, 'yes') && ~strcmp(params.signal_type, 'Running speed')
    signal = zscore(signal);
end

%% Select SWRs for analysis

% If sleep session, get indicies of different SWR-types
if sleep_idx == 1
    select_swr_idx = sort([sData.ephysdata.NREM_spindle_uncoupled_swr, sData.ephysdata.spindle_coupled_swr]);
end

% Get the indicies of user specified SWR types
event_idx = get_swr_idx(params.swr_for_analysis, sData, select_swr_idx, params);


n_ripples = length(event_idx);

% If sleep session, get indicies of different SWR-types
% if strcmp(params.beh_state, 'sleep')
% if sleep_idx == 1
%     unclassified_SWR_idx = ismember(event_idx, sData.ephysdata.unclassified_swr);
%     single_swr_idx       = ismember(event_idx, sData.ephysdata.NREM_spindle_uncoupled_swr);
%     spindle_coupled_idx  = ismember(event_idx, sData.ephysdata.spindle_coupled_swr);
% 
%     event_idx = sort(single_swr_idx, spindle_coupled_idx);
% end
%% Run main analysis

% Preallocate
SWR_aligned_data = deal( zeros(n_ripples, win_length));


% Loop over nr of SWRs
for swr_nr = 1:n_ripples
    
    % Extract peri-SWR ephys data
    SWR_aligned_data = extract_avg_activity(event_idx, signal, SWR_aligned_data, nr_of_seconds, params);
end
 
% If sleep session, split SWR-aligned data into subcategories.
% if strcmp(params.beh_state, 'sleep')
%     SWR_ephys = cell(2,3);
% 
%     SWR_ephys{1,1} = SWR_aligned_data(unclassified_SWR_idx,:);
%     SWR_ephys{1,2} = SWR_aligned_data(single_swr_idx,:);
%     SWR_ephys{1,3} = SWR_aligned_data(spindle_coupled_idx,:);
% else
    SWR_ephys = cell(2,1);
    SWR_ephys = SWR_aligned_data;
% end

%% Select signal to plot
% 
% if length(varargin) > 1
%     if sig_to_plot == 1
%         awake_sig_to_plot   = signal_swr_awake;
%         SE_awake   = std(awake_sig_to_plot, 'omitnan') ./ sqrt(size(awake_sig_to_plot,1));
% 
%     elseif sig_to_plot == 2
%         awake_sig_to_plot   = signal_swr_awake_zscore;
%         SE_awake   = std(awake_sig_to_plot, 'omitnan') ./ sqrt(size(awake_sig_to_plot,1));
%     end
% end
% %% Plot results
% if length(varargin) > 1
%     figure,
%     subplot(3,2,1:4),
%     x1       = [time(1) time(end)];
%     y1_awake = [1, size(awake_sig_to_plot,1)];
%     imagesc(x1, y1_awake, awake_sig_to_plot)
%     colorbar
% %     if sig_to_plot == 2
% %         caxis([-3 3])
% %     end
%     ylabel('SWR #')
% 
%     subplot(3,2, 5:6)
%     shadedErrorBar(time, mean(awake_sig_to_plot,'omitnan'),SE_awake,'lineprops', 'r');
%     xlabel('Time from ripple peak (sec)')
% end
