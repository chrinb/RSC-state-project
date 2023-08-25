function SWR_ephys = swr_aligned_ephys(varargin)

% Written by Christoffer Berge | Vervaeke Lab

% Calculate mean SWR-aligned ephys activity. This could for example be
% either ECoG or EMG signals. Works with awake and sleep sessions. 

sData  = varargin{1,1};
params = varargin{1,2};

% Load and assign variables
nr_of_seconds  = 3;
win_length     = (nr_of_seconds*2500*2)+1;
time           = (-(nr_of_seconds*2500):(nr_of_seconds*2500))./2500;
sessionID      = sData.sessionInfo.sessionID;
select_swr_idx = sData.ephysdata.absRipIdx; 

% Select signal
if strcmp(params.signal_type, 'ECoG')
    signal = sData.ephysdata2.lfp;
elseif strcmp(params.signal_type, 'EMG')
    signal = sData.ephysdata3.lfp;
elseif strcmp(params.signal_type,'Running speed')
    signal = sData.daqdata.runSpeed;
end

% Z-score signal (not running speed)
if strcmp(params.signal_type,'Running speed')
    signal_zscore = signal;
else
    signal_zscore = zscore(signal);
end

%% Select SWRs for analysis
% Get the indicies of user specified SWR types
RippleIdx = get_swr_idx(params.swr_for_analysis, sData, select_swr_idx, params);
n_ripples = length(RippleIdx);

% If sleep session, get indicies of different SWR-types
if strcmp(params.beh_state, 'sleep')
    unclassified_SWR_idx = ismember(RippleIdx, sData.ephysdata.unclassified_swr);
    single_swr_idx       = ismember(RippleIdx, sData.ephysdata.NREM_spindle_uncoupled_swr);
    spindle_coupled_idx  = ismember(RippleIdx, sData.ephysdata.spindle_coupled_swr);
end
%% Run main analysis

% Preallocate
[SWR_ephys_signals, SWR_ephys_signals_zscore] = deal( zeros(n_ripples, win_length));


% Loop over nr of SWRs
for swr_nr = 1:n_ripples
    
    % Extract peri-SWR ephys data
    [SWR_ephys_signals,SWR_ephys_signals_zscore] = extract_avg_activity(RippleIdx, signal, signal_zscore,...
         SWR_ephys_signals, SWR_ephys_signals_zscore, nr_of_seconds, params);
end
 
% If sleep session, split SWR-aligned data into subcategories.
if strcmp(params.beh_state, 'sleep')
    SWR_ephys = cell(2,3);

    SWR_ephys{1,1} = SWR_ephys_signals(unclassified_SWR_idx,:);
    SWR_ephys{2,1} = SWR_ephys_signals_zscore(unclassified_SWR_idx,:);
    SWR_ephys{1,2} = SWR_ephys_signals(single_swr_idx,:);
    SWR_ephys{2,2} = SWR_ephys_signals_zscore(single_swr_idx,:);
    SWR_ephys{1,3} = SWR_ephys_signals(spindle_coupled_idx,:);
    SWR_ephys{2,3} = SWR_ephys_signals_zscore(spindle_coupled_idx,:);
else
    SWR_ephys = cell(2,1);
    SWR_ephys{1,1} = SWR_ephys_signals;
    SWR_ephys{2,1} = SWR_ephys_signals_zscore;
end

%% Select signal to plot

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
