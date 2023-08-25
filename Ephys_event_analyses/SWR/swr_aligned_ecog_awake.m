function mean_signal = swr_aligned_ecog_awake(varargin)

% Written by Christoffer Berge | Vervaeke Lab

% Calculate mean SWR-aligned ephys activity. This could for example be
% either ECoG or EMG signals. Works with awake and sleep sessions. 

sData = varargin{1,1};

% Load and assign variables
nr_of_seconds  = 3;
ECoG           = sData.ephysdata2.lfp;
EMG            = sData.ephysdata3.lfp;
ECoG_zscore    = zscore(ECoG);
EMG_zscore     = zscore(EMG);
win_length     = (nr_of_seconds*2500*2)+1;
time           = (-(nr_of_seconds*2500):(nr_of_seconds*2500))./2500;
sessionID      = sData.sessionInfo.sessionID;
select_swr_idx = sData.ephysdata.absRipIdx; 



% Select signal
if checkParameter(params.signal_type, 1 , 'principal cells')
    celltype = 'PC';
elseif checkParameter(opts.split_rois, 2 , 'inhibitory cells')
    celltype = 'IN';
else
%% Select SWRs for analysis
% Get the indicies of user specified SWR types
RippleIdx = get_swr_idx(opts.swr_for_analysis, sData, select_swr_idx, opts);
n_ripples = length(RippleIdx);


%% Run main analysis

% Preallocate
[SWR_ECoG_signals, SWR_ECoG_signals_zscore] = deal( zeros(n_ripples, win_length));

% Loop over nr of SWRs
for swr_nr = 1:n_ripples
    
    [SWR_ECoG_signals,SWR_ECoG_signals_zscore] = extract_avg_activity(RippleIdx, ECoG, ECoG_zscore,...
         nr_of_seconds, ECoG, opts.baseSub);
%     
%     % find beginning and end of ECoG snippet centered on SWR
%     ECoG_signal_start = RippleIdx(swr_nr) - (nr_of_seconds * 2500);
%     ECoG_signal_end = RippleIdx(swr_nr) + (nr_of_seconds * 2500);
% 
%     % if ECoG snippet is starting before or outlasts the actual ECoG signal
%     % from which it is taken (because SWR occurs at the very beginning or
%     % end), exclude that ECoG snippet
%     if ECoG_signal_start > 1 && ECoG_signal_end < length(ECoG)
%         SWR_ECoG_signals(swr_nr,:) = ECoG(ECoG_signal_start:ECoG_signal_end);
%     end
% 
%     % Baseline subtraction
%     if baseSub == 1
%         % Baseline = mean activity in -1 to -0.8 sec before ripple peak
%         baseline_sig = nanmean(SWR_ECoG_signals(swr_nr,1:200));
%         SWR_ECoG_signals(swr_nr,:) = SWR_ECoG_signals(swr_nr,:)-baseline_sig;
%     end

end
 
signal_swr_awake        = SWR_ECoG_signals;
signal_swr_awake_zscore = zscore( signal_swr_awake, 0,2);


% Save results as cell 
mean_signal = {signal_swr_awake, signal_swr_awake_zscore};

%% Select signal to plot

if length(varargin) > 1
    if sig_to_plot == 1
        awake_sig_to_plot   = signal_swr_awake;
        SE_awake   = std(awake_sig_to_plot, 'omitnan') ./ sqrt(size(awake_sig_to_plot,1));
       
    elseif sig_to_plot == 2
        awake_sig_to_plot   = signal_swr_awake_zscore;
        SE_awake   = std(awake_sig_to_plot, 'omitnan') ./ sqrt(size(awake_sig_to_plot,1));
    end
end
%% Plot results
if length(varargin) > 1
    figure,
    subplot(3,2,1:4),
    x1       = [time(1) time(end)];
    y1_awake = [1, size(awake_sig_to_plot,1)];
    imagesc(x1, y1_awake, awake_sig_to_plot)
    colorbar
%     if sig_to_plot == 2
%         caxis([-3 3])
%     end
    ylabel('SWR #')

    subplot(3,2, 5:6)
    shadedErrorBar(time, mean(awake_sig_to_plot,'omitnan'),SE_awake,'lineprops', 'r');
    xlabel('Time from ripple peak (sec)')
end
