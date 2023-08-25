function output = swr_aligned_ecog(sData, z_score, spindle_select, swr_select, baseSub, plotR)

% Written by Christoffer Berge | Vervaeke Lab

% Compute mean SWR-aligned ECoG activity for different SWR-types. 

nr_of_seconds = 1;

% if length(varargin) > 1
%     prompt = sprintf('Plot raw signal (1) or z-score signal (2)');
%     sig_to_plot = input(prompt);
% end
sig_to_plot = z_score;

%% Select spindle freq

% BY DEFAULT SET TO 8-18 Hz. Uncomment code below to analyze detected
% 10-16Hz spindles instead. 

% prompt = sprintf('8-18Hz (1) or 10-16hz (2)? ');
% spindle_band_select = input(prompt);
% 
% if spindle_band_select == 1
%     spindle_select = [];
% elseif spindle_band_select == 2
%     spindle_select = '1016';
% end
% spindle_select                 = [];
unclassified_swr_str           = strcat('unclassified_swr', spindle_select);
NREM_spindle_uncoupled_swr_str = strcat('NREM_spindle_uncoupled_swr', spindle_select);
spindle_coupled_swr_str        = strcat('spindle_coupled_swr', spindle_select);
spin_idx_str                   = strcat('NREMAbsSpindleIdx', spindle_select);
nr_of_spindles_in_analysis     = length(sData.ephysdata2.(spin_idx_str) );

%% Select SWRs for analysis
% prompt = sprintf('All SWRs (1) | Rem. locom. SWRs (2) | Rem. clust. SWRs (3) | Rem- locom. & clust. SWRs (4) ');
% swr_for_analysis = input(prompt);
swr_for_analysis = swr_select;
RippleIdx = sData.ephysdata.absRipIdx;

if swr_for_analysis == 2
    [RippleIdx,~] = riprun2(sData, RippleIdx);
elseif swr_for_analysis == 3
    [RippleIdx] = RemoveRip(RippleIdx, [], nr_of_seconds);
elseif swr_for_analysis == 4
    [RippleIdx]   = RemoveRip(RippleIdx, [], nr_of_seconds);
    [RippleIdx,~] = riprun2(sData, RippleIdx);
end

% Find indicies of different SWR types
unclassified_SWR_idx = ismember(RippleIdx, sData.ephysdata.(unclassified_swr_str) );
single_swr_idx       = ismember(RippleIdx, sData.ephysdata.(NREM_spindle_uncoupled_swr_str));
spindle_coupled_idx  = ismember(RippleIdx, sData.ephysdata.(spindle_coupled_swr_str));

% Create variables for analysis
n_ripples        = length(RippleIdx);
RSC_ECoG         = sData.ephysdata2.lfp;
win_length       = (nr_of_seconds*2500*2)+1;
time             = (-(nr_of_seconds*2500):(nr_of_seconds*2500))./2500;
sessionID        = sData.sessionInfo.sessionID;

SWR_ECoG_signals = zeros(n_ripples, win_length);


% Specify whether to do baseline subtraction
% prompt  = sprintf('Do baseline subtraction? (1) | everything else = no) ');
% baseSub = input(prompt);
% if baseSub ~= 1
%     baseSub = [];
% end

%% Run main analysis

% Loop over nr of SWRs
for swr_nr = 1:n_ripples

    % find beginning and end of ECoG snippet centered on SWR
    ECoG_signal_start = RippleIdx(swr_nr) - (nr_of_seconds * 2500);
    ECoG_signal_end   = RippleIdx(swr_nr) + (nr_of_seconds * 2500);

    % skip SWRs occurring too early or late (i.e. the SWR-window
    % extends before or beyond recording length) 
    if ECoG_signal_start > 1 && ECoG_signal_end < length(RSC_ECoG)
        SWR_ECoG_signals(swr_nr,:) = RSC_ECoG(ECoG_signal_start:ECoG_signal_end);
    end

    % Baseline subtraction
    if baseSub == 1
        % Baseline = mean activity in -1 to -0.8 sec before ripple peak
        baseline_vm                = nanmean(SWR_ECoG_signals(swr_nr,1:200));
        SWR_ECoG_signals(swr_nr,:) = SWR_ECoG_signals(swr_nr,:)-baseline_vm;
    end
end
 
% Index results into different SWR-type matrices
signal_swr_awake           = SWR_ECoG_signals(unclassified_SWR_idx,:);
signal_swr_spindle_single  = SWR_ECoG_signals(single_swr_idx,:);
signal_swr_spindle_coupled = SWR_ECoG_signals(spindle_coupled_idx,:);

% Z-score normalize
signal_swr_awake_zscore   = zscore( signal_swr_awake, 0,2);
signal_swr_single_zscore  = zscore( signal_swr_spindle_single, 0,2);
signal_swr_coupled_zscore = zscore( signal_swr_spindle_coupled, 0,2);

% Save results in cell array 
output = {signal_swr_awake,          signal_swr_awake_zscore, ...
         signal_swr_spindle_single,  signal_swr_single_zscore, ...
         signal_swr_spindle_coupled, signal_swr_coupled_zscore};

%% Plot results
if ~isempty(plotR)
    plot_swr_aligned_ecog(sessionID, time ,sig_to_plot, output, baseSub)
end