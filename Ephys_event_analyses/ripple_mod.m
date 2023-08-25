function output = ripple_mod(varargin)

% Written by Christoffer Berge | Vervaeke Lab

% Function that finds SWR-modulated ROIs during either (1) all SWRs, SWRs
% classified as (2) awake, (3) spindle-coupled, or (4) spindle-uncoupled.
% User specifies which of the four options the shuffle procedure operates on
% to find SWR-modulated ROIs. 

sData                   = varargin{1,1};
opts.swr_for_analysis   = varargin{1,2};
opts.swr_cluster_thresh = varargin{1,3};
%% Select signal for shuffle analysis
prompt = sprintf('Which signal? (1 = deconvolved | 2 = DF/F) ');
signalSelection = input(prompt); 

% select signal for analysis, and additional signal to plot just for
% visualization
if signalSelection == 1
    signal_for_an  = sData.imdata.roiSignals(2).ciaDeconvolved;
    signal_plot    = sData.imdata.roiSignals(2).newdff;
    label1  = 'Mean deconvolved dF/F';
    label2  = 'Mean dF/F';
elseif signalSelection == 2
    signal_for_an = sData.imdata.roiSignals(2).newdff;
    signal_plot   = sData.imdata.roiSignals(2).ciaDeconvolved;
    label1   = 'Mean dF/F';
    label2  = 'Mean deconvolved dF/F';
end


% Add roisignals folder to path to get roi array
try
    folder_name = dir;
    addpath( [folder_name(1).folder '\roisignals'])
catch
end


% Length of the window before/after SWR peak
nr_of_seconds = 3;

frames    = sData.daqdata.frame_onset_reference_frame;
sessionID = sData.sessionInfo.sessionID;

% nr of shuffle iterations
nr_of_shuffles = 1000; 


%% Select spindle freq
% prompt = sprintf('8-18Hz (1) or 10-16hz (2)? ');
% spindle_band_select = input(prompt);
% 
% if spindle_band_select == 1
%     spindle_select = [];
% elseif spindle_band_select == 2
%     spindle_select = '1016';
% end

% awake_swr                  = strcat('awake_swr', spindle_select);
% NREM_spindle_uncoupled_swr = strcat('NREM_spindle_uncoupled_swr', spindle_select);
% spindle_coupled_swr        = strcat('spindle_coupled_swr', spindle_select);
%% Select SWRs for analysis

select_swr_idx = sData.ephysdata.absRipIdx; 

% Get the indicies of user specified SWR types
RippleIdx = get_swr_idx(opts.swr_for_analysis ,sData,select_swr_idx, opts);
% prompt = sprintf('Select SWRs to analyze: 1 = All | 2 = Awake | 3 Spindle-coupled | 4 Spindle-Uncoupled ');
% swrs_to_analyze = input(prompt);
% 
% if swrs_to_analyze == 1 
%     swr_idx = sData.ephysdata.absRipIdx;
% elseif swrs_to_analyze == 2
%     swr_idx = sData.ephysdata.(awake_swr);
% elseif swrs_to_analyze == 3 
%     swr_idx = sData.ephysdata.(NREM_spindle_uncoupled_swr);
% elseif swrs_to_analyze == 4 
%     swr_idx = sData.ephysdata.(spindle_coupled_swr);
% end
% 
% prompt = sprintf('Remove locomotion SWR? (1 = yes) ');
% rem_run_swr = input(prompt);
%     
% prompt = sprintf('Remove temporally close SWR? (1 = yes) ');
% rem_cluster_swr = input(prompt);
%     
% if ~isempty(rem_run_swr) && isempty(rem_cluster_swr)
%     [swr_idx,~] = riprun2(sData, swr_idx); 
% elseif isempty(rem_run_swr) && ~isempty(rem_cluster_swr)
%     swr_idx = RemoveRip(swr_idx,[]);
% elseif ~isempty(rem_run_swr) && ~isempty(rem_cluster_swr)
%    [swr_idx,~] = riprun2(sData, swr_idx); 
%    swr_idx = RemoveRip(swr_idx,[]);
% end

%% Specify whether to threshold deconvolved dF/F or not
% prompt = sprintf('Threshold deconvolved dF/F? (y = yes | everything else = no) ');
% dothreshold = input(prompt, 's');
% if strcmp(dothreshold, 'y')
%     threshold = true(1,1);
% else
%     threshold = false(1,1);
% end

%% Initialize various variables
win_length           = (nr_of_seconds*31)*2+1;
All_ROIs             = zeros(nr_of_shuffles,1);
mean_swr_activity    = zeros( size(signal_for_an,1),win_length);
% mean_activity_all_shuffle_iterations = zeros(nr_of_shuffles, win_length);
swr_activity         = zeros( length(RippleIdx), win_length);
% threshold_percentile = zeros( size(signal_for_an,1),1);
time                 = linspace(-nr_of_seconds,nr_of_seconds,win_length);
swr_idx              = frames(round(RippleIdx));
time_win             = (62:108); % -1s to +500ms after SWR peak
min_nr_events        = 5;
%% Shuffle analysis

% check if any of the SWR categories are empty
% awake_check     = ~isempty(awakeSWRidx);
% spindleUC_check = ~isempty(NREMspindleUncoupledSWRidx);
% spindleC_check  = ~isempty(NREMspindleCoupledSWRidx);
% 
% n_swr_cat      = awake_check+spindleUC_check+spindleC_check;
% test_cell      = {awakeSWRidx, NREMspindleUncoupledSWRidx,NREMspindleCoupledSWRidx};
% log_idx        = [awake_check, spindleUC_check, spindleC_check];
% swr_to_analyse = test_cell(log_idx);

% Loop over different SWR categories
% for i = 1:n_swr_cat
%     swr_cell  = swr_to_analyse(i);
%     swr_idx   = frames(round(swr_cell{1,1}));
cleaned_rois = remove_cells();

% Loop over nr of ROIs
for roinr = cleaned_rois

    % Select the activity for a given ROI
    ROI_signal = signal_for_an(roinr, :);

    % Loop over SWRs and extract 2P segments
    for swr_nr = 1:length(swr_idx)
        swr_window_start = swr_idx(swr_nr) - (nr_of_seconds*31); 
        swr_window_end   = swr_idx(swr_nr) + (nr_of_seconds*31); 

        % skip SWRs at the beginning or end with a window shorter than
        % length specified in seconds by user above
        if swr_window_start > 1 && swr_window_end < length(signal_for_an)

%             if length(swr_window_start:swr_window_end) == win_length
            swr_activity(swr_nr, :) = ...
                ROI_signal(swr_window_start:swr_window_end); 
        end
    end

    % Threshold signal using K-means clustering to separate "noise" vs
    % true deconvolved events.
    if signalSelection == 1 && threshold == 1 
        [threshold_percentile, ~] = threshold_deconvolved(signal_for_an, ROI_signal, roinr);
    % Remove deconvolved values below threshold
    ROI_idx = reshape(swr_activity, 1,[]);
    ROI_idx(ROI_idx < threshold_percentile(roinr)) = 0;
    swr_activity = reshape(ROI_idx, size(swr_activity));
    end

    mean_swr_activity(roinr,:) = nanmean(swr_activity); 
    mean_zscore_swr_activity   = zscore(mean_swr_activity);
    DeconvWeightedAvg(roinr)   = {swr_activity};

    % shuffle spindle activity
    All_ROIs = shuffle_analysis(win_length,...
        nr_of_shuffles, swr_activity, roinr, All_ROIs, mean_swr_activity,...
        time_win, min_nr_events);
end
    
%     swr_mod_ROIs{i} = nonzeros(All_ROIs);
%     if length(swr_idx) == length(awakeSWRidx)
%         awake_swr_all_ROIs = nonzeros(All_ROIs);
%     elseif length(swr_idx) == length(NREMspindleCoupledSWRidx)
%         spindleC_swr_all_ROIs = nonzeros(All_ROIs);
%     elseif length(swr_idx) == length(NREMspindleUncoupledSWRidx)
%         spindleUC_swr_all_ROIs = nonzeros(All_ROIs);
%     end
% end

%% Find activated and suppressed ROIs

% swr_mod_ROIs = {awake_swr_all_ROIs,spindleC_swr_all_ROIs,spindleUC_swr_all_ROIs };
% mod_ROI_to_analyse = swr_mod_ROIs(log_idx);
% zscore_swr_mod_neurons = mean_zscore_swr_activity(awake_swr_mod_ROI_nr,:);
% for i = 1:n_swr_cat
%     swr_cell  = mod_ROI_to_analyse(i);

swr_modulated_ROI_nr = nonzeros(All_ROIs);

ROIs_activated = [];
ROIs_suppressed = [];
ROIs_unclassified = [];
m = 1;
n = 1;
p = 1;
% loop over modulated ROIs and classify as activated or suppressed based on
% comparing mean activity in a baseline window (-3s to -2s) vs a test
% period (-1s to 0s before SWR peak)
baseline_window = (1:31);
test_window     = (86:101);

for mod_neurons = 1:length( swr_modulated_ROI_nr )
    neuron_n = swr_modulated_ROI_nr(mod_neurons);
    if mean( mean_swr_activity(neuron_n,baseline_window)) > mean( mean_swr_activity(neuron_n,test_window))
        ROIs_suppressed(m,:) = mean_swr_activity( neuron_n,:);
        m = m+1;
    elseif mean( mean_swr_activity(neuron_n,baseline_window)) < mean( mean_swr_activity(neuron_n, test_window))
        ROIs_activated(n,:) = mean_swr_activity( neuron_n,:);
        n = n+1;
    else
        ROIs_unclassified(p,:) = mean_swr_activity( neuron_n,:);
        p = p+1;
    end
end
    
%% Find the indicies of SWR modulated ROIs
ROIs_activated_idx = zeros(size(ROIs_activated,1),1);
% loop over nr of activated ROIs
for j = 1:size(ROIs_activated,1)
    % find the indicies of the ROIs that match the activated ROI signal.
    temp_var      = sum( mean_swr_activity == ROIs_activated(j, :),2);
    ROI_idx       = temp_var == win_length;
    ROI_list      = (1:roinr)'; 
    ROIs_activated_idx(j) = ROI_list(ROI_idx);
end

%list the identity of SWR-inhibited ROIs
ROIs_suppressed_idx = zeros(size(ROIs_suppressed,1),1);
for j = 1:size(ROIs_suppressed,1)
    temp_var      = sum( mean_swr_activity == ROIs_suppressed(j, :),2);
    ROI_idx       = temp_var == win_length;
    ROI_list      = (1:roinr)'; 
    ROIs_suppressed_idx(j) = ROI_list(ROI_idx);
end

%list the identity of unclassified ROIs
ROIs_unclassified_idx = zeros(size(ROIs_unclassified,1),1);
for j = 1:size(ROIs_unclassified,1)
    temp_var      = sum( mean_swr_activity == ROIs_unclassified(j, :),2);
    ROI_idx       = temp_var == win_length;
    ROI_list      = (1:roinr)'; 
    ROIs_unclassified_idx(j) = ROI_list(ROI_idx);
end
%     activity_mod_ROIs{i} = {suppressed_ROIs, activated_ROIs, unclassified_ROIs};
%     idx_mod_ROIs{i} = {idx_activated_ROIs,idx_suppressed_ROIs};
%     clearvars suppressed_ROIs activated_ROIs unclassified_ROIs n m p...
%         idx_activated_ROIs idx_suppressed_ROIs
% end
%% Find corresponding dF/F / deconvolved dF/F signals for plotting

% first find the swr segments for activated ROIs
activated_sig_to_plot = zeros(length(ROIs_activated_idx), win_length);
if ~isempty(ROIs_activated_idx) 
    % same code to extract 2P swr segments as above
    for i = 1:length(ROIs_activated_idx)
        mod_ROI = signal_plot(ROIs_activated_idx(i),:);
        for swr_nr = 1:length(swr_idx)
            swr_window_start = swr_idx(swr_nr) - (nr_of_seconds*31); 
            swr_window_end = swr_idx(swr_nr) + (nr_of_seconds*31); 
            if swr_window_start > 1 && swr_window_end < length(signal_for_an)
                swr_activity(swr_nr, :) = ...
                    mod_ROI(swr_window_start:swr_window_end); 
            end
        end
        activated_sig_to_plot(i,:) = nanmean(swr_activity);
    end
end

% ... and now for the spindle suppressed ROis
suppressed_sig_to_plot = zeros(length(ROIs_suppressed_idx), win_length);
if ~isempty(ROIs_suppressed_idx) 
    for i = 1:length(ROIs_suppressed_idx)
        mod_ROI = signal_plot(ROIs_suppressed_idx(i),:);
        for swr_nr = 1:length(swr_idx)
            swr_window_start = swr_idx(swr_nr) - (nr_of_seconds*31); 
            swr_window_end = swr_idx(swr_nr) + (nr_of_seconds*31); 
            if swr_window_start > 1 && swr_window_end < length(signal_for_an)
                swr_activity(swr_nr, :) = ...
                    mod_ROI(swr_window_start:swr_window_end); 
            end
        end
        suppressed_sig_to_plot(i,:) = nanmean(swr_activity);
    end
end

output = {ROIs_activated, ROIs_activated_idx, ROIs_suppressed,ROIs_suppressed_idx, ...
    activated_sig_to_plot, suppressed_sig_to_plot, ROIs_unclassified, ...
    ROIs_unclassified_idx};
%% Plotting

figure, 
if ~isempty(ROIs_activated) && size(ROIs_activated,1) > 1
    subplot(241)
    plot(time,ROIs_activated' + (size(ROIs_activated,1):-1:1)*.005)
    set(gca, 'xlim', [time(1), time(end)])

    subplot(245)
    shadedErrorBar(time, nanmean(ROIs_activated), ...
        std(ROIs_activated)./sqrt(size(ROIs_activated,1)) ,'lineprops', 'b');
    xlabel('Time from SWR peak (sec)')
    ylabel(label1)
    set(gca, 'xlim',[min(time) max(time)])

    subplot(242)
    plot(time,activated_sig_to_plot' + (size(activated_sig_to_plot,1):-1:1)*.03)
    set(gca, 'xlim', [time(1), time(end)])
    
    subplot(246)
    shadedErrorBar(time, nanmean(activated_sig_to_plot), ...
        std(activated_sig_to_plot)./sqrt(size(activated_sig_to_plot,1)) ,'lineprops', 'b');
    xlabel('Time from SWR peak (sec)')
    ylabel(label1)
    set(gca, 'xlim',[min(time) max(time)])
end

if ~isempty(ROIs_suppressed) && size(ROIs_suppressed,1) > 1
    subplot(243)
    plot(time,ROIs_suppressed' + (size(ROIs_suppressed,1):-1:1)*.005)
    set(gca, 'xlim', [time(1), time(end)])

    subplot(244)
    plot(time,suppressed_sig_to_plot' + (size(suppressed_sig_to_plot,1):-1:1)*.03)
    set(gca, 'xlim', [time(1), time(end)])
 
    subplot(247)
    shadedErrorBar(time, nanmean(ROIs_suppressed), ...
        std(ROIs_suppressed)./sqrt(size(ROIs_suppressed,1)) ,'lineprops', 'b');
    xlabel('Time from SWR peak (sec)')
    ylabel(label1)
    set(gca, 'xlim',[min(time) max(time)])

    subplot(248)
    shadedErrorBar(time, nanmean(suppressed_sig_to_plot), ...
        std(suppressed_sig_to_plot)./sqrt(size(suppressed_sig_to_plot,1)) ,'lineprops', 'b');
    xlabel('Time from SWR peak (sec)')
    ylabel(label2)
    set(gca, 'xlim',[min(time) max(time)])
end





