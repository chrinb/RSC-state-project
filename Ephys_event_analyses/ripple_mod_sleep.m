function [activated_ROIs, idx_activated_ROIs, suppressed_ROIs,idx_suppressed_ROIs, ...
    activated_sig_to_plot, suppressed_sig_to_plot] = ripple_mod_sleep(sData)

% Written by Christoffer Berge | Vervaeke Lab

% Function that finds SWR-modulated ROIs

%% Select signal for shuffle analysis
prompt = sprintf('Which signal? (1 = deconvolved | 2 = DF/F ');
signalSelection = input(prompt); 

% select signal for analysis, and additional signal to plot just for
% visualization
if signalSelection == 1
    signal_for_an  = sData.imdata.roiSignals(2).ciaDeconvolved;
    signal_plot    = sData.imdata.roiSignals(2).newdff;
    label1   = 'Mean deconvolved dF/F';
    label2  = 'Mean dF/F';
elseif signalSelection == 2
    signal_for_an = sData.imdata.roiSignals(2).newdff;
    signal_plot   = sData.imdata.roiSignals(2).ciaDeconvolved;
    label1   = 'Mean dF/F';
    label2  = 'Mean deconvolved dF/F';
end

% Length of the window before/after SWR peak
nr_of_seconds = 3;

frames      = sData.daqdata.frame_onset_reference_frame;
sessionID   = sData.sessionInfo.sessionID;


spindleIdxStart = frames(round(sData.ephysdata2.spindleStartEnd(:,1)));

% nr of shuffle iterations
nr_of_shuffles = 1000; 

%% Select SWRs for analysis
prompt = sprintf('All ripples? (y = yes | everything else = no) ');
allrip = input(prompt,'s');

if strcmp(allrip,'y') %keep all ripples
    awakeSWRidx = sData.ephysdata.awake_swr;
    NREMspindleUncoupledSWRidx = sData.ephysdata.NREM_spindle_uncoupled_swr;
    NREMspindleCoupledSWRidx = sData.ephysdata.spindle_coupled_swr;
    
else
    prompt = sprintf('Remove locomotion SWR? (y = yes | everything else = no) ');
    riprun = input(prompt, 's');
    
    prompt = sprintf('Remove temporally close SWR? (y = yes | everything else = no) ');
    removerip = input(prompt, 's');
    
    % if remove locomotion SWR but not temporally close SWR
    if strcmp(riprun, 'y') && ~strcmp(removerip, 'y')
        [awakeSWRidx, NREMspindleUncoupledSWRidx, NREMspindleCoupledSWRidx ] = ripRunAn(sData); 
    % if remove temporally close SWR but not locomotion SWR
    elseif strcmp(removerip, 'y') && ~strcmp(riprun, 'y')
        [awakeSWRidx, NREMspindleUncoupledSWRidx, NREMspindleCoupledSWRidx] = removeCloseRip(sData);

    % if remove both temporally close and locotion SWR
    elseif strcmp(removerip, 'y') && strcmp(riprun, 'y')
        [awakeSWRidx, NREMspindleUncoupledSWRidx, NREMspindleCoupledSWRidx ] = ripRunAn(sData,1);
    end
end

prompt = sprintf('Which SWR type? 1 = awake | 2 spindle-coupled | 3 = spindle-uncoupled ');
swr_type_select = input(prompt);
if swr_type_select == 1
    swr_idx = awakeSWRidx;
elseif swr_type_select == 2
    swr_idx = NREMspindleCoupledSWRidx;
elseif swr_type_select == 3
    swr_idx = NREMspindleUncoupledSWRidx;
end


%% Specify whether to threshold deconvolved dF/F or not
prompt = sprintf('Threshold deconvolved dF/F? (y = yes | everything else = no) ');
dothreshold = input(prompt, 's');
if strcmp(dothreshold, 'y')
    threshold = true(1,1);
else
    threshold = false(1,1);
end

%% Initialize various variables
win_length           = (nr_of_seconds*31)*2+1;
All_ROIs             = zeros(nr_of_shuffles,1);
mean_swr_activity    = zeros( size(signal_for_an,1),win_length);
% mean_activity_all_shuffle_iterations = zeros(nr_of_shuffles, win_length);
swr_activity         = zeros( length(swr_idx), win_length);
threshold_percentile = zeros( size(signal_for_an,1),1);
time                 = linspace(-nr_of_seconds,nr_of_seconds,win_length);
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
    
    % Loop over nr of ROIs
for roinr = 1:size(signal_for_an,1)

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
    end

    % Remove deconvolved values below threshold
    ROI_idx = reshape(swr_activity, 1,[]);
    ROI_idx(ROI_idx < threshold_percentile(roinr)) = 0;
    swr_activity = reshape(ROI_idx, size(swr_activity));

    mean_swr_activity(roinr,:) = nanmean(swr_activity); 
    mean_zscore_swr_activity   = zscore(mean_swr_activity);
    DeconvWeightedAvg(roinr)   = {swr_activity};

    % shuffle spindle activity
    All_ROIs = shuffle_analysis(win_length,...
        nr_of_shuffles, swr_activity, roinr, All_ROIs, mean_swr_activity);
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

activated_ROIs = [];
suppressed_ROIs = [];
m = 1;
n = 1;
p = 1;
% loop over modulated ROIs and classify as activated or suppressed based on
% comparing mean activity in a baseline window (-3s to -2s) vs a test
% period (-1s to 0s before SWR peak)
for mod_neurons = 1:length( swr_modulated_ROI_nr )
    neuron_n = swr_modulated_ROI_nr(mod_neurons);
    if mean( mean_swr_activity(neuron_n,1:31)) > mean( mean_swr_activity(neuron_n,62:93))
        suppressed_ROIs(m,:) = mean_swr_activity( neuron_n,:);
        m = m+1;
    elseif mean( mean_swr_activity(neuron_n,1:31)) < mean( mean_swr_activity(neuron_n, 62:93))
        activated_ROIs(n,:) = mean_swr_activity( neuron_n,:);
        n = n+1;
    else
        unclassified_ROIs(p,:) = mean_swr_activity( neuron_n,:);
        p = p+1;
    end
end
    
%% Find the indicies of SWR modulated ROIs
idx_activated_ROIs = zeros(size(activated_ROIs,1),1);
% loop over nr of activated ROIs
for j = 1:size(activated_ROIs,1)
    % find the indicies of the ROIs that match the activated ROI signal.
    temp_var      = sum( mean_swr_activity == activated_ROIs(j, :),2);
    ROI_idx       = temp_var == win_length;
    ROI_list      = (1:roinr)'; 
    idx_activated_ROIs(j) = ROI_list(ROI_idx);
end

%list the identity of SWR-inhibited ROIs
idx_suppressed_ROIs = zeros(size(suppressed_ROIs,1),1);
for j = 1:size(suppressed_ROIs,1)
    temp_var      = sum( mean_swr_activity == suppressed_ROIs(j, :),2);
    ROI_idx       = temp_var == win_length;
    ROI_list      = (1:roinr)'; 
    idx_suppressed_ROIs(j) = ROI_list(ROI_idx);
end
%     activity_mod_ROIs{i} = {suppressed_ROIs, activated_ROIs, unclassified_ROIs};
%     idx_mod_ROIs{i} = {idx_activated_ROIs,idx_suppressed_ROIs};
%     clearvars suppressed_ROIs activated_ROIs unclassified_ROIs n m p...
%         idx_activated_ROIs idx_suppressed_ROIs
% end
%% Find corresponding dF/F / deconvolved dF/F signals for plotting

% first find the swr segments for activated ROIs
activated_sig_to_plot = zeros(length(idx_activated_ROIs), win_length);
if ~isempty(idx_activated_ROIs) 
    % same code to extract 2P swr segments as above
    for i = 1:length(idx_activated_ROIs)
        mod_ROI = signal_plot(idx_activated_ROIs(i),:);
        for swr_nr = 1:length(swr_idx)
            swr_window_start = swr_idx(swr_nr) - (nr_of_seconds*31); 
            swr_window_end = swr_idx(swr_nr) + (nr_of_seconds*31); 
            if length(swr_window_start:swr_window_end) == win_length
                swr_activity(swr_nr, :) = ...
                    mod_ROI(swr_window_start:swr_window_end); 
            end
        end
        activated_sig_to_plot(i,:) = nanmean(swr_activity);
    end
end

% ... and now for the spindle suppressed ROis
suppressed_sig_to_plot = zeros(length(idx_suppressed_ROIs), win_length);
if ~isempty(idx_suppressed_ROIs) 
    for i = 1:length(idx_suppressed_ROIs)
        mod_ROI = signal_plot(idx_suppressed_ROIs(i),:);
        for swr_nr = 1:length(swr_idx)
            swr_window_start = swr_idx(swr_nr) - (nr_of_seconds*31); 
            swr_window_end = swr_idx(swr_nr) + (nr_of_seconds*31); 
            if length(swr_window_start:swr_window_end) == win_length
                swr_activity(swr_nr, :) = ...
                    mod_ROI(swr_window_start:swr_window_end); 
            end
        end
        suppressed_sig_to_plot(i,:) = nanmean(swr_activity);
    end
end

%% Plotting

figure, 
subplot(241)
plot(time,activated_ROIs' + (size(activated_ROIs,1):-1:1)*.005)
set(gca, 'xlim', [time(1), time(end)])

subplot(242)
plot(time,activated_sig_to_plot' + (size(activated_sig_to_plot,1):-1:1)*.005)
set(gca, 'xlim', [time(1), time(end)])
plot(time,suppressed_ROIs' + (size(suppressed_ROIs,1):-1:1)*.005)
set(gca, 'xlim', [time(1), time(end)])

subplot(243)
plot(time,suppressed_ROIs' + (size(suppressed_ROIs,1):-1:1)*.005)
set(gca, 'xlim', [time(1), time(end)])

subplot(244)
plot(time,suppressed_sig_to_plot' + (size(suppressed_sig_to_plot,1):-1:1)*.005)
set(gca, 'xlim', [time(1), time(end)])

subplot(245)
shadedErrorBar(time, nanmean(activated_ROIs), ...
    std(activated_ROIs)./sqrt(size(activated_ROIs,1)) ,'lineprops', 'b');
xlabel('Time from SWR peak (sec)')
ylabel(label1)
set(gca, 'xlim',[min(time) max(time)])

subplot(246)
shadedErrorBar(time, nanmean(activated_sig_to_plot), ...
    std(activated_sig_to_plot)./sqrt(size(activated_sig_to_plot,1)) ,'lineprops', 'b');
xlabel('Time from SWR peak (sec)')
ylabel(label2)
set(gca, 'xlim',[min(time) max(time)])

subplot(247)
shadedErrorBar(time, nanmean(suppressed_ROIs), ...
    std(suppressed_ROIs)./sqrt(size(suppressed_ROIs,1)) ,'lineprops', 'b');
xlabel('Time from SWR peak (sec)')
ylabel(label1)
set(gca, 'xlim',[min(time) max(time)])

subplot(248)
shadedErrorBar(time, nanmean(suppressed_sig_to_plot), ...
    std(suppressed_sig_to_plot)./sqrt(size(suppressed_sig_to_plot,1)) ,'lineprops', 'b');
xlabel('Time from SWR peak (sec)')
ylabel(label2)
set(gca, 'xlim',[min(time) max(time)])





