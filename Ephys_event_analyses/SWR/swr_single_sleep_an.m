function [varargout] = swr_single_sleep_an(varargin)

%Written by Christoffer Berge | Vervaeke Lab

% Calculate mean/median activity (e.g., dF/F or deconvolved dF/F) around SWRs for 
% a population of cells. Most up-to-date code. 

sData = varargin{1,1};


%% Select signal 
prompt = sprintf('Which signal? (1 = deconvolved dF/F | 2 = dF/F | 3 = other) ');
signalSelection = input(prompt); 

if signalSelection == 1
    signal = sData.imdata.roiSignals(2).ciaDeconvolved;
    text = 'Mean deconvolved dF/F';
    text2 = 'Mean z-score deconvolved dF/F';
elseif signalSelection == 2
    signal = sData.imdata.roiSignals(2).newdff;
    text = 'mean dF/F';
    text2 = 'mean zscored dF/F';
elseif signalSelection == 3
    prompt = sprintf('Type structname: '); 
    signal = input(prompt);
    text = 'Mean neuropil?';
    text2 = 'Mean z-score neuropil?';
end

% Seconds before/after ripple peak
nr_of_seconds = 3;
sessionID = sData.sessionInfo.sessionID;
nr_of_frames = (nr_of_seconds*31*2)+1;
time = (-(31*nr_of_seconds):(31*nr_of_seconds))./31;
%% Select spindle freq
prompt = sprintf('8-18Hz (1) or 10-16hz (2)? ');
spindle_band_select = input(prompt);

if spindle_band_select == 1
    spindle_select = [];
elseif spindle_band_select == 2
    spindle_select = '1016';
end

unclassified_swr_str           = strcat('unclassified_swr', spindle_select);
NREM_spindle_uncoupled_swr_str = strcat('NREM_spindle_uncoupled_swr', spindle_select);
spindle_coupled_swr_str        = strcat('spindle_coupled_swr', spindle_select);
spin_idx_str                   = strcat('NREMAbsSpindleIdx', spindle_select);
nr_of_spindles_in_analysis     = length(sData.ephysdata2.(spin_idx_str) );

%% Select SWRs for analysis
prompt = sprintf('All SWRs (1) | Rem. locom. SWRs (2) | Rem. clust. SWRs (3) | Rem- locom. & clust. SWRs (4) ');
swr_for_analysis = input(prompt);

frames  = sData.daqdata.frame_onset_reference_frame;
swr_idx = sData.ephysdata.absRipIdx;

if swr_for_analysis == 1
    RippleIdx            = sData.ephysdata.absRipIdx;
    unclassified_SWR_idx = ismember(RippleIdx, sData.ephysdata.(unclassified_swr_str) );
    single_swr_idx       = ismember(RippleIdx, sData.ephysdata.(NREM_spindle_uncoupled_swr_str));
    spindle_coupled_idx  = ismember(RippleIdx, sData.ephysdata.(spindle_coupled_swr_str));
elseif swr_for_analysis == 2 % Remove movement-SWRs
    [RippleIdx,~] = riprun2(sData, swr_idx);
elseif swr_for_analysis == 3 % Remove close SWRs 
    [RippleIdx] = RemoveRip(swr_idx);
elseif swr_for_analysis == 4 % Remove movemement and close SWRs
    [RippleIdx,~] = riprun2(sData, swr_idx);
    [RippleIdx]   = RemoveRip(RippleIdx);
end

if swr_for_analysis ~= 1 
    unclassified_SWR_idx = ismember(RippleIdx, sData.ephysdata.(unclassified_swr_str) );
    single_swr_idx       = ismember(RippleIdx, sData.ephysdata.(NREM_spindle_uncoupled_swr_str));
    spindle_coupled_idx  = ismember(RippleIdx, sData.ephysdata.(spindle_coupled_swr_str));
      
    % Preallocate
    signal_swr_activity_unclassified = zeros(size(signal,1), nr_of_frames);
    signal_swr_activity_single       = zeros(size(signal,1), nr_of_frames);
    signal_swr_activity_coupled      = zeros(size(signal,1), nr_of_frames);

end

% Specify whether to threshold deconvolved dF/F or not
if signalSelection == 1
    prompt = sprintf('Threshold deconvolved dF/F? (1) ');
    dothreshold = input(prompt);
    if dothreshold == 1
        threshold = true(1,1);
    else
        threshold = false(1,1);
    end
end


nr_of_frames = (nr_of_seconds*31*2)+1;
time = (-(31*nr_of_seconds):(31*nr_of_seconds))./31;


% convert to imaging frames
RippleIdx = frames(round(RippleIdx));
% preallocate 
signal_swr_activity = zeros(length(RippleIdx), nr_of_frames);
%% Run analysis
%loop over ROIs
for roinr = 1:size(signal,1) %runs through the nr of coloumns (corresponding to nr of rois)
    roi_signal = signal(roinr, :); %creates vector of the frame signals for a particualr roi

    for swr_nr = 1:length(RippleIdx)
        % find SWR window start end
        swr_window_start = RippleIdx(swr_nr) - (nr_of_seconds*31); 
        swr_window_end   = RippleIdx(swr_nr) + (nr_of_seconds*31); 

        % skip SWRs at the beginning or end with a window shorter than
        % length specified in seconds by user above
        if swr_window_start > 1 && swr_window_end < size(signal,2)
            signal_swr_activity(swr_nr, :) = ...
                roi_signal(swr_window_start:swr_window_end); 
        end
    end
    % Threshold signal using K-means clustering to separate "noise" vs
    % true deconvolved events.
    if signalSelection == 1 && threshold == 1 
        [threshold_percentile, ~] = threshold_deconvolved(signal, roi_signal, roinr);

        % Remove deconvolved values below threshold
        ROI_idx = reshape(signal_swr_activity, 1,[]);
        ROI_idx(ROI_idx < threshold_percentile(roinr)) = 0;
        signal_swr_activity = reshape(ROI_idx, size(signal_swr_activity));
    end
        

    signal_swr_unclassified    = signal_swr_activity(unclassified_SWR_idx,:);
    signal_swr_spindle_single  = signal_swr_activity(single_swr_idx,:);
    signal_swr_spindle_coupled = signal_swr_activity(spindle_coupled_idx,:);

    signal_swr_activity_unclassified(roinr,:) = nanmean(signal_swr_unclassified);
    signal_swr_activity_single(roinr,:)       = nanmean(signal_swr_spindle_single);
    signal_swr_activity_coupled(roinr,:)      = nanmean(signal_swr_spindle_coupled);
end

% Z-score data
signal_swr_activity_unclassifiedZ   = zscore( signal_swr_activity_unclassified, 0,2);
signal_swr_activity_singleZ         = zscore( signal_swr_activity_single, 0,2);
signal_swr_activity_coupledZ        = zscore( signal_swr_activity_coupled, 0,2);

% Store output in cell array
mean_signal = {signal_swr_activity_unclassified, signal_swr_activity_unclassifiedZ,...
    signal_swr_activity_single, signal_swr_activity_singleZ, ...
    signal_swr_activity_coupled, signal_swr_activity_coupledZ};

varargout{1} = mean_signal;
%% Calculate SE of average activity
SE_spindle_A  = std(signal_swr_unclassified,'omitnan') ./ sqrt(size(signal_swr_unclassified,1));
SE_spindle_AZ = std(signal_swr_activity_unclassifiedZ,'omitnan') ./ sqrt(size(signal_swr_activity_unclassifiedZ,1));

SE_spindle_C  = std(signal_swr_spindle_single,'omitnan') ./ sqrt(size(signal_swr_spindle_single,1));
SE_spindle_CZ = std(signal_swr_activity_coupledZ,'omitnan') ./ sqrt(size(signal_swr_activity_coupledZ,1));

SE_spindle_U  = std(signal_swr_spindle_coupled,'omitnan') ./ sqrt(size(signal_swr_spindle_coupled,1));
SE_spindle_UZ = std(signal_swr_activity_singleZ,'omitnan') ./ sqrt(size(signal_swr_activity_singleZ,1));

%% Plot results
if length(varargin) > 1
    figure,
    sgtitle(sessionID),

    x1 = [time(1) time(end)];
    y_unclassified = [1 size(signal_swr_activity_unclassified,1) ];
    y_spindleC = [1 size(signal_swr_activity_unclassified,1) ];
    y_spindleU = [1 size(signal_swr_activity_unclassified,1) ];

    %% unclassified SWRs
    if length(nanmean(signal_swr_unclassified)) > 1
    subplot(231)
    imagesc(x1,y_unclassified, signal_swr_activity_unclassifiedZ) %  colorplot of average individual roi activity during ripples
    ylabel('ROI #')
    xlabel('Time from ripple peak (sec)')
    title(['unclassified SWR (n = ',  num2str(sum(unclassified_SWR_idx)), ')'])

    subplot(234)
    shadedErrorBar(time, nanmean(signal_swr_activity_unclassifiedZ),SE_spindle_AZ,'lineprops', 'r');
    xlabel('Time from ripple peak (sec)')
    ylabel(text2)
    %title('unclassified SWR')
    end

    %% Spindle-uncoupled SWRs
    subplot(232)
    imagesc(x1, y_spindleU, signal_swr_activity_singleZ) %  colorplot of average individual roi activity during ripples
    ylabel('ROI #')
    xlabel('Time from ripple peak (sec)')
    title(['Spindle-uncoupled SWR (n = ', num2str(sum(single_swr_idx)), ')'])


    subplot(235)
    shadedErrorBar(time, nanmean(signal_swr_activity_singleZ),SE_spindle_UZ,'lineprops', 'b');
    xlabel('Time from ripple peak (sec)')
    ylabel(text2)
    %title('Spindle-uncoupled SWR')

    %% Spindle-coupled SWRs
    subplot(233)
    imagesc(x1,y_spindleC, signal_swr_activity_coupledZ) %  colorplot of average individual roi activity during ripples
    ylabel('ROI #')
    xlabel('Time from ripple peak (sec)')
    title(['Spindle-coupled SWR (n = ', num2str(sum(spindle_coupled_idx)),...
           ' spindle n = ',  num2str(nr_of_spindles_in_analysis), ')'])

    subplot(236)
    shadedErrorBar(time, nanmean(signal_swr_activity_coupledZ),SE_spindle_CZ,'lineprops', 'k');
    xlabel('Time from ripple peak (sec)')
    ylabel(text2)
end

