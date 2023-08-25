function [varargout] = spin_bulk_an(varargin)

sData = varargin{1,1};

%Written by Christoffer Berge | Vervaeke Lab

% Calculate mean activity (e.g., mean dF/F or deconvolved dF/F) during all
% NREM sleep spindles for bulk imaging sessions. 

%% Select signal and temporal window around SWRs
prompt = sprintf('Which signal? (1 = deconvolved | 2 = DF/F | 3 = other) ');
signalSelection = input(prompt); 

if signalSelection == 1
    signal = sData.imdata.roiSignals(2).ciaDeconvolved;
    text = 'Mean deconvolved Ca2+';
elseif signalSelection == 2
    signal = sData.imdata.roiSignals(2).newdff;
    text = 'mean DF/F';
    text2 = 'mean zscored DF/F';
elseif signalSelection == 3
    prompt = sprintf('Type structname: '); 
    signal = input(prompt);
end

nr_of_seconds = 6;

% Choose spindle 
prompt = sprintf('8-18Hz (1) or 10-16hz (2)? ');
spindle_band_select = input(prompt);

if spindle_band_select == 1
    spindle_select = [];
elseif spindle_band_select == 2
    spindle_select = '1016';
end

spin_idx_str = strcat('NREMAbsSpindleIdx', spindle_select);
spin_start_end_str = strcat('NREMspindleStartEnd', spindle_select);

frames = sData.daqdata.frame_onset_reference_frame;
spindleIdx = frames(round(sData.ephysdata2.(spin_idx_str)));

spindleIdxStart = frames(round(sData.ephysdata2.(spin_start_end_str)(:,1)));

sessionID = sData.sessionInfo.sessionID;

prompt = sprintf('Do baseline subtraction? (1) | everything else = no) ');
baseSub = input(prompt);

%% Calculate average imaging activity 
nrOfFrames = (nr_of_seconds*31*2)+1;
time = (-(31*nr_of_seconds):(31*nr_of_seconds))./31;

% Find average DF/F signal across all ROIs
mean_DFF = mean(signal);

% preallocate various matrices
spindleActivity = zeros(length(spindleIdx),nrOfFrames); 
signal_spindle_activity = zeros(length(spindleIdxStart),nrOfFrames); 

%% Calculate mean DF/F from spindle start

% loop over spindle onset times
for spindle_nr = 1:length(spindleIdxStart)
    spindle_window_start = spindleIdxStart(spindle_nr) - (nr_of_seconds*31); 
    spindle_window_end = spindleIdxStart(spindle_nr) + (nr_of_seconds*31); 

    % skip spindles at the beginning or end with a window shorter than
    % length specified in seconds by user above
    if spindle_window_start > 1 && spindle_window_end < length(mean_DFF)
        signal_spindle_activity(spindle_nr, :) = ...
            mean_DFF(spindle_window_start:spindle_window_end); 
    end

    % Baseline subtraction
    if baseSub == 1
    baselineDFF = nanmean(signal_spindle_activity(spindle_nr,1:31));
    signal_spindle_activity(spindle_nr,:) = signal_spindle_activity(spindle_nr,:)-baselineDFF;
    end
end

% average all awake SWR-aligned snippets for a given ROI
signal_spindle_activityZ = zscore(signal_spindle_activity,0,2);
mean_signal              = {signal_spindle_activity, signal_spindle_activityZ};

varargout{1} = signal_spindle_activity;
varargout{2} = signal_spindle_activityZ;
varargout{3} = mean_signal;
%% Sort ROIs according to activity 
% [~, pt] = sort(sum(avg_spindleActivity, 2));
% sortedAvgAwake = avg_spindleActivity(pt, :);

%% Calculate SE of average activity
SE_spindleStart        = std(rmmissing(signal_spindle_activity))./sqrt(numel(rmmissing(signal_spindle_activity(:, 1))));
SE_zscore_spindleStart     = std(rmmissing(signal_spindle_activityZ))./sqrt(numel(rmmissing(signal_spindle_activityZ(:, 1))));

%% Plot DF/F at spindle start
if length(varargin) > 1
    x1 = [time(1) time(end)];
    y1 = [1 size(signal_spindle_activity,1) ];

    figure,
    sgtitle(sessionID) 
    
    subplot(3,2,[1,3]);
    imagesc(x1, y1, signal_spindle_activity) %  colorplot of average individual roi activity during ripples
    ylabel('Spindle #')
    xlabel('Time from spindle onset (s)')
    title(['Spindle n = ', num2str(length(spindleIdxStart))])
    colorbar

    subplot(325)
    shadedErrorBar(time, nanmean(signal_spindle_activity),SE_spindleStart,'lineprops', 'r');
    xlabel('Time from spindle onset (s)')
    ylabel(text)
    set(gca, 'xlim',[min(time) max(time)])
    legend('DF/F')

    subplot(3,2,[2,4]);
    imagesc(x1, y1, signal_spindle_activityZ) %  colorplot of average individual roi activity during ripples
    ylabel('Spindle #')
    xlabel('Time from spindle onset (s)')
    colorbar
    caxis([-3 3])
    
    subplot(326)
    shadedErrorBar(time, nanmean(signal_spindle_activityZ),SE_zscore_spindleStart,'lineprops', 'b');
    xlabel('Time from spindle onset (s)')
    ylabel(text2)
    set(gca, 'xlim',[min(time) max(time)])
    legend('z-score DF/F')
end

