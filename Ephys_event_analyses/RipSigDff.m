function [varargout] = RipSigDff(sData)

%Written by Christoffer Berge | Vervaeke Lab

sessionID = sData.sessionInfo.sessionID;

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

prompt = sprintf('How many seconds before/after ripple peak? ');
nrOfSeconds = input(prompt);


%% Select SWRs for analysis

frames = sData.daqdata.frame_onset_reference_frame;

prompt = sprintf('All ripples? (y = yes | everything else = no) ');
allrip = input(prompt,'s');

if strcmp(allrip,'y')
    RippleIdx = sData.ephysdata.frameRipIdx; %keep all ripples
else
    prompt = sprintf('Only no movement ripples? (y = yes | everything else = no) ');
    riprun = input(prompt, 's');
    prompt = sprintf('Remove temporally close ripples? (y = yes | everything else = no) ');
    removerip = input(prompt, 's');
    if strcmp(riprun, 'y') && ~strcmp(removerip, 'y')
        RippleIdx = riprun2(sData); %removes all ripples if there is movement 2s before or after ripple peak
        RippleIdx = frames(RippleIdx);

    elseif strcmp(removerip, 'y') && ~strcmp(riprun, 'y')
        RippleIdx = RemoveRip(sData);
        RippleIdx = frames(RippleIdx);

    elseif strcmp(removerip, 'y') && strcmp(riprun, 'y')
        sData.ephysdata.absRipIdx = riprun2(sData);
        RippleIdx = RemoveRip(sData);
        RippleIdx = frames(RippleIdx);
    end
end

prompt = sprintf('Do baseline subtraction? (y = yes | everything else = no) ');
baseSub = input(prompt,'s');

%% Calculate average imaging activity 
nrOfFrames = (nrOfSeconds*31*2)+1;
time = (-(31*nrOfSeconds):(31*nrOfSeconds))./31;

% preallocate various matrices
awakeSWRactivity = zeros(length(RippleIdx),nrOfFrames); 

% Find average DF/F signal across all ROIs
meanDFF = mean(signal);

% loop over awake SWRs
for awake_SWRnr = 1:length(RippleIdx) % creates vector of nr of ripples during recording
    x = RippleIdx(awake_SWRnr) - (nrOfSeconds*31); %gives frame nr X  before ripple peak 
        
        % check if first time point in SWR-aligned time series begins
        % before first imaging frame, or ends after last imaging frame in
        % in session. If so, set time point to 1 or last imaging frame,
        % respectively. 
        if x <= 0
            x = 1;
        end
        y = RippleIdx(awake_SWRnr) + (nrOfSeconds*31); %gives frame nr X after ripple peak
        if y > length(signal)
            y = length(signal);
        end
        
        if length(x:y) == (nrOfSeconds*31)*2+1 %if the number of frames from 2s prior to 2s post ripple occurence equals 125 frames, do this
            awakeSWRactivity(awake_SWRnr, :) = meanDFF(x:y); %index every individual ripple frame segment for the mean DF/F )
        elseif length(x:y) < (nrOfSeconds*31)*2+1 && length(x:y) > (nrOfSeconds*31)*2+1 %skip ripples that occurred so early/late in recording that there isn't 2 full seconds before/after ripple peak 
            awakeSWRactivity(awake_SWRnr, :) = NaN;
        end
        
    % Baseline subtraction
    if strcmp(baseSub, 'y')
        baselineDFF = nanmean(awakeSWRactivity(awake_SWRnr,1:31));
        awakeSWRactivity(awake_SWRnr,:) = awakeSWRactivity(awake_SWRnr,:)-baselineDFF;
        clear baselineDFF
    end
end
    
    % average all awake SWR-aligned snippets for a given ROI
    zscore_awakeSWRactivity = zscore(awakeSWRactivity,0,2);
    varargout{1} = awakeSWRactivity;
    varargout{2} = zscore_awakeSWRactivity;

%% Sort ROIs according to activity 

% [~, pt] = sort(sum(avg_awakeSWRactivity, 2));
% sortedAvgAwake = avg_awakeSWRactivity(pt, :);
% [~,t] = sort(sum(zscored_avg_awakeSWRActivity,2));
% sortedZscoresAvgAwake = zscored_avg_awakeSWRActivity(t, :);
% 
% [~, pt] = sort(sum(avg_NREMspindleUncoupledSWRactivity, 2));
% sortedAvgUncoupled = avg_NREMspindleUncoupledSWRactivity(pt, :);
% [~,t] = sort(sum(zscored_avg_NREMspindleUncoupledSWRactivity,2));
% sortedZscoresAvgUncoupled = zscored_avg_NREMspindleUncoupledSWRactivity(t, :);
% 
% [~, pt] = sort(sum(avg_NREMspindleCoupledSWRactivity, 2));
% sortedAvgCoupled = avg_NREMspindleCoupledSWRactivity(pt, :);
% [~,t] = sort(sum(zscored_avg_NREMspindleCoupledSWRactivity,2));
% sortedZscoresAvgCoupled = zscored_avg_NREMspindleCoupledSWRactivity(t, :);

%% Calculate SE of average activity
SE_avg_awake        = std(rmmissing(awakeSWRactivity))./sqrt(numel(rmmissing(awakeSWRactivity(:, 1))));
SE_zscore_awake     = std(rmmissing(zscore_awakeSWRactivity))./sqrt(numel(rmmissing(zscore_awakeSWRactivity(:, 1))));

%% Plot results
figure,
sgtitle(sessionID),
subplot(221)
imagesc(awakeSWRactivity) %  colorplot of average individual roi activity during ripples
xticks(1:31:((nrOfSeconds*31)*2+1))
set(gca, 'XTickLabel', (-(nrOfSeconds):1:nrOfSeconds))
ylabel('SWR #')
xlabel('Time from ripple peak (sec)')

subplot(223)
shadedErrorBar(time, nanmean(awakeSWRactivity),SE_avg_awake,'lineprops', 'r');
set(gca, 'xlim',[time(1), time(end)])
xlabel('Time from ripple peak (sec)')
ylabel(text)

%% Plot z-scored results

subplot(222)
imagesc(zscore_awakeSWRactivity) %  colorplot of average individual roi activity during ripples
xticks(1:31:((nrOfSeconds*31)*2+1))
set(gca, 'XTickLabel', (-(nrOfSeconds):1:nrOfSeconds))
ylabel('SWR #')
xlabel('Time from ripple peak (sec)')

subplot(224)
shadedErrorBar(time, nanmean(zscore_awakeSWRactivity),SE_zscore_awake,'lineprops', 'r');
set(gca, 'xlim',[time(1), time(end)])
xlabel('Time from ripple peak (sec)')
ylabel(text2)

