function [varargout] = RipSigAn(sData)

%Written by Christoffer Berge | Vervaeke Lab

% DO NOT USE!! Something is wrong with SWR type indexing here... use
% "rip_sig_an" for mean population activity or "RipSigAn2" for mean bulk
% activity instead (note: 12.11.2021)


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

%% Select spindle freq
prompt = sprintf('8-18Hz (1) or 10-16hz (2)? ');
spindle_band_select = input(prompt);

if spindle_band_select == 1
    spindle_select = [];
elseif spindle_band_select == 2
    spindle_select = '1016';
end

awake_swr_str                  = strcat('awake_swr', spindle_select);
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
%     RippleIdx = sData.ephysdata.absRipIdx;
    awake_SWR_idx = frames(sData.ephysdata.(awake_swr_str));
    NREM_spinU_SWR_idx = frames(sData.ephysdata.(NREM_spindle_uncoupled_swr_str));
    NREM_spinC_SWR_idx = frames(sData.ephysdata.(spindle_coupled_swr_str));
elseif swr_for_analysis == 2
    [RippleIdx,~] = riprun2(sData, swr_idx);
elseif swr_for_analysis == 3
    [RippleIdx] = RemoveRip(swr_idx);
elseif swr_for_analysis == 4
    [RippleIdx,~] = riprun2(sData, swr_idx);
    [RippleIdx]   = RemoveRip(RippleIdx);
end

if swr_for_analysis ~= 1 
    temp_var1 = ismember(sData.ephysdata.(awake_swr_str), RippleIdx );
    awake_SWR_idx = frames(sData.ephysdata.absRipIdx(temp_var1));

    temp_var2 = ismember( sData.ephysdata.(NREM_spindle_uncoupled_swr_str),RippleIdx);
    NREM_spinU_SWR_idx = frames(sData.ephysdata.absRipIdx(temp_var2));

    temp_var3 = ismember(sData.ephysdata.(spindle_coupled_swr_str),RippleIdx);
    NREM_spinC_SWR_idx = frames(sData.ephysdata.absRipIdx(temp_var3));
end

% Specify whether to threshold deconvolved dF/F or not
prompt = sprintf('Threshold deconvolved dF/F? (y = yes | everything else = no) ');
dothreshold = input(prompt, 's');
if strcmp(dothreshold, 'y')
    threshold = true(1,1);
else
    threshold = false(1,1);
end

%% Calculate average imaging activity 
nr_of_frames = (nr_of_seconds*31*2)+1;
time = (-(31*nr_of_seconds):(31*nr_of_seconds))./31;

% preallocate various matrices
awakeSWRactivity = zeros(length(awake_SWR_idx),nr_of_frames); 
avg_awakeSWRactivity = zeros(size(signal,1),nr_of_frames);
zscored_avg_awakeSWRActivity = zeros(size(signal,1),nr_of_frames);

NREMspindleUncoupledSWRactivity = zeros(length(NREM_spinU_SWR_idx),nr_of_frames); 
avg_NREMspindleUncoupledSWRactivity = zeros(size(signal,1),nr_of_frames);
zscored_avg_NREMspindleUncoupledSWRactivity = zeros(size(signal,1),nr_of_frames);

NREMspindleCoupledSWRactivity = zeros(length(NREM_spinC_SWR_idx),nr_of_frames); 
avg_NREMspindleCoupledSWRactivity = zeros(size(signal,1),nr_of_frames);
zscored_avg_NREMspindleCoupledSWRactivity = zeros(size(signal,1),nr_of_frames);

%loop over ROIs
for roinr = 1:size(signal,1) %runs through the nr of coloumns (corresponding to nr of rois)
    roi_signal = signal(roinr, :); %creates vector of the frame signals for a particualr roi
    
   
%     % loop over awake SWRs
    for awake_SWRnr = 1:length(awake_SWR_idx) % creates vector of nr of ripples during recording
        x = awake_SWR_idx(awake_SWRnr) - (nr_of_seconds*31); %gives frame nr X  before ripple peak 
        
        % check if first time point in SWR-aligned time series begins
        % before first imaging frame, or ends after last imaging frame in
        % in session. If so, set time point to 1 or last imaging frame,
        % respectively. 
        if x <= 0
            x = 1;
        end
        y = awake_SWR_idx(awake_SWRnr) + (nr_of_seconds*31); %gives frame nr X after ripple peak
        if y > length(signal)
            y = length(signal);
        end
        
        if length(x:y) == (nr_of_seconds*31)*2+1 %if the number of frames from 2s prior to 2s post ripple occurence equals 125 frames, do this
            awakeSWRactivity(awake_SWRnr, :) = roi_signal(x:y); %index every individual ripple frame segment into a new row for that particular roi(i) )
        elseif length(x:y) < (nr_of_seconds*31)*2+1 && length(x:y) > (nr_of_seconds*31)*2+1 %skip ripples that occurred so early/late in recording that there isn't 2 full seconds before/after ripple peak 
            awakeSWRactivity(awake_SWRnr, :) = NaN;
        end
    end
    
        % loop over NREM spindle-uncoupled SWRs
    for NREM_spindle_uncoupled_SWRnr = 1:length(NREM_spinU_SWR_idx) 
        x1 = NREM_spinU_SWR_idx(NREM_spindle_uncoupled_SWRnr) - (nr_of_seconds*31); 
 
        if x1 <= 0
            x1 = 1;
        end
        y1 = NREM_spinU_SWR_idx(NREM_spindle_uncoupled_SWRnr) + (nr_of_seconds*31); 
        if y1 > length(signal)
            y1 = length(signal);
        end
        
        if length(x1:y1) == (nr_of_seconds*31)*2+1 
            NREMspindleUncoupledSWRactivity(NREM_spindle_uncoupled_SWRnr, :) = roi_signal(x1:y1); 
        elseif length(x1:y1) < (nr_of_seconds*31)*2+1 && length(x1:y1) > (nr_of_seconds*31)*2+1 
            NREMspindleUncoupledSWRactivity(NREM_spindle_uncoupled_SWRnr, :) = NaN;
        end
    end
    
        % loop over NREM spindle-coupled SWRs
    for NREM_spindle_coupled_SWRnr = 1:length(NREM_spinC_SWR_idx)
        x3 = NREM_spinC_SWR_idx(NREM_spindle_coupled_SWRnr) - (nr_of_seconds*31); 

        if x3 <= 0
            x3 = 1;
        end
        
        y3 = NREM_spinC_SWR_idx(NREM_spindle_coupled_SWRnr) + (nr_of_seconds*31); 
        if y3 > length(signal)
            y3 = length(signal);
        end
        
        if length(x3:y3) == (nr_of_seconds*31)*2+1 
            NREMspindleCoupledSWRactivity(NREM_spindle_coupled_SWRnr, :) = roi_signal(x3:y3); 
        elseif length(x3:y3) < (nr_of_seconds*31)*2+1 && length(x3:y3) > (nr_of_seconds*31)*2+1
            NREMspindleCoupledSWRactivity(NREM_spindle_coupled_SWRnr, :) = NaN;
        end
    end
    
    % average all awake SWR-aligned snippets for a given ROI
    avg_awakeSWRactivity(roinr, :) = nanmean(awakeSWRactivity); 
    zscored_avg_awakeSWRActivity(roinr, :) = zscore(avg_awakeSWRactivity(roinr, :)); 
    awakeWeightedAvg(roinr) = {awakeSWRactivity};
    varargout{1} = avg_awakeSWRactivity;
    varargout{2} = zscored_avg_awakeSWRActivity;
    varargout{3} = awakeWeightedAvg;
    
    % average all NREM spindle-uncoupled SWR-aligned snippets for a given ROI
    avg_NREMspindleUncoupledSWRactivity(roinr, :) = nanmean(NREMspindleUncoupledSWRactivity);    
    zscored_avg_NREMspindleUncoupledSWRactivity(roinr, :) = zscore(avg_NREMspindleUncoupledSWRactivity(roinr, :)); 
    uncoupledWeightedAvg(roinr) = {NREMspindleUncoupledSWRactivity};
    varargout{4} = avg_NREMspindleUncoupledSWRactivity;
    varargout{5} = zscored_avg_NREMspindleUncoupledSWRactivity;
    varargout{6} = uncoupledWeightedAvg;
   
    % average all NREM spindle-coupled SWR-aligned snippets for a given ROI
    avg_NREMspindleCoupledSWRactivity(roinr, :) = nanmean(NREMspindleCoupledSWRactivity); 
    zscored_avg_NREMspindleCoupledSWRactivity(roinr, :) = zscore(avg_NREMspindleCoupledSWRactivity(roinr, :)); 
    coupledWeightedAvg(roinr) = {NREMspindleUncoupledSWRactivity};
    varargout{7} = avg_NREMspindleCoupledSWRactivity;
    varargout{8} = zscored_avg_NREMspindleCoupledSWRactivity;
    varargout{9} = coupledWeightedAvg;
end


% signal_swr_activity_awakeZ     = zscore( signal_swr_activity_awake, 0,2);
% signal_swr_activity_coupledZ   = zscore( signal_swr_activity_coupled, 0,2);
% signal_swr_activity_uncoupledZ = zscore( signal_swr_activity_uncoupled, 0,2);

%% Calculate SE of average activity
SE_avg_awake        = std(rmmissing(avg_awakeSWRactivity))./sqrt(numel(rmmissing(avg_awakeSWRactivity(:, 1))));
SE_zscore_awake     = std(rmmissing(zscored_avg_awakeSWRActivity))./sqrt(numel(rmmissing(zscored_avg_awakeSWRActivity(:, 1))));


SE_uncoupled        = std(rmmissing(avg_NREMspindleUncoupledSWRactivity))./sqrt(numel(rmmissing(avg_NREMspindleUncoupledSWRactivity(:, 1))));
SE_zscore_uncoupled = std(rmmissing(zscored_avg_NREMspindleUncoupledSWRactivity))./sqrt(numel(rmmissing(zscored_avg_NREMspindleUncoupledSWRactivity(:, 1))));

SE_spindle_coupled  = std(rmmissing(avg_NREMspindleCoupledSWRactivity))./sqrt(numel(rmmissing(avg_NREMspindleCoupledSWRactivity(:, 1))));
SE_zscore_coupled   = std(rmmissing(zscored_avg_NREMspindleCoupledSWRactivity))./sqrt(numel(rmmissing(zscored_avg_NREMspindleCoupledSWRactivity(:, 1))));


% %% Plot results
figure,
sgtitle(sessionID),
subplot(231)
imagesc(avg_awakeSWRactivity) %  colorplot of average individual roi activity during ripples
xticks(1:31:((nr_of_seconds*31)*2+1))
set(gca, 'XTickLabel', (-(nr_of_seconds):1:nr_of_seconds))
ylabel('ROI #')
xlabel('Time from ripple peak (sec)')
title(['Awake SWR (n = ', num2str(awake_SWRnr), ')'])

subplot(234)
shadedErrorBar(time, nanmean(avg_awakeSWRactivity),SE_avg_awake,'lineprops', 'r');
xlabel('Time from ripple peak (sec)')
ylabel(text)
%title('Awake SWR')

subplot(232)
imagesc(avg_NREMspindleUncoupledSWRactivity) %  colorplot of average individual roi activity during ripples
xticks(1:31:((nr_of_seconds*31)*2+1))
set(gca, 'XTickLabel', (-(nr_of_seconds):1:nr_of_seconds))
ylabel('ROI #')
xlabel('Time from ripple peak (sec)')
title(['Spindle-uncoupled SWR (n = ', num2str(NREM_spindle_uncoupled_SWRnr), ')'])


subplot(235)
shadedErrorBar(time, nanmean(avg_NREMspindleUncoupledSWRactivity),SE_uncoupled,'lineprops', 'b');
xlabel('Time from ripple peak (sec)')
ylabel(text)
%title('Spindle-uncoupled SWR')

subplot(233)
imagesc(avg_NREMspindleCoupledSWRactivity) %  colorplot of average individual roi activity during ripples
xticks(1:31:((nr_of_seconds*31)*2+1))
set(gca, 'XTickLabel', (-(nr_of_seconds):1:nr_of_seconds))
ylabel('ROI #')
xlabel('Time from ripple peak (sec)')
title(['Spindle-coupled SWR (SWR n = ', num2str(NREM_spindle_coupled_SWRnr),...
       ' spindle n = ', num2str(nr_of_spindles_in_analysis), ')'])

subplot(236)
shadedErrorBar(time, nanmean(avg_NREMspindleCoupledSWRactivity),SE_spindle_coupled,'lineprops', 'k');
xlabel('Time from ripple peak (sec)')
ylabel(text)
%title('Spindle-coupled SWR')

%% Plot z-scored results
figure,
sgtitle(sessionID),
subplot(231)
imagesc(zscored_avg_awakeSWRActivity) %  colorplot of average individual roi activity during ripples
xticks(1:31:((nr_of_seconds*31)*2+1))
set(gca, 'XTickLabel', (-(nr_of_seconds):1:nr_of_seconds))
ylabel('ROI #')
xlabel('Time from ripple peak (sec)')
title(['Awake SWR (n = ', num2str(awake_SWRnr), ')'])

subplot(234)
shadedErrorBar(time, nanmean(zscored_avg_awakeSWRActivity),SE_zscore_awake,'lineprops', 'r');
xlabel('Time from ripple peak (sec)')
ylabel(text2)
%title('Awake SWR')

subplot(232)
imagesc(zscored_avg_NREMspindleUncoupledSWRactivity) %  colorplot of average individual roi activity during ripples
xticks(1:31:((nr_of_seconds*31)*2+1))
set(gca, 'XTickLabel', (-(nr_of_seconds):1:nr_of_seconds))
ylabel('ROI #')
xlabel('Time from ripple peak (sec)')
title(['Spindle-uncoupled SWR (n = ', num2str(NREM_spindle_uncoupled_SWRnr), ')'])


subplot(235)
shadedErrorBar(time, nanmean(zscored_avg_NREMspindleUncoupledSWRactivity),SE_zscore_uncoupled,'lineprops', 'b');
xlabel('Time from ripple peak (sec)')
ylabel(text2)
%title('Spindle-uncoupled SWR')

subplot(233)
imagesc(zscored_avg_NREMspindleCoupledSWRactivity) %  colorplot of average individual roi activity during ripples
xticks(1:31:((nr_of_seconds*31)*2+1))
set(gca, 'XTickLabel', (-(nr_of_seconds):1:nr_of_seconds))
ylabel('ROI #')
xlabel('Time from ripple peak (sec)')
title(['Spindle-coupled SWR (n = ', num2str(NREM_spindle_coupled_SWRnr),...
       ' spindle n = ', num2str(nr_of_spindles_in_analysis), ')'])

subplot(236)
shadedErrorBar(time, nanmean(zscored_avg_NREMspindleCoupledSWRactivity),SE_zscore_coupled,'lineprops', 'k');
xlabel('Time from ripple peak (sec)')
ylabel(text2)
title('Spindle-coupled SWR')

