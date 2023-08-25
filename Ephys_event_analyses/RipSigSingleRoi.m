function plot_single_roi_activity(sData, roinr)

%Written by Christoffer Berge | Vervaeke Lab

prompt = sprintf('How many seconds before/after ripple peak? ');
nrOfSeconds = input(prompt);

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
        RippleIdx = RipRun(sData); %removes all ripples if there is movement 2s before or after ripple peak
    elseif strcmp(removerip, 'y') && ~strcmp(riprun, 'y')
        RippleIdx = RemoveRip(sData);
    elseif strcmp(removerip, 'y') && strcmp(riprun, 'y')
        sData.ephysdata.frameRipIdx = RipRun(sData);
        RippleIdx = RemoveRip(sData);
    end
end

threshold = true(1,1);

% compute ripple-aligned activity for neuropil, DF/F, and deconvolved data
% in sequential order
for jjj = 1:3
    if jjj == 1
        v = sData.imdata.roiSignals(2).npilMediF;
    elseif jjj == 2
        v = sData.imdata.roiSignals(2).newdff;
    elseif jjj == 3
        v = sData.imdata.roiSignals(2).ciaDeconvolved;
    end

rippleactivity = zeros(length(RippleIdx),((nrOfSeconds*31)*2+1)); % 31 = frames per second, +1 frame for ripple peak
signalActivity = zeros(length(RippleIdx),((nrOfSeconds*31)*2+1));
zscoredRippleActivity = zeros(length(RippleIdx),((nrOfSeconds*31)*2+1));

    V = v(roinr, :); %creates vector of the frame signals for a particular roi
    for ripplenr = 1:length(RippleIdx) % creates vector of nr of ripples during recording
        x = RippleIdx(ripplenr) - (nrOfSeconds*31); %gives frame nr X  before ripple peak 
        if x < 0
            x = 1;
        end
        y = RippleIdx(ripplenr) + (nrOfSeconds*31); %gives frame nr X after ripple peak
        if y > length(v)
            y = length(v);
        end
        if length(x:y) == (nrOfSeconds*31)*2+1 %if the number of frames from 2s prior to 2s post ripple occurence equals 125 frames, do this
            rippleactivity(ripplenr, :) = V(x:y); %index every individual ripple frame segment into a new row for that particular roi(i) )
       %skip ripples that occurred so early/late in recording that there isn't 2 full seconds before/after ripple peak 
        elseif length(x:y) < (nrOfSeconds*31)*2+1 && length(x:y) > (nrOfSeconds*31)*2+1 
            % do nothing
        end
        
    signalActivity(ripplenr, :) = rippleactivity(ripplenr, :); %averages all ripple activity for a roi and assigns it into a new row
    zscoredRippleActivity(ripplenr, :) = zscore(signalActivity(ripplenr, :)); %create zscores for each frame value of the roi 
    end
    
        % with an SNR of 0.99, the deconvolution can detect small signal but
    % also includes false positives, i.e. noise. To remove this, the
    % following thresholding procedure is applied. First, very small
    % deconvolved values are removed ( < 0.0000001). Then the remaining
    % deconvolved values are split into two clusters using k-means. Usually
    % the first cluster captures the noise while the second cluster
    % captures deconvolved signal reflecting actual Ca2+ transients. Next,
    % a linearly spaced vector from the first centroid in cluster 1 to the
    % second centroid in cluster 2 is created. Finally, the threshold is
    % set to the 15th percentile of this vector. (Additionally, extremely
    % large deconvolved values (above 1) are set to 1).
    
    if jjj == 3 && threshold == 1
    unfilteredRippleActivity = rippleactivity;
    tempSignal = v(roinr,:);
    %check if there are (1) no deconvolved events (clustering is not
    %possible) or (2) all values are NaNs
    if sum(tempSignal) == 0 || sum(isnan(tempSignal)) == size(v,2) 
        % do nothing
    else
        index1 = tempSignal > 0.0000001;
        index2 = tempSignal(index1);
        [~, C] = kmeans(index2', 2, 'Replicates', 100);
        Vector = linspace(C(1), C(2));
        ThresholdPercentile = prctile(Vector, 15);
    end
            for rippleNR = 1:size(rippleactivity,1)
                for ii = 1:size(rippleactivity,2)
                    if rippleactivity(rippleNR, ii) < ThresholdPercentile 
                        rippleactivity(rippleNR, ii) = 0; 
                    elseif rippleactivity(rippleNR, ii) > 1
                        rippleactivity(rippleNR, ii) = 1; 
                    else
                        % do nothing
                    end
                end
            end   
    end
    
if jjj == 1
figure,
subplot(3,3,1)
imagesc(signalActivity)
xticks(1:31:((nrOfSeconds*31)*2+1))
set(gca, 'XTickLabel', (-(nrOfSeconds):1:nrOfSeconds))
xlabel('Time from ripple peak (sec)')
ylabel('SWR #')
colorbar
title('Neuropil')

subplot(3,3,7)
imagesc(zscoredRippleActivity)
xticks(1:31:((nrOfSeconds*31)*2+1))
set(gca, 'XTickLabel', (-(nrOfSeconds):1:nrOfSeconds))
xlabel('Time from ripple peak (sec)')
ylabel('SWR #')
colorbar
title('z-score Neuropil')

subplot(3,3,4)
shadedErrorBar([], nanmean(signalActivity), std(signalActivity)/sqrt(numel(signalActivity(:, 1))) ,'lineprops', 'b');
xticks(1:31:((nrOfSeconds*31)*2+1))
set(gca, 'XTickLabel', (-(nrOfSeconds):1:nrOfSeconds))
xlabel('Time from ripple peak (sec)')
ylabel('Mean fluorescence')
title('Mean neuropil fluorescence')

elseif jjj == 2
subplot(3,3,2)
imagesc(signalActivity)
xticks(1:31:((nrOfSeconds*31)*2+1))
set(gca, 'XTickLabel', (-(nrOfSeconds):1:nrOfSeconds))
xlabel('Time from ripple peak (sec)')
ylabel('SWR #')
colorbar
title('DF/F subtracted npil')

subplot(3,3,8)
imagesc(zscoredRippleActivity)
xticks(1:31:((nrOfSeconds*31)*2+1))
set(gca, 'XTickLabel', (-(nrOfSeconds):1:nrOfSeconds))
xlabel('Time from ripple peak (sec)')
ylabel('SWR #')
colorbar
title('z-score DF/F subtracted npil')

subplot(3,3,5)
shadedErrorBar([], nanmean(signalActivity), std(signalActivity)/sqrt(numel(signalActivity(:, 1))) ,'lineprops', 'b');
xticks(1:31:((nrOfSeconds*31)*2+1))
set(gca, 'XTickLabel', (-(nrOfSeconds):1:nrOfSeconds))
xlabel('Time from ripple peak (sec)')
ylabel('Mean DF/F')
title('Mean DF/F subtracted npil')

elseif jjj == 3
subplot(3,3,3)
imagesc(rippleactivity)
xticks(1:31:((nrOfSeconds*31)*2+1))
set(gca, 'XTickLabel', (-(nrOfSeconds):1:nrOfSeconds))
xlabel('Time from ripple peak (sec)')
ylabel('SWR #')
colorbar
title('deconvolved Ca2+ signal')

subplot(3,3,9)
imagesc(signalActivity)
xticks(1:31:((nrOfSeconds*31)*2+1))
set(gca, 'XTickLabel', (-(nrOfSeconds):1:nrOfSeconds))
xlabel('Time from ripple peak (sec)')
ylabel('SWR #')
colorbar
title('deconvolved Ca2+ signal (without threshold)')

subplot(3,3,6)
shadedErrorBar([], nanmean(rippleactivity), std(rippleactivity)/sqrt(numel(rippleactivity(:, 1))), 'lineprops', 'b');
xticks(1:31:((nrOfSeconds*31)*2+1))
set(gca, 'XTickLabel', (-(nrOfSeconds):1:nrOfSeconds))
xlabel('Time from ripple peak (sec)')
ylabel('Mean deconvolved Ca2+')
title('Mean deconvolved Ca2+ signal')

x = roinr;
txt = ['ROI #', num2str(x)];
sgtitle(txt)
end
end


%figure
%subplot(2,1,1)
%imagesc(deconvSWRmodulation)
%xticks(1:31:125)
%xticklabels({'-2', '-1', '0', '1', '2'})
%xlabel('Time from ripple peak (sec)')
%ylabel('SWR #')
%colorbar
%title('deconvolved Ca2+ signal')

%subplot(2,1,2)
%errorbar(mean(deconvSWRmodulation), std(deconvSWRmodulation)/sqrt(numel(deconvSWRmodulation(:, 1))));
%xticks(1:31:125)
%xticklabels({'-2', '-1', '0', '1', '2'})
%xlabel('Time from ripple peak (sec)')
%ylabel('Mean deconvolved Ca2+')
%title('Mean deconvolved Ca2+ signal')