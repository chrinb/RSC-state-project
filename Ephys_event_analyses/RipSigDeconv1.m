function [zscoredRippleActivity, deconvAverage, DeconvWeightedAvg] = RipSigDeconv1(sData)

%Written by Christoffer Berge | Vervaeke Lab

prompt = sprintf('All ripples? (y = yes | everything else = no) ');
allrip = input(prompt,'s');
if strcmp(allrip,'y')
    RippleIdx = sData.ephysdata.frameRipIdx; %keep all ripples
else
    prompt = sprintf('Only no movement ripples? (y = yes | everything else = no) ');
    riprun = input(prompt, 's');
    prompt = sprintf('Remove within 4s of another? (y = yes | everything else = no) ');
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

threshold = false(1,1);
prompt = sprintf('Threshold deconvolved signal? (y = yes | everything else = no) ');
dothreshold = input(prompt, 's');
if strcmp(dothreshold, 'y')
    threshold = true(1,1);
end

rippleactivity = zeros(length(RippleIdx),125); % 125 reflects total nr of frames surrounding each ripple (2s), will be different if looking for longer time periods
deconvAverage = zeros(size(sData.imdata.roiSignals(2).deconv,1),125);
zscoredRippleActivity = zeros(size(sData.imdata.roiSignals(2).deconv,1),125);
for i = 1:size(sData.imdata.roiSignals(2).deconv,1) %runs through the nr of coloumns (corresponding to nr of rois)
    roinr = i;
    v = sData.imdata.roiSignals(2).deconv(roinr, :); %creates vector of the frame signals for a particualr roi
    for j = 1:length(RippleIdx) % creates vector of nr of ripples during recording
        ripplenr = j;
        x = RippleIdx(ripplenr) - 62; %gives frame nr 62 frames before ripple peak (2s)
        if x < 0
            x = 1;
        end
        y = RippleIdx(ripplenr) + 62; %gives frame nr 62 frames after ripple peak (2s)
        if y > length(sData.imdata.roiSignals(2).deconv)
            y = length(sData.imdata.roiSignals(2).deconv);
        end
        if length(x:y) == 125 %if the number of frames from 2s prior to 2s post ripple occurence equals 125 frames, do this
            rippleactivity(j, :) = v(x:y); %index every individual ripple frame segment into a new row for that particular roi(i) )
        elseif length(x:y) < 125 && length(x:y) > 125 %skip ripples that occurred so early/late in recording that there isn't 2 full seconds before/after ripple peak 
            % do nothing
        end
        
        % if nr of input arguments is 4, do a simple thresholding of the
        % deconvolved signal to remove
        % events that are most likely noise, i.e. errors in the deconvolution
        % due to noise in the DF/F trace. Using the
        % pipe_extractSignalsFromROIsS function, a threshold of 0.05 keeps
        % most deconvolved signals reflecting actual Ca2+ transients while
        % removing the smaller deconvolved signals likely reflecting noise.

        if threshold == 1
            for rippleNR = 1:size(rippleactivity,1)
                for ii = 1:size(rippleactivity,2)
                    if rippleactivity(rippleNR, ii) > 0.05
                        % do nothing
                    else
                        rippleactivity(rippleNR, ii) = 0; 
                    end
                end
            end
        else
            % do nothing
        end
        
    end
    
    deconvAverage(roinr, :) = nanmean(rippleactivity); %averages all ripple activity for a roi and assigns it into a new row
    zscoredRippleActivity(roinr, :) = zscore(deconvAverage(roinr, :)); %create zscores for each frame value of the roi 
    DeconvWeightedAvg(roinr) = {[rippleactivity]};
end

ripplenr

[~, pt] = sort(sum(deconvAverage, 2));
sorteddeconvAverage = deconvAverage(pt, :);
[~,t] = sort(sum(zscoredRippleActivity,2));
sortedZscores = zscoredRippleActivity(t, :);

figure,
subplot(2,2,1)
imagesc(deconvAverage) %  colorplot of average individual roi activity during ripples
xticks(1:31:125)
xticklabels({'-2', '-1', '0', '1', '2'})
ylabel('ROI #')
colorbar
xlabel('Time from ripple peak (sec)')
title('Individual ROI mean deconvolved Ca2+')

subplot(2, 2, 2);
imagesc(zscoredRippleActivity) % colorplot of zscored average individual roi activity during ripples
colormap('jet')
xticks(1:31:125)
xticklabels({'-2', '-1', '0', '1', '2'})
xlabel('Time from ripple peak (sec)')
ylabel('ROI #')
title('Individual ROI z-scored deconvolved Ca2+')
colorbar

MeandeconvAverage = nanmean(deconvAverage);
ZscoreGrandMeanAllROIs = zscore(MeandeconvAverage);
AverageZscore = nanmean(sortedZscores);

subplot(2,2,3)
shadedErrorBar([],MeandeconvAverage, std(rmmissing(deconvAverage))./sqrt(numel(rmmissing(deconvAverage(:, 1)))),'lineprops', 'b');
xticks(1:31:125)
xticklabels({'-2', '-1', '0', '1', '2'})
xlabel('Time from ripple peak (sec)')
ylabel('deconvolved')
title('Grand mean deconvolved Ca2+')

subplot(2, 2, 4);
shadedErrorBar([],AverageZscore, std(rmmissing(sortedZscores))./sqrt(numel(rmmissing(sortedZscores(:, 1)))),'lineprops', 'b');
xticks(1:31:125)
xticklabels({'-2', '-1', '0', '1', '2'})
xlabel('Time from ripple peak (sec)')
ylabel('zScore deconvolved')
title('Mean z-scored deconvolved Ca2+')
