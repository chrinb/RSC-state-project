function [zscoredRippleActivity, dffSubtractedNpilAvg, DffWeightedAvg] = RipSigDffs(sData, NoRunRipples, ~)

%Written by Christoffer Berge | Vervaeke Lab

% nr and type of inputs determines which ripples are to be included in the
% analysis. "sData" is always input argument nr 1. 
if nargin == 2 && NoRunRipples == 1 
    RippleIdx = RipRun(sData); %removes all ripples if there is movement 2s before or after ripple peak
elseif nargin == 2 && NoRunRipples == 0
    RippleIdx = RemoveRip(sData); %removes ripples occuring within 4s of another
elseif nargin == 3 
    sData.ephysdata.frameRipIdx = RipRun(sData); %first removes movement-ripples, then ripples occuring within 4s of another
    RippleIdx = RemoveRip(sData);
else
    RippleIdx = sData.ephysdata.frameRipIdx; % do analysis with all ripples
end

rippleactivity = zeros(length(RippleIdx),125); % 125 reflects total nr of frames surrounding each ripple (2s), will be different if looking for longer time periods
dffSubtractedNpilAvg = zeros(size(sData.imdata.roiSignals(2).newdff,1),125);
zscoredRippleActivity = zeros(size(sData.imdata.roiSignals(2).newdff,1),125);

for i = 1:size(sData.imdata.roiSignals(2).newdff,1) %runs through the nr of coloumns (corresponding to nr of rois)
    roinr = i;
    v = sData.imdata.roiSignals(2).newdff(roinr, :); %creates vector of the frame signals for a particualr roi
    for j = 1:length(RippleIdx) % creates vector of nr of ripples during recording
        ripplenr = j;
        x = RippleIdx(ripplenr) - 62; %gives frame nr 62 frames before ripple peak (2s)
        if x < 0
            x = 1;
        end
        y = RippleIdx(ripplenr) + 62; %gives frame nr 62 frames after ripple peak (2s)
        if y > length(sData.imdata.roiSignals(2).newdff)
            y = length(sData.imdata.roiSignals(2).newdff);
        end
        if length(x:y) == 125 %if the number of frames from 2s prior to 2s post ripple occurence equals 125 frames, do this
            rippleactivity(j, :) = v(x:y); %index every individual ripple frame segment into a new row for that particular roi(i) )
        elseif length(x:y) < 125 
            % do nothing
        end
    end
    
    dffSubtractedNpilAvg(roinr, :) = nanmean(rippleactivity); %averages all ripple activity for a roi and assigns it into a new row
    zscoredRippleActivity(roinr, :) = zscore(dffSubtractedNpilAvg(roinr, :)); %create zscores for each frame value of the roi 
    DffWeightedAvg(roinr) = {[rippleactivity]};

end

ripplenr

[~, index1] = sort(sum(dffSubtractedNpilAvg, 2));
sortedDffsAverage = dffSubtractedNpilAvg(index1, :);
[~,index2] = sort(sum(zscoredRippleActivity,2));
sortedZscores = zscoredRippleActivity(index2, :);

figure,
subplot(2,2,1)
imagesc(dffSubtractedNpilAvg) % sorted colorplot of average individual roi activity during ripples
xticks(1:31:125)
xticklabels({'-2', '-1', '0', '1', '2',})
xlabel('Time from ripple peak (sec)')
ylabel('ROI #')
colorbar
title('mean DFF')

subplot(2,2,2)
imagesc(zscoredRippleActivity) 
xticks(1:31:125)
xticklabels({'-2', '-1', '0', '1', '2',})
caxis([-3 3])
xlabel('Time from ripple peak (sec)')
ylabel('ROI #')
colorbar
title('z-score mean DFF')

MeanDffsAverage = nanmean(sortedDffsAverage);
ZscoreGrandMeanAllROIs = zscore(MeanDffsAverage);
AverageZscore = nanmean(sortedZscores);

subplot(2,2,3)
SEM = std(rmmissing(sortedDffsAverage))./sqrt(numel(rmmissing(sortedDffsAverage(:, 1))));
shadedErrorBar([],MeanDffsAverage, SEM,'lineprops', 'b');
xticks(1:31:125)
xticklabels({'-2', '-1', '0', '1', '2'})
xlabel('Time from ripple peak (sec)')
ylabel('DFF')
title('Grand mean DFF')

subplot(2,2,4);
SEM = std(rmmissing(sortedZscores))./sqrt(numel(rmmissing(sortedZscores(:,1))));
shadedErrorBar([],AverageZscore, SEM,'lineprops', 'b');
xticks(1:31:125)
xticklabels({'-2', '-1', '0', '1', '2'})
xlabel('Time from ripple peak (sec)')
ylabel('zScore DFF')
title('Grand z-score mean DFF')