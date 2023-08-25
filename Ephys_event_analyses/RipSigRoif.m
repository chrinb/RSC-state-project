function [RoifAverage, zscoredRippleActivity] = RipSigRoif(sData)

rippleactivity = zeros(length(sData.ephysdata.frameRipIdx),125); % 125 reflects total nr of frames surrounding each ripple (2s), will be different if looking for longer time periods
RoifAverage = zeros(size(sData.imdata.roiSignals(2).roif,1),125);
zscoredRippleActivity = zeros(size(sData.imdata.roiSignals(2).roif,1),125);
for i = 1:size(sData.imdata.roiSignals(2).roif,1) %runs through the nr of coloumns (corresponding to nr of rois)
    roinr = i;
    v = sData.imdata.roiSignals(2).roif(roinr, :); %creates vector of the frame signals for a particualr roi
    for j = 1:length(sData.ephysdata.frameRipIdx) % creates vector of nr of ripples during recording
        ripplenr = j;
        x = sData.ephysdata.frameRipIdx(ripplenr) - 62; %gives frame nr 62 frames before ripple peak (2s)
        if x < 0
            x = 1;
        end
        y = sData.ephysdata.frameRipIdx(ripplenr) + 62; %gives frame nr 62 frames after ripple peak (2s)
        if y > length(sData.imdata.roiSignals(2).roif)
            y = length(sData.imdata.roiSignals(2).roif);
        end
        if length(x:y) == 125 %if the number of frames from 2s prior to 2s post ripple occurence equals 125 frames, do this
            rippleactivity(j, :) = v(x:y); %index every individual ripple frame segment into a new row for that particular roi(i) )
        elseif length(x:y) < 125 && length(x:y) > 125
            % do nothing
        end
    end
    
    RoifAverage(roinr, :) = nansum(rippleactivity)/length(sData.ephysdata.frameRipIdx); %averages all ripple activity for a roi and assigns it into a new row
    zscoredRippleActivity(roinr, :) = zscore(RoifAverage(roinr, :)); %create zscores for each frame value of the roi 

end

size(RoifAverage);
[~, pt] = sort(sum(RoifAverage, 2));
sortedRoifAverage = RoifAverage(pt, :);
[~,t] = sort(sum(zscoredRippleActivity,2));
sortedZscores = zscoredRippleActivity(t, :);

%figure, imagesc(roifAverage) % colorplot of average individual roi activity during ripples
%xticks(1:31:125)

figure, 
subplot(3,2,1)
imagesc(sortedRoifAverage) % sorted colorplot of average individual roi activity during ripples
xticks(1:31:125)
xticklabels({'-2', '-1', '0', '1', '2'})
xlabel('Time (sec)')
ylabel('Nr of ROIs')
colorbar
title('Individual mean ROI fluorescence')

%subplot(3,2,2)
%imagesc(sortedZscores) % sorted colorplot of zscored average  individual roi activity during ripples
%xticks(1:31:125)
%xticklabels({'-2', '-1', '0', '1', '2'})
%caxis([-3 3])
%xlabel('Time (sec)')
%ylabel('Nr of ROIs')
%colorbar
%title('z-score individual ROI fluorescence')

subplot(3,2,2)
imagesc(zscoredRippleActivity) % sorted colorplot of zscored average  individual roi activity during ripples
xticks(1:31:125)
xticklabels({'-2', '-1', '0', '1', '2'})
caxis([-3 3])
xlabel('Time (sec)')
ylabel('Nr of ROIs')
colorbar
title('z-score individual ROI fluorescence')

MeanRoifAverage = mean(sortedRoifAverage);
ZscoreGrandMeanAllROIs = zscore(MeanRoifAverage);
AverageZscore = mean(sortedZscores);

subplot(3,2,3)
SEM = std(sortedRoifAverage)./sqrt(numel(sortedRoifAverage(:, 1)));
errorbar(MeanRoifAverage, SEM);
xticks(1:31:125)
xticklabels({'-2', '-1', '0', '1', '2'})
xlabel('Time (sec)')
ylabel('Fluorescence')
title('Grand mean ROI fluorescence')

subplot(3,2,4);
SEM = std(sortedZscores)./sqrt(numel(sortedZscores(:,1)));
errorbar(AverageZscore, SEM);
xticks(1:31:125)
xticklabels({'-2', '-1', '0', '1', '2'})
xlabel('Time (sec)')
ylabel('zScore fluorescence')
title('Mean z-score fluorescence')

subplot(3,2,5);
plot(ZscoreGrandMeanAllROIs)
xticks(1:31:125)
xticklabels({'-2', '-1', '0', '1', '2'})
xlabel('Time (sec)')
ylabel('zScore fluorescence')
title('z-score grand mean ROI fluorescence')
