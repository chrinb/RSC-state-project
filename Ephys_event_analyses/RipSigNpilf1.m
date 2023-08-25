function [sortedNpilfAverage, MeanNpilfAverage, ZscoreGrandMeanAllROIs]= RipSigNpilf1(sData)

rippleactivity = zeros(length(sData.ephysdata.frameRipIdx),125); % 125 reflects total nr of frames surrounding each ripple (2s), will be different if looking for longer time periods
npilfAverage = zeros(size(sData.imdata.roiSignals(2).npilf,1),125);
zscoredRippleActivity = zeros(size(sData.imdata.roiSignals(2).npilf,1),125);
for i = 1:size(sData.imdata.roiSignals(2).npilf,1) %runs through the nr of coloumns (corresponding to nr of rois)
    roinr = i;
    v = sData.imdata.roiSignals(2).npilf(roinr, :); %creates vector of the frame signals for a particualr roi
    for j = 1:length(sData.ephysdata.frameRipIdx) % creates vector of nr of ripples during recording
        ripplenr = j;
        x = sData.ephysdata.frameRipIdx(ripplenr) - 62; %gives frame nr 62 frames before ripple peak (2s)
        if x < 0
            x = 1;
        end
        y = sData.ephysdata.frameRipIdx(ripplenr) + 62; %gives frame nr 62 frames after ripple peak (2s)
        if y > length(sData.imdata.roiSignals(2).npilf)
            y = length(sData.imdata.roiSignals(2).npilf);
        end
        if length(x:y) == 125 %if the number of frames from 2s prior to 2s post ripple occurence equals 125 frames, do this
            rippleactivity(j, :) = v(x:y); %index every individual ripple frame segment into a new row for that particular roi(i) )
        elseif length(x:y) < 125 && length(x:y) > 125
            % do nothing
        end
    end
    
    npilfAverage(roinr, :) = nansum(rippleactivity)/length(sData.ephysdata.frameRipIdx); %averages all ripple activity for a roi and assigns it into a new row
    zscoredRippleActivity(roinr, :) = zscore(npilfAverage(roinr, :)); %create zscores for each frame value of the roi 

end

size(npilfAverage)
[~, pt] = sort(sum(npilfAverage, 2));
sortedAverage = npilfAverage(pt, :);
[~,t] = sort(sum(zscoredRippleActivity,2));
sortedZscores = zscoredRippleActivity(t, :);

%figure, imagesc(npilfAverage) % colorplot of average individual roi activity during ripples
%xticks(1:31:125)

subplot(2,2,1)
imagesc(sortedAverage) % sorted colorplot of average individual roi activity during ripples
xticks(1:31:125)
xticklabels({'-2', '-1', '0', '1', '2'})
xlabel('Time (sec)')
ylabel('Nr of ROIs')
colorbar
title('Individual ROI Mean Npil')

subplot(2,2,2)
imagesc(sortedZscores) % sorted colorplot of zscored average  individual roi activity during ripples
xticks(1:31:125)
xticklabels({'-2', '-1', '0', '1', '2'})
caxis([-3 3])
xlabel('Time (sec)')
ylabel('Nr of ROIs')
colorbar
title('zScore Npil')

Mean_of_AverageROI = mean(sortedAverage);
ZscoreAverage = zscore(Mean_of_AverageROI);

subplot(2,2,3)
Errorbars = std(sortedAverage);
errorbar(Mean_of_AverageROI, Errorbars);
xticks(1:31:125)
xticklabels({'-2', '-1', '0', '1', '2'})
xlabel('Time (sec)')
ylabel('Npil')
title('All ROI Mean Npil')

subplot(2,2,4);
Errorbars1 = std(sortedZscores);
errorbar(ZscoreAverage, Errorbars1)
xticks(1:31:125)
xticklabels({'-2', '-1', '0', '1', '2'})
xlabel('Time (sec)')
ylabel('zScore Npil')
title('Mean zScore Npil')