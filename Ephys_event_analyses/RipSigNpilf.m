function [sortednpilfAverage, sortedZscores]= RipSigNpilf(sData, NoRunRipples)

NoRunRipples = 1; 

if nargin < 2
    NoRunRipples = false;
end

if NoRunRipples
RippleIdx = RipRun(sData);
else
    RippleIdx = removerip(sData); %remove ripples occuring within 4s of anothe
end

rippleactivity = zeros(length(RippleIdx),125); % 125 reflects total nr of frames surrounding each ripple (2s), will be different if looking for longer time periods
npilfAverage = zeros(size(sData.imdata.roiSignals(2).npilMediF,1),125);
zscoredRippleActivity = zeros(size(sData.imdata.roiSignals(2).npilMediF,1),125);
for i = 1:size(sData.imdata.roiSignals(2).npilMediF,1) %runs through the nr of coloumns (corresponding to nr of rois)
    roinr = i;
    v = sData.imdata.roiSignals(2).npilMediF(roinr, :); %creates vector of the frame signals for a particualr roi
    for j = 1:length(RippleIdx) % creates vector of nr of ripples during recording
        ripplenr = j;
        x = RippleIdx(ripplenr) - 62; %gives frame nr 62 frames before ripple peak (2s)
        if x < 0
            x = 1;
        end
        y = RippleIdx(ripplenr) + 62; %gives frame nr 62 frames after ripple peak (2s)
        if y > length(sData.imdata.roiSignals(2).npilMediF)
            y = length(sData.imdata.roiSignals(2).npilMediF);
        end
        if length(x:y) == 125 %if the number of frames from 2s prior to 2s post ripple occurence equals 125 frames, do this
            rippleactivity(j, :) = v(x:y); %index every individual ripple frame segment into a new row for that particular roi(i) )
        elseif length(x:y) < 125 && length(x:y) > 125
            % do nothing
        end
    end
    
    npilfAverage(roinr, :) = nanmean(rippleactivity); %averages all ripple activity for a roi and assigns it into a new row
    zscoredRippleActivity(roinr, :) = zscore(npilfAverage(roinr, :)); %create zscores for each frame value of the roi 

end

ripplenr

[~, pt] = sort(sum(npilfAverage, 2));
sortednpilfAverage = npilfAverage(pt, :);
[~,t] = sort(sum(zscoredRippleActivity,2));
sortedZscores = zscoredRippleActivity(t, :);

%figure, imagesc(npilfAverage) % colorplot of average individual roi activity during ripples
%xticks(1:31:125)

figure, 
subplot(2,2,1)
imagesc(npilfAverage) % sorted colorplot of average individual roi activity during ripples
xticks(1:31:125)
xticklabels({'-2', '-1', '0', '1', '2'})
xlabel('Time from ripple peak (sec)')
ylabel('ROI #')
colorbar
title('Individual ROI Mean npilf')

%subplot(3,2,2)
%imagesc(sortedZscores) % sorted colorplot of zscored average  individual roi activity during ripples
%xticks(1:31:125)
%xticklabels({'-2', '-1', '0', '1', '2'})
%caxis([-3 3])
%xlabel('Time from ripple peak (sec)')
%ylabel('ROI #')
%colorbar
%title('z-score individual ROI npilf')

subplot(2,2,2)
imagesc(zscoredRippleActivity) % sorted colorplot of zscored average  individual roi activity during ripples
xticks(1:31:125)
xticklabels({'-2', '-1', '0', '1', '2'})
caxis([-3 3])
xlabel('Time from ripple peak (sec)')
ylabel('ROI #')
colorbar
title('Individual ROI z-scored npilf')

MeannpilfAverage = nanmean(sortednpilfAverage);
ZscoreGrandMeanAllROIs = zscore(MeannpilfAverage);
AverageZscore = nanmean(sortedZscores);

subplot(2,2,3)
SEM = std(rmmissing(sortednpilfAverage))./sqrt(numel(rmmissing(sortednpilfAverage(:, 1))));
shadedErrorBar([],MeannpilfAverage, SEM,'lineprops', 'b');
xticks(1:31:125)
xticklabels({'-2', '-1', '0', '1', '2'})
xlabel('Time from ripple peak (sec)')
ylabel('npilf')
title('Grand Mean ROI npilf')

subplot(2,2,4);
SEM = std(rmmissing(sortedZscores))./sqrt(numel(rmmissing(sortedZscores(:,1))));
shadedErrorBar([],AverageZscore, SEM,'lineprops', 'b');
xticks(1:31:125)
xticklabels({'-2', '-1', '0', '1', '2'})
xlabel('Time from ripple peak (sec)')
ylabel('zScore npilf')
title('Mean z-scored npilf')

%subplot(3,2,5);
%plot(ZscoreGrandMeanAllROIs)
%xticks(1:31:125)
%xticklabels({'-2', '-1', '0', '1', '2'})
%xlabel('Time from ripple peak (sec)')
%ylabel('zScore npilf')
%title('z-score grand mean ROI npilf')