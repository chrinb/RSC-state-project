function [zscoredRippleActivity, ciaDeconvolvedAverage, ciaDeconvWeightedAvg] = RipSigciaDeconv(sData, NoRunRipples)

%If NoRunRipples is set to 1 when calling the function, this function will 
%only include ripple segments that has 0 running speed (determined by the function RipRun)
% If sData is the only input argument, all ripple segments will be analyzed
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
ciaDeconvolvedAverage = zeros(size(sData.imdata.roiSignals(2).ciaDeconvolved,1),125);
zscoredRippleActivity = zeros(size(sData.imdata.roiSignals(2).ciaDeconvolved,1),125);
for i = 1:size(sData.imdata.roiSignals(2).ciaDeconvolved,1) %runs through the nr of coloumns (corresponding to nr of rois)
    roinr = i;
    v = sData.imdata.roiSignals(2).ciaDeconvolved(roinr, :); %creates vector of the frame signals for a particualr roi
    for j = 1:length(RippleIdx) % creates vector of nr of ripples during recording
        ripplenr = j;
        x = RippleIdx(ripplenr) - 62; %gives frame nr 62 frames before ripple peak (2s)
        if x < 0
            x = 1;
        end
        y = RippleIdx(ripplenr) + 62; %gives frame nr 62 frames after ripple peak (2s)
        if y > length(sData.imdata.roiSignals(2).ciaDeconvolved)
            y = length(sData.imdata.roiSignals(2).ciaDeconvolved);
        end
        if length(x:y) == 125 %if the number of frames from 2s prior to 2s post ripple occurence equals 125 frames, do this
            rippleactivity(j, :) = v(x:y); %index every individual ripple frame segment into a new row for that particular roi(i) )
        elseif length(x:y) < 125 && length(x:y) > 125
            % do nothing
        end
    end
    
    ciaDeconvolvedAverage(roinr, :) = nanmean(rippleactivity); %averages all ripple activity for a roi and assigns it into a new row
    zscoredRippleActivity(roinr, :) = zscore(ciaDeconvolvedAverage(roinr, :)); %create zscores for each frame value of the roi 
    
    ciaDeconvWeightedAvg(roinr) = {[rippleactivity]};
end

ripplenr

[~, pt] = sort(sum(ciaDeconvolvedAverage, 2));
sortedciaDeconvolvedAverage = ciaDeconvolvedAverage(pt, :);
[~,t] = sort(sum(zscoredRippleActivity,2));
sortedZscores = zscoredRippleActivity(t, :);

figure,
subplot(2,2,1)
imagesc(ciaDeconvolvedAverage)
xticks(1:31:125)
xticklabels({'-2', '-1', '0', '1', '2'})
ylabel('ROI #')
colorbar
xlabel('Time from ripple peak (sec)')
title('Mean deconvolved Ca2+')

%subplot(3, 2, 2);
%imagesc(sortedZscores) % sorted colorplot of zscored average  individual roi activity during ripples
%xticks(1:31:125)
%xticklabels({'-2', '-1', '0', '1', '2'})
%caxis([-3 3])
%xlabel('Time from ripple peak (sec)')
%ylabel('ROI #')
%title('z-score ciaDeconvolved')
%colorbar

subplot(2, 2, 2);
imagesc(zscoredRippleActivity) 
colormap('jet')
xticks(1:31:125)
xticklabels({'-2', '-1', '0', '1', '2'})
xlabel('Time from ripple peak (sec)')
ylabel('ROI #')
title('z-scored mean deconvolved Ca2+')
colorbar

%RippleAligned = sortedZscores(sortedZscores(:, 62) >2 | sortedZscores(:, 63) >2 |sortedZscores(:, 64) >2 | sortedZscores(:, 65) >2 | sortedZscores(:, 66) >2| sortedZscores(:, 67) >2| sortedZscores(:, 67) >2,:);

%subplot(3,2,5);
%imagesc(RippleAligned)
%xticks(1:31:125)
%xticklabels({'-2', '-1', '0', '1', '2'})
%xlabel('Time from ripple peak (sec)')
%ylabel('ROI #')
%title('z-score Active ciaDeconvolved')
%caxis([-3 3])
%colorbar

MeanciaDeconvolvedAverage = nanmean(ciaDeconvolvedAverage);
ZscoreGrandMeanAllROIs = zscore(MeanciaDeconvolvedAverage);
AverageZscore = nanmean(sortedZscores);

subplot(2,2,3)
shadedErrorBar([],MeanciaDeconvolvedAverage, std(rmmissing(ciaDeconvolvedAverage))./sqrt(numel(rmmissing(ciaDeconvolvedAverage(:, 1)))),'lineprops', 'b');
xticks(1:31:125)
xticklabels({'-2', '-1', '0', '1', '2'})
xlabel('Time from ripple peak (sec)')
ylabel('Deconvolved Ca2+')
title('Grand mean deconvolved Ca2+')

subplot(2, 2, 4);
shadedErrorBar([], AverageZscore, std(rmmissing(sortedZscores))./sqrt(numel(rmmissing(sortedZscores(:, 1)))),'lineprops', 'b');
xticks(1:31:125)
xticklabels({'-2', '-1', '0', '1', '2'})
xlabel('Time from ripple peak (sec)')
ylabel('z-scored Deconvolved Ca2+')
title('Grand mean z-scored Deconvolved Ca2+')

%subplot(3, 2, 5);
%plot(ZscoreGrandMeanAllROIs);
%xticks(1:31:125)
%xticklabels({'-2', '-1', '0', '1', '2'})
%xlabel('Time from ripple peak (sec)')
%ylabel('zScore ciaDeconvolved')
%title('z-score grand mean ciaDeconvolved')