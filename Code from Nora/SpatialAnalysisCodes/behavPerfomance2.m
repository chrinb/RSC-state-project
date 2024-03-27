function sData = behavPerfomance2(sData)

%%% additional calculations, such as first lick, lick median if I discard first 20% of track, and a combined lick score

% get data
nBins = sData.behavior.meta.nBins;
% nTrials = size(sData.behavior.binning.lickBinned,1);
% nBinsPlus = nBins + 10;
% LapLength = sData.behavior.stats.LapLengthCm; % wheel circumference
%RewardZoneSize = 6; % cm
% BinSize = sData.behavior.meta.binSize;
% Settings:
% DeleteLickBeforeCm = 0.2*round(LapLength); % delete licks before 20% length of the track, because sometimes animals licks the reward later than RZ, and it would contaminate lick calculation

% Set Saving - no figures yet
%mkdir(strcat(sData.sessionInfo.savePath,'\Behavior'),'BehavPerformance');
%savePath = strcat(sData.sessionInfo.savePath,'\Behavior\BehavPerformance');

%%% generate binned lick data for every trial until reward is given in each trials
% waterGivenMatrix = sData.behavior.binning.waterGivenBinnedExtended; % matrix, rows trials, columns bins (at the end 10 extra bin is given from the beginning of the next trial)
% waterGivenMatrixGOL = waterGivenMatrix;
% waterGivenMatrixGOL(:,1:10) = 0; % set zero the beginning of trials, I want to see the reward only at the end (+10 bin is given at the end)
% LickMatrixBeforeRew = sData.behavior.binning.lickBinnedExtended;
% waterGivenFirst = NaN(nTrials,1); % sometimes I gave more than one drop of reward, now I search for the first reward drop in each trial's end
% % GOL task
% if max(mean(waterGivenMatrixGOL))> 0.5 % gol task
%     for i = 1:1:nTrials
%         waterGivenFirst(i) = find(waterGivenMatrixGOL(i,:) >= 1,1,'first'); % search for first reward, sometimes more than one reward is given
%         LickMatrixBeforeRew(i,waterGivenFirst(i):end) = 0; % delete licks from that bin when reward given until end of trial, 
%     end
%     % delete licks before DeleteLickBeforeCm (31.4 cm), because sometimes animals licks the reward later than RZ, and it would contaminate lick calculation
%     % then the last 125 cm remain to analyse in each trial (80% of track)
%     LickCalculationStartBin = round(DeleteLickBeforeCm/sData.behavior.meta.binSize) + 1; % 
%     %ExpectedLickBeforeRZStartBin = round(ExpectedLickFromCm/sData.behavior.meta.binSize); % expert animal should lick mostly before reward, like mostly at the last 10% of the track, from 141 cm
%     LickMatrixBeforeRew(:,1:LickCalculationStartBin-1) = 0; % delete licks before that bin from which I want to calculate them
% else % rf task, delete licks after a reward (at same length as in GOL task 'DeleteLickBeforeCm' ) 
%     LickMatrixBeforeRew(waterGivenMatrix > 0) = -1; % delete licks from that bin when reward given until end of trial, 
%     DeletedLicksBins = round(DeleteLickBeforeCm/sData.behavior.meta.binSize);
%     for i = nTrials:-1:1
%         for j = nBinsPlus:-1:1
%             if LickMatrixBeforeRew(i,j) == -1
%                LickMatrixBeforeRew(i,j:(j+DeletedLicksBins-1)) = 0; % discard licks in same amount of bins as in GOL task
%             end
%         end
%     end
%     LickMatrixBeforeRew(2:nTrials,1:DeletedLicksBins) = LickMatrixBeforeRew(1:nTrials-1,nBins+1:nBins+DeletedLicksBins);
%     LickCalculationStartBin = 1;
% end
% 
% %%% Calculate lick performance
% performance2lick = struct;
% 
% % position of the 25, 50, 75% quantiles of the lick distibution
% performance2lick.noteQuantiles = 'Lick were deleted in the first 20% of track, 31,4 cm';
% performance2lick.quantilesBins(nTrials,3) = NaN;
% performance2lick.quantilesCm(nTrials,3) = NaN;
% for i = 1:1:nTrials
%     tempLicks = LickMatrixBeforeRew(i,LickCalculationStartBin:nBins);
%     LicksPositionArray = NaN(1,sum(tempLicks));
%     counter = 1;
%     for j=1:1:nBins-LickCalculationStartBin+1
%         if tempLicks(j)>0
%             for k = 1:1:tempLicks(j)
%                 LicksPositionArray(counter) = j + LickCalculationStartBin -1;
%                 counter = counter + 1;
%             end
%         end
%     end
%     performance2lick.quantilesBins(i,:) = round(quantile(LicksPositionArray,[0.25 0.50 0.75])); %/nBins*100
%     performance2lick.quantilesCm(i,:) = quantile(LicksPositionArray,[0.25 0.50 0.75])/nBins*LapLength; %/nBins*100
% end
% 
% % calculation of mean of 25-50-75 percentiles of lick-position for a session 
% performance2lick.meanQ25 = nanmean(performance2lick.quantilesCm(:,1));
% performance2lick.meanQ50 = nanmean(performance2lick.quantilesCm(:,2));
% performance2lick.meanQ75 = nanmean(performance2lick.quantilesCm(:,3));
% 
% % calculation of first lick in gol task, discard fist 57 cm, so only last 100 cm is checked
% performance2lick.FirstLickNote = 'Licks were deleted before 57 cm, so only last 100 cm was considered';
% %if max(mean(waterGivenMatrixGOL))> 0.5 % gol task
%     FirstLickBin = NaN(nTrials,1);
%     LickCalculationStartBin2 = round(57/sData.behavior.meta.binSize); % 
%     LickMatrixBeforeRew(:,1:LickCalculationStartBin2-1) = 0; 
%     for i = 1:1:nTrials
%        if sum(LickMatrixBeforeRew(i,:))==0
%            continue
%        else
%            FirstLickBin(i) = find(LickMatrixBeforeRew(i,:),1,'first'); 
%        end
%     end
%     performance2lick.FirstLickCm = FirstLickBin * sData.behavior.meta.binSize;
%     performance2lick.MeanFirstLickCm = nanmean(performance2lick.FirstLickCm);
%     performance2lick.FirstLickBin = FirstLickBin;
%     performance2lick.MeanFirstLickBin = nanmean(performance2lick.FirstLickBin);
% %end

% combined lick score: expert behavior is when the animal does not lick at all in the middle of the trial
% e.g. at 50 cm, licks a lot before reward (150-157 cm), and start to lick about 30 cm before reward (100-120 cm)
% I wanted to sample 10 cm-s around 50 cm (50-60cm),  120 cm (120-130cm), 150 cm (147-157)
LickMatrix         = sData.behavior.binning.lickBinned;
BinsToAdd          = round(20/sData.behavior.meta.binSize) - 1; % 10 cm equals with that amount of bins
BinBeginning       = round(40/sData.behavior.meta.binSize); 
% BinMiddle          = round(120/sData.behavior.meta.binSize); 
BinEnd             = nBins - BinsToAdd; 
LickBeginningPerCm = mean( mean( LickMatrix(:, BinBeginning:BinBeginning + BinsToAdd))) / ((BinsToAdd+1)*sData.behavior.meta.binSize);
% LickMiddlePerCm    = mean( mean( LickMatrix(:, BinMiddle:BinMiddle + BinsToAdd))) / ((BinsToAdd+1)*sData.behavior.meta.binSize);
LickEndPerCm       = mean( mean( LickMatrix(:, BinEnd:BinEnd + BinsToAdd))) / ((BinsToAdd+1)*sData.behavior.meta.binSize);

performance2lick.LickBeginningPerCm = LickBeginningPerCm;
% performance2lick.LickMiddlePerCm    = LickMiddlePerCm;
performance2lick.LickEndPerCm       = LickEndPerCm;
% LickScores
performance2lick.LickScore1Note = 'LickEndPerCm / LickBeginningPerCm';
performance2lick.LickScore1 = LickEndPerCm / LickBeginningPerCm; 
% performance2lick.LickScore2Note = 'LickEndPerCm / LickMiddlePerCm';
% performance2lick.LickScore2 = LickEndPerCm / LickMiddlePerCm; 
% performance2lick.LickScore3Note = 'LickEndPerCm / LickMiddlePerCm * (1-LickBeginningPerCm)';
% performance2lick.LickScore3 = LickEndPerCm / LickMiddlePerCm * (1-LickBeginningPerCm); 


% Save file 
sData.behavior.performance2lick = performance2lick;
% save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');
%save(fullfile(savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');

end