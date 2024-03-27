function sData = behavPerfomance(sData)

% BAsed on animals behavior it seemed that they accomodate to a new task quickly usually within 8-10 laps, maximum 15 laps. 
% So I called the behavioral change within the early 1-15 laps accomodation , and later changes development (e.g. between laps 20-30 and 90-100).

% the animal  performs well the task, if:
% 50% of the licks should be after BinStartLickBeforeReward=70 bin. It means, that the animal licks 10 times more/bin in bin 70-80 compared to bin 20-70
% check when it happens in 3-5-10 consecutive trials

% the velocity in the first bin of the RZ should be less than 25 percent of max speed, if min speed in the trial is considered 0
% check when it happens in 3-5-10 consecutive trials


% get data
nBins = sData.behavior.meta.nBins;
[nTrials, nBinsPlus] = size(sData.behavior.binning.lickBinned);
LapLength = sData.behavior.stats.LapLengthCm; % wheel circumference
%RewardZoneSize = 6; % cm
% BinSize = sData.behavior.meta.binSize;
% DeleteLickBeforeCm = 10 * BinSize; %It is 19.6350 for 80 bins, I want to discard licks before it
% ExpectedLickFromCm = 73 * BinSize; % It is 143.3352 for 80 bins, I expect at least 50% of the licks here after learning
% Settings:
AccomodationLaps = 15; % in the beginning of GOL task animal behavior changes (accomodate to the task), I want to calculate steady performance after this period
DevelopmentTestLaps = 10; % After the accomodation I will check the development during the session in every 10 laps
DeleteLickBeforeCm = 15.7; % delete licks before 10% length of the track, because sometimes animals licks the reward later than RZ, and it would contaminate lick calculation
ExpectedLickFromCm = 141.3; % I expect most of the lick in the last 10% of the track (157 x 0.9 = 141.3) 
 % 25 percent velocity in speed range. Speed range= maxspeed - minspeed in a trial


% Set Saving - no figures yet
%mkdir(strcat(sData.sessionInfo.savePath,'\Behavior'),'BehavPerformance');
%savePath = strcat(sData.sessionInfo.savePath,'\Behavior\BehavPerformance');

%%% generate binned lick data for every trial until reward is given in each trials
waterGivenMatrix = sData.behavior.binning.waterGivenBinned; % matrix, rows trials, columns bins (at the end 10 extra bin is given from the beginning of the next trial)
waterGivenMatrix(:,1:10) = 0; % set zero the beginning of trials, I want to see the reward only at the end (+10 bin is given at the end)
LickMatrixBeforeRew = sData.behavior.binning.lickBinned;
waterGivenFirst = NaN(nTrials,1); % sometimes I gave more than one drop of reward, now I search for the first reward drop in each trial's end
for i = 1:1:nTrials
    if max(mean(waterGivenMatrix))> 0.5 % gol task
        waterGivenFirst(i) = find(waterGivenMatrix(i,:) >= 1,1,'first');
    else
        waterGivenFirst(i) = 81; % in RF trials lets set reward position at the end of the track
    end
    LickMatrixBeforeRew(i,waterGivenFirst(i):end) = 0; % delete licks from that bin when reward given until end of trial
end
% delete licks before 15.7 cm, because sometimes animals licks the reward later than RZ, and it would contaminate lick calculation
% then the last 141 cm remain to analyse in each trial (90% of track)
LickCalculationStartBin = round(DeleteLickBeforeCm/sData.behavior.meta.binSize) + 1; % I set 17 because then I will calculate the licks between 157-17 cm = 140 cm
ExpectedLickBeforeRZStartBin = round(ExpectedLickFromCm/sData.behavior.meta.binSize); % expert animal should lick mostly before reward, like mostly at the last 10% of the track, from 141 cm
LickMatrixBeforeRew(:,1:LickCalculationStartBin-1) = 0; % delete licks before that bin from which I want to calculate them
%nBinsFirstHalf = ExpectedLickBeforeRZStartBin - 1 - LickCalculationStartBin;
%nBinsSecondHalf = waterGivenFirst - ExpectedLickBeforeRZStartBin;

%%% Calculate lick performance
performance = struct;
performance.lick = struct;

% position of the 25, 50, 75% quantiles of the lick distibution
performance.lick.quantilesBins(nTrials,3) = NaN;
performance.lick.quantilesCm(nTrials,3) = NaN;
for i = 1:1:nTrials
    tempLicks = LickMatrixBeforeRew(i,LickCalculationStartBin:nBins);
    LicksPositionArray = NaN(1,sum(tempLicks));
    counter = 1;
    for j=1:1:nBins-LickCalculationStartBin+1
        if tempLicks(j)>0
            for k = 1:1:tempLicks(j)
                LicksPositionArray(counter) = j + LickCalculationStartBin -1;
                counter = counter + 1;
            end
        end
    end
    performance.lick.quantilesBins(i,:) = round(quantile(LicksPositionArray,[0.25 0.50 0.75])); %/nBins*100
    performance.lick.quantilesCm(i,:) = quantile(LicksPositionArray,[0.25 0.50 0.75])/nBins*LapLength; %/nBins*100
end
% additional calculation from 2021.03.25
% calculation of 25-50-75 percentiles of mean of meadian position of licks for a session 
performance.lick.meanQ25 = nanmean(performance.lick.quantilesCm(:,1));
performance.lick.meanQ50 = nanmean(performance.lick.quantilesCm(:,2));
performance.lick.meanQ75 = nanmean(performance.lick.quantilesCm(:,3));
% 50% (the median) of the licks should be after the bin what I set , now it is last 10% of track (bin 72). It means approximately that the animal licks 10 times more/bin in this bins compared to first part 
performance.lick.Q50.note = 'Expected licks: Q50 of licks is in the last 10% of the track, licks are deleted between 0-15.7 cm';

% Which trial is the 3/5/10th when the animal 50 quantile is at >= as expected bin?
performance.lick.Q50.TrialsLickExpectedBin = find(performance.lick.quantilesBins(:,2) >= ExpectedLickBeforeRZStartBin);
performance.lick.Q50.TrialsLickExpectedBin(length(performance.lick.Q50.TrialsLickExpectedBin)+1:nTrials) = NaN;
performance.lick.Q50.details.NotConsTrialsExp3 = performance.lick.Q50.TrialsLickExpectedBin(3);
performance.lick.Q50.details.NotConsTrialsExp5 = performance.lick.Q50.TrialsLickExpectedBin(5);
performance.lick.Q50.details.NotConsTrialsExp10 = performance.lick.Q50.TrialsLickExpectedBin(10);
% Find the trial when it is true that in 3/5/7 consecutive trials 50Q is larger= than expected bin
counterArray = NaN(10,1);
counter = 0;
performance.lick.Q50.details.ConsTrialsExp3 = NaN;
performance.lick.Q50.details.ConsTrialsExp5 = NaN;
performance.lick.Q50.details.ConsTrialsExp10 = NaN;
for i = 1:1:nTrials
    if performance.lick.quantilesBins(i,2) >= ExpectedLickBeforeRZStartBin
       counter = counter + 1;
       counterArray(counter,1) = i;
       if counter == 3 && isnan(performance.lick.Q50.details.ConsTrialsExp3) 
           performance.lick.Q50.details.ConsTrialsExp3 = i;
       elseif counter == 5 && isnan(performance.lick.Q50.details.ConsTrialsExp5)
           performance.lick.Q50.details.ConsTrialsExp5 = i;
       elseif counter == 10 && isnan(performance.lick.Q50.details.ConsTrialsExp10)
           performance.lick.Q50.details.ConsTrialsExp10 = i;
       end
    else 
       counterArray = NaN(10,1);
       counter = 0; 
    end
end

% percentage of trials in which Q50 is above expected bin
performance.lick.Q50.PercentageOfTrialsLickExp = sum(~isnan(performance.lick.Q50.TrialsLickExpectedBin))/nTrials;
% as first 8-15 laps usualy needed for accomodation to the new task , calculate percentage of Q50 from trial 16
TrialsLickExpectedAfterLap15 = performance.lick.Q50.TrialsLickExpectedBin(performance.lick.Q50.TrialsLickExpectedBin>AccomodationLaps);
performance.lick.Q50.PercentageOfTrialsLickExpAfterT15 = length(TrialsLickExpectedAfterLap15)/(nTrials-AccomodationLaps);
% calculate the change in lap 1-15 : Accomodation (3-3 laps)
performance.lick.Q50.details.Accomodation_meanQ50Lap1_3 = nanmean(performance.lick.quantilesBins(1:3,2));
performance.lick.Q50.details.Accomodation_meanQ50Lap12_15 = nanmean(performance.lick.quantilesBins(AccomodationLaps-3:AccomodationLaps,2));
performance.lick.Q50.details.Accomodation_Lap1215_subtr_Lap13 = performance.lick.Q50.details.Accomodation_meanQ50Lap12_15 - performance.lick.Q50.details.Accomodation_meanQ50Lap1_3;
performance.lick.Q50.details.Accomodation_Lap1215_per_Lap13 = performance.lick.Q50.details.Accomodation_meanQ50Lap12_15 / performance.lick.Q50.details.Accomodation_meanQ50Lap1_3;
% calculate the change after lap 15 : Development (10-10 laps) DevelopmentTestLaps
performance.lick.Q50.details.Development_meanQ50Lap16_25 = nanmean(performance.lick.quantilesBins(AccomodationLaps+1:AccomodationLaps+DevelopmentTestLaps,2));
performance.lick.Q50.details.Development_meanQ50Last10Lap = nanmean(performance.lick.quantilesBins(nTrials-DevelopmentTestLaps:nTrials,2));
performance.lick.Q50.details.Development_Last10Lap_subtr_Lap16_25 = performance.lick.Q50.details.Development_meanQ50Last10Lap - performance.lick.Q50.details.Development_meanQ50Lap16_25;
performance.lick.Q50.details.Development_Last10Lap_per_Lap16_25 = performance.lick.Q50.details.Development_meanQ50Last10Lap / performance.lick.Q50.details.Development_meanQ50Lap16_25;


%%% Velocitiy performance
performance.velo = struct;

% calcuate the max and min speed and its position in each trial, discard the first part of the trials, as in the Lick performance calculation (first 19 cm. This belongs to the previous trial)
[performance.velo.veloMax, performance.velo.veloMaxPosBin]= max(sData.behavior.binning.veloBinned(1:nTrials,LickCalculationStartBin:nBins),[],2);
[performance.velo.veloMin, performance.velo.veloMinPosBin] = min(sData.behavior.binning.veloBinned(1:nTrials,LickCalculationStartBin:nBinsPlus),[],2);

% calculate velocity at Reward Zone(RZ) first bin (bin1 or 81) on a max-min scale (where 0 = the min velo is at RZ and and after 10 bins, 1= the max velo in this trial)
performance.velo.BeforeRZspeed.veloBeforeRZperc = NaN(nTrials,1); 
performance.velo.BeforeRZspeed.veloBeforeRZCms(1:nTrials,1) = sData.behavior.binning.veloBinned(:,nBins);
performance.velo.BeforeRZspeed.veloBeforeRZperc(1:nTrials,1) = (sData.behavior.binning.veloBinned(:,nBins) - performance.velo.veloMin) ./ (performance.velo.veloMax-performance.velo.veloMin)*100;

% What is the percentge when the speed goes below 10/25/50% of speed range before the reward?
performance.velo.BeforeRZspeed.veloLessThan50percTrials = NaN(10,1);
TempArray = find(performance.velo.BeforeRZspeed.veloBeforeRZperc(:,1) < 50);
performance.velo.BeforeRZspeed.veloLessThan50percTrials(1:length(TempArray)) = find(performance.velo.BeforeRZspeed.veloBeforeRZperc(:,1) < 50);
performance.velo.BeforeRZspeed.veloLessThan50percPerc = length(performance.velo.BeforeRZspeed.veloLessThan50percTrials)/nTrials;

performance.velo.BeforeRZspeed.veloLessThan25percTrials = NaN(10,1);
TempArray = find(performance.velo.BeforeRZspeed.veloBeforeRZperc(:,1) < 25);
performance.velo.BeforeRZspeed.veloLessThan25percTrials(1:length(TempArray)) = find(performance.velo.BeforeRZspeed.veloBeforeRZperc(:,1) < 25);
performance.velo.BeforeRZspeed.veloLessThan25percPerc = length(performance.velo.BeforeRZspeed.veloLessThan25percTrials)/nTrials;

performance.velo.BeforeRZspeed.veloLessThan10percTrials = NaN(10,1);
TempArray = find(performance.velo.BeforeRZspeed.veloBeforeRZperc(:,1) < 10);
performance.velo.BeforeRZspeed.veloLessThan10percTrials(1:length(TempArray)) = find(performance.velo.BeforeRZspeed.veloBeforeRZperc(:,1) < 10);
performance.velo.BeforeRZspeed.veloLessThan10percPerc = length(TempArray)/nTrials;

% only after lap 15
TrialsSpeedLess50AfterLap15 = performance.velo.BeforeRZspeed.veloLessThan50percTrials(performance.velo.BeforeRZspeed.veloLessThan50percTrials>AccomodationLaps);
TrialsSpeedLess25AfterLap15 = performance.velo.BeforeRZspeed.veloLessThan25percTrials(performance.velo.BeforeRZspeed.veloLessThan25percTrials>AccomodationLaps);
TrialsSpeedLess10AfterLap15 = performance.velo.BeforeRZspeed.veloLessThan10percTrials(performance.velo.BeforeRZspeed.veloLessThan10percTrials>AccomodationLaps);
performance.velo.veloLessThan50percPercAfterT15 = length(TrialsSpeedLess50AfterLap15)/(nTrials-AccomodationLaps)*100;
performance.velo.veloLessThan25percPercAfterT15 = length(TrialsSpeedLess25AfterLap15)/(nTrials-AccomodationLaps)*100;
performance.velo.veloLessThan10percPercAfterT15 = length(TrialsSpeedLess10AfterLap15)/(nTrials-AccomodationLaps)*100;

% Which trial is the 3/5/10th when at the RZ the animal speed is below 10/25/50 percent of max speed?
performance.velo.BeforeRZspeed.NotConsTrials_P50_3th = performance.velo.BeforeRZspeed.veloLessThan50percTrials(3);
performance.velo.BeforeRZspeed.NotConsTrials_P50_5th = performance.velo.BeforeRZspeed.veloLessThan50percTrials(5);
performance.velo.BeforeRZspeed.NotConsTrials_P50_10th = performance.velo.BeforeRZspeed.veloLessThan50percTrials(10);

performance.velo.BeforeRZspeed.NotConsTrials_P25_3th = performance.velo.BeforeRZspeed.veloLessThan25percTrials(3);
performance.velo.BeforeRZspeed.NotConsTrials_P25_5th = performance.velo.BeforeRZspeed.veloLessThan25percTrials(5);
performance.velo.BeforeRZspeed.NotConsTrials_P25_10th = performance.velo.BeforeRZspeed.veloLessThan25percTrials(10);

performance.velo.BeforeRZspeed.NotConsTrials_P10_3th = performance.velo.BeforeRZspeed.veloLessThan10percTrials(3);
performance.velo.BeforeRZspeed.NotConsTrials_P10_5th = performance.velo.BeforeRZspeed.veloLessThan10percTrials(5);
performance.velo.BeforeRZspeed.NotConsTrials_P10_10th = performance.velo.BeforeRZspeed.veloLessThan10percTrials(10);

% Find the trial when it is true that in 3/5/10 consecutive trials the speed is below 10/25/50 percent of max-min before RZ
% if speed is below 50% at RZ
counterArray = NaN(10,1);
counter = 0;
performance.velo.BeforeRZspeed.ConsTrials_P50_3th = NaN;
performance.velo.BeforeRZspeed.ConsTrials_P50_5th = NaN;
performance.velo.BeforeRZspeed.ConsTrials_P50_10th = NaN;
for i = 1:1:nTrials
    if performance.velo.BeforeRZspeed.veloBeforeRZperc(i,1) < 50
       counter = counter + 1;
       counterArray(counter,1) = i;
       if counter == 3 && isnan(performance.velo.BeforeRZspeed.ConsTrials_P50_3th)
           performance.velo.BeforeRZspeed.ConsTrials_P50_3th = i;
       elseif counter == 5 && isnan(performance.velo.BeforeRZspeed.ConsTrials_P50_5th)
           performance.velo.BeforeRZspeed.ConsTrials_P50_5th = i;
       elseif counter == 10 && isnan(performance.velo.BeforeRZspeed.ConsTrials_P50_10th)
           performance.velo.BeforeRZspeed.ConsTrials_P50_10th = i;
       end
    else 
       counterArray = NaN(10,1);
       counter = 0; 
    end
end
% if speed is below 25% at RZ
counterArray = NaN(10,1);
counter = 0;
performance.velo.BeforeRZspeed.ConsTrials_P25_3th = NaN;
performance.velo.BeforeRZspeed.ConsTrials_P25_5th = NaN;
performance.velo.BeforeRZspeed.ConsTrials_P25_10th = NaN;
for i = 1:1:nTrials
    if performance.velo.BeforeRZspeed.veloBeforeRZperc(i,1) < 25
       counter = counter + 1;
       counterArray(counter,1) = i;
       if counter == 3 && isnan(performance.velo.BeforeRZspeed.ConsTrials_P25_3th)
           performance.velo.BeforeRZspeed.ConsTrials_P25_3th = i;
       elseif counter == 5 && isnan(performance.velo.BeforeRZspeed.ConsTrials_P25_5th)
           performance.velo.BeforeRZspeed.ConsTrials_P25_5th = i;
       elseif counter == 10 && isnan(performance.velo.BeforeRZspeed.ConsTrials_P25_10th)
           performance.velo.BeforeRZspeed.ConsTrials_P25_10th = i;
       end
    else 
       counterArray = NaN(10,1);
       counter = 0; 
    end
end
% if speed is below 10% at RZ
counterArray = NaN(10,1);
counter = 0;
performance.velo.BeforeRZspeed.ConsTrials_P10_3th = NaN;
performance.velo.BeforeRZspeed.ConsTrials_P10_5th = NaN;
performance.velo.BeforeRZspeed.ConsTrials_P10_10th = NaN;
for i = 1:1:nTrials
    if performance.velo.BeforeRZspeed.veloBeforeRZperc(i,1) < 10
       counter = counter + 1;
       counterArray(counter,1) = i;
       if counter == 3 && isnan(performance.velo.BeforeRZspeed.ConsTrials_P10_3th)
           performance.velo.BeforeRZspeed.ConsTrials_P10_3th = i;
       elseif counter == 5 && isnan(performance.velo.BeforeRZspeed.ConsTrials_P10_5th)
           performance.velo.BeforeRZspeed.ConsTrials_P10_5th = i;
       elseif counter == 10 && isnan(performance.velo.BeforeRZspeed.ConsTrials_P10_10th)
           performance.velo.BeforeRZspeed.ConsTrials_P10_10th = i;
       end
    else 
       counterArray = NaN(10,1);
       counter = 0; 
    end
end

% claculate the change in lap 1-15 : Accomodation (3-3 laps)
performance.velo.BeforeRZspeed.Accomodation_VeloLap1_3 = nanmean(performance.velo.BeforeRZspeed.veloBeforeRZperc(1:3,1));
performance.velo.BeforeRZspeed.Accomodation_VeloLap12_15 = nanmean(performance.velo.BeforeRZspeed.veloBeforeRZperc(AccomodationLaps-3:AccomodationLaps,1));
performance.velo.BeforeRZspeed.Accomodation_VeloLap1215_subtr_Lap13 = nanmean(performance.velo.BeforeRZspeed.veloBeforeRZperc(1:3,1)) - nanmean(performance.velo.BeforeRZspeed.veloBeforeRZperc(AccomodationLaps-3:AccomodationLaps,1));
performance.velo.BeforeRZspeed.Accomodation_VeloLap1215_div_Lap13 = nanmean(performance.velo.BeforeRZspeed.veloBeforeRZperc(1:3,1)) ./ nanmean(performance.velo.BeforeRZspeed.veloBeforeRZperc(AccomodationLaps-3:AccomodationLaps,1));

% claculate the change after lap 15 : Development (10-10 laps) DevelopmentTestLaps
performance.velo.BeforeRZspeed.Development_VeloLap16_25 = nanmean(performance.velo.BeforeRZspeed.veloBeforeRZperc(AccomodationLaps+1:AccomodationLaps+DevelopmentTestLaps,1));
performance.velo.BeforeRZspeed.Development_VeloLast10Lap = nanmean(performance.velo.BeforeRZspeed.veloBeforeRZperc(nTrials-DevelopmentTestLaps:nTrials,1));
performance.velo.BeforeRZspeed.Development_VeloLast10_subtr_Lap1625 = performance.velo.BeforeRZspeed.Development_VeloLast10Lap - performance.velo.BeforeRZspeed.Development_VeloLap16_25;
performance.velo.BeforeRZspeed.Development_VeloLast10_div_Lap1625 = performance.velo.BeforeRZspeed.Development_VeloLast10Lap / performance.velo.BeforeRZspeed.Development_VeloLap16_25;


%%% Speed decelaration calculation: at which bin does the animal reach the half of the max speed ? (max-min/2)
meanSpeed = mean(sData.behavior.binning.veloBinned,1); % mean speed profile of session
% discard speed in the very beginning of the trial as in licks
[meanSpeedMax, meanSpeedMaxPos] = max(meanSpeed(1,LickCalculationStartBin:nBins),[],2);
[meanSpeedMin, ~]  = min(meanSpeed(1,LickCalculationStartBin:nBinsPlus),[],2);
meanSpeed(1:meanSpeedMaxPos-1) = NaN; % delete speed data before max
meanSpeed(nBins:nBinsPlus) = NaN; % delete speed data after reward
halfspeed = (meanSpeedMax - meanSpeedMin)/2 + meanSpeedMin;
performance.velo.HalfspeedBin = NaN;
[~,closestIndex] = min(abs(meanSpeed-halfspeed)); % search the bin in which the velo is the most similar to the half speed

performance.velo.HalfspeedCm = closestIndex*sData.behavior.meta.binSize;

% Save file 
sData.behavior.performance = performance;
% save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');
%save(fullfile(savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');

end