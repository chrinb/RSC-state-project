function sData = SA_CalcBehav(sData,nBins) 

% script to calculate behavior data for the session in sData.behavior struct

% SET:
% nBins = 80; number of bins for spatial binning
WheelDiameterCm    = 48.7; 
MovMean            = 7 ; % size of moving window for smoothing velocity curve
SizeOfRewardZoneCm = 6; % in cm

sData.behavior.meta.nBins   = nBins;
sData.behavior.meta.binSize = WheelDiameterCm*pi/nBins;        %Bin size in Cm when binning data spatially. 

% Calculate Trial number (TRNu) and RewardStart array (array of zeros and ones, one when lap starts)
% Note: having 3 khz sampling, if the mouse run 100 cm/s, the photodiode signal (2cm black part) takes 20 ms, it is 60 samples. 
PhotoDiode = sData.daqdataB.wheelDiode; % original analog data from tdms file
if any(PhotoDiode(:) > 2.0)
    PhotoDiode(PhotoDiode < 1.5)  = 0; % binarized to zeros and ones
    PhotoDiode(PhotoDiode >= 1.5) = 1;
    PdZeroOne                     = diff(PhotoDiode) == 1; % search when PD signal change from zero to one - potential lap starts (and arteficial back-and-forth as well)
else
    PdZeroOne = PhotoDiode;
end

PotRewardStartIndex = find(PdZeroOne); % potential lap starts
RewardStartIndexPre = NaN(numel(PotRewardStartIndex),1); % real lap-starts
LapCounter          = 0;

for i = 2:1:numel(PotRewardStartIndex)  % the first 0-1 PD signal assumed to be real (but does not really matter if not)
    if sData.daqdataB.distanceCm(PotRewardStartIndex(i)) - sData.daqdataB.distanceCm(PotRewardStartIndex(i-1)) > 150 % lapstart if photodiode signal is 1 and animal runned at least 150 cm before previous reward
        LapCounter                      = LapCounter + 1; % sum up already runned laps
        RewardStartIndexPre(LapCounter) = PotRewardStartIndex(i);
    end
end
sData.behavior.details.rewardStartIndices = RewardStartIndexPre(~isnan(RewardStartIndexPre));
sData.behavior.wheelLap                   = LapCounter;


% Calculate number of imaging frames based on frame signal recorded by LabView
% generate fake frame signal (needed if there was no imaging but I want to downsample data the same way)
FrameStart = diff(sData.daqdataB.frameSignal)== 1; % array with ones if the frame starts
if sum(FrameStart) == 0
    msgbox('No frame signal');
    FrameSignalDur                                  = 98; %Samples
    nSamples                                        = size(sData.daqdataB.distanceCm,1);
    Samples                                         = 1:1:nSamples;
    FakeFrameSignal                                 = zeros(nSamples,1);
    FakeFrameSignal(rem(Samples,FrameSignalDur)==0) = 1;
    FrameStart                                      = diff(FakeFrameSignal)== 1;
    FakeFrameSignal2                                = 'yes';
end

sData.behavior.details.nFrames = sum(FrameStart); % frame number = sampleNumber usually behav.SampleNu = behav.FRNu;

% calculate number of frame signals and imaging frequency from frame signals  
FrameStartIndex                          = find(FrameStart); % array, indices when frame Starts , returns with the indices when FrameStart is not zero
sData.behavior.details.frameStartIndices = FrameStartIndex;
sData.behavior.meta.daqSamplingRate      = sData.daqdataB.meta.fs;
sData.behavior.meta.imagingSamplingRate  = 1/((mean(diff(FrameStartIndex)))/sData.behavior.meta.daqSamplingRate); % frame frequency. Calculate mean sample number between frames, take the reciproc to get Hz

% discard recording before first reward and before frame signal starts
Distance = sData.daqdataB.distanceCm;
behavStart = min(sData.behavior.details.rewardStartIndices(sData.behavior.details.rewardStartIndices > FrameStartIndex(1)));
Distance(1:behavStart-1) = NaN;

% Aim: Calculate Absolute Distance (position on the wheel, zero at reward point) from cumulative distance data, downsample to imaging sampling rate
% Calculate when rewards start, 
DistSubtract                       = Distance(sData.behavior.details.rewardStartIndices); % cumulative distances at lap starts (to be subtracted from cumDist data)
sData.behavior.stats.LapLengthCm   = mean(diff(DistSubtract), 'omitnan'); % mean length of laps, should be equal to Circumference
sData.behavior.stats.SDLapLengthCm = std(diff(DistSubtract), 'omitnan');

% inform the user
msgbox(sprintf('Full trial number is: %d. Mean LapLength: %g, Stdev LapLength: %g. Number of frames catched by LV frame signal: %g. Imaging sampling frequency is %g Hz.',...
    sData.behavior.wheelLap, sData.behavior.stats.LapLengthCm, sData.behavior.stats.SDLapLengthCm, sData.behavior.details.nFrames, sData.behavior.meta.imagingSamplingRate));

% Calculate position on the wheel, zero in each lap is photodiode signal.
sData.behavior.wheelPos = NaN(numel(sData.daqdataB.distanceCm),1);
LapCounter = 1;
for i = sData.behavior.details.rewardStartIndices(1):1:numel(sData.daqdataB.distanceCm)
    if i > sData.behavior.details.rewardStartIndices(LapCounter+1)
        LapCounter = LapCounter + 1;
        if LapCounter == numel(sData.behavior.details.rewardStartIndices)
            break
        end
    end
    sData.behavior.wheelPos(i) = Distance(i) - DistSubtract(LapCounter);
end
% Downsample distance data to imaging frequency, consider only full laps with frame signal
sData.behavior.wheelPosDs = sData.behavior.wheelPos(FrameStartIndex);
% make it monotonic incr - in some cases,especially close to the reward, animals grab the wheel and move it back and forward, while their body is stationary. 
% Therefore in the position signal sometimes there are false backward movements. 
% I wanted to compensate for this and make the position signal monotonically increasing (if position moves backward, I gave the previous position value)
sData.behavior.wheelPosDsMonIncr = sData.behavior.wheelPosDs; % make it monotonically increasing 
for i = 1:1:numel(sData.behavior.wheelPosDsMonIncr)-1
   if sData.behavior.wheelPosDsMonIncr(i)>150 && sData.behavior.wheelPosDsMonIncr(i+1)<5 % neglect lap start region
      continue
   elseif  sData.behavior.wheelPosDsMonIncr(i)> sData.behavior.wheelPosDsMonIncr(i+1) 
       sData.behavior.wheelPosDsMonIncr(i+1) = sData.behavior.wheelPosDsMonIncr(i);
   elseif  sData.behavior.wheelPosDsMonIncr(i) < 5   && sData.behavior.wheelPosDsMonIncr(i+1) >150 
       sData.behavior.wheelPosDsMonIncr(i+1) = sData.behavior.wheelPosDsMonIncr(i);
   end
end

% if the maximum backtracking is larger then 5 cm, discard the session
sData.behavior.wheelPosBackTrack    = sData.behavior.wheelPosDsMonIncr - sData.behavior.wheelPosDs;
sData.behavior.details.MaxBackTRack = max(sData.behavior.wheelPosBackTrack);
msgbox(sprintf('The maximum backtracing was (cm): %g ',sData.behavior.details.MaxBackTRack));


% Velocity calculation
Circumfer                                 = pi*WheelDiameterCm; % Wheel circumference
sData.behavior.stats.TheoreticLapLengthCm = Circumfer;
%MovMean = 7 ; % moving window for smoothing
sData.behavior.runSpeedPre = movmean(diff(sData.daqdataB.distanceCm),(sData.daqdataB.meta.fs/5)); % smoothing for 3000/5=600 samples, 200 ms;
sData.behavior.runSpeed    = sData.behavior.runSpeedPre / (1/sData.daqdataB.meta.fs);

test = sData.behavior.runSpeed(FrameStartIndex);

VelSubtr   = diff(sData.behavior.wheelPosDsMonIncr); % calculate distance elapsed between frames
% VelSubtr   = [NaN; diff(sData.behavior.wheelPosDsMonIncr)]; % calculate distance elapsed between frames

VelSubtr(1) = sData.behavior.wheelPosDsMonIncr(2); % first should have been calculated manually
% at lap starts have to adjust do not have negative valuse (157 -> 0)
for i = 1:length(VelSubtr)
   if VelSubtr(i) < 0 && Circumfer > sData.behavior.wheelPosDsMonIncr(i)
       VelSubtr(i) = Circumfer - sData.behavior.wheelPosDsMonIncr(i) + sData.behavior.wheelPosDsMonIncr(i+1);
       % VelSubtr(i) = Circumfer - sData.behavior.wheelPosDsMonIncr(i-1) + sData.behavior.wheelPosDsMonIncr(i);

   elseif VelSubtr(i) < 0 % sometimes abs distance is calculated bigger than Circumference
       VelSubtr(i) = sData.behavior.wheelPosDsMonIncr(i+1);
   elseif VelSubtr(i) == 0
       VelSubtr(i) = 0.00001; % zeros is not good for later processing. 
   end
end
% the last value cannot be computed (now it is NAN), set it to the previous value, and also set a few more data to be able to calculate movmean
i = sum(~isnan(VelSubtr),1);  % search the last value which is non Nan
for j = 1:ceil(MovMean/2)
    VelSubtr(i+j) = VelSubtr(i);
end

VelPre                         = VelSubtr./(1/(sData.behavior.meta.imagingSamplingRate)); 
VelSmooth                      = movmean(VelPre,MovMean); % calculate movmean
NaNReplace                     = find(isnan(VelSmooth)); % in the beginning and end for movemean data, I have to substitute with original data
VelSmooth(NaNReplace,1)        = VelPre(NaNReplace,1);
VelSmooth(VelSmooth<0)         = 0.00001; % set negative values to almost zero
sData.behavior.runSpeedDs      = VelSmooth;
sData.behavior.stats.VelMaxCmS = max(sData.behavior.runSpeedDs);

% generate downsampled lick data and downsampled water-given data
sData.behavior.lickDs        = NaN(sData.behavior.details.nFrames-1,1); 
sData.behavior.waterRewardDs = NaN(sData.behavior.details.nFrames-1,1);
for i = 1:1:sData.behavior.details.nFrames-1
    TempLick  = NaN(max(diff(FrameStartIndex))+1,1); % temporary array for lick signal during a frame (scan)
    TempWater = NaN(max(diff(FrameStartIndex))+1,1); % temporary array for water signal during a frame (scan)
    TempLick  = sData.daqdataB.lickSignal(FrameStartIndex(i):FrameStartIndex(i+1));
    TempWater = sData.daqdataB.waterValve(FrameStartIndex(i):FrameStartIndex(i+1));
    if sum(TempLick)> max(diff(FrameStartIndex))/2 % if there is lick within this frame, put 1 into DS lick array (if lick signal is one more than half time during frame scanning (one frame is 320 samples, one lick is usually 500 samples, sampling is 10 kHz)
        sData.behavior.lickDs(i) = 1;
    else
        sData.behavior.lickDs(i) = 0;
    end
    if sum(TempWater) > max(diff(FrameStartIndex))/2 % if there is water within this frame, put 1 into DS lick array (if lick signal is one more than half time during frame scanning (one frame is 320 samples, one lick is usually 500 samples, sampling is 10 kHz)
        sData.behavior.waterRewardDs(i) = 1;
    else
        sData.behavior.waterRewardDs(i) = 0;
    end
end
sData.behavior.stats.LicksDSNu     = sum(sData.behavior.lickDs);
sData.behavior.stats.WaterRewardNu = sum(sData.behavior.waterRewardDs);

% Check licking duarion
sData.behavior.stats.OneFrameDurInSample = 98;
LickStartTemp                            = find(diff(sData.daqdataB.lickSignal)==1); % search lick signal start
LickEndTemp                              = find(diff(sData.daqdataB.lickSignal)==-1); % lick signal ends
if LickEndTemp(1)> LickStartTemp(1)   % do not start during a lick
    LickEnd = LickEndTemp; % lick signal ends
else
    LickEnd = LickEndTemp(2:end);
end
if LickEnd(end)> LickStartTemp(end)   % do not end during a lick
    LickStart = LickStartTemp; % lick signal ends
else
    LickStart = LickStartTemp(1:end-1);
end

LickLengthSmpl                              = LickEnd - LickStart; % duration of lick signal
LickLengthSmpl                              = LickLengthSmpl(LickLengthSmpl> sData.behavior.stats.OneFrameDurInSample/5); % many single ones, the bottom of distribution is about 100
sData.behavior.stats.medianLickLengthSample = median(LickLengthSmpl);
sData.behavior.stats.SDLickLengthSample     = std(LickLengthSmpl);

% Time axis
sData.behavior.timeInSec                  = ((0:length(sData.behavior.wheelPosDsMonIncr)-1)/sData.behavior.meta.imagingSamplingRate)';
sData.behavior.stats.RecordingDurationSec = max(sData.behavior.timeInSec);
sData.behavior.stats.RecordingDurationMin = max(sData.behavior.timeInSec)/60;

%% BINNING
sData                          = SA_binning(sData,nBins);
sData.behavior.wheelLapImaging = sum(sData.behavior.binning.enterIntoBinIndex(:,1)>0)-1;

%% MAKING PLOTS
% save behavior plots into BEHAVIOR subfolder
% mkdir(sData.sessionInfo.savePath,'Behavior');
savePath = sData.sessionInfo.savePath;

% plot binned velocity heatplot 
[sData.behavior.binning.veloBinnedExtended, MeanVeloBinExtended, ...
    sData.behavior.WheelLapImagingExtended]= SA_plotHeatBinVelo(sData,round(sData.behavior.stats.VelMaxCmS,-1)); % input= sData.behavior.runSpeedDs, Ymax for velo
FileName = strcat('VeloHeatBin-',sData.sessionInfo.sessionID);

savefig(fullfile(savePath,FileName));
saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

% plot binned licks heatplot
[sData.behavior.binning.lickBinnedExtended, MeanLickBinExtended, ...
    sData.behavior.binning.lickPerCmBinnedExtended,MeanLickCmExtended] = SA_plotHeatBinLicks(sData); 
FileName = strcat('LickHeatBin-',sData.sessionInfo.sessionID);

savefig(fullfile(savePath,FileName));
saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

% plot WaterGiven plot / hit-miss visible
[sData.behavior.binning.waterGivenBinnedExtended, sData.behavior.details.hitTrials, ...
    sData.behavior.details.hitRate] = SA_plotHeatBinWaterGiven(sData,SizeOfRewardZoneCm);
FileName = strcat('WaterGivenHeatBin-',sData.sessionInfo.sessionID);

savefig(fullfile(savePath,FileName));
saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

% write data into sData
if exist('FakeFrameSignal2')
    sData.behavior.meta.fakeframesignal = FakeFrameSignal2;
end

sData.behavior.binning.veloBinned          = sData.behavior.binning.veloBinnedExtended(1:sData.behavior.WheelLapImagingExtended,1:sData.behavior.meta.nBins);
sData.behavior.binning.lickBinned          = sData.behavior.binning.lickBinnedExtended(1:sData.behavior.WheelLapImagingExtended,1:sData.behavior.meta.nBins);
sData.behavior.binning.lickPerCmBinned     = sData.behavior.binning.lickPerCmBinnedExtended(1:sData.behavior.WheelLapImagingExtended,1:sData.behavior.meta.nBins);
sData.behavior.binning.waterGivenBinned    = sData.behavior.binning.waterGivenBinnedExtended(1:sData.behavior.WheelLapImagingExtended,1:sData.behavior.meta.nBins);
sData.behavior.binning.meanVeloBinned      = MeanVeloBinExtended(1:sData.behavior.meta.nBins);
sData.behavior.binning.meanLickBinned      = MeanLickBinExtended(1:sData.behavior.meta.nBins);
sData.behavior.binning.meanLickPerCmBinned = MeanLickCmExtended(1:sData.behavior.meta.nBins);

%% PLOT average lick / trials
% Plot settings
Xstep    = sData.behavior.meta.binSize;
XaxisEnd = nBins * sData.behavior.meta.binSize;

%% PLOT average of licks throughout the session
figure('Color','white'); 
Xaxis = Xstep/2:Xstep:XaxisEnd-Xstep/2;
plot(Xaxis,MeanLickCmExtended(1:sData.behavior.meta.nBins)); hold on
xlabel('Position on wheel (cm)');
ylabel('Licks/cm');
xlim([0 XaxisEnd])
ylim([0 2])
FileName = strcat('LickSum-',sData.sessionInfo.sessionID);
savefig(fullfile(savePath,FileName));
saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));  

%% PLOT average of speed throughout the session
figure('Color','white'); 
Xaxis = Xstep/2:Xstep:XaxisEnd-Xstep/2;
plot(Xaxis,MeanVeloBinExtended(1:nBins)); hold on
xlabel('Position on wheel (cm)');
ylabel('Velocity (cm/s)');
xlim([0 XaxisEnd])
ylim([0 70])
FileName = strcat('VeloSum-',sData.sessionInfo.sessionID);
savefig(fullfile(savePath,FileName));
saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));


%%% calculate behavioral performnace
if sData.behavior.wheelLap > 20
   sData = behavPerfomance(sData);
   sData = behavPerfomance2(sData);
end

% Save file to same path where LV files can be found 
% save(fullfile(sData.sessionInfo.savePath(1:54),...
%     strcat(sData.sessionInfo.sessionID(1:14), sData.sessionInfo.sessionID(25:28),'.mat')),'sData');

save(fullfile(sData.sessionInfo.savePath(1:54),...
    strcat(sData.sessionInfo.sessionID,'.mat')),'sData');
end