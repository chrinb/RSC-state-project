

clear

% Select and load sData file
[sessionID,filePath,~] = uigetfile('*.mat','Select sData File','C:\Users\Mate Neubrandt\Documents\RECORDINGS','MultiSelect','off' );
load(fullfile(filePath,sessionID));
sessionID = strsplit(sessionID,'.mat');
sessionID = sessionID{1};

%sData = vr.downsampleBehavior(sData);
%% Prepare data and create trial matrices for all ROIs

%samplingRate = sData.daqdata.meta.samplingRate;
binNumber = sData.behavior.trialMatrices.meta.binNumber;
binSize = sData.behavior.trialMatrices.meta.binSize;
homeBoxLength = abs(sData.stats.sessionAvs(1).plotXAxis(1) - binSize);
% rewardZone = sData.behavior.rewardZone;
corridorLength = binNumber*binSize;
viewDistance = sData.behavior.viewDistance;

nHomeBoxBins = homeBoxLength/binSize;
% binnedPosition = discretize(sData.daqdata.unityPosition,-homeBoxLength:binSize:-homeBoxLength+corridorLength); 
corridorPosition = sData.behavior.signals(1).corridorPosition;
binnedPosition = discretize(corridorPosition,floor(min(corridorPosition))+binSize:binSize:ceil(max(corridorPosition))); 
  
trialStartIndexes = find(diff(sData.daqdata.unityPosition)< -10)+1;
allTrials = numel(trialStartIndexes)-1;
%velocity = sData.behavior.signals.velocity;


pupilDiameter = sData.behavior.pupildata.diameter;
pupilCenter = sData.behavior.pupildata.center;

% pupilDiameter = RadiusCorr*2; pupilCenter = CenterRotatedCorr;


%%

%Xaxis = sData.daqdata.unityPosition(find(diff(sData.daqdata.frameIndex)==1));
frameIndexes = find(diff(sData.daqdata.pupilCamFrameIndex)==1); % Sampleindexes of image frames in original data 
%frameIndexes = fixFrameIndexes(frameSignal, samplePerFrame);

frameSignalInd(1:numel(frameIndexes)) = zeros; % This modified frame signal contains the index in the frame array ..001..0002.. 003

for i = 1:1:numel(frameIndexes)

    frameSignalInd(frameIndexes(i)) = i;
    
end




%binNumber = binNumber + homeBoxLength/binSize;

%binnedRoisDff = nan(allTrials,binNumber,nROIs);
%binnedRoisDeconv = nan(allTrials,binNumber,nROIs);
%binnedRoisDeconvRate = nan(allTrials,binNumber,nROIs);

binnedPupilDiameter = nan(allTrials,binNumber);
binnedPupilCenterX = nan(allTrials,binNumber);
binnedPupilCenterY = nan(allTrials,binNumber);

binnedPupilFrameInd = nan(allTrials,binNumber);

pupilCamFrameIndexCorr = sData.daqdata.pupilCamFrameIndexCorr;

% Create roi matrices 
for r = 1:1:allTrials
    tempPos = binnedPosition(trialStartIndexes(r):(trialStartIndexes(r+1)-1));
    tempFrameSInd = frameSignalInd(trialStartIndexes(r):(trialStartIndexes(r+1)-1));
    %tempVel = velocity(trialStartIndexes(r):(trialStartIndexes(r+1)-1));
    %tempFrameSInd(find(tempVel < velThreshold)) = 0;
    tempBinPos = tempPos(find(tempFrameSInd > 1));
    tempFrameSIndCorr = pupilCamFrameIndexCorr(trialStartIndexes(r):(trialStartIndexes(r+1)-1));
    
   % for roi = 1:1:nROIs
        tempSignalDiam = pupilDiameter(tempFrameSInd(tempFrameSInd>0));
        tempSignalX = pupilCenter(tempFrameSInd(tempFrameSInd>0),1);
        tempSignalY = pupilCenter(tempFrameSInd(tempFrameSInd>0),2);
        for c = 1:1:binNumber
            binnedPupilDiameter(r,c) = mean(tempSignalDiam(find(tempBinPos == c)));
            binnedPupilCenterX(r,c) = mean(tempSignalX(find(tempBinPos == c)));  
            binnedPupilCenterY(r,c) = mean(tempSignalY(find(tempBinPos == c))); 
            
            binnedPupilFrameInd(r,c) = nanmean(tempFrameSIndCorr(find(tempPos == c))); 
        end
        clear('tempSignalDiam','tempSignalX','tempSignalY','tempFrameSIndCorr');
   % end
    clear('tempPos','tempFrameSInd','tempBinPos');
end


%% Save matrices to sData

sData.behavior.trialMatrices.binnedPupilDiameter = binnedPupilDiameter;
sData.behavior.trialMatrices.binnedPupilCenterX = binnedPupilCenterX;
sData.behavior.trialMatrices.binnedPupilCenterY = binnedPupilCenterY;

% save
save(fullfile(filePath,sessionID),'sData');

%%

Xax = sData.stats.sessionAvs(1).plotXAxis;

%{
figure;
hold on
plot(Xax,nanmean(binnedPupilCenterX)-nanmean(nanmean(binnedPupilCenterX)))
plot(Xax,nanmean(binnedPupilCenterY)-nanmean(nanmean(binnedPupilCenterY)))
title(sData.sessionInfo.sessionID(1:5))
%plot(Xax,nanmean(binnedPupilDiameter))
%}


figure('Color','white','Position',[0 0 300 600])

dataRange = max(nanmean(binnedPupilDiameter(:,:))) - min(nanmean(binnedPupilDiameter(:,:)));
xlims = [0 max(Xax)];
ylims = [floor((min(nanmean(binnedPupilDiameter(:,:)))-(0.5*dataRange))/5)*5, ceil((max(nanmean(binnedPupilDiameter(:,:)))+(2.5*dataRange))/5)*5];



subplot(2,1,1)
imagesc(Xax,1:size(binnedPupilDiameter,1),binnedPupilDiameter)
title([sData.sessionInfo.sessionID(1:17) ' - Pupil diameter'])
%caxis([20 50])
xlim(xlims)
ylabel('Trials')
xlabel('Track position (cm)')

subplot(2,1,2)
hold on
for t = 1:1:length(sData.trials.contextsMeta)
    trials = sData.trials.contextsMeta(t).trials;
    plot(Xax,nanmean(binnedPupilDiameter(trials,:)))
end

legend({sData.trials.contextsMeta.name}); %{sData.trials.contextsMeta(1).name,sData.trials.contextsMeta(2).name,sData.trials.contextsMeta(3).name,sData.trials.contextsMeta(4).name})
xlim(xlims)
ylim(ylims)
ylabel('Pupil diameter (px)')
xlabel('Track position (cm)')


%% Save Fig

if ~exist([filePath 'pupilFigs'], 'dir')
   mkdir([filePath 'pupilFigs']);
end

saveas(gcf,strcat(fullfile(filePath,'pupilFigs',[sData.sessionInfo.sessionID(1:17) '_pupilDiameter']),'.png'));
%close(gcf)





%%

figure('Color','white','Position',[0 0 350 350])
hold on
for i = 2:2:4
rectangle('Position',[sData.trials.contextsMeta(i).blockStart, ylims(1), sData.trials.contextsMeta(i).nTrials, ylims(2)], 'FaceColor', [0.95 0.95 0.95], 'EdgeColor', 'none')
end
plot(1:allTrials,nanmean(binnedPupilDiameter'),'Linewidth',1.5)
ylim(ylims)

xlabel('Trials')
ylabel('Pupil diameter (px)')

saveas(gcf,strcat(fullfile(filePath,'pupilFigs',[sData.sessionInfo.sessionID(1:17) '_pupilDiameterOverTrials']),'.png'));
%close(gcf)


%%

% imagesc(Xax,1:size(binnedPupilDiameter,1),binnedPupilCenterX)

%{
vels = binVel(:,71);
[~, order] = sort(vels);


figure
imagesc(Xax,1:size(binnedPupilDiameter,1),binVel(order,:))

figure
imagesc(Xax,1:size(binnedPupilDiameter,1),binnedPupilDiameter(order,:))


figure
plot(sData.behavior.pupildata.lostFramesArray)

figure
plot(sData.daqdata.pupilCamFrameIndex)

%}


%% Settings

nTrialTypes = size(sData.trials.trialTypesMeta,2);
if nTrialTypes >1
for i = 1:1:nTrialTypes-1
stimulusType = i+1; %Use it to derermine which stimulus prorocol to plot in addition to the control

switch stimulusType
    case 2
        plotLegend = {'no stim'; ''; 'Full trial'; ''};
    case 3
            plotLegend = {'no stim'; ''; 'Home-box'; ''};
    case 4
            plotLegend = {'no stim'; ''; 'On-the-path'; ''};
end

smoothSpan = 5;      
% plotLegend = {'no stim'; ''; 'Full trial'; ''};
rewardZone = sData.behavior.rewardZone;
rewardZoneWidth = 10;
corridorLength = rewardZone + rewardZoneWidth;
corridorStart = -140;

% mymap = lines(nTrialTypes);
mymap = lines;
Xax = sData.behavior.trialMatrices.plotXAxis;
binVel = sData.behavior.trialMatrices.binVel;
lickFreq = sData.behavior.trialMatrices.lickFreqInBin;
optStimMatrix = sData.behavior.trialMatrices.optStimMatrix;
binnedPupilDiameter = sData.behavior.trialMatrices.binnedPupilDiameter;
binnedPupilCenterX = sData.behavior.trialMatrices.binnedPupilCenterX;
binnedPupilCenterY = sData.behavior.trialMatrices.binnedPupilCenterY;

%lickDistrCtrl = sData.trials.trialTypesMeta(1).lickQuartiles(2)*rewardZone/100;
%lickDistrStim = sData.trials.trialTypesMeta(nTrialTypes).lickQuartiles(2)*rewardZone/100;
firstLickCtrl = sData.trials.trialTypesMeta(1).firstLicksCm(2);
firstLickStim = sData.trials.trialTypesMeta(stimulusType).firstLicksCm(2);

% for plotting stimulus location
stimCurve = nanmean(optStimMatrix(sData.trials.trialTypesMeta(stimulusType).trials,:));
stimMax = max(stimCurve);
stimMin = min(stimCurve); 
stimLocation = stimCurve;
stimLocation(stimCurve < (stimMax - stimMin)/3) = 0;
stimLocation(stimCurve >= (stimMax - stimMin)/3) = 1;
stimStart = Xax(min(find(stimLocation == 1)));
stimEnd = Xax(max(find(stimLocation == 1)));






figure('Color','white','Position',[0 0 850 850]); %pos of figure [left bottom width height]

suptitle([sData.sessionInfo.sessionID newline])

subplot(2,2,1); % Velocity
hold on

rectangle('Position',[corridorStart+1,0.5,-corridorStart,80],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
rectangle('Position',[rewardZone,0.5,10,80],'FaceColor',[0.9 0.97 0.92],'EdgeColor','none'); % RZ
rectangle('Position',[stimStart,76,(stimEnd-stimStart),80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % Stimulus
%rectangle('Position',[(rewardZone-viewDistance),0.5,0.5,80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

for t = 1:stimulusType-1:stimulusType
    vr.trialType;
    color = mymap(t,:);
    trialMatrix = binVel(trials,:);
    %plot(Xax,smoothdata(sessionAvs(t).avBinVel,'gaussian',smoothSpan),'Color',map(t,:));
    plotQuartiles(Xax,trialMatrix,color,smoothSpan)
end


axis([corridorStart corridorLength 0 80]);
xlabel('Position in unity (virtual cm)'); %  num2str(binSize) ' cm/bin)'
ax = gca;
ax.TickDir = 'out';
ylabel('Velocity (cm/s)');



subplot(2,2,2); % Lick freq
hold on

rectangle('Position',[corridorStart+1,0.05,-corridorStart,80],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
rectangle('Position',[rewardZone,0.05,10,80],'FaceColor',[0.9 0.97 0.92],'EdgeColor','none'); % RZ
rectangle('Position',[stimStart,7.6,(stimEnd-stimStart),80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % Stimulus
%rectangle('Position',[(rewardZone-viewDistance),0.5,0.5,80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

for t = 1:stimulusType-1:stimulusType
    vr.trialType;
    color = mymap(t,:);
    trialMatrix = lickFreq(trials,:);
    
    plotQuartiles(Xax,trialMatrix,color,smoothSpan)
    
end


%rectangle('Position',[lickDistrCtrl-1, 4.0025, 2, 10],'FaceColor',mymap(1,:),'EdgeColor','none'); % control lick distr
%rectangle('Position',[lickDistrStim-1, 4.0025, 2, 10],'FaceColor',mymap(nTrialTypes,:),'EdgeColor','none'); % stim lick distr

rectangle('Position',[firstLickCtrl-0.5, 0.0025, 1, 4],'FaceColor',mymap(1,:),'EdgeColor','none'); % control lick distr
rectangle('Position',[firstLickStim-0.5, 0.0025, 1, 4],'FaceColor',mymap(stimulusType,:),'EdgeColor','none'); % stim lick distr


axis([corridorStart corridorLength 0 8]);
xlabel('Position in unity (virtual cm)'); %  num2str(binSize) ' cm/bin)'
ax = gca;
ax.TickDir = 'out';

%ylabel('Lick before reward (lick/cm)');
ylabel('Lick rate (Hz)');
legend(plotLegend,'Location','northwest')



subplot(2,2,3); % Pupil position
hold on

rectangle('Position',[corridorStart+1,-6.95,-corridorStart,80],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
rectangle('Position',[rewardZone,-6.95,10,80],'FaceColor',[0.9 0.97 0.92],'EdgeColor','none'); % RZ
rectangle('Position',[stimStart,6.5,(stimEnd-stimStart),80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % Stimulus
%rectangle('Position',[(rewardZone-viewDistance),0.5,0.5,80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

for t = 1:stimulusType-1:stimulusType
    vr.trialType;
    color = mymap(t,:);
    % trials = sData.trials.trialTypesMeta(t).trials;
    trialMatrix = binnedPupilCenterX(trials,:);
    trialMatrix = trialMatrix-nanmean(nanmean(trialMatrix))+2;
    %plot(Xax,smoothdata(sessionAvs(t).avBinVel,'gaussian',smoothSpan),'Color',map(t,:));
    plotQuartiles(Xax,trialMatrix,color,smoothSpan)
    
    text(-120,4,'X axis')
    %{
    color = mymap(t+1,:);
    trialMatrix = binnedPupilCenterY(:,:);
    trialMatrix = (trialMatrix-nanmean(nanmean(trialMatrix)))*-1;
    %plot(Xax,smoothdata(sessionAvs(t).avBinVel,'gaussian',smoothSpan),'Color',map(t,:));
    plotQuartiles(Xax,trialMatrix,color,smoothSpan)
    %}
end
for t = 1:stimulusType-1:stimulusType
    vr.trialType;
    color = mymap(t,:);
    % trials = sData.trials.trialTypesMeta(t).trials;
    trialMatrix = binnedPupilCenterY(trials,:);
    trialMatrix = trialMatrix-nanmean(nanmean(trialMatrix))+2;
    trialMatrix = trialMatrix *-1; %invert axis to plot upward movement up
    %plot(Xax,smoothdata(sessionAvs(t).avBinVel,'gaussian',smoothSpan),'Color',map(t,:));
    plotQuartiles(Xax,trialMatrix,color,smoothSpan)
    
    text(-120,-4,'Y axis')
    %{
    color = mymap(t+1,:);
    trialMatrix = binnedPupilCenterY(:,:);
    trialMatrix = (trialMatrix-nanmean(nanmean(trialMatrix)))*-1;
    %plot(Xax,smoothdata(sessionAvs(t).avBinVel,'gaussian',smoothSpan),'Color',map(t,:));
    plotQuartiles(Xax,trialMatrix,color,smoothSpan)
    %}
end

axis([corridorStart corridorLength -7 7]);
xlabel('Position in unity (virtual cm)'); %  num2str(binSize) ' cm/bin)'
ax = gca;
ax.TickDir = 'out';
ylabel('Pupile position (px)');




subplot(2,2,4); % Pupil diameter
hold on

rectangle('Position',[corridorStart+1,0.5,-corridorStart,80],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
rectangle('Position',[rewardZone,0.5,10,80],'FaceColor',[0.9 0.97 0.92],'EdgeColor','none'); % RZ
rectangle('Position',[stimStart,76,(stimEnd-stimStart),80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % Stimulus
%rectangle('Position',[(rewardZone-viewDistance),0.5,0.5,80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

for t = 1:stimulusType-1:stimulusType
    vr.trialType;
    color = mymap(t,:);
    trialMatrix = binnedPupilDiameter(trials,:);
    %plot(Xax,smoothdata(sessionAvs(t).avBinVel,'gaussian',smoothSpan),'Color',map(t,:));
    plotQuartiles(Xax,trialMatrix,color,smoothSpan)
end


axis([corridorStart corridorLength 0 80]);
xlabel('Position in unity (virtual cm)'); %  num2str(binSize) ' cm/bin)'
ax = gca;
ax.TickDir = 'out';
ylabel('Pupil diameter (px)');







%% Save Fig


if ~exist([filePath 'pupilFigs'], 'dir')
   mkdir([filePath 'pupilFigs']);
end

saveas(gcf,strcat(fullfile(filePath,'pupilFigs',[sData.sessionInfo.sessionID '_pupilFigs_' num2str(stimulusType)]),'.png'));
%close(gcf)


end
end

% msgbox(['Lost frames: ' num2str(sData.behavior.pupildata.nLostFrames)])





%% Plot not opto data



if nTrialTypes == 1
% settings

smoothSpan = 5;      
% plotLegend = {'no stim'; ''; 'Full trial'; ''};
rewardZone = sData.behavior.rewardZone;
rewardZoneWidth = 10;
%corridorLength = rewardZone + rewardZoneWidth;
%corridorStart = -140;

% mymap = lines(nTrialTypes);
mymap = lines;
Xax = sData.behavior.trialMatrices.plotXAxis;
binVel = sData.behavior.trialMatrices.binVel;
lickFreq = sData.behavior.trialMatrices.lickFreqInBin;
optStimMatrix = sData.behavior.trialMatrices.optStimMatrix;
binnedPupilDiameter = sData.behavior.trialMatrices.binnedPupilDiameter;
binnedPupilCenterX = sData.behavior.trialMatrices.binnedPupilCenterX;
binnedPupilCenterY = sData.behavior.trialMatrices.binnedPupilCenterY;

%lickDistrCtrl = sData.trials.trialTypesMeta(1).lickQuartiles(2)*rewardZone/100;
%lickDistrStim = sData.trials.trialTypesMeta(nTrialTypes).lickQuartiles(2)*rewardZone/100;
firstLickCtrl = sData.trials.trialTypesMeta(1).firstLicksCm(2);

    
    
    
figure('Color','white','Position',[0 0 850 850]); %pos of figure [left bottom width height]

suptitle([sData.sessionInfo.sessionID newline])

subplot(2,2,1); % Velocity
hold on

rectangle('Position',[min(Xax)+1,0.5,-min(Xax),80],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
rectangle('Position',[max(Xax),0.5,10,80],'FaceColor',[0.9 0.97 0.92],'EdgeColor','none'); % RZ
% rectangle('Position',[stimStart,76,(stimEnd-stimStart),80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % Stimulus
%rectangle('Position',[(rewardZone-viewDistance),0.5,0.5,80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

t=1;
%for t = 1:stimulusType-1:stimulusType
    %vr.trialType;
    color = mymap(t,:);
    trialMatrix = binVel(sData.trials.trialTypesMeta(t).trials,:);
    %plot(Xax,smoothdata(sessionAvs(t).avBinVel,'gaussian',smoothSpan),'Color',map(t,:));
    plotQuartiles(Xax,trialMatrix,color,smoothSpan)
%end


axis([min(Xax) max(Xax) 0 80]);
xlabel('Position in unity (virtual cm)'); %  num2str(binSize) ' cm/bin)'
ax = gca;
ax.TickDir = 'out';
ylabel('Velocity (cm/s)');



subplot(2,2,2); % Lick freq
hold on

rectangle('Position',[min(Xax)+1,0.05,-min(Xax),80],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
rectangle('Position',[max(Xax),0.05,10,80],'FaceColor',[0.9 0.97 0.92],'EdgeColor','none'); % RZ
% rectangle('Position',[stimStart,7.6,(stimEnd-stimStart),80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % Stimulus
%rectangle('Position',[(rewardZone-viewDistance),0.5,0.5,80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

t = 1;
%for t = 1:nTrialTypes-1:nTrialTypes
    
    color = mymap(t,:);
    trialMatrix = lickFreq(sData.trials.trialTypesMeta(t).trials,:);
    
    plotQuartiles(Xax,trialMatrix,color,smoothSpan)
    
%end

%{
rectangle('Position',[lickDistrCtrl-1, 4.0025, 2, 10],'FaceColor',mymap(1,:),'EdgeColor','none'); % control lick distr
rectangle('Position',[lickDistrStim-1, 4.0025, 2, 10],'FaceColor',mymap(nTrialTypes,:),'EdgeColor','none'); % stim lick distr

rectangle('Position',[firstLickCtrl-0.5, 0.0025, 1, 4],'FaceColor',mymap(1,:),'EdgeColor','none'); % control lick distr
rectangle('Position',[firstLickStim-0.5, 0.0025, 1, 4],'FaceColor',mymap(nTrialTypes,:),'EdgeColor','none'); % stim lick distr
%}

axis([min(Xax) max(Xax) 0 8]); 
xlabel('Position in unity (virtual cm)'); %  num2str(binSize) ' cm/bin)'
ax = gca;
ax.TickDir = 'out';

%ylabel('Lick before reward (lick/cm)');
ylabel('Lick rate (Hz)');
% legend(plotLegend,'Location','northwest')



subplot(2,2,3); % Pupil position
hold on

rectangle('Position',[min(Xax)+1,-6.95,-min(Xax),80],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
rectangle('Position',[max(Xax),-6.95,10,80],'FaceColor',[0.9 0.97 0.92],'EdgeColor','none'); % RZ
% rectangle('Position',[stimStart,66.5,(stimEnd-stimStart),80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % Stimulus
%rectangle('Position',[(rewardZone-viewDistance),0.5,0.5,80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

t=1;
%for t = 1:stimulusType-1:stimulusType
    %vr.trialType;
    color = mymap(t,:);
    trials = sData.trials.trialTypesMeta(t).trials;
    trialMatrix = binnedPupilCenterX(:,:);
    trialMatrix = (trialMatrix-nanmean(nanmean(trialMatrix)))+2;
    %plot(Xax,smoothdata(sessionAvs(t).avBinVel,'gaussian',smoothSpan),'Color',map(t,:));
    plotQuartiles(Xax,trialMatrix,color,smoothSpan)
        
%    color = mymap(t+1,:);
    trialMatrix = binnedPupilCenterY(:,:);
    trialMatrix = ((trialMatrix-nanmean(nanmean(trialMatrix)))+2)*-1;
    %plot(Xax,smoothdata(sessionAvs(t).avBinVel,'gaussian',smoothSpan),'Color',map(t,:));
    plotQuartiles(Xax,trialMatrix,color,smoothSpan)
    
%end

text(-120,4,'X axis')
text(-120,-4,'Y axis')

axis([min(Xax) max(Xax) -7 7]);
xlabel('Position in unity (virtual cm)'); %  num2str(binSize) ' cm/bin)'
ax = gca;
ax.TickDir = 'out';
ylabel('Pupile position (px)');




subplot(2,2,4); % Pupil diameter
hold on

rectangle('Position',[min(Xax)+1,0.5,-min(Xax),80],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
rectangle('Position',[max(Xax),0.5,10,80],'FaceColor',[0.9 0.97 0.92],'EdgeColor','none'); % RZ
% rectangle('Position',[stimStart,66.5,(stimEnd-stimStart),80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % Stimulus
%rectangle('Position',[(rewardZone-viewDistance),0.5,0.5,80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

t=1;
%for t = 1:stimulusType-1:stimulusType
    %vr.trialType;
    color = mymap(t,:);
    trialMatrix = binnedPupilDiameter(:,:);
    %plot(Xax,smoothdata(sessionAvs(t).avBinVel,'gaussian',smoothSpan),'Color',map(t,:));
    plotQuartiles(Xax,trialMatrix,color,smoothSpan)
%end


axis([min(Xax) max(Xax) 0 80]);
xlabel('Position in unity (virtual cm)'); %  num2str(binSize) ' cm/bin)'
ax = gca;
ax.TickDir = 'out';
ylabel('Pupil diameter (px)');




%% Save Fig


if ~exist([filePath 'pupilFigs'], 'dir')
   mkdir([filePath 'pupilFigs']);
end

saveas(gcf,strcat(fullfile(filePath,'pupilFigs',[sData.sessionInfo.sessionID '_pupilFigs']),'.png'));
%close(gcf)


end

msgbox(['Lost frames: ' num2str(sData.behavior.pupildata.nLostFrames)])





%% Plot pupil diameter from different sessions %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{

clear
sDataFiles = vr.loadData;

mymap = lines(2);

figure
hold on

x = sDataFiles{1, 1}.behavior.trialMatrices.plotXAxis;
trialMatrix = sDataFiles{1, 1}.behavior.trialMatrices.binnedPupilDiameter; 
color = mymap(1,:);
plotQuartiles(x,trialMatrix,color,5)

trialMatrix = sDataFiles{1, 2}.behavior.trialMatrices.binnedPupilDiameter; 
color = mymap(2,:);
plotQuartiles(x,trialMatrix,color,5)

title('Pupil diameter: masking light vs no-masling light')



trialMatrix = sDataFiles{1, 1}.behavior.trialMatrices.binnedPupilDiameter; 
trialMatrix = sDataFiles{1, 1}.behavior.trialMatrices.binnedPupilCenterX; 

figure
hold on


x = sDataFiles{1, 1}.behavior.trialMatrices.plotXAxis;

color = mymap(1,:);
trials = sDataFiles{1, 1}.trials.trialTypesMeta(1).trials;
plotQuartiles(x,trialMatrix(trials,:),color,5)


color = mymap(2,:);
trials = sDataFiles{1, 1}.trials.trialTypesMeta(4).trials;
plotQuartiles(x,trialMatrix(trials,:),color,5)

title('Pupil diameter: masking light vs no-masling light')






figure
hold on

rectangle('Position',[corridorStart+1,0.5,-corridorStart,80],'FaceColor',[0.95 0.95 0.95],'EdgeColor','none'); % home-box
rectangle('Position',[rewardZone,0.5,10,80],'FaceColor',[0.9 0.97 0.92],'EdgeColor','none'); % RZ
% rectangle('Position',[stimStart,66.5,(stimEnd-stimStart),80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % Stimulus
%rectangle('Position',[(rewardZone-viewDistance),0.5,0.5,80],'FaceColor',[1 0.5 0.5],'EdgeColor','none'); % line indicating the view distance from reward zone

t=1;
%for t = 1:nTrialTypes-1:nTrialTypes
    %vr.trialType;
    color = mymap(t,:);
    trialMatrix = binnedPupilDiameter(:,:);
    %plot(Xax,smoothdata(sessionAvs(t).avBinVel,'gaussian',smoothSpan),'Color',map(t,:));
    plotQuartiles(Xax,trialMatrix,color,smoothSpan)
%end


axis([corridorStart corridorLength 0 80]);
xlabel('Position in unity (virtual cm)'); %  num2str(binSize) ' cm/bin)'
ax = gca;
ax.TickDir = 'out';
ylabel('Pupil diameter (px)');





%}

