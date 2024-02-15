%% correct pupil data for missing frames 
%inputs:
%sessionID
%filePath
%sData
%output:
%sData

function sData = processPupilData(sData,sessionID,filePath)

% read saved frame times
timeArray = vr.importPupilFrameTimes(fullfile(filePath,['-' sessionID '_frametimes.txt']));
diffArray = milliseconds(diff(timeArray));


% calculate frame times from frame index signal
frameIndexes = find(diff(sData.daqdata.pupilCamFrameIndex) == 1)';
frameIndexIntervalls = diff(frameIndexes)/(sData.daqdata.meta.fs/1000); % frame sample interwall in ms
frameIndexIntervalls = frameIndexIntervalls(2:numel(frameIndexIntervalls));


% Calculate the cummulative difference between the two signal
cummDiff(min(numel(frameIndexIntervalls)-1, numel(diffArray)),1) = nan; 
for i = 1:1:min(numel(frameIndexIntervalls)-1, numel(diffArray))

    cummDiff(i,1) =  sum(frameIndexIntervalls(1:i)) - sum(diffArray(1:i));

end

% baseline substraction
for i=1:1:floor(numel(cummDiff)/1000) %fit linear to all 1000 sample fragment then take the median steapness
   p(i,:) = polyfit(1:1000, cummDiff(1+(i-1)*1000:i*1000)', 1);   
end

m = median(p(:,1));
fitY = polyval([m,1], 1:1:numel(cummDiff))';


fs = 1000/mean(frameIndexIntervalls);

lostFramesArray = round((cummDiff-fitY)/(1000/fs) * -1);
[~,ind] = min(lostFramesArray);
lostFramesArray(2:numel(cummDiff)+1) = lostFramesArray - min(lostFramesArray);
lostFramesArray(1:ind) = 0;
% lostFramesArray(2:numel(cummDiff)+1) = lostFramesArray - lostFramesArray(1); 
nLostFrames = lostFramesArray(end);



%{
figure
hold on
%plot(cummDiff)
plot(lostFramesArray)
plot(fitY)
%}

%% correct pupil frame signal

pupilCamFrameIndexCorr = sData.daqdata.pupilCamFrameIndex;

if nLostFrames > 0
    
    lostFramesLocations = find(diff(lostFramesArray)>0)' + 1 ;
    lostFramesCounts = lostFramesArray(lostFramesLocations)';
    
    lostFramesLocations = lostFramesLocations(diff([0 lostFramesCounts])>0);
    lostFramesCounts = lostFramesCounts(diff([0 lostFramesCounts])>0);
    lostFramesLocations = ((lostFramesLocations + lostFramesCounts)-lostFramesCounts(1))';
    lostFramesCounts = diff([0 lostFramesCounts])';
    
    deleteIndexes = zeros(numel(lostFramesArray)+ nLostFrames,1);
    
    for i = 1:1:numel(lostFramesLocations)
        
        deleteIndexes(lostFramesLocations(i):lostFramesLocations(i) + lostFramesCounts(i)) = 1;
        
    end
    
    delInd = find(deleteIndexes);
    
    
    for i = 1:1:numel(delInd)
        pupilCamFrameIndexCorr(frameIndexes(delInd(i)+1):frameIndexes(delInd(i)+ 1)) = 0;
    end
    
end

sData.daqdata.pupilCamFrameIndexCorr = pupilCamFrameIndexCorr;
sData.behavior.signals.pupilCamFrameIndexCorr = pupilCamFrameIndexCorr;

%% correct data for lost frames
try 
    load(fullfile(filePath,['-' sessionID '_pupilData']));
catch
    msgbox(['ERROR:' newline ['-' sessionID '_pupilData'] newline 'is not available in the defaule folder!'])
end

diffLostFramesArray = diff(lostFramesArray);
breakingPointIndexes = find(diffLostFramesArray);
framesToShift = lostFramesArray(breakingPointIndexes+1);
% framesToShift = diffLostFramesArray(breakingPointIndexes);
breakingPointIndexes(numel(breakingPointIndexes)+1) = numel(Radius); %add the last data point to the end 

RadiusCorr = nan(size(Radius,1)+nLostFrames,1) ;
CenterCorr = nan(size(Center,1)+nLostFrames,2);
CenterRotatedCorr = nan(size(CenterRotated,1) + nLostFrames,2);

if nLostFrames <= 0
    RadiusCorr = Radius;
    CenterCorr = Center;
    CenterRotatedCorr = CenterRotated;
else
    
    for i = 0:1:numel(framesToShift)
        
        if i == 0 
            j = 1;
            k = breakingPointIndexes(i+1);
            l = j;
            m = k;
        else
            j = breakingPointIndexes(i) + 1;
            k = breakingPointIndexes(i+1);
            l = j + framesToShift(i);
            m = k + framesToShift(i);
        end
        
        RadiusCorr(l:m) = Radius(j:k);
        CenterCorr(l:m,:) = Center(j:k,:);
        CenterRotatedCorr(l:m,:) = CenterRotated(j:k,:);
        
    end
end


%% Filter noise (Method 2, see below)
th1 = 10;
th2 = 5;

[~, score] = pca(CenterCorr);
pcaXz = abs(normalize(score(:,1)));
pcaYz = abs(normalize(score(:,2)));

errors = union(find(pcaXz>th2),find(pcaYz>th2));

if max(pcaXz)>th1 || max(pcaYz)>th1
    RadiusCorr(errors) = NaN;
    CenterCorr(errors,:) = NaN;
    CenterRotatedCorr(errors,:) = NaN;
end   



%{
figure
hold on
plot(lostFramesArray)
plot(RadiusCorr*(nanmean(lostFramesArray)/nanmean(RadiusCorr)))

figure
plot(isnan(Radius))
%}


sData.behavior.signals.pupilCamFrameIndexCorr = pupilCamFrameIndexCorr;

sData.behavior.pupildata.diameter = RadiusCorr*2;
sData.behavior.pupildata.center = CenterRotatedCorr;
sData.behavior.pupildata.centerNotRotated = CenterCorr;

sData.behavior.pupildata.diameterRaw = Radius*2;
sData.behavior.pupildata.centerRaw = CenterRotated;
sData.behavior.pupildata.centerNotRotatedRaw = Center;

sData.behavior.pupildata.lostFramesArray = lostFramesArray;
sData.behavior.pupildata.nLostFrames = nLostFrames;

nLostFrames


end


%{ 

%   CenterRotatedCorr = CenterRotated; RadiusCorr = Radius;

% Method 1

CenterXz = abs(normalize(CenterRotatedCorr(:,1)));
%{
figure
hold on
plot(CenterRotatedCorr(:,1),CenterRotatedCorr(:,2),'.')
plot(CenterRotatedCorr(CenterXz>5,1),CenterRotatedCorr(CenterXz>5,2),'.')
%}


% Method 2

[~, score] = pca(CenterRotatedCorr);
x= score(:,1);
y= score(:,2);

pcaXz = abs(normalize(x));
pcaYz = abs(normalize(y));

% figure; histogram(pcaYz)

%{
figure
hold on
plot(CenterRotatedCorr(:,1),CenterRotatedCorr(:,2),'.')
plot(CenterRotatedCorr(pcaXz>5,1),CenterRotatedCorr(pcaXz>5,2),'.')
plot(CenterRotatedCorr(pcaYz>5,1),CenterRotatedCorr(pcaYz>5,2),'.')
%}


% Method 3

diffRadCorr = diff(RadiusCorr);
threshold = 5; 
outliers = find(abs(diffRadCorr)>threshold)+1;
diffOutliers = diff(outliers);
errors = false(size(RadiusCorr));

for i = 1:1:numel(diffOutliers)
    
    if diffOutliers(i)< 20
        j = outliers(i);
        k = outliers(i+1)-1;
        
        errors(j:k) = true;
    end
    
end

%{
figure
hold on
plot(CenterRotatedCorr(:,1),CenterRotatedCorr(:,2),'.')
plot(CenterRotatedCorr(errors,1),CenterRotatedCorr(errors,2),'.')
%}


%figure; histogram(diffRadCorr)



figure
subplot(1,3,1)
hold on
plot(CenterRotatedCorr(:,1),CenterRotatedCorr(:,2),'.')
plot(CenterRotatedCorr(CenterXz>5,1),CenterRotatedCorr(CenterXz>5,2),'.')


subplot(1,3,2)
hold on
plot(CenterRotatedCorr(:,1),CenterRotatedCorr(:,2),'.')
plot(CenterRotatedCorr(pcaXz>5,1),CenterRotatedCorr(pcaXz>5,2),'.')
plot(CenterRotatedCorr(pcaYz>10,1),CenterRotatedCorr(pcaYz>10,2),'.')

subplot(1,3,3)
hold on
plot(CenterRotatedCorr(:,1),CenterRotatedCorr(:,2),'.')
plot(CenterRotatedCorr(errors,1),CenterRotatedCorr(errors,2),'.')



%}













