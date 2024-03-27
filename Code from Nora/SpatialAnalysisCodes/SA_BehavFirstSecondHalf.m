function sData = SA_BehavFirstSecondHalf(sData,LapsTested) % LapsTested = 30; % how many laps take to analyse in the beginning and end

sData.behavior.FirstSecondHalf = struct;
%savePath = 'C:\MATLAB\SAVE';
mkdir(sData.sessionInfo.savePath,'FirstSecondHalf');
savePath = strcat(sData.sessionInfo.savePath,'\FirstSecondHalf');

% set which laps belong to the first half and second half of the session 
LastLap = sData.behavior.wheelLapImaging-1; 
%LapsTested = 10; 
FirstStartLap  = 1;
FirstEndLap    = FirstStartLap + LapsTested - 1;

% Correct for last lap count being higher than the binned matrix
if LastLap > size(sData.behavior.binning.lickPerCmBinned,1)
    LastLap = size(sData.behavior.binning.lickPerCmBinned,1);
end
SecondStartLap = LastLap - LapsTested + 1;
SecondEndLap   = LastLap;

   
%%% PLOT mean Licks
sData.behavior.FirstSecondHalf.LickPerCmBinned.First  = sData.behavior.binning.lickPerCmBinned(FirstStartLap:FirstEndLap,:);
sData.behavior.FirstSecondHalf.LickPerCmBinned.Second = sData.behavior.binning.lickPerCmBinned(SecondStartLap:SecondEndLap,:);

% mean Lick First-Second half 
RewardZoneLength = 6;
figure('Color','white'); 
Xaxis = sData.behavior.meta.binSize:sData.behavior.meta.binSize:sData.behavior.meta.binSize*(sData.behavior.meta.nBins);
%Ymax = ceil(max(mean(sData.behavior.FirstSecondHalf.LickPerCmBinned.First)));
plot(Xaxis,mean(sData.behavior.FirstSecondHalf.LickPerCmBinned.First)); hold on;
plot(Xaxis,mean(sData.behavior.FirstSecondHalf.LickPerCmBinned.Second)); hold on;
%line([157 157],[0 Ymax],'Color','black','LineStyle','--'); hold on;
%line([157+RewardZoneLength 157+RewardZoneLength],[0 Ymax],'Color','black','LineStyle','--');
xlabel('Position on wheel (cm)');
ylabel('Licks/cm');
legend(strcat('Lap:',num2str(FirstStartLap),'-',num2str(FirstEndLap)),strcat('Lap:',num2str(SecondStartLap),'-',num2str(SecondEndLap)),'Location','North');
FileName = strcat(sData.sessionInfo.sessionID,'LickFirstSecondHalf-lap-',num2str(FirstStartLap),'-',num2str(FirstEndLap),'-',num2str(SecondStartLap),'-',num2str(SecondEndLap));
savefig(fullfile(savePath,FileName));
saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

%%% PLOT mean Velo
sData.behavior.FirstSecondHalf.VeloBinned.First = sData.behavior.binning.veloBinned(FirstStartLap:FirstEndLap,:);
sData.behavior.FirstSecondHalf.VeloBinned.Second = sData.behavior.binning.veloBinned(SecondStartLap:SecondEndLap,:);

% mean Velo First-Second half  %(sData.behavior.FirstSecondHalf.VeloBinned.OptoOffFirst)
figure('Color','white'); 
Xaxis = sData.behavior.meta.binSize:sData.behavior.meta.binSize:sData.behavior.meta.binSize*(sData.behavior.meta.nBins);
%Ymax = ceil(max(mean(sData.behavior.FirstSecondHalf.VeloBinned.First)));
plot(Xaxis,mean(sData.behavior.FirstSecondHalf.VeloBinned.First)); hold on;
plot(Xaxis,mean(sData.behavior.FirstSecondHalf.VeloBinned.Second)); hold on;
%line([157 157],[0 Ymax],'Color','black','LineStyle','--'); hold on;
%line([157+RewardZoneLength 157+RewardZoneLength],[0 Ymax],'Color','black','LineStyle','--');
xlabel('Position on wheel (cm)');
ylabel('Speed (cm/s)');
legend(strcat('Lap:',num2str(FirstStartLap),'-',num2str(FirstEndLap)),strcat('Lap:',num2str(SecondStartLap),'-',num2str(SecondEndLap)),'Location','South');
FileName = strcat(sData.sessionInfo.sessionID,'SpeedFirstSecondHalf-lap-',num2str(FirstStartLap),'-',num2str(FirstEndLap),'-',num2str(SecondStartLap),'-',num2str(SecondEndLap));
savefig(fullfile(savePath,FileName));
saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));     


% Save file to same path where other files can be found 
% save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.sessionID,'_sData.mat')),'sData');

end

