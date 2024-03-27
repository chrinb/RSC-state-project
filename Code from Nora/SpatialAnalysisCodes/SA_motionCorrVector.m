function sData = SA_motionCorrVector(sData,filePath,savePath) % filePath: raw data recording folder, containing frameShifts.mat

%%% calculate motion correction error matrix (sum of squares)

% open motion correction / imreg_variables / _frameShifts_ch1.mat
str1 = filePath;
str2 = regexprep(str1,'roisignals(\w*)','','ignorecase');
%str21 = str2(1:end-1);
str3 = strcat(str2,'imreg_variables\'); % str21 instead of str2
List = dir(fullfile(str3,'*_frameShifts_ch1.mat'));
load(fullfile(str3,List.name)); %#ok<LOAD>

% set nBins, nTrials
nBin = sData.behavior.meta.nBins;
nTrials = sData.behavior.wheelLapImaging;
MaxSamplesinBin = max(sData.behavior.binning.samplesSpentInBin(:)); 
nSamples = sData.imdata.nSamples;

% calculate squared root error of x and y shifts:
Motion = rssq(frameShifts,2);
MotionMatrix = NaN(nTrials,nBin); 
  
% load matrices with data in trial and bins
for j = 1:1:nTrials-1
    for k = 1:1:nBin
        SampleInd = sData.behavior.binning.enterIntoBinIndex(j,k); % get sample index when enter into a bin
        if isnan(SampleInd) % if recording ends in LV stop calculation
           break 
        end
        SpentInBin = sData.behavior.binning.leaveBinIndex(j,k) - sData.behavior.binning.enterIntoBinIndex(j,k) + 1; % get how many samples spent in that bin
        ErrorInBin = NaN(MaxSamplesinBin,1); % temporary array to calcuate mean Ca-value in a bin during time spent in a bin (> velo lim)
        for m = 1:1:SpentInBin
            if SampleInd+m > nSamples % if recording ends in SciScan stop calculation (sometimes not the same size data in SciScan and LV)
                break
            end
            %if SampleInBinLim(SampleInd+m-1,k) > 0 % limited values were set to -1, I do not want to contain them
              ErrorInBin(m) = Motion(SampleInd+m-1); % Error
            %end
        end
        MotionMatrix(j,k) = nanmean(ErrorInBin); % mean of data within a bin witihn a trial
    end
end

MeanMotionMatrix = nanmean(MotionMatrix,1);

sData.imdata.binned.MotCorrBinned = MotionMatrix;
sData.imdata.binned.MeanMotCorrBinned = MeanMotionMatrix;


% plot
% cue positions: %{
C1A = 23; C1B = 29; % velcro
C2A = 43; C2B = 49; % hot glue
C3A = 63; C3B = 69; % hot glue
C4A = 83; C4B = 89; % velcro

%savePath = 'C:\MATLAB\SAVE';

%figure();
nTrials = sData.behavior.wheelLapImaging;
SA_plotHeatBinCa(sData.imdata.binned.MotCorrBinned,sData.sessionInfo.fileID,1,'MotionCorrVector',sData.behavior.meta.binSize,sData.behavior.meta.nBins,nTrials,'on'); %plotHeatBinCa(data,fileID,roi,ylab,BinSize,nTrials,figure visibility)
caxis([0 inf]); hold on;
line([C1A C1A],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on; line([C1B C1B],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on;
line([C2A C2A],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C2B C2B],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
line([C3A C3A],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C3B C3B],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
line([C4A C4A],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on; line([C4B C4B],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on;
axis([0 160 0 nTrials]); % ceil(Ymax)
title(strcat(sData.sessionInfo.fileID,' binned motion artefact (calc. from aligment vectors)'));
xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
ylabel('Binned rssq value');
fname = strcat(sData.sessionInfo.fileID,'-BinnedMotionArtHeat');
savefig(fullfile(savePath,fname));
saveas(gcf,(fullfile(savePath,[fname '.jpg'])));

figure();
Xaxis = sData.behavior.meta.binSize:sData.behavior.meta.binSize:sData.behavior.meta.binSize*nBin;
Ymax = (max(MeanMotionMatrix(1,:)))*1.1;
plot(Xaxis,MeanMotionMatrix(1,1:nBin),'LineWidth',2)
line([C1A C1A],[0 Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on; line([C1B C1B],[0 Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on;
line([C2A C2A],[0 Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C2B C2B],[0 Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
line([C3A C3A],[0 Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C3B C3B],[0 Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
line([C4A C4A],[0 Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on; line([C4B C4B],[0 Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on;
axis([0 160 0 Ymax]); % ceil(Ymax)
title(strcat(sData.sessionInfo.fileID,' binned motion artefact (calc. from aligment vectors)'));
xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
ylabel('Mean binned rssq value');
fname = strcat(sData.sessionInfo.fileID,'-meanBinnedMotionArt');
savefig(fullfile(savePath,fname));
saveas(gcf,(fullfile(savePath,[fname '.jpg'])));

end