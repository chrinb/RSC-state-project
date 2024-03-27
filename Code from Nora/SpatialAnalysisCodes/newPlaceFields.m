function newPlaceFields(sData)

% Note: Grienbergers definition for new place field:
% First she looks for a Ca2+ transient that is more than 3 standard deviations (SD) about the noise (me: baseline). 
% Then, in the 6 laps after that, there needs to be a transient that is at least 3 SD about the noise in 3 out of 6 laps.
% My modified version: after induction there needs to be 5 transients at the same bin(s) in 5 out of 15 laps (30%)


mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'NewPC');
savePath = strcat(sData.sessionInfo.savePath,'\Imaging\NewPC');

% set which data set you wanna use, parameters:
VelMin = 0.1;
FigVisible = 'off';
TestTrials = 6; % how many laps to test after induction
Reliability = 50; % percentage
TestBins = 20; % how many bins test after indution for the shift
% set
nBin = sData.behavior.meta.nBins;
nTrials = sData.behavior.wheelLapImaging-1; 
nROIs = sData.imdata.nROIs;
%signal = sData.imdata.MaoPC_deconv_IO3_RI01; % Note: I set the deconvolution parameters such, that it detects spikes if it is above 3 SNR. So every 'spike' in the deconvolution data are considered as a potential spike
signalPre = sData.imdata.roiSignals(2).dff;
signal = NaN(size(signalPre)); 
% search for transients which are 3x higher than SD of data
SD = NaN(nROIs,1);
for i = 1:1:nROIs
    SD(i) = std(signalPre(i,:));
    signalTemp = signalPre(i,:);
    signalTemp(signalTemp < 3*SD(i)) = 0;
    signal(i,:) = signalTemp;
end

%%% plot signal in heatplots
signalBinned = signalHeatPlot(signal,sData,VelMin,FigVisible);

% palce cells find based on certain ctiteria (I used similar to Mao, expcet that Reliability needs to be > 50%)
SignalZeroOne = signalBinned;
for roi = 1:1:nROIs
    SignalZeroOne{1,roi}(SignalZeroOne{1,roi} > 0) = 1; % logical transformation of the spike-matrix. 0: no activity, 1: activity
end
NewPlaceCellCollection = struct;
NewPlaceCellCollection.PF = NaN(1,8);
NewPlaceCellCollection.legend(1,1) = convertCharsToStrings('roiID'); NewPlaceCellCollection.legend(1,2) = convertCharsToStrings('PFindTrial'); NewPlaceCellCollection.legend(1,3) = convertCharsToStrings('PFindBin');
NewPlaceCellCollection.legend(1,4) = convertCharsToStrings('PFindLength'); NewPlaceCellCollection.legend(1,5) = convertCharsToStrings('PFindBinCM'); NewPlaceCellCollection.legend(1,6) = convertCharsToStrings('PFshiftBinCM'); 
NewPlaceCellCollection.legend(1,7) = convertCharsToStrings('veloInd'); NewPlaceCellCollection.legend(1,8) = convertCharsToStrings('PFestabLength');
nNewPC = 0; % number of new place fields
for roi = 1:1:nROIs
    SignalBefore = circshift(SignalZeroOne{1,roi},1,1);
    SignalBefore(1,:) = 0;
    SignalAfter = circshift(SignalZeroOne{1,roi},-1,1);
    SignalAfter(end,:) = 0;
    TempBinnedSignal = horzcat(SignalBefore,SignalZeroOne{1,roi},SignalAfter); % circularize data
    % transform data: 0: no activity, 1: start of transient, 2: ongoing transient
    for p = 1:1:nTrials
        for q = 2:1:3*nBin
            if TempBinnedSignal(p,q)==1 && TempBinnedSignal(p,q-1)==1 || TempBinnedSignal(p,q)==1 && TempBinnedSignal(p,q-1)==2
               TempBinnedSignal(p,q) = 2;
            end
        end
    end    
%%% looking for new place fields    
    for i = 1:1:nTrials 
       for j = nBin+1:1:2*nBin
          if TempBinnedSignal(i,j) == 1 % looks for the starting point of transients, there should not be activity in the previous 6 laps 
              %%% test if activity on the previous 6 laps was zero in this and neighbouring bins 
              PreActivity = zeros(nTrials,3); 
               if i > 6
                  ii = i-6;
               else   
                  ii = 1; % if activity if before 6th laps, I can test fewer laps before
               end
               for r = ii:1:i-1 % e.g. test from lap 10 to 16
                  PreActivity(r,1:3) = TempBinnedSignal(r,j-1:j+1); 
               end
               if max(max(PreActivity))> 0 % if there was actity in the previous 6 laps, go to the next candidate
                   continue
               end
               %%% continue if it passes the criteria
               PFStartTrial = i;
               PFStartBin = j-nBin;
               PFlength = 1; 
               for r = j+1:1:(4*nBin-j)
                    if TempBinnedSignal(i,r) > 0
                        PFlength = PFlength + 1;
                    else
                        break
                    end
               end
               PFBinCMStart = floor(PFStartBin + (PFlength / 2)); % center of mass of induction PF
               if PFStartTrial + TestTrials > nTrials
                   break
               end
               % generate an array showing in which laps after induction there was activity in the same bins (+1, -1 bins) as during induction
               AfterInductionArray = TempBinnedSignal(PFStartTrial + 1:PFStartTrial + TestTrials,PFStartBin+nBin-1:PFStartBin+nBin+1);
               AfterInductionArray(AfterInductionArray == 2) = 1;
               AfterInductionArraySum = max(AfterInductionArray,[],2); 
               if sum(AfterInductionArraySum) >= TestTrials*(Reliability/100) % considered as newly formed place field, if 1/2 of the induction following trials have transients
                  nNewPC = nNewPC + 1;
                  NewPlaceCellCollection.PF(nNewPC,1) = roi; % write the ROI number into column one is it is a place cell
                  NewPlaceCellCollection.PF(nNewPC,2) = PFStartTrial; % induction lap
                  NewPlaceCellCollection.PF(nNewPC,3) = PFStartBin; % induction bin
                  NewPlaceCellCollection.PF(nNewPC,4) = PFlength; % induction transient length (bin) 
                  NewPlaceCellCollection.PF(nNewPC,5) = PFBinCMStart; % Center of Mass of induction lap transient
                  % calculate center of mass and length of established place field after the induction laps
                  PFestabLengthArray = NaN(TestTrials,1);
                  ActTempMeanPos = NaN(TestTrials,1);
                  for n = 1:1:TestTrials % check all test trials after induction where there was activity in place-bins (where in induction lap)
                      if  AfterInductionArraySum(n) > 0
                          ActivityTemp = TempBinnedSignal(PFStartTrial+n,nBin+PFBinCMStart-TestBins:nBin+PFBinCMStart+TestBins);
                          % get rid of activation which happens searately from induction position
                          % if there is a peak at induction bin, keep it, and discard other activity in more distant bins
                          
                          if ActivityTemp(TestBins+1) > 0
                              for s = TestBins:-1:1
                                  if ActivityTemp(s)==0
                                     ActivityTemp(1:s)=0; 
                                  end
                              end
                              for s = TestBins+2:1:2*TestBins+1
                                  if ActivityTemp(s)==0
                                     ActivityTemp(s:2*TestBins+1)=0; 
                                  end
                              end
                          else % if no actiivty in induction bin, but there is in negihboring bins
                          %}
                              [~,peaklocation,peakWidth,~] = findpeaks(ActivityTemp); % find largest/widest peak
                              ActivityTemp(1:peaklocation(find(max(peakWidth)))-2)=0;
                              ActivityTemp(peaklocation(find(max(peakWidth)))+floor(peakWidth(find(max(peakWidth)))):2*TestBins+1)=0;
                          end                              
                          for o = 1:1:2*TestBins+1
                              if ActivityTemp(o) > 0
                                 ActivityTemp(o) = o - TestBins - 1; % zero will be the same position as on CM in induction lap
                              elseif ActivityTemp(o) == 0
                                  ActivityTemp(o) = NaN;
                              end
                          end
                           % calculate length of activation (in bins)
                          PFestabLengthArray(n) = sum(~isnan(ActivityTemp));
                          ActTempMeanPos(n) = nanmean(ActivityTemp);
                      end
                  end
                  NewPlaceCellCollection.PF(nNewPC,6) = nanmean(ActTempMeanPos);
                  NewPlaceCellCollection.PF(nNewPC,7) = sData.behavior.binning.veloBinned(i,j-nBin);
                  NewPlaceCellCollection.PF(nNewPC,8) = nanmean(PFestabLengthArray);
               end
           end
       end
    end
end

%%% figures

% induction velo vs PF width 
figure('Color','white');
scatter(NewPlaceCellCollection.PF(:,7),NewPlaceCellCollection.PF(:,8)*sData.behavior.meta.binSize)
xlabel('Velocity during induction (cm/s)');
ylabel('Place Field Width (cm)');
title(strcat(sData.sessionInfo.fileID,'-Velocity vs established place field width'));
FileName = strcat(sData.sessionInfo.fileID,'-VeloVsPFwidth');
savefig(fullfile(savePath,FileName));
saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

% PFwidth shift 
figure('Color','white');
histogram(NewPlaceCellCollection.PF(:,6)*sData.behavior.meta.binSize)
xlabel('Shift in place tuning after induction (cm)');
ylabel('Count (place fields)');
title(strcat(sData.sessionInfo.fileID,'-Shift In Place Tuning'));
FileName = strcat(sData.sessionInfo.fileID,'-ShiftInPlaceTuning');
savefig(fullfile(savePath,FileName));
saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

% Cumulative distribution plot for place field formation
% convert nTrials into a percentage where 100% is the end of session
MaxTrials = nTrials - TestTrials;
CumDistr = NewPlaceCellCollection.PF(:,2)/MaxTrials*100;

figure('Color','white');
cdfplot(CumDistr) 
xlabel('Lap of place field induction');
ylabel('Cumulative fraction of place cells');
title(strcat(sData.sessionInfo.fileID,'-Cumulative distribution of lap of place field induction'));
FileName = strcat(sData.sessionInfo.fileID,'-CumDistPF');
savefig(fullfile(savePath,FileName));
saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

% Save file to same path where LV files can be found 
save(fullfile(savePath,strcat(sData.sessionInfo.fileID,'_NewPlaceCellCollection.mat')),'NewPlaceCellCollection');

end


