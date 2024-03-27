function sData = SA_binning(sData,nBins) 

% BINNING (nBins: number of spatial bins)
% calculating which sample belong to which bin
nSamples     = sData.behavior.details.nFrames;
BinSize      = (48.7*pi)/nBins;
trial_nr     = sData.behavior.wheelLap;
SampleInBin  = NaN(nSamples,nBins);    

% first search all data for the specified distance bin (indexes for 
% entering a new bin can be calculated based on AbsDistDS data, j is distance cm, bin on the wheel (e.g. 0-1cm)

% Loop over nr of bins - 1 (why -1?
for j = 1:(nBins-1)

    % Loop over nr of samples (frames)
    for i = 1:nSamples 
       % scan through the whole dataset and sort into bin j, between trials, there are NaNs on the column

       % For each position sample (sData.behavior.wheelPosDsMonIncr(i)),
       % check if it is bigger than previous bin j and smaller than current
       % bin j. If true, put sample index into the bin column. 
       if  sData.behavior.wheelPosDsMonIncr(i) >= (j-1)*BinSize && sData.behavior.wheelPosDsMonIncr(i) < j*BinSize  %check ALL distanceDS datapoint if it belongs to the specified bin, if yes, continue
           SampleInBin(i,j) = i; 
       end
    end 
end
if j == nBins-1 % last Bin, which is not a full bin or bigger, collect together what is before lapstart into last bin
   for i = 1:1:nSamples
      if  sData.behavior.wheelPosDsMonIncr(i) > j*BinSize 
          SampleInBin(i,j+1) = i;          
      end
   end  
end

%% Correction for too fast trials, when one sample must be assign to two/more bins
% circularize data:
HeightPlus10Bin                                         = size(SampleInBin,1)-1;
SampleInBinPlusBins                                     = NaN(HeightPlus10Bin,nBins+10); % Why add 10 bins???
SampleInBinPlusBins(1:HeightPlus10Bin,1:nBins)          = SampleInBin(1:HeightPlus10Bin,1:nBins);

% Below, the first 10 columns of the original SampleInBin matrix is
% inserted into the added 10 columns of SameInBinPlusBins matrix
SampleInBinPlusBins(1:HeightPlus10Bin,nBins+1:nBins+10) = SampleInBin(1:HeightPlus10Bin,1:10);

% Correction

% Loop over nr of samples -1
for j = 1:nSamples-1 
    % Loop over nr of bins + 10 - 1
    for i = 1:(nBins+10-1) 
       
       % If current sample is not NaN and next sample is NaN
       if  ~isnan(SampleInBinPlusBins(j,i)) && isnan(SampleInBinPlusBins(j+1,i))

           % Then check if the sample in the next column (bin) is NaN AND
           % that the subsequent sample in that column is NaN. If true,
           % correct. 
           if isnan(SampleInBinPlusBins(j,i+1)) && isnan(SampleInBinPlusBins(j+1,i+1))
               SampleInBinPlusBins(j,i+1) = SampleInBinPlusBins(j,i); % if in a given bin no sample assign (too fast running), copy the previous bin's data there
           end
       end
    end
end
SampleInBinCorrect(1:HeightPlus10Bin,1:nBins) = SampleInBinPlusBins(1:HeightPlus10Bin,1:nBins);
SampleInBinCorrect(1:HeightPlus10Bin,1:5)     = SampleInBinPlusBins(1:HeightPlus10Bin,nBins+1:nBins+5);


%% Make 3 matrices:
% (1) EnterIntoBinSampleInd/ (2) LeaveBinSampleInd: when the animal enter/leave(in the next sample) a bin (sample ind, row: trial, col: bin); 
% (3) SampleSpentInBin: how many samples spend the aimal in this bin. At the end of recording there are NaNs in the matrices 

% Fix: For one mouse the 

SampleInBinIsNaN1 = isnan(SampleInBinCorrect);

% Fix for m6162 session RF day 2...
% SampleInBinIsNaN1(723, 1:5) = true;
% SampleInBinIsNaN1(1031, 1:5) = true;

SampleInBinIsNaN2 = diff(SampleInBinIsNaN1,1) ==-1; % matrix shows first sample in a bin. It is 1 and other is 0. (row: trial, column: bin)
if SampleInBinIsNaN1(1,1) == 0 % rare case when recording start in bin 1, it is needed to be changed in order to detect
    SampleInBinIsNaN2(1,1) = 1; 
end
SampleInBinIsNaN3                           = diff(SampleInBinIsNaN1,1)==1; % matrix when last sample spent in a given bin, set to 1 and other is 0.
%LastSampleFix = find(SampleInBinIsNaN1(Samples,:) == 0);
LastSampleFix                               = SampleInBinIsNaN1(nSamples-1,:) == 0;  % I have to drop the last sample in order the code to function prperly
SampleInBinIsNaN3(nSamples-1,LastSampleFix) = 1 ; % drop last sample 
EnterIntoBinSampleInd                       = NaN(trial_nr+1,nBins); % sample index when animal enter into given bin given trial
LeaveBinSampleInd                           = NaN(trial_nr+1,nBins); % sample index when animal leaves a given bin

for i = 1:nBins % I need to do it in a complicated way using temporary arrays because if the last trial did not go to the end it gave error (mismatch in TRNu and bin-start at later bins)
    
%     TempArray1 = zeros(trial_nr+1,1);
    TempArray1 = find(SampleInBinIsNaN2(:,i)==1)+1; % first sample spent in a given bin, given trial
    if numel(TempArray1)<trial_nr+1
        TempArray1(trial_nr+1) = NaN;
    end
    EnterIntoBinSampleInd(:,i) = TempArray1;
    
%     TempArray1 = zeros(trial_nr+1,1);
    TempArray2 = find(SampleInBinIsNaN3(:,i)==1); % last sample spent in a given bin, given trial
    if numel(TempArray2)<trial_nr+1
        TempArray2(trial_nr+1) = NaN;
    end
    LeaveBinSampleInd(:,i) = TempArray2; % last sample spent in a given bin, given trial
end

SampleSpentInBin = LeaveBinSampleInd - EnterIntoBinSampleInd + 1; % spent in bin = last sample - first sample +1
% correction for bins in which the animal spent half bin time
% circularize matrix for calculation
HeightPlus10Bin = size(EnterIntoBinSampleInd,1)-1;
EnterIntoBinSampleIndPlusBins = NaN(HeightPlus10Bin,nBins+10);
EnterIntoBinSampleIndPlusBins(1:HeightPlus10Bin,1:nBins) = EnterIntoBinSampleInd(1:HeightPlus10Bin,:);
EnterIntoBinSampleIndPlusBins(1:HeightPlus10Bin,nBins+1:nBins+10) = EnterIntoBinSampleInd(2:HeightPlus10Bin+1,1:10);
SampleSpentInBinPlusBins = NaN(HeightPlus10Bin,nBins+10);
SampleSpentInBinPlusBins(1:HeightPlus10Bin,1:nBins) = SampleSpentInBin(1:HeightPlus10Bin,:);
SampleSpentInBinPlusBins(1:HeightPlus10Bin,nBins+1:nBins+10) = SampleSpentInBin(2:HeightPlus10Bin+1,1:10);


for i = 1:1:trial_nr
    for j = 1:1:nBins+10-1 %10 extra bin was added for circularization
        counter = 1;
        for k = 0:1:nBins+10-j-1
            if EnterIntoBinSampleIndPlusBins(i,j+k) == EnterIntoBinSampleIndPlusBins(i,j+k+1)
               counter = counter + 1;
            else
                break
            end
        end
        if counter > 1
            FrameSpent = SampleSpentInBinPlusBins(i,j)/counter;
            for k = 0:1:counter-1
                SampleSpentInBinPlusBins(i,j+k) = FrameSpent; % share the time spent there
            end
        end
    end
end
SampleSpentInBinCorrect(1:HeightPlus10Bin,1:nBins) = SampleSpentInBinPlusBins(1:HeightPlus10Bin,1:nBins);
SampleSpentInBinCorrect(2:HeightPlus10Bin,1:5)     = SampleSpentInBinPlusBins(1:HeightPlus10Bin-1,nBins+1:nBins+5);

sData.behavior.binning.enterIntoBinIndex            = EnterIntoBinSampleInd;
sData.behavior.binning.leaveBinIndex                = LeaveBinSampleInd;
sData.behavior.binning.samplesSpentInBin            = SampleSpentInBinCorrect;
sData.behavior.binning.enterIntoBinIndexExtended    = EnterIntoBinSampleIndPlusBins;
sData.behavior.binning.SampleSpentInBinExtendedBins = SampleSpentInBinPlusBins; 
sData.behavior.binning.samplesInBinIndex            = SampleInBinCorrect;

% control for the monotonic transfer, what was the maximum difference
% MaxPositionDiff = max(behav.wheelPosDsMonIncr - behav.wheelPosDs);

end
