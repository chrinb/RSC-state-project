function sData = SA_placeCellMaosData(sData, datatype) 
% Calculating which ROI is considered as a place cell based on Mao et al., 2017, criteria#1
% Note: Mao used 1.5 cm bins, I use 2 cm bins and minimum velocity of 1 cm/s (set in previous scripts). 
% Mao summated the activity within a bin and normalized by the occupancy (the time spent) in each bin. I used mean activity within a bin, which is the same (just have to divide all data with SampleTime= 32.3 msec)

multi_session_idx              = sData.imdata.roi_classification == 1;

%% PARAMATERS TO SET:
% type = 0; 0: use dff data, 1: use deconvolved data, -1: use dff, where slow transients are removed 
% InOutRatio = 3; Reliability = 0.34; regarding Mao, for VIP cells IOR = 2, Rel=0.2
% ROIs: array containing ROI IDs to analyse. If it is 0 : analyse all ROIs

% CRITERIA 1 the maximum activity (peak) in the potential place field must have been at least (1+activityTreshold) -fold larger (Mao used 1.3x for deconv signal) compared to minimum activity in pos tuning curve

if datatype == 0 || datatype == -1 % for dff data, I set new paramterers on 2021.08.15.
    activityTreshold  = 0.5; % I set place cell tuning around peak to be 50% as larger as minimum, and place field act should be 2x larger than out of field with 20% reliability
    MinPlaceFieldSize = 2; % cm minimum size of place field Mao:15
    MaxPlaceFieldSize = 100; % cm maximum size of place field Mao:120
    
    % CRITERIA 2 the mean in-field activity must be at least three times larger than the mean out-of-field activity;
    % Note: setting it to a low value (0.3) allow larger jitter in peak to fit into place cell category 
    InOutRatio = 1.5; 
    
    % CRITERIA 3 more than one third of the trials must have peak position-mapped activity fall within the potential place field
    ReliabiliyIndex = 0.25; 
    % tests on a session (m8058-20210527-01) showed the best parameters to find place cells (which looks like place cells based on visual inspection of spatially binned Calcium data)
    % activityTreshold = 0.7; InOutRatio = 2; ReliabiliyIndex = 0.2; OR activityTreshold = 0.5; InOutRatio = 2; ReliabiliyIndex = 0.3; OR 

else % deconvolved signal
    activityTreshold  = 0.3; % Mao = 0.3
    GaussianFilter    = round(10/sData.behavior.meta.binSize); % I set 10 cms to smooth (5 bins for 2 cm Binsize)
    MinPlaceFieldSize = 1; % cm minimum size of place field Mao used 15 for deconvolved and little bit smoothed data
    MaxPlaceFieldSize = 60; % cm maximum size of place field Mao used 120 for deconvolved and little bit smoothed data
    
    %CRITERIA 2 the mean in-field activity must be at least three times larger than the mean out-of-field activity;
    InOutRatio = 2; % Mao = 3
    
    %CRITERIA 3 more than one third of the trials must have peak position-mapped activity fall within the potential place field
    ReliabiliyIndex = 0.25; % Mao = 0.34
end

%% INITIALIZE PARAMETERS 
% Load parameters:
Mao                          = struct; % collect place cell data
Mao.params.activityTreshold  = activityTreshold;
Mao.params.InOutRatio        = InOutRatio;
Mao.params.ReliabiliyIndex   = ReliabiliyIndex;
Mao.params.MinPlaceFieldSize = MinPlaceFieldSize;
Mao.params.MaxPlaceFieldSize = MaxPlaceFieldSize;
bin_nr                       = sData.behavior.meta.nBins;
trial_nr                     = sData.behavior.wheelLapImaging;
BinSize                      = sData.behavior.meta.binSize;
% roi_nr                       = sum(sData.imdata.roi_classification==1);
roi_nr                       = sData.imdata.nROIs;

% Create arrays for processed data, set up saving path
dataset = struct;
if datatype == 0 
    mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'MaoPlaceCells-dFF');
    SavePath     = strcat(sData.sessionInfo.savePath,'\Imaging\MaoPlaceCells-dFF');
    dataset      = sData.imdata.binned.RoidFF;
    Mao.datatype = 'dFF';

elseif datatype == -1
    mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'MaoPlaceCells-dFF-SR');
    SavePath     = strcat(sData.sessionInfo.savePath,'\Imaging\MaoPlaceCells-dFF-SR');
    dataset      = sData.imdata.binned.RoidFFSR;
    Mao.datatype = 'dFF-slowremoved';

elseif datatype == 1
    mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'MaoPlaceCells-deconv');
    SavePath     = strcat(sData.sessionInfo.savePath,'\Imaging\MaoPlaceCells-deconv');
    dataset      = sData.imdata.binned.RoiDeconvolved;
    Mao.datatype = 'deconvolved';
end

% calculate position tuning curve for each ROI (mean activity at spatial bins throughout the session)
Mao.PosTuning = NaN(roi_nr,bin_nr);
for i = 1:roi_nr
    if multi_session_idx(i) == 1
        Mao.PosTuning(i,1:bin_nr) = mean(dataset{1,i},1, 'omitnan');  
    end
end

%%% use 3 bin gaussian smoothing for deconvolved data.
if datatype == 1
    GaussianSmoothed = NaN(roi_nr,bin_nr);

    for i = 1:roi_nr
        if multi_session_idx(i) == 1
            GaussianSmoothed(i,:) = smoothdata(Mao.PosTuning(i,:),'gaussian',GaussianFilter);
        end
    end
    Mao.PosTuning = GaussianSmoothed;
end


% Normalize mean activity for visualization
Mao.PosTuningNorm = NaN(roi_nr, bin_nr);
MaxInROI          = max( Mao.PosTuning,[],2); % search for maximum value (in a bin) in each ROI
Mao.PosTuningNorm = Mao.PosTuning ./ MaxInROI; % normalize

%% CRITERIA #1: place fields must be a continuous region between MinPlaceFieldSize and MaxPlaceFieldSize (cm) in width, within which the activity magnitude 
%  must be above a certain percentage of the difference between the maximum and minimum activity in the position tuning curve;
MinInROI               = min(Mao.PosTuning,[],2);
MinInROI(MinInROI<0)   = 0;
Treshold               = MinInROI + (MaxInROI - MinInROI)*activityTreshold;   %%% activity must be larger in place field than this value (e.g. 130% of minimum)
LargerThanTresholdBins = zeros(roi_nr,bin_nr);

for i = 1:roi_nr
    if multi_session_idx(i) == 1

        SignalTemp3 = Mao.PosTuning(i,:);
        SignalTemp3(SignalTemp3<Treshold(i)) = NaN; % put NaN to bins where activity is smaller than treshold
        LargerThanTresholdBins(i,:) = SignalTemp3; % matrix tresholded activity 
    end
end

LargerThanTresholdBins2 = ~isnan(LargerThanTresholdBins); % ones when the activity is above treshold. Count ones next to each other
ExtraZeros              = zeros(roi_nr,1);
LargerThanTresholdBins3 = [ExtraZeros LargerThanTresholdBins2 LargerThanTresholdBins2]; % concatenate three matrices to lengten trials 
LargerThanTresholdBins4 = NaN(roi_nr,2*bin_nr+1); % calculate potential place field size. 
kmax                    = MaxPlaceFieldSize/BinSize; % do not check biger size then what was set as maximum place field size 

for i = 1:roi_nr
    if multi_session_idx(i) == 1

        for j = 1:2*bin_nr % - ceil(kmax) %
          if LargerThanTresholdBins3(i,j+1) == 1 % find the first / next one, in this matrix data are shifted with one column ->(j+1)
              for k = 1:bin_nr+kmax-j % count how many ones are after it
                  if LargerThanTresholdBins3(i,j+1+k) == 0 % if there is a zero (low activity), write the length of place field and jump to the next ones           
                      LargerThanTresholdBins4(i,j+1) = k; % how many bin long place field is it , I put an extra column in the beginning
                      break
                  end
              end
          end
        end
    end
end

PlaceFieldSize = NaN(roi_nr,bin_nr);
for i = 1:roi_nr
    if multi_session_idx(i) == 1

        SignalTemp5 = LargerThanTresholdBins4(i,1:bin_nr+1);
        if ~isnan(SignalTemp5(2)) && ~isnan(SignalTemp5(bin_nr+1)) % if there is activity at the end and beginning of track as well
            % delete values from the beginning of the track if there is already high activity at the end
            SignalTemp5(2:SignalTemp5(2)+1) = 0;
        end
        for j = 1:bin_nr
            if isnan(SignalTemp5(j)) && ~isnan(SignalTemp5(j+1)) 
                PlaceFieldSize(i,j) = SignalTemp5(j+1); % j+1, because it was sifted in LargerThanTresholdBins matrix to the right
            end
        end
    end
end
% Collect ROIs having at least one potential place field
Mao.Criteria1Passed = NaN(roi_nr,1);
FieldMinBin = round(MinPlaceFieldSize/BinSize);
FieldMaxBin = round(MaxPlaceFieldSize/BinSize);
for i = 1:roi_nr
    if multi_session_idx(i) == 1

        SignalTemp6 = PlaceFieldSize(i,:);
        if any(SignalTemp6<FieldMaxBin) && any(SignalTemp6>FieldMinBin)  
            Mao.Criteria1Passed(i)=i;
        end
    end
end

% How many potential place fields do they have and in which bin?
PotPlaceFieldNu       = NaN(roi_nr,1);
PotPlaceFieldStartBin = NaN(roi_nr,20);
PotPlaceFieldLength   = NaN(roi_nr,20);
for i = 1:roi_nr
    if multi_session_idx(i) == 1

        counter = 0;
        for j = 1:bin_nr
            if PlaceFieldSize(i,j) < FieldMaxBin && PlaceFieldSize(i,j) > FieldMinBin % if the size of the potential place field fits to what was set up as min-max size
                counter = counter + 1;
                PotPlaceFieldStartBin(i,counter) = j;
                PotPlaceFieldLength(i,counter) = PlaceFieldSize(i,j);
            end
        end
        PotPlaceFieldNu(i) = counter; 
    end
end
Mao.PotPlaceFieldStartBin = PotPlaceFieldStartBin;
Mao.PotPlaceFieldLength   = PotPlaceFieldLength;


%% CRITERIA 2: the mean in-field activity must be at least X times larger than the mean out-of-field activity; 
% version from 2021.08.20. for background dFF level calculation: mask out the bins which are potential place fields
ActAllTrials     = [Mao.PosTuning Mao.PosTuning]; % concatenate two matrices to have circular data 
ActInPlaceField  = NaN(roi_nr,size(PotPlaceFieldStartBin,2)); 
ActOutPlaceField = NaN(roi_nr,size(PotPlaceFieldStartBin,2));

for i = 1:roi_nr 
    if multi_session_idx(i) == 1

        for j = 1:PotPlaceFieldNu(i) % if ROI has more than one potential place fields, check all of them
            ActInPlaceField(i,j) = mean(ActAllTrials(i,PotPlaceFieldStartBin(i,j):(PotPlaceFieldStartBin(i,j)+PotPlaceFieldLength(i,j)-1)));  
            % If tested if place cell recognition gets better if I delete values of all potential place fields for IO ratio calculation. I thought landmark cells are easily masked out having two peaks.
            % But did not give better place cell recognition, indeed more wide-place-fields (more unselectivity) popped up, so I use the original version
            %{
            ActOutPlaceFieldTemp = Mao.PosTuning(i,1:BinNu);
            for k = 1:1:PotPlaceFieldNu(i)
                ActOutPlaceFieldTemp(PotPlaceFieldStartBin(i,k):(PotPlaceFieldStartBin(i,k)+PotPlaceFieldLength(i,k)-1)) = NaN; 
            end
            ActOutPlaceField(i,j) = nanmean(ActOutPlaceFieldTemp);
            %}
            ActOutPlaceField(i,j) = mean(ActAllTrials(i,(PotPlaceFieldStartBin(i,j)+PotPlaceFieldLength(i,j)):(PotPlaceFieldStartBin(i,j)+bin_nr-1)));
        end
    end
end
ActRatioInOutPlaceField = ActInPlaceField ./ ActOutPlaceField;
Mao.Criteria12Passed    = Mao.Criteria1Passed;

for i = 1:roi_nr
    if multi_session_idx(i) == 1

        if ~any(ActRatioInOutPlaceField(i,:) > InOutRatio) % if non-of the potential place fields pass through criteria 2 , delete from potential place cells (overwrite to NaN)
           Mao.Criteria12Passed(i) = NaN;
        end
    end
end

% update potential place field matrix 
PotPlaceFieldStartBin2                                     = PotPlaceFieldStartBin;
PotPlaceFieldStartBin2(ActRatioInOutPlaceField<InOutRatio) = NaN;
PotPlaceFieldLength2                                       = PotPlaceFieldLength;
PotPlaceFieldLength2(ActRatioInOutPlaceField<InOutRatio)   = NaN;
PotPlaceFieldNu2                                           = sum(~isnan(PotPlaceFieldStartBin2),2); 

%% CRITERIA 3: more than X percentage of the trials must have a peak positioned within the potential place field
Reliability = NaN(roi_nr,1);
PFStartBin = NaN(roi_nr,1); % collect in which bin the place field starts
PFLength = NaN(roi_nr,1);
for i = 1:roi_nr
    if multi_session_idx(i) == 1

        if any(~isnan(PotPlaceFieldStartBin2(i,:))) % if no potential place field, jump to next roi
            counter             = 0; % counter for reliable trials
            PosTuningPeakSearch = NaN(1,bin_nr); % in position tuning curve mask out those parts where no acceptable peak was found  
            PotBins             = PotPlaceFieldStartBin2(i,:); % look for potential PF positions after CRITERIA 2 is tested
            PotBins             = PotBins(~isnan(PotBins));
            PotLength           = PotPlaceFieldLength2(i,:);
            PotLength           = PotLength(~isnan(PotLength));
            %PeakMax = max(Mao.PosTuning(i,PotPlaceFieldStartBin2(i,~isnan(PotPlaceFieldStartBin2(i,:))))); % take the larger peak, find the position
            for j = 1:PotPlaceFieldNu2(i)
                PosTuningPeakSearch(PotBins(j):(PotBins(j)+PotLength(j)-1)) = ActAllTrials(i,(PotBins(j):PotBins(j)+PotLength(j)-1)); % get data from circularized Position Tuning dataset
            end
    
            PeakMax                       = max(PosTuningPeakSearch); % if more peaks pass the criteria, choose the larger peak
            PeakMaxLoc                    = find(PosTuningPeakSearch==PeakMax,1); % find the position (bin) 
            PotBins(PotBins > PeakMaxLoc) = NaN;
            PotBins                       = PotBins(~isnan(PotBins));
            PFStartBin(i)                 = max(PotBins); % search to which potential place field the max peak belongs to 
            PFLength(i)                   = PotPlaceFieldLength2(i,find(PotPlaceFieldStartBin2(i,:)==PFStartBin(i)));
            PFEnd                         = PFStartBin(i) + PFLength(i)-1; % place field end
    
            for j = 1:trial_nr-1
                TrialData = dataset{1,i}(j,:); % put each trial data into a container
                if any(isnan(TrialData)) % after the end  of recording there are NaNs. Go to next ROI
                   break 
                end
                if mean(TrialData) == 0 % sometimes there are only zeros in a row. it makes matlab crazy, jump to next trial
                    continue
                end
                MaxBinFirst = find(TrialData == max(TrialData(:)),1); % find maximum (which bin), sometimes the same max data is in two neighboring position, take the first
                if PFEnd <= bin_nr
                    if MaxBinFirst >= PFStartBin(i) && MaxBinFirst <= PFEnd % peak should be between place field start and end. Place field might continued in to next bin
                        counter = counter + 1;
                    end
                elseif PFEnd > bin_nr % if place field goes  through trial start
                    if MaxBinFirst >= PFStartBin(i) || MaxBinFirst <= rem(PFEnd,bin_nr)
                        counter = counter + 1;
                    end
                end
            end
            Reliability(i,1) = counter / (trial_nr-1);  % reliability data, 1 is the 100%. 
        end
    end
end
Mao.Criteria123Passed = Mao.Criteria12Passed;
for i = 1:roi_nr
    if multi_session_idx(i) == 1

        if Reliability(i,1) < ReliabiliyIndex
            Mao.Criteria123Passed(i) = NaN;
        end
    end
end
Mao.PlaceCells =  Mao.Criteria123Passed(~isnan(Mao.Criteria123Passed));

% update potential place field matrix 
PotPlaceFieldStartBin3                                = PotPlaceFieldStartBin2;
PotPlaceFieldStartBin3(Reliability < ReliabiliyIndex) = NaN;
PotPlaceFieldLength3                                  = PotPlaceFieldLength2;
PotPlaceFieldLength3(Reliability < ReliabiliyIndex)   = NaN;
PotPlaceFieldNu3                                      = sum(~isnan(PotPlaceFieldStartBin3),2); 
Mao.PlaceFieldStartBin                                = PotPlaceFieldStartBin3;
Mao.PlaceFieldBinLength                               = PotPlaceFieldLength3;
Mao.PlaceFieldNu                                      = PotPlaceFieldNu3;
Mao.ReliabilityTreshold                               = Reliability;

% PLOT place cell ROI averaged (all trial activity)
Mao.placeROINu         = length(Mao.PlaceCells);  
Mao.placeROINormActBin = Mao.PosTuningNorm(Mao.PlaceCells,:); % collect only place cell data


if datatype == 0 
    sData.imdata.MaoPC_dff = struct;
    sData.imdata.MaoPC_dff = Mao;
    fn                     = 'dff';
    dataForSortPC          = sData.imdata.binned.MeanRoiAct(Mao.PlaceCells,:);
    dataForSort            = sData.imdata.binned.MeanRoiAct;

elseif datatype == -1 
    sData.imdata.MaoPC_dffSR = struct;
    sData.imdata.MaoPC_dffSR = Mao;
    fn                       = 'dffSR';
    dataForSortPC            = sData.imdata.binned.MeanRoiActSR(Mao.PlaceCells,:);
    dataForSort              = sData.imdata.binned.MeanRoiActSR;

elseif datatype == 1 
    sData.imdata.MaoPC_deconv = struct;
    sData.imdata.MaoPC_deconv = Mao;
    fn                        = 'deconv';
    dataForSortPC             = sData.imdata.binned.MeanRoiAct_Deconv(Mao.PlaceCells,:);
    dataForSort               = sData.imdata.binned.MeanRoiAct_Deconv;
end

% plot position tuning curves with potential place fields
% plot only place cells (in light off trials):
% for j = 1:length(Mao.PlaceCells)
%     i = Mao.PlaceCells(j);
%     f1 = figure(4); 
%     f1.Visible = 'on';
%     Xaxis = BinSize:BinSize:BinSize*bin_nr;
%     Ymax = (max(Mao.PosTuning(i,:)))+0.00000000001;
%     for k = 1:size(Mao.PlaceFieldStartBin,2)
%         if ~isnan(Mao.PlaceFieldStartBin(i,k))
%             if Mao.PlaceFieldStartBin(i,k) + Mao.PlaceFieldBinLength(i,k)-1 < bin_nr +1
%                 rectangle('Position',[Mao.PlaceFieldStartBin(i,k)*BinSize Ymax/500 (Mao.PlaceFieldBinLength(i,k)-1)*BinSize Ymax],'FaceColor',[1 0.95 1],'EdgeColor','none'); hold on;
%                 line([Mao.PlaceFieldStartBin(i,k)*BinSize Mao.PlaceFieldStartBin(i,k)*BinSize],[0 Ymax],'Color',[1 0.3 0.2],'LineStyle','--','LineWidth',1); hold on; line([(Mao.PlaceFieldStartBin(i,k)+Mao.PlaceFieldBinLength(i,k)-1)*BinSize (Mao.PlaceFieldStartBin(i,k)+Mao.PlaceFieldBinLength(i,k)-1)*BinSize],[0 Ymax],'Color',[1 0.3 0.2],'LineStyle','--','LineWidth',1); hold on;
%             else
%                 rectangle('Position',[Mao.PlaceFieldStartBin(i,k)*BinSize Ymax/500 (bin_nr-Mao.PlaceFieldStartBin(i,k))*BinSize Ymax],'FaceColor',[1 0.95 1],'EdgeColor','none'); hold on;
%                 rectangle('Position',[BinSize Ymax/500 (Mao.PlaceFieldStartBin(i,k)+Mao.PlaceFieldBinLength(i,k)-bin_nr-1)*BinSize Ymax],'FaceColor',[1 0.95 1],'EdgeColor','none'); hold on;
%                 line([Mao.PlaceFieldStartBin(i,k)*BinSize Mao.PlaceFieldStartBin(i,k)*BinSize],[0 Ymax],'Color',[1 0.3 0.2],'LineStyle','--','LineWidth',1); hold on; 
%                 line([(Mao.PlaceFieldStartBin(i,k)+Mao.PlaceFieldBinLength(i,k)-bin_nr-1)*BinSize (Mao.PlaceFieldStartBin(i,k)+Mao.PlaceFieldBinLength(i,k)-bin_nr-1)*BinSize],[0 Ymax],'Color',[1 0.3 0.2],'LineStyle','--','LineWidth',1); hold on;
%             end
%         end
%     end
%     plot(Xaxis,Mao.PosTuning(i,1:bin_nr),'LineWidth',2); hold on; % pos tuning curve
%     line([0 BinSize*bin_nr],[Treshold(i) Treshold(i)],'Color','red','LineWidth',1); hold on; 
%     axis([0 160 0 Ymax]); % ceil(Ymax)
%     title(strcat(sData.sessionInfo.sessionID,' ROI #',num2str(i),'-',fn));
%     xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
%     ylabel('Position tuning of activity');
%     fname = strcat(sData.sessionInfo.sessionID,'-roi',num2str(i),fn);
% 
% %     % Save fig
%     savefig(fullfile(SavePath,fname));
%     saveas(gcf,(fullfile(SavePath,[fname '.jpg'])));
%     close(4)
%     % close all;
% end

% sort place cell activity based on peak activity
pc_nr           = length(Mao.PlaceCells);
MaxData         = max(dataForSortPC, [], 2); % search the max in each row
MaxDataBin      = NaN(pc_nr,2);
MaxDataBin(:,2) = Mao.PlaceCells; % search Max position in bins
for i = 1:pc_nr
    if find(dataForSortPC(i,:)== MaxData(i)) % needed if the max number is twice in the dataset
        MaxDataBin(i,1) = find(dataForSortPC(i,:)== MaxData(i));   % search Max position in bins
        continue
    end
end
SortingOrder                 = sortrows(MaxDataBin, 1); % second column is the sorted ROI order. plot in these order
Mao.PlaceCellActSortingOrder = SortingOrder(:,2);
%SortedData = NaN(PCNu,BinNu);
SortedData = dataForSort(SortingOrder(:,2), :); % new matrix containing sorted data.
% Normalize
MaxData2 = max(SortedData, [], 2);
SortedDataMax = SortedData ./ MaxData2;

%% PLOT FIGURE WITH NORMALIZED DATA

x_axis = 1:sData.behavior.meta.binSize:sData.behavior.meta.binSize*sData.behavior.meta.nBins;
y_axis = 1:pc_nr;

figure('Color','white'); % smoothdata(SortedDataMax,2,'gaussian',5)

imagesc(x_axis ,y_axis, smoothdata(SortedDataMax,2,'gaussian',5)) %(1:number of bins;1:number of trials)
c = colorbar;
colormap(viridis);
c.Label.String   = fn;
c.Label.FontSize = 11;
c.TickDirection  = 'out'; 
%caxis([CADATA.VelMin 50]); %set limits for color plot, below 1st black, above 2nd white
line([sData.imdata.cues.C1A sData.imdata.cues.C1A],[0 pc_nr],'Color','black','LineStyle','--','LineWidth',2); hold on; line([sData.imdata.cues.C1B sData.imdata.cues.C1B],[0 pc_nr],'Color','black','LineStyle','--','LineWidth',2); hold on;
line([sData.imdata.cues.C2A sData.imdata.cues.C2A],[0 pc_nr],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([sData.imdata.cues.C2B sData.imdata.cues.C2B],[0 pc_nr],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
line([sData.imdata.cues.C3A sData.imdata.cues.C3A],[0 pc_nr],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([sData.imdata.cues.C3B sData.imdata.cues.C3B],[0 pc_nr],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
line([sData.imdata.cues.C4A sData.imdata.cues.C4A],[0 pc_nr],'Color','black','LineStyle','--','LineWidth',2); hold on; line([sData.imdata.cues.C4B sData.imdata.cues.C4B],[0 pc_nr],'Color','black','LineStyle','--','LineWidth',2); hold on;
xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticklabels = 0:20:160; %cannot set the axis to show 157 as the end%
Xticks = linspace(1, 160, numel(xticklabels)); set(gca, 'XTick', Xticks, 'XTickLabel', xticklabels)
ylabel('ROIs'); yticklabels = 0:10:roi_nr; yticks = linspace(1, roi_nr, numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
title(strcat(sData.sessionInfo.sessionID,'Place-cell activity sorted based on peak activity')); 

% Save fig
savefig(fullfile(SavePath,strcat(sData.sessionInfo.sessionID,'SortedMeanPlaceROI-',fn,'.fig')));
saveas(gcf,(fullfile(SavePath,strcat(sData.sessionInfo.sessionID,'SortedMeanPlaceROI-',fn,'.jpg'))));

%% PLOT FIGURE NOT NORMALIZED
figure('Color','white'); 
imagesc(x_axis , y_axis, smoothdata(SortedData,2,'gaussian',5)) %(1:number of bins;1:number of trials)
c = colorbar;
colormap(viridis);
c.Label.String   = fn;
c.Label.FontSize = 11;
c.TickDirection  = 'out';
caxis([0 0.5]);
%caxis([CADATA.VelMin 50]); %set limits for color plot, below 1st black, above 2nd white
line([sData.imdata.cues.C1A sData.imdata.cues.C1A],[0 pc_nr],'Color','black','LineStyle','--','LineWidth',2); hold on; line([sData.imdata.cues.C1B sData.imdata.cues.C1B],[0 pc_nr],'Color','black','LineStyle','--','LineWidth',2); hold on;
line([sData.imdata.cues.C2A sData.imdata.cues.C2A],[0 pc_nr],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([sData.imdata.cues.C2B sData.imdata.cues.C2B],[0 pc_nr],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
line([sData.imdata.cues.C3A sData.imdata.cues.C3A],[0 pc_nr],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([sData.imdata.cues.C3B sData.imdata.cues.C3B],[0 pc_nr],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
line([sData.imdata.cues.C4A sData.imdata.cues.C4A],[0 pc_nr],'Color','black','LineStyle','--','LineWidth',2); hold on; line([sData.imdata.cues.C4B sData.imdata.cues.C4B],[0 pc_nr],'Color','black','LineStyle','--','LineWidth',2); hold on;
xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticklabels = 0:20:160; %cannot set the axis to show 157 as the end%
Xticks = linspace(1, 160, numel(xticklabels)); set(gca, 'XTick', Xticks, 'XTickLabel', xticklabels)
ylabel('ROIs'); yticklabels = 0:10:roi_nr; yticks = linspace(1, roi_nr, numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
title(strcat(sData.sessionInfo.sessionID,'Place-cell activity sorted based on peak activity')); 

% save fig
% savefig(fullfile(SavePath,strcat(sData.sessionInfo.sessionID,'SortedMeanPlaceROINotNorm.fig')));
% saveas(gcf,(fullfile(SavePath,strcat(sData.sessionInfo.sessionID,'SortedMeanPlaceROINotNorm.jpg'))));

%% PLOT FIGURE ZSCORED
figure('Color','white'); 
imagesc(x_axis , y_axis, zscore(smoothdata(SortedData,2,'gaussian',5)')') %(1:number of bins;1:number of trials)
c = colorbar;
colormap(viridis);
c.Label.String   = strcat(fn,'(ZScored)');
c.Label.FontSize = 11;
c.TickDirection  = 'out'; 
%caxis([CADATA.VelMin 50]); %set limits for color plot, below 1st black, above 2nd white
line([sData.imdata.cues.C1A sData.imdata.cues.C1A],[0 pc_nr],'Color','black','LineStyle','--','LineWidth',2); hold on; line([sData.imdata.cues.C1B sData.imdata.cues.C1B],[0 pc_nr],'Color','black','LineStyle','--','LineWidth',2); hold on;
line([sData.imdata.cues.C2A sData.imdata.cues.C2A],[0 pc_nr],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([sData.imdata.cues.C2B sData.imdata.cues.C2B],[0 pc_nr],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
line([sData.imdata.cues.C3A sData.imdata.cues.C3A],[0 pc_nr],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([sData.imdata.cues.C3B sData.imdata.cues.C3B],[0 pc_nr],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
line([sData.imdata.cues.C4A sData.imdata.cues.C4A],[0 pc_nr],'Color','black','LineStyle','--','LineWidth',2); hold on; line([sData.imdata.cues.C4B sData.imdata.cues.C4B],[0 pc_nr],'Color','black','LineStyle','--','LineWidth',2); hold on;
xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticklabels = 0:20:160; %cannot set the axis to show 157 as the end%
Xticks = linspace(1, 160, numel(xticklabels)); set(gca, 'XTick', Xticks, 'XTickLabel', xticklabels)
ylabel('ROIs'); yticklabels = 0:10:roi_nr; yticks = linspace(1, roi_nr, numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
title(strcat(sData.sessionInfo.sessionID,' Place-cell activity sorted based on peak activity')); 

% Save fig
savefig(fullfile(SavePath,strcat(sData.sessionInfo.sessionID,'SortedMeanPlaceROINotNormZScore.fig')));
saveas(gcf,(fullfile(SavePath,strcat(sData.sessionInfo.sessionID,'SortedMeanPlaceROINotNormZScore.jpg'))));

% Save file to same path where other files can be found 
if datatype == 0 
    sData.imdata.MaoPC_dff = Mao; 
elseif datatype == -1 
    sData.imdata.MaoPC_dffSR = Mao;
elseif datatype == 1 
    sData.imdata.MaoPC_deconv = Mao;    
end

% save(fullfile(sData.sessionInfo.savePath(1:54),...
%     strcat(sData.sessionInfo.sessionID,'.mat')),'sData');

% save(fullfile(sData.sessionInfo.savePath(1:54),...
%     strcat(sData.sessionInfo.sessionID(1:14), sData.sessionInfo.sessionID(25:28),'.mat')),'sData');

end
