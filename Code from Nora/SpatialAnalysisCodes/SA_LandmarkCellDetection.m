function sData = SA_LandmarkCellDetection(sData)

%%% CHECK FOR LANDMARK CELLS: 
% there are two set of landmarks on the wheel at certain distance: 
% velcro at 26 and 86 cm and hot glue spikes at 46 and 66 cm. All tactile cues span for 2 cms.
% The distance between HGS are 20 cm and velcros are 60 cm. I look for double peaks with these distances. 
% Since in the dF/F data the second peak is often masked by the first peak's dacay, I check these double peaks in the dconvolved (unsmoothed) data

bin_size = sData.behavior.meta.binSize; 

% 6 cm , search for landmark cells where the landmark was presented and a little addition
LandmarkArea = floor(8/sData.behavior.meta.binSize); 

VelcroStart1 = floor(sData.imdata.cues.C1A  /bin_size);
VelcroEnd1 = VelcroStart1 + LandmarkArea;

VelcroStart2 = floor(sData.imdata.cues.C4A  /bin_size);
VelcroEnd2 = VelcroStart2 + LandmarkArea;

HGSStart1 = floor(sData.imdata.cues.C2A  /bin_size);
HGSEnd1 = HGSStart1 + LandmarkArea;

HGSStart2 = floor(sData.imdata.cues.C3A  /bin_size);
HGSEnd2 = HGSStart2 + LandmarkArea; % delay in Ca transient due to previous Ca transient evoked by first HGS


% DiffVelcro = floor(60/sData.behavior.meta.binSize); %30,55 bin
% DiffHGS = floor(20/sData.behavior.meta.binSize); %10,2 bin
DiffVelcro = VelcroStart2- VelcroStart1; %30,55 bin
DiffHGS = HGSStart2 - HGSStart1; %10,2 bin

AllowedJitter = 2; % bins

PCdff = sData.imdata.MaoPC_dff; % struct containing which cells are place cells based on dF/F data
% datasetDeconv = sData.imdata.binned.RoiDeconvolved;
MinPeakDetection = 0.0001; % spike detection treshold in the deconvolved dataset, 0.0001

LandmarkCells = struct;
LandmarkCells.LandmarkCells = NaN(sData.imdata.nROIs,1); 
LandmarkCells.VelcroCells = NaN(sData.imdata.nROIs,1);
LandmarkCells.HGSCells = NaN(sData.imdata.nROIs,1);
LandmarkCells.Peak1PosBin = NaN(sData.imdata.nROIs,1);
LandmarkCells.Peak2PosBin = NaN(sData.imdata.nROIs,1);

% Quick fix until deconv. data is computed (CNB, 06.02.23)
binned_data = sData.imdata.binned.RoidFF;
n_bins      = sData.behavior.meta.nBins;
bin_vector  = 1:n_bins;

% Area +- cue position to search for activity (in cm)
cue_zone = 8;

% Cue positions (in cm)
x_soft1 = median( [sData.imdata.cues.C1A,  sData.imdata.cues.C1B]) ;
x_soft2 = median( [sData.imdata.cues.C4A,  sData.imdata.cues.C4B]) ;
x_hard1 = median( [sData.imdata.cues.C2A,  sData.imdata.cues.C2B]) ;
x_hard2 = median( [sData.imdata.cues.C3A,  sData.imdata.cues.C3B]) ;

% Area around each cue to search for activity
cue_soft1_bins = round( [x_soft1-cue_zone, x_soft1+cue_zone] /bin_size);
cue_soft2_bins = round( [x_soft2-cue_zone, x_soft2+cue_zone]/bin_size);
cue_hard1_bins = round( [x_hard1-cue_zone, x_hard1+cue_zone]/bin_size);
cue_hard2_bins = round ( [x_hard2-cue_zone, x_hard2+cue_zone]/bin_size);

all_cue_zones = [cue_soft1_bins; cue_soft2_bins; cue_hard1_bins; cue_hard2_bins];


%%
n_rois = sData.imdata.nROIs;

cue_cells = zeros(n_rois, 1);
for roi_n = 1:n_rois % comparing peak positions (for a ROI) if the distance of peaks equals to distance between landmarks 
    % if PCdff.Criteria123Passed(i) > 0 % if the cell was a place cell

        activity_map = binned_data{1,roi_n};
        activity_map_mean = mean( activity_map, 'omitnan');

        % check if it can be a velcro or a HGS cell, then check if there is a peak at the other similar landmark (one bin jitter allowed for Velcro, 2 bin for HGS)
        if PeakAmplitudeLocation >= VelcroStart1 &&  PeakAmplitudeLocation <= VelcroEnd1 % Velcro1
            DiffVelcroArray = PeakAmplitudeLocation + DiffVelcro - AllowedJitter : PeakAmplitudeLocation + DiffVelcro + AllowedJitter;
            DiffVelcroArray(DiffVelcroArray < VelcroStart2) = NaN;
            DiffVelcroArray(DiffVelcroArray > VelcroEnd2) = NaN;
            if  any(ismember(locs,DiffVelcroArray)>0)
                LandmarkCells.LandmarkCells(roi_n) = roi_n; 
                LandmarkCells.VelcroCells(roi_n) = roi_n;
                LandmarkCells.Peak1PosBin(roi_n) = PeakAmplitudeLocation;
                LandmarkCells.Peak2PosBin(roi_n) = min(locs(ismember(locs,DiffVelcroArray)>0));
            end
        elseif PeakAmplitudeLocation >= VelcroStart2 &&  PeakAmplitudeLocation <= VelcroEnd2 %Velcro2
            DiffVelcroArray = PeakAmplitudeLocation - DiffVelcro - AllowedJitter : PeakAmplitudeLocation - DiffVelcro + AllowedJitter;
            DiffVelcroArray(DiffVelcroArray < VelcroStart1) = NaN;
            DiffVelcroArray(DiffVelcroArray > VelcroEnd1) = NaN;
            if any(ismember(locs,DiffVelcroArray)>0) 
                LandmarkCells.LandmarkCells(roi_n) = roi_n; 
                LandmarkCells.VelcroCells(roi_n) = roi_n;
                LandmarkCells.Peak1PosBin(roi_n) = PeakAmplitudeLocation;
                LandmarkCells.Peak2PosBin(roi_n) = min(locs(ismember(locs,DiffVelcroArray)>0));
            end
        elseif PeakAmplitudeLocation >= HGSStart1 &&  PeakAmplitudeLocation <= HGSEnd1
            DiffHGSArray = PeakAmplitudeLocation + DiffHGS - AllowedJitter : PeakAmplitudeLocation + DiffHGS + AllowedJitter + 1; % an extra bin is added for jitter because of the Ca transient of the previous HGS decays here and mask new transient
            DiffHGSArray(DiffHGSArray < HGSStart2) = NaN;
            DiffHGSArray(DiffHGSArray > HGSEnd2) = NaN;
            if any(ismember(locs,DiffHGSArray)>0) 
                LandmarkCells.LandmarkCells(roi_n) = roi_n; 
                LandmarkCells.HGSCells(roi_n) = roi_n;
                LandmarkCells.Peak1PosBin(roi_n) = PeakAmplitudeLocation;
                LandmarkCells.Peak2PosBin(roi_n) = min(locs(ismember(locs,DiffHGSArray)>0));
            end
        elseif PeakAmplitudeLocation >= HGSStart2 &&  PeakAmplitudeLocation <= HGSEnd2
            DiffHGSArray = PeakAmplitudeLocation - DiffHGS - AllowedJitter -1 : PeakAmplitudeLocation - DiffHGS + AllowedJitter; % an extra bin is added for jitter because of the Ca transient of the previous HGS decays here and mask new transient
            DiffHGSArray(DiffHGSArray < HGSStart1) = NaN;
            DiffHGSArray(DiffHGSArray > HGSEnd1) = NaN;
            if any(ismember(locs,DiffHGSArray)>0)
                LandmarkCells.LandmarkCells(roi_n) = roi_n; 
                LandmarkCells.HGSCells(roi_n) = roi_n;
                LandmarkCells.Peak1PosBin(roi_n) = PeakAmplitudeLocation;
                LandmarkCells.Peak2PosBin(roi_n) = min(locs(ismember(locs,DiffHGSArray)>0));
            end
        end
    end
end      

LandmarkCells.Params.LandmarkArea = LandmarkArea; 
LandmarkCells.Params.VelcroStart1 = VelcroStart1;
LandmarkCells.Params.VelcroEnd1 = VelcroEnd1;
LandmarkCells.Params.VelcroStart2 = VelcroStart2;
LandmarkCells.Params.VelcroEnd2 = VelcroEnd2;
LandmarkCells.Params.HGSStart1 = HGSStart1;
LandmarkCells.Params.HGSEnd1 = HGSEnd1;
LandmarkCells.Params.HGSStart2 = HGSStart2;
LandmarkCells.Params.HGSEnd2 = HGSEnd2;
LandmarkCells.Params.DiffVelcro = DiffVelcro;
LandmarkCells.Params.DiffHGS = DiffHGS;
LandmarkCells.Params.AllowedJitter = AllowedJitter; 
LandmarkCells.Params.MinPeakDetection = MinPeakDetection;
LandmarkCells.LandmarkCellsList = LandmarkCells.LandmarkCells(~isnan(LandmarkCells.LandmarkCells)); 
LandmarkCells.VelcroCellsList = LandmarkCells.VelcroCells(~isnan(LandmarkCells.VelcroCells)); 
LandmarkCells.HGSCellsList = LandmarkCells.HGSCells(~isnan(LandmarkCells.HGSCells));

sData.imdata.LandmarkCells = LandmarkCells;

%% Detect place cells with only one place field
datasetDFF = sData.imdata.binned.RoidFF; % use dFF data

PlaceCellsWithOnePF = struct;
PlaceCellsWithOnePF.PlaceCellsWithOnePF = NaN(sData.imdata.nROIs,1);
PlaceCellsWithOnePF.PeakPosdFFBin = NaN(sData.imdata.nROIs,1);
PlaceCellsWithOnePF.PeakPosDeconvBin = NaN(sData.imdata.nROIs,1);
for roi_n = 1:1:sData.imdata.nROIs % search for place cells with one peak 
    if sData.imdata.MaoPC_dff.Criteria123Passed(roi_n) > 0 && isnan(sData.imdata.LandmarkCells.LandmarkCells(roi_n))  % if the cell was a place cell but not a landmark cell
        %DataTemp = smoothdata(nanmean(dataset{1,i},1),'Gaussian',3);  % use little bit of smoothing on deconvolved data, because sometimes larger transients cause double peaks in deconvolution
        %[pks,locs] = findpeaks(DataTemp,'MinPeakProminence',MinPeakDetection); % peaks and indices at which the peaks occur.
        DataTemp2 = [nanmean(datasetDFF{1,roi_n},1) nanmean(datasetDFF{1,roi_n},1)];  % use little bit of smoothing on deconvolved data, because sometimes larger transients cause double peaks in deconvolution
        DataTempMax = max(DataTemp2);
        MinPeakDetection = DataTempMax/5;
        [pks,locs] = findpeaks(DataTemp2,'MinPeakProminence',MinPeakDetection); % peaks and indices at which the peaks occur.
        PeakAmplitude = max(pks);
        PeakAmplitudeLocation = locs(find(pks == max(pks),1)); % search for the peak amplitude's location (which bin)
        %if length(PeakAmplitudeLocation) < 1
        %    continue
        %end
        pksMod = pks;
        pksMod(pksMod == PeakAmplitude) = 0;
        if any((pksMod)> PeakAmplitude/3) % if there is another big peak in the dFF data (1/3 of biggest), it is not considered as a single peak place cell
            continue
        else
            PlaceCellsWithOnePF.PlaceCellsWithOnePF(roi_n) = roi_n;
            PlaceCellsWithOnePF.PeakPosdFFBin(roi_n) = PeakAmplitudeLocation;
            if PeakAmplitudeLocation > sData.behavior.meta.nBins
                PlaceCellsWithOnePF.PeakPosdFFBin(roi_n) = PeakAmplitudeLocation - sData.behavior.meta.nBins;
            end
            % find the peak location in the deconvolved signal
            DataTemp3 = nanmean(datasetDeconv{1,roi_n},1);
            if sum(DataTemp3)==0
                continue
            end
            % DataTempMax = max(DataTemp3);
            % MinPeakDetection = DataTempMax/4;
            % [~,locs] = findpeaks(DataTemp3,'MinPeakProminence',MinPeakDetection); % peaks and indices at which the peaks occur.
            % locs2 = locs(locs <= PeakAmplitudeLocation);
            % if isempty(locs2)
            %    locs2 = locs;
            %    PeakAmplitudeLocation = PeakAmplitudeLocation + sData.behavior.meta.nBins; % if the dff peak is in the beginning of trial and deconv peak is at the end of previous trial
            % end
            % locs3 = PeakAmplitudeLocation - locs2;
            % PlaceCellsWithOnePF.PeakPosDeconvBin(i) = locs2(locs3 == min(locs3));
        end
    end
end
PlaceCellsWithOnePF.PlaceCellsWithOnePFList = PlaceCellsWithOnePF.PlaceCellsWithOnePF(~isnan(PlaceCellsWithOnePF.PlaceCellsWithOnePF));

sData.imdata.PlaceCellsWithOnePF = PlaceCellsWithOnePF;

% saving
% save(fullfile(sData.sessionInfo.savePath(1:54),...
%     strcat(sData.sessionInfo.sessionID(1:14), sData.sessionInfo.sessionID(25:28),'.mat')),'sData');

end