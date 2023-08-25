function plot_roi_swr_timestamps(varargin)

% Written by Christoffer Berge || Vervaeke lab

% Plot DF/F & Deconvolved DF/F for a single roi and timestamps of detected
% SWRs for close-up inspection of interesting ROIs.

sData    = varargin{1,1};
roinr    = varargin{1,2};

% If a third input argument called extract wanted cell type
if nargin > 2
    % Find indicies of principal cells
    [pc_rois, in_rois] = remove_cells(sData.imdata.roi_arr);
    if contains(varargin{1,3},'axon')
%             ROImtx                     = sData.imdata.roiSignals(2).newdff(pc_rois,:);
%             [dff,deconv, roiClustIDs, n_unchanged_rois] = hierClust_axons(ROImtx,1, sData, pc_rois);
        dff    = sData.imdata.roiSignals(2).mergedAxonsDff;
        deconv = sData.imdata.roiSignals(2).mergedAxonsDec; 
    elseif contains(varargin{1,3},'pc')
        dff    = sData.imdata.roiSignals(2).newdff(pc_rois,:);
        deconv = sData.imdata.roiSignals(2).ciaDeconvolved(pc_rois,:);
    elseif contains(varargin{1,3},'in')
        dff    = sData.imdata.roiSignals(2).newdff(in_rois,:);
        deconv = sData.imdata.roiSignals(2).ciaDeconvolved(in_rois,:);

    end
else
    dff   = sData.imdata.roiSignals(2).newdff;
    deconv = sData.imdata.roiSignals(2).ciaDeconvolved;
end

dff_roi    = dff(roinr,:);
deconv_roi = deconv(roinr,:);

dff_max    = max(dff_roi);
deconv_max = max(deconv_roi);


% Get SWR timestamps
swr_timestamps1 = zeros(1, length(sData.imdata.roiSignals(2).newdff));
swr_timestamps1(sData.ephysdata.frameRipIdx) = dff_max;

swr_timestamps2 = zeros(1, length(sData.imdata.roiSignals(2).newdff));
swr_timestamps2(sData.ephysdata.frameRipIdx) = deconv_max;

% Create time vector for plotting
time_vec = (0:length(sData.imdata.roiSignals(2).newdff)-1) / 31;

RippleIdx = length(sData.ephysdata.absRipIdx);
sessionID = sData.sessionInfo.sessionID;

%% Plot
figure,

sgtitle([ sessionID, ', SWR n = ' num2str(RippleIdx)], 'Interpreter', 'none')

fAx1 = subplot(211);
plot(time_vec, dff_roi), hold on
plot(time_vec, swr_timestamps1);
set(gca,'xlim',[0, time_vec(end)])
ylabel('DF/F')

fAx2 = subplot(212);
plot(time_vec, deconv_roi), hold on
plot(time_vec, swr_timestamps2);
set(gca,'xlim',[0, time_vec(end)])
xlabel('Time (sec)')
ylabel('Deconvolved')

linkaxes([fAx1, fAx2], 'x');