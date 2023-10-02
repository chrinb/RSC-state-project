function [pc_rois, in_rois] = remove_cells_longitudinal(sData)

% Written by Christoffer Berge | Vervaeke lab

% Check if recording is multi-day and if so, exclude ROIs not present in
% current session. Returns "corrected" ROI indices for excitatory and
% inhibitory cells

%% Get exctitatory and inhibitory indices (THIS NEEDS ITS OWN FUNCTION!!)
% Find excitatory/inhibitory indicies 
[pc_rois, in_rois] = remove_cells(sData);

% If session is part of multi-day recordings, remove ROIs not present in
% current session
long_roi_classification = [];

if isfield(sData.imdata, 'roi_classification')
    long_roi_classification = sData.imdata.roi_classification;
end
    % Add check for grid ROI classification if exists
%     try grid_roi_classification = sData.imdata.ch2_grid_classficiation;
%         tmp = [long_roi_classification, grid_roi_classification];
%         long_roi_classification = sum(tmp,2) ==2;
%     catch
%     end
grid_roi_classification =[];

if isfield(sData.imdata, 'ch2_grid_classficiation')
    grid_roi_classification = sData.imdata.ch2_grid_classficiation;
end
        tmp = [long_roi_classification, grid_roi_classification];
        long_roi_classification = sum(tmp,2) ==2;

    % Out of ALL multi-day ROIs (roi_classification), index the different
    % cell types to find which of them are principal or inhibitory cells 
    pc_roi_idx = long_roi_classification(pc_rois); 
    in_roi_idx = long_roi_classification(in_rois); 

    % Next, determine which of these ROIs are present in current session
    log_idx_pc = pc_roi_idx == 1;
    log_idx_in = in_roi_idx == 1;

    % Use index of specific cell type present in current session to select  
    % out of all principal or inhibitory cells
    pc_rois    = pc_rois(log_idx_pc);
    in_rois    = in_rois(log_idx_in);
end
