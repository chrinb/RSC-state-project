function fraction = calc_in_fraction(sData, params)

% Written by Christoffer Berge | Vervaeke lab

% Compute the fraction of inhibitory neurons per session

% Load ROI array
roi_arr = sData.imdata.roi_arr;

% Find excitatory/inhibitory indicies 
[pc_rois, in_rois] = remove_cells(sData, params);

% If session is part of multi-day recordings, remove ROIs not present in
% current session
% if isfield(sData.imdata, 'roi_classification')
% 
%     roi_classification = sData.imdata.roi_classification;
% 
%     pc_roi_idx = roi_classification(pc_rois); % of all ROIs, index out cell type
%     in_roi_idx = roi_classification(in_rois); % of all ROIs, index out cell type
%     log_idx_pc = pc_roi_idx == 1;
%     log_idx_in = in_roi_idx == 1;
%     pc_rois    = pc_rois(log_idx_pc);
%     in_rois    = in_rois(log_idx_in);
% end

% Check if IN roi array is empty. This indicates an indexing/detection
% error in the code
if isempty(in_rois)
    warndlg(['No inhibitory cells found for session ', sData.sessionInfo.sessionID] );
end

% Compute ratio
fraction = numel(in_rois)/ ( numel( pc_rois) + numel(in_rois) );