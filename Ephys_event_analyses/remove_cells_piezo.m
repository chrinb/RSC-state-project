function [plane_idx_pc, plane_idx_in, plane_idx_axon, all_idx_pc, all_idx_in, all_idx_axon] = remove_cells_piezo(sData)

% Written by Christoffer Berge || Vervaeke lab

% Similar to "remove_cells" function, this function finds the indices of 
% principal cells, inhibitory cells and axons in multi-plane recordings.

% Loop over planes 
roi_arr   = cell(1,4);
all_rois = [];
for plane_nr = 1:4

    % Get ROI arrays for each plane and put in cell
    roi_arr{:, plane_nr} = sData.imdata.roi_arr(plane_nr, 1).roiArray;
    
    % In this separate variable, concatenate all ROI arrays
    all_rois = [all_rois, roi_arr{1, plane_nr}];
end


% Find the true indicies of different cell types
all_idx_pc   = find( cellfun(@(c) strcmp('', c), {all_rois.celltype}, 'uni', 1) );
all_idx_in   = find( cellfun(@(c) strcmp('inhibitory', c),  {all_rois.celltype} , 'uni', 1) );
all_idx_axon = find( cellfun(@(c) strcmp('axon', c),  {all_rois.celltype} , 'uni', 1) );


% Loop over ROI arrays and find cell-type indices
plane_idx_pc   = cell(1,4);
plane_idx_in   = cell(1,4);
plane_idx_axon = cell(1,4);

for plane_nr = 1:4
    
    % Get ROI array from current plane
    temp_roi_arr = roi_arr{1, plane_nr}  ;
    
    % Get cell-type indices from current plane
    plane_idx_pc{:, plane_nr}   = find( cellfun(@(c) strcmp('', c), {temp_roi_arr.celltype}, 'uni', 1) );
    plane_idx_in{:, plane_nr}   = find( cellfun(@(c) strcmp('inhibitory', c),  {temp_roi_arr.celltype} , 'uni', 1) );
    plane_idx_axon{:, plane_nr} = find( cellfun(@(c) strcmp('axon', c),  {temp_roi_arr.celltype} , 'uni', 1) );
end


