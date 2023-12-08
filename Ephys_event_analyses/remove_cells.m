function [pc_rois, in_rois] = remove_cells(sData, params)

% Written by Christoffer Berge || Vervaeke lab

% Split ROIs into principal cells (PC) or inhibitory cells (IN) based on
% whether the INs in roi_arr are tagged as "overfilled" in roimanager. Also
% discards ROIs tagged as "gliopil/AG".

% opts.roi_list = x;
% sData   = varargin{1,1};
% params  = varargin{1,2};

roi_arr = sData.imdata.roi_arr;


% Find indicies of principal cells and inhibitory cells
if size(roi_arr,2) == 2
    roi_arr = roi_arr{1, 2};
end
   
cell_idx_pc = find( cellfun(@(c) ~strcmp('inhibitory', c), {roi_arr.celltype}, 'uni', 1) );
cell_idx_in = find( cellfun(@(c) strcmp('inhibitory', c),  {roi_arr.celltype} , 'uni', 1) );

% Remove channel 2 ROIs within channel 1 grid ROIs that fail correlation
% criteria
total_nr_rois = 1:length(roi_arr);

switch params.use_roi_classification
    case 'grid'
          total_nr_rois = total_nr_rois(logical(sData.imdata.ch2_grid_classficiation) );
    case 'ch2_across_session'
        total_nr_rois = total_nr_rois(logical(sData.imdata.ch2_acrossSession_classification) );
end


pc_rois = intersect(total_nr_rois, cell_idx_pc);
in_rois = intersect(total_nr_rois, cell_idx_in);