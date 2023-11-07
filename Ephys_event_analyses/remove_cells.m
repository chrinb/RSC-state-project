function [pc_rois, in_rois] = remove_cells(varargin)

% Written by Christoffer Berge || Vervaeke lab

% Split ROIs into principal cells (PC) or inhibitory cells (IN) based on
% whether the INs in roi_arr are tagged as "overfilled" in roimanager. Also
% discards ROIs tagged as "gliopil/AG".

% opts.roi_list = x;
sData   = varargin{1,1};
roi_arr = sData.imdata.roi_arr;


% Find indicies of principal cells and inhibitory cells
if size(roi_arr,2) == 2
    roi_arr = roi_arr{1, 2};
end

total_nr_rois = 1:length(roi_arr);

cell_idx_pc            = find( cellfun(@(c) ~strcmp('inhibitory', c), {roi_arr.celltype}, 'uni', 1) );
cell_idx_in            = find( cellfun(@(c) strcmp('inhibitory', c),  {roi_arr.celltype} , 'uni', 1) );

try
    total_nr_rois = total_nr_rois(logical(sData.imdata.ch2_grid_classficiation) );
catch
    msg = 'No grid ROI classification available!';
end

pc_rois = intersect(total_nr_rois, cell_idx_pc);
in_rois = intersect(total_nr_rois, cell_idx_in);