function [pc_rois, in_rois] = remove_cells(varargin)

% Written by Christoffer Berge || Vervaeke lab

% Split ROIs into principal cells (PC) or inhibitory cells (IN) based on
% whether the INs in roi_arr are tagged as "overfilled" in roimanager. Also
% discards ROIs tagged as "gliopil/AG".

% opts.roi_list = x;
sData   = varargin{1,1};
roi_arr = sData.imdata.roi_arr;

% if nargin == 0 || isempty(roi_arr) 
%     folder_name   = dir;
%     temp_var      = folder_name(1).folder;
% %     temp_var2     = sessionObject.DataLocation.Subfolders
%     roiFolder     = [temp_var '\roisignals'];
%     filesInFolder = dir(roiFolder);
%     
%     for x = 1:length(filesInFolder)
%         if ~isempty(strfind(filesInFolder(x).name,'rois.mat')) == 1 % if the file is a signals.mat file
%             try 
%                 load(filesInFolder(x).name)
%             catch
%                 warning('Add "roisignals" to path')
%             end          
%         end    
%     end
% end

% Check if user has specified subset of ROIs to analyze (e.g., for multiple
% session with overlapping ROIs)
% if sum(opts.roi_list) > 1
%     roi_arr = roi_arr(opts.roi_list);
% end

% Find indicies of principal cells and inhibitory cells
if size(roi_arr,2) == 2
    roi_arr = roi_arr{1, 2};
end

% If longitudinal classification, re
% if isfield(sData.imdata, 'roi_classification')
%     total_nr_rois = 1:length(roi_arr(sData.imdata.roi_classification == 1));
%     roi_arr       = roi_arr(sData.imdata.roi_classification == 1);
% else
total_nr_rois = 1:length(roi_arr);
% end

% New sessions (analyzed after dec 2022) or old sessions reanalyzed in NANSEN have
% different inhibitory interneuron tag compared to "old" datasets. 
% if isfield(sData.imdata, 'roi_classification')
cell_idx_pc            = find( cellfun(@(c) ~strcmp('inhibitory', c), {roi_arr.celltype}, 'uni', 1) );
cell_idx_in            = find( cellfun(@(c) strcmp('inhibitory', c),  {roi_arr.celltype} , 'uni', 1) );
% else
%     cell_idx_pc            = find( cellfun(@(c) contains('Overfilled', c), {roi_arr.tags}, 'uni', 1) );
%     cell_idx_in            = find( cellfun(@(c) ~contains('Overfilled', c), {roi_arr.tags}, 'uni', 1) );
% end




% Also remove ROIs tagged by 'AG' = gliopill in roimanager. These are ROIs
% that was detected in one session and not found in subsequent sessions
% when trying to match imaging FOV. 
% missing_cells             = find( cellfun(@(c) contains('AG', c), {roi_arr.tags}, 'uni', 1) ); 
% all_cells1                = logical(total_nr_rois);
% all_cells2                = logical(total_nr_rois);
% % all_cells1(missing_cells) = false;
% % all_cells2(missing_cells) = false;
% all_cells2(cell_idx_in)   = false;
% all_cells1(cell_idx_pc)   = false;

% Discard ROIs 
try
    total_nr_rois = total_nr_rois(logical(sData.imdata.ch2_grid_classficiation) );
catch
    msg = ['No grid ROI classification available!'];
end
% 
% pc_rois = total_nr_rois(cell_idx_pc);
% in_rois = total_nr_rois(cell_idx_in);

pc_rois = intersect(total_nr_rois, cell_idx_pc);
in_rois = intersect(total_nr_rois, cell_idx_in);