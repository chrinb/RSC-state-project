function [roi_mean_unclassified, roi_mean_single, roi_mean_coupled] = join_arr_sleep(varargin)

% Written by Christoffer Berge | Vervaeke Lab

% User inputs output cell from "swr_sleep_an" function for two(?) different 
% sessions with overlapping/matching ROIs. The code below finds the cell
% containing the ROI x SWR x time PETH for each cell. 

% Note: nr of ROIs have to match across the cell arrays

for i = 1:size(varargin,2)
    cell_arr{1,i} = varargin{1,i};
end

roi_mean_unclassified = zeros( size(cell_arr{1}{1,7},1), size(cell_arr{1}{1,9},3));
roi_mean_single       = zeros( size(cell_arr{1}{1,7},1), size(cell_arr{1}{1,9},3));
roi_mean_coupled      =  zeros( size(cell_arr{1}{1,7},1), size(cell_arr{1}{1,9},3));

for roinr = 1:size(cell_arr{1}{1,7},1)

    combined_roi_dat_unclassified = [];
    combined_roi_dat_single       = [];
    combined_roi_dat_coupled      = [];

    for session_nr = 1:size(cell_arr,2)
           combined_roi_dat_unclassified = [combined_roi_dat_unclassified; squeeze(cell_arr{1,session_nr}{1,9}(roinr,:,:))];
           combined_roi_dat_single       = [combined_roi_dat_single; squeeze(cell_arr{1,session_nr}{1,10}(roinr,:,:))]; 
           combined_roi_dat_coupled      = [combined_roi_dat_coupled; squeeze(cell_arr{1,session_nr}{1,11}(roinr,:,:))];
    end

    roi_mean_unclassified(roinr,:) = mean(combined_roi_dat_unclassified, 'omitnan');
    roi_mean_single(roinr,:)       = mean(combined_roi_dat_single, 'omitnan');
    roi_mean_coupled(roinr,:)      = mean(combined_roi_dat_coupled, 'omitnan');

end







