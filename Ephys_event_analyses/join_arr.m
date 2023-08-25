function roi_mean = join_arr(varargin)

% Written by Christoffer Berge | Vervaeke Lab

% User inputs output cell from "swr_awake_an" function for two(?) different 
% sessions with overlapping/matching ROIs. The code below finds the cell
% containing the ROI x SWR x time PETH for each cell. 

% Note: nr of ROIs have to match across the cell arrays

for i = size(varargin)
    cell_arr{i} = varargin{1,i};
end

roi_mean = zeros( size(cell_arr{1}{1,8},1), size(cell_arr{1}{1,8},3));

for roinr = 1:size(cell_arr{1}{1,8},1)
    combined_roi_data = [];

    for session_nr = 1:size(cell_arr,2)
           combined_roi_data = [combined_roi_data; squeeze(cell_arr{1,session_nr}{1,8}(roinr,:,:))];
    end

%     roi_data1 = squeeze(cell_arr{1,8}(roinr,:,:));
%     roi_data2 = squeeze(cell_arr{1,8}(roinr,:,:));

    roi_mean(roinr,:) = mean(combined_roi_data, 'omitnan');
    end


end





