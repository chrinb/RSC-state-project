function roi_list = match_roi_arr(roi_arr1, roi_arr2)

% Written by Christoffer Berge || Vervaeke lab

% User inputs to roi arrays: one from the "original" session containing all
% ROIs, and one from the other session were ROIs matching those of the
% first session have been kept. To find the corresponding ROI indicies
% after deleting non-overlapping ROIs this function loops over ROI unique
% identifiers (of the smallest ROI array) and searches for the
% corresponding ROI uid in the larger array. That way one can match the
% ROIs across sessions. 

roi_list = [];
for i = 1:length(roi_arr1)
    roi_id = roi_arr1(1,i).uid;
        for j = 1:length(roi_arr2)
            roi_id_to_match = roi_arr2(1,j).uid;
 
            if strcmp(roi_id, roi_id_to_match)
                roi_list = [roi_list, j];
                break % exit inner for loop
            end
        end
end
