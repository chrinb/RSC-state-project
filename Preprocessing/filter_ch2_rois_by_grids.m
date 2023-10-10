function filter_ch2_rois_by_grids(sessionObjects, roi_array, grids_to_keep_corr)

% Written by Christoffer Berge | Vervaeke lab

%{
Use list of channel 1 grid ROIs that pass image stack correlation analysis
(to check for z-drift during rec.) to sort channel 2 ROIs. The function
finds the center coordinate for each channel 2 ROI and check if its inside
one of the channel 1 grid ROIs.
%}

% roi_array = sessionObjects.loadData('RoiArray');
sData = sessionObjects.loadData('sData');

n_ch2_rois = size( roi_array(2).roiArray, 2);
n_ch1_rois = size( roi_array(1).roiArray, 2);

list_of_rois = zeros(n_ch2_rois, 1);

% Lamda func to check if ch2 ROI center coordinates is inside ch1 grid ROI
% boundaries
is_center_inside_boundaries = @(x, y, z) x(1) > y(1) && x(1) < z(1) && x(2) > y(2) && x(2) < z(2);

for roi_nr_ch2 = 1:n_ch2_rois

    ch2_roi_center = roi_array(2).roiArray(1,roi_nr_ch2).center;

    for roi_nr_ch1 = 1:n_ch1_rois

        if strcmp(roi_array(1).roiArray(1, roi_nr_ch1).tags, 'keep')
%         if grids_to_keep_idx(roi_nr_ch1) == 1
            ch1_roi_min_coordinates = min(roi_array(1).roiArray(1, roi_nr_ch1).coordinates);
            ch1_roi_max_coordinates = max(roi_array(1).roiArray(1, roi_nr_ch1).coordinates); 

            if is_center_inside_boundaries(ch2_roi_center, ch1_roi_min_coordinates, ch1_roi_max_coordinates)
%                 ch2_roi_center(1) > ch1_roi_min_coordinates(1) && ch2_roi_center(1) < ch1_roi_max_coordinates(1) && ch2_roi_center(2) > ch1_roi_min_coordinates(2) && ch2_roi_center(2) < ch1_roi_max_coordinates(2)

                list_of_rois(roi_nr_ch2, 1) = 1;
                break % Break out of loop if ch2 ROI is within ch1 grid boundaries
            end
        end
    end
end

%% Visualize remaning ROIs
roi_masks_to_keep = zeros( roi_array(1, 2).FovImageSize );
for i = 1:n_ch2_rois
    
    if list_of_rois(i) == 1
        roi_masks_to_keep = roi_masks_to_keep + roi_array(2).roiArray(1, i).mask ; 
    end
end
figure, imagesc(roi_masks_to_keep)

%% Store indices of ROIs to keep in sData
sData.imdata.ch2_grid_classficiation = list_of_rois;
sData.imdata.grids_to_keep_idx       = grids_to_keep_corr;
sessionObjects.saveData('sData', sData);
sessionObjects.Data.resetCache('sData');