function [roi_classification, roi_corr] = roi_triplet_image_classification(sessionObject, params)

%{
Calculate 
%}

import nansen.twophoton.roi.getRoiAppData

% Load images
imageStack = sessionObject.loadData( params.ImageStackVariableName );

% Load rois
roiGroup = sessionObject.loadData(params.RoiArrayVariableName);

imageStack.CurrentChannel = params.ChannelNumber;
% imageStack.CurrentChannel = 1;

planeNumber = 1;

roiArray      = roiGroup(planeNumber, params.ChannelNumber).roiArray;
roiArray_2use = roiArray;

sData = sessionObject.Data.sData;

roi_image_type = 'ActivityWeightedMean';

% Select epochs in recording
params = select_roi_im_epochs(sData, params);

 % Adjust ROI image size for inhibitory cells in axon recordings
if strcmp(sessionObject.Exp_type, 'Axons') && strcmp(params.cell_type, 'in')
    [~, in_idx]         = remove_cells(sData);
    roiArray_2use            = roiArray_2use(in_idx);
    params.RoiImageSize = [40 40];
end

roiImageTriplets = deal( cell(1, numel(params.EpochFrameStart)));

for i = 1:numel(params.EpochFrameStart)
    epochStart     = params.EpochFrameStart(i);
    epochStop      = params.EpochFrameStart(i)-1 + params.EpochLength;
    imArray        = imageStack.getFrameSet(epochStart:epochStop);
    imArray        = squeeze(imArray);
    [roiImages, ~] = getRoiAppData(imArray, roiArray_2use, params); % Imported function

    roiImageTriplets{i} = roiImages;
end

% Concatenate ROI images before-during-after REM
im_before = cat(3, roiImageTriplets{1,1}(:).(roi_image_type));
im_during = cat(3, roiImageTriplets{1,2}(:).(roi_image_type));

if strcmp(params.rem, 'yes')
    im_after  = cat(3, roiImageTriplets{1,3}(:).(roi_image_type));
    all_images = cat(2, im_before, im_during, im_after);
elseif strcmp(params.rem, 'no')
    all_images = cat(2, im_before, im_during);
end

%% Sort ROIs based on correlation between image 1-2 and 1-3
n_rois = size(im_before,3);

roi_corr = zeros(n_rois, 2);
for roi_nr = 1:n_rois
    
    tmp1 = squeeze(im_before(:, :, roi_nr));
    tmp2 = squeeze(im_during(:, :, roi_nr));

    tmp_corr_val1 = corr2(tmp1, tmp2);

    roi_corr(roi_nr, 1) = tmp_corr_val1;

    if strcmp(params.rem, 'yes')
                tmp3               = squeeze(im_after(:, :, roi_nr));
                tmp_corr_val2       = corr2(tmp1, tmp3);
                roi_corr(roi_nr, 2) = tmp_corr_val2;
    end
end

%% Threshold

roi_classification = zeros(1, size(roiArray,2) );

if strcmp(sessionObject.Exp_type, 'Axons') && strcmp(params.cell_type, 'in')
    roi_classification(in_idx) = roi_corr(:,1) > 0.85;
else
    roi_classification = roi_corr(:,1) > 0.85;
end

% [sorted_roi_corr, sorted_idx] = sortrows(roi_corr, 1, 'descend');