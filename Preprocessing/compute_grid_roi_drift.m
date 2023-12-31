function [grids_to_keep_corr, grids_to_avoid_corr, corr_1st_frame] = compute_grid_roi_drift(sessionObjects, params)

% Written by Christoffer Berge | Vervaeke lab

%{ 
This function loops over each grid ROI and creates an image stack for a
downsampled two photon series stack for that ROI.
%}

image_stack = sessionObjects.loadData('TwoPhotonSeries_Downsampled500');
roi_array   = sessionObjects.loadData('RoiArray');

stack_dim = image_stack.Data.StackSize;

n_rois     = size(roi_array(1, 1).roiArray, 2 );
channel_nr = 1;

%% Loop over ROIS and compute image correlations across stack per ROI 
[grid_stacks, all_roi_corr]  = deal( cell(n_rois, 1));
% all_roi_corr = zeros(n_rois, stack_dim(4));

for roi_nr = 1:n_rois

    min_coordinates = min(roi_array(1, 1).roiArray(1, roi_nr).coordinates);
    max_coordinates = max(roi_array(1, 1).roiArray(1, roi_nr).coordinates); 

    length_dim1 = min_coordinates(1):max_coordinates(1);
    length_dim2 = min_coordinates(2):max_coordinates(2);

    grid_stacks{roi_nr} = squeeze(image_stack.Data(length_dim2, length_dim1  ,channel_nr,:)); 

    results = ophys.twophoton.analysis.computeDriftSummary(grid_stacks{roi_nr});

    all_roi_corr{roi_nr} = results.ImageCorrelations;
end


cat_1st_corrs = @(x) vertcat(x(1,:));
corr_1st_frame = cell2mat( cellfun(cat_1st_corrs, all_roi_corr, 'UniformOutput', false));

%% Plot correlation matrix for all ROIs first frame vs rest
% figure, 
% imagesc(corr_1st_frame)
% c = colorbar;
% caxis([0.5 1])
% c.Label.String = 'Correlation coeff.';
% c.FontSize =14;
% ylabel('ROIs', FontSize=16)
% xlabel('Frames (downsampled x 500)', FontSize=16)
%% Find rows (ROIs) where corr. coefs never dips below threshold
threshold = 0.8;
mat_log = corr_1st_frame >threshold;

row_without_zero = @(x) sum(x,2) == size(x,2);

grids_to_keep_idx = row_without_zero(mat_log);   

grid_list = 1:n_rois;

grids_to_keep_corr  = grid_list(grids_to_keep_idx);
grids_to_avoid_corr = grid_list(~grids_to_keep_idx);

% mean_corr_1st_frame = mean(corr_1st_frame(:,2:end),2);
% med_corr_1st_frame  = median(corr_1st_frame(:, 2:end),2);


% figure, 
% % sgtitle(['Threshold = ', num2str(threshold)],FontSize=16)
% subplot(131)
% imagesc(corr_1st_frame(grids_to_keep_corr,:))
% ylabel('ROIs', FontSize=16)
% colorbar
% caxis([0.7 1])
% axis square
% title('Rois to keep')
% xlabel('Frames (downsampled x 500)', FontSize=16)
% 
% subplot(132)
% imagesc(corr_1st_frame(grids_to_avoid_corr,:))
% colorbar
% caxis([0.7 1])
% title('Rois to exclude')
% xlabel('Frames (downsampled x 500)', FontSize=16)
% 
% axis square
% 
% subplot(133)
% plot( mean(corr_1st_frame(grids_to_keep_corr,:)), 'LineWidth',1), hold on
% plot( mean(corr_1st_frame(grids_to_avoid_corr,:)), 'LineWidth',1)
% legend({'Keep'; 'Exclude'})
% set(gca, 'xlim',[1 37])
% ylabel('Mean corr.', FontSize=16)
% xlabel('Frames (downsampled x 500)', FontSize=16)
% axis square

%% Tag ch1 grid ROIs to keep

if strcmp(params.tag_ch1_grid_rois, 'yes')

    % Load ROI group
    roi_group = load(sessionObjects.getDataFilePath('RoiArray'));
    plane_nr  = 1;
    
    % Covert to ROI array obj
    roi_array_ch1 = roimanager.utilities.struct2roiarray(roi_group.roiArray{plane_nr, channel_nr});
    
    % Tag ch1 grid ROIS to keep
    for i = 1:n_rois
    
        if grids_to_keep_idx(i) ==  1
            roi_array_ch1(1,i) = roi_array_ch1(1,i).addTag('keep');
        end
    end
    
    % Convert back to struct and save
    roiStruct = roimanager.utilities.roiarray2struct(roi_array_ch1);
    
    roi_group.roiArray{plane_nr, channel_nr} = roiStruct;
    
    save( sessionObjects.getDataFilePath('RoiArray'), '-struct', 'roi_group' );
    sessionObjects.Data.resetCache('RoiArray')

end
%% Filter channel 2 ROIs by channel 1 grid ROIS to keep
if strcmp(params.filter_ch2_rois, 'yes')
    filter_ch2_rois_by_grids(sessionObjects, roi_array, grids_to_keep_corr)
end

%% PLot correlation for first vs remaining frames for a single ROI
% roi = 43;
% figure, 
% plot(corr_1st_frame(roi,2:end), 'LineWidth',1), yline(med_corr_1st_frame(roi))
% set(gca, 'ylim',[.6 1])
% xlabel('Frames (downsampled x 500)', FontSize=16)
% ylabel('Corr. coefs', FontSize=16)
% title(['ROI ', num2str(roi)],FontSize=14)

%% Plot population mean F
% F = sessionObjects.loadData('RoiSignals_MeanF');
% F = table2array(F)';
% 
% zf = zscore(F, 0, 2);
% zftrend = detrend(zf);
% 
% figure, 
% 
% subplot(131)
% imagesc(zf(rois_to_keep_corr,:))
% colorbar
% caxis([-2 2])
% title('Rois to keep mean z-score F')
% xlabel('Frames', FontSize=16)
% ylabel('ROIs', Fontsize=16)
% axis square
% 
% subplot(132)
% imagesc(zf(rois_to_avoid_corr,:))
% colorbar
% caxis([-2 2])
% title('Rois to exclude mean z-score F')
% xlabel('Frames', FontSize=16)
% axis square
% 
% subplot(133)
% plot(mean(zf(rois_to_avoid_corr, :)), 'LineWidth',.5), hold on
% plot(mean(zf(rois_to_keep_corr, :)), 'LineWidth',.5)
% legend({'Exclude'; 'Keep'})
% xlabel('Frames', FontSize=16)
% ylabel('Mean z-score fluorescence', Fontsize=16)
% set(gca, 'ylim', [-3 3 ])
% axis square
