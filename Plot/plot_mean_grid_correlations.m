function plot_mean_grid_correlations(roi_list_cat)

% Written by Christoffer Berge | Vervaeke lab

%{
Compute mean correlations over all channel 1 grid ROIs to be included vs
mean correlations for grid ROIs to exclude. 
%}

[mean_corr_rois_include, mean_corr_rois_exclude] = deal([]);
for i = 1:num_sessions
    tmp_idx_rois_include = roi_list_cat{i, 1};
    tmp_idx_rois_exclude = roi_list_cat{i, 2};

    tmp_roi_i            = roi_list_cat{i,3}(tmp_idx_rois_include, :);
    tmp_roi_e            = roi_list_cat{i,3}(tmp_idx_rois_exclude, :);
    mean_corr_rois_include = [mean_corr_rois_include; tmp_roi_i];
    mean_corr_rois_exclude = [mean_corr_rois_exclude; tmp_roi_e];
end

figure, 
plot(mean(mean_corr_rois_include), 'linew',1), hold on
plot(mean(mean_corr_rois_exclude), 'linew',1)
xlabel('Frames (downsampled x 500)', FontSize=16)
ylabel('Corr. coef', FontSize=16)
legend({'Include'; 'Exclude'})
title('Mean corr.', FontSize=14)