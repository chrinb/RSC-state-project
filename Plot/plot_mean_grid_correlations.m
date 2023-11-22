function plot_mean_grid_correlations(roi_list_cat, num_sessions)

% Written by Christoffer Berge | Vervaeke lab

%{
Compute mean correlations over all channel 1 grid ROIs to be included vs
mean correlations for grid ROIs to exclude. 
%}

[mean_corr_rois_include, mean_corr_rois_exclude] = deal([]);
for i = 1:num_sessions
    tmp_idx_rois_include = roi_list_cat{i, 1};
    tmp_idx_rois_exclude = roi_list_cat{i, 2};
    all_roi_corrs        = roi_list_cat{i,3};
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

figure, hold on
for i = 1:size(all_roi_corrs,1)

    if ismember(i, tmp_idx_rois_include)
        col = 'blue';
    elseif ismember(i, tmp_idx_rois_exclude)
        col = 'magenta';
    end

    plot(all_roi_corrs(i,:),'Color', col, 'linew',1)
    xlabel('Frames (downsampled x 500)', FontSize=16)
    ylabel('Corr. coef', FontSize=16)
    title('Corr. with 1st frame per ROI', FontSize=14)
end
yline(0.8, 'r--', 'LineWidth',2)


min_roi_corr_val         = min(all_roi_corrs,[],2);
min_roi_corr_val_include = min(tmp_roi_i, [],2);
min_roi_corr_val_exclude = min(tmp_roi_e, [],2);

%% Histogram minimum corr vals per ROI
bin_edges             = 0:.05:1;
[counts, bin_centers] = histcounts(min_roi_corr_val, bin_edges);
x_bins =(bin_centers(1:end-1)+bin_centers(2:end))/2;


figure;
b = bar(x_bins, counts, 1, 'FaceColor',[.5 .5 .5],'EdgeColor',[0 0 0]);


% figure,
% hist = histogram(min_roi_corr_val,20);
% histogram(min_roi_corr_val_include, 20)
% histogram(min_roi_corr_val_exclude, 20)
ylabel('N grid ROIs', FontSize=16)
xlabel('Minimum ROI corr. values', FontSize=16)
threshold_val = .8;  % Example x-positions
xline(threshold_val, 'r--','LineWidth',3)
% 
% % Define the new color for the bars at the specified x-positions
% new_color = 'r';  % Red color, you can use any valid color specification
% 
% above_threshold = hist.BinEdges > threshold_val;
% hist.FaceColor = new_color;

