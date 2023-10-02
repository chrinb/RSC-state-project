function plot_grid_heatmap(num_sessions, roi_list_cat)  

% Written by Christoffer Berge | Vervaeke lab

%{
Function that plots heatmap showing which channel 1 grid ROIs to keep after 
z-drift correlation analysis. If a single session was provided as input, the 
plot simply shows which grid ROIs were excluded. If data from multiple sessions 
are provided the plot shows proportion of times a ROI was excluded across sessions
%}

grid_list_e = [];
for i = 1:num_sessions
    tmp_dat_e =  roi_list_cat{i,2};
    tmp_dat_i = roi_list_cat{i,1};
    grid_list_e = [grid_list_e, tmp_dat_e];
end

c = unique(grid_list_e);

 for i = 1:length(c)
   counts(i,1) = sum(grid_list_e==c(i)); % number of times each unique value is repeated
 end

 lol_e      = [c', counts];
 lol_e(:,2) = lol_e(:,2)./num_sessions;

 [n_rois_vec1, n_rois_vec2] = deal(1:64);
 n_rois_vec1(lol_e(:,1)) = [];
 n_rois_vec2(n_rois_vec1) = 0;
 n_rois_vec2(lol_e(:,1)) = lol_e(:,2);

 grid_mat = reshape(n_rois_vec2, 8, 8);
 grid_mat = grid_mat';

figure,
imagesc(grid_mat)
xlabel('Column ROIs', FontSize=16)
ylabel('Row ROIs', FontSize=16)
colorbar
title('Proportion of excluded ROIs in FOV across sessions', FontSize=12)