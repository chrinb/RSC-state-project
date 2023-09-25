function plot_ch1_grid_meanF_std(F, roi_nr)

std_vec = std(F,0,2);

figure, 

plot(F(roi_nr,:)),
% plot(F(roi_nr,:)-mean(F(roi_nr,:))), 

% yline(std_vec(roi_nr)), 
% yline(-std_vec(roi_nr))