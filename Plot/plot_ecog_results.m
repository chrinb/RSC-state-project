function plot_ecog_results(data1, data2, data3)

% Written by Christoffer Berge || Vervaeke lab

% Plot colormaps and mean/median SWR-aligned ECoG data. 

% Compute standard error of the mean
SE_data1 = std(data1, 'omitnan')./ sqrt(size(data1,1));
SE_data2 = std(data2, 'omitnan')./ sqrt(size(data2,1));
SE_data3 = std(data3, 'omitnan')./ sqrt(size(data3,1));

%% Determine max-min values of Y limits in average plots

max_data1 = max( mean(data1, 'omitnan'));
max_data2       = max( mean(data2, 'omitnan'));
max_data3      = max( mean(data3, 'omitnan'));

min_data1 = min( mean(data1, 'omitnan'));
min_data2 = min( mean(data2, 'omitnan'));
min_data3 = min( mean(data3, 'omitnan'));

Y_lim_max = max([max_data1, max_data2,max_data3]) + 0.01;
Y_lim_min = min([min_data1, min_data2,min_data3]) - 0.01;

%% Plot result
y_name   = 'SWR #';
x_label  = 'Time from SWR peak (sec)'; 
cLim     = [-.4 .4];
time_vec = (-1:1/2500:1);

% Set up X and Y limits of imagesc function
x1 = [time_vec(1) time_vec(end)];
data1_y = [1, size(data1,1)]; 
data2_y = [1, size(data2,1)]; 
data3_y = [1, size(data3,1)]; 

% Plot data
figure,
subplot(3, 3, [1 4])
imagesc(x1, data1_y, data1)
ylabel(y_name)
colorbar
caxis(cLim)
title(['SWR n = ', num2str(size(data1,1)) ])


subplot(3,3,7)
shadedErrorBar(time_vec, mean(data1, 'omitnan'),SE_data1,'lineprops', 'r');
xlabel(x_label)
set(gca, 'ylim', [Y_lim_min, Y_lim_max])
legend('Unclassified')

subplot(3, 3, [2 5])
imagesc(x1, data2_y, data2)
ylabel(y_name)
colorbar
caxis(cLim)
title(['SWR n = ', num2str(size(data2,1)) ])

subplot(3,3,8)
shadedErrorBar(time_vec, mean(data2, 'omitnan'),SE_data2,'lineprops', 'b');
xlabel(x_label)
set(gca, 'ylim', [Y_lim_min, Y_lim_max])
legend('Spindle-uncoupled')

subplot(3, 3, [3 6])
imagesc(x1, data3_y, data3)
ylabel(y_name)
colorbar
caxis(cLim)
title(['SWR n = ', num2str(size(data3,1)) ])

subplot(3,3,9)
shadedErrorBar(time_vec, mean(data3, 'omitnan'),SE_data3,'lineprops', 'k'), hold on
% plot(time_vec, median(data1, 1, 'omitnan'), 'linew',1)
xlabel(x_label)
set(gca, 'ylim', [Y_lim_min, Y_lim_max])
legend('Spindle-coupled')

