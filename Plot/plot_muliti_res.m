%% Data
data_awake_all = vertcat(m6134_in{2},m6137_in{2},m6139_in{2},m6140_in{2});
data_sleep_all = vertcat(m6134_in{1},m6137_in{1},m6139_in{1},m6140_in{1});


%% Variables
time           = (-(31*3):(31*3))./31;
n_rois_sleep = [1 size( rmmissing(data_sleep_all),1)];
n_rois_awake = [1 size( rmmissing(data_awake_all),1)];

data_sleep_all = rmmissing(data_sleep_all);
data_awake_all = rmmissing(data_awake_all);


frameshift    = round(.5/(1/31));
interval_mean = mean( data_sleep_all(:, ( median(1:187)-frameshift:median(1:187)+frameshift) ),2);
[max_val,~]   = max(interval_mean,[],2);
[~, sort_idx_sleep] = sortrows(max_val);

interval_mean = mean( data_awake_all(:, ( median(1:187)-frameshift:median(1:187)+frameshift) ),2);
[max_val,~]   = max(interval_mean,[],2);
[~, sort_idx_awake] = sortrows(max_val);
%% Plot 
figure, 
subplot(3,2,[1, 3])
imagesc(time, n_rois_sleep, data_sleep_all)
% colorbar
caxis([-.3 .3])
ylabel('Neuron #', 'FontSize', 16)
title('Sleep SWRs', 'FontSize',16)

subplot(3,2,5)
SE_sleep = std(data_sleep_all, 'omitnan')./ sqrt( size( data_sleep_all,1));
shadedErrorBar(time, mean(data_sleep_all, 'omitnan'), SE_sleep)
set(gca, 'xlim', [time(1) time(end)])
xline(0, 'r--', 'linew', 1) % Mark SWR peak time
ylabel('Mean DF/F (z-score)', 'FontSize', 16)

subplot(3,2,[2, 4])
imagesc(time, n_rois_awake, data_awake_all)
caxis([-.3 .3])
C = colorbar;
C.Position(1) = 0.92;
title('Awake SWRs', 'FontSize',16)

subplot(3,2,6)
SE_awake = std(data_awake_all, 'omitnan')./ sqrt( size(data_awake_all,1));
shadedErrorBar(time, mean(data_awake_all, 'omitnan'), SE_awake)
set(gca, 'xlim', [time(1) time(end)])
xline(0, 'r--', 'linew', 1) % Mark SWR peak time
xlabel('Time from SWR peak (s)', 'FontSize', 16)
