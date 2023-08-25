function plot_swr_cleaning(varargin)

% Written by Christoffer Berge | Vervaeke Lab


% Function to plot SWRs after removing movement-related, close-SWRs or both
% to see how well the "SWR cleaning" using "riprun2", "RemoveRip" or both 
% function works.

sData = varargin{1,1};

RippleIdx  = sData.ephysdata.absRipIdx;
color_plot = [0 0.4470 0.7410];

[RippleIdx2,~, loc_vec] = remove_movement_swr(sData, RippleIdx);
color_plot2 = 'r';

[RippleIdx3]   = remove_clustered_swr(RippleIdx,[]);
color_plot3 = 'k';

[RippleIdx41]   = remove_clustered_swr(RippleIdx,[]);
[RippleIdx42,~] = remove_movement_swr(sData, RippleIdx41);
color_plot4 = 'm';

%% Find SWR onset/offset
[swr_start_stop,~, ~] = mark_ripple_onset_offset(sData);

% Find time stamps and duration of SWR clusters
swr_idx = ismember(sData.ephysdata.absRipIdx,RippleIdx);
swr_idx2 = ismember(sData.ephysdata.absRipIdx,RippleIdx2);
swr_idx3 = ismember(sData.ephysdata.absRipIdx,RippleIdx3);
swr_idx4 = ismember(sData.ephysdata.absRipIdx,RippleIdx42);

swr_dur   = swr_start_stop(swr_idx, :);
swr_dur2   = swr_start_stop(swr_idx2, :);
swr_dur3   = swr_start_stop(swr_idx3, :);
swr_dur4   = swr_start_stop(swr_idx4, :);

%% Plot

low_patchV   = .2;
high_patchV  = .25;
low_patchV2  = .25;
high_patchV2 = .3;
low_patchV3  = .3;
high_patchV3 = .35;
low_patchV4  = .35;
high_patchV4 = .40;

ylim_low   = mean(sData.ephysdata.lfp)-1;
ylim_high  = mean(sData.ephysdata.lfp)+1; 

figure,
ax(1) = subplot(211);
time_vec = (0:length(sData.ephysdata.lfp)-1)./2500;
plot(time_vec, sData.ephysdata.lfp), hold on

% Plot all SWRs
for nr_swr = 1:size(swr_dur,1)
    x = [ swr_dur(nr_swr,1) swr_dur(nr_swr,1) swr_dur(nr_swr,2) swr_dur(nr_swr,2)]/2500;
    y = [low_patchV high_patchV high_patchV low_patchV];
    patch(x, y, color_plot, 'edgecolor', 'none', 'FaceAlpha', .7,'LineWidth',2);
end



% Plot SWRs - movement-related SWRs
for nr_swr = 1:size(swr_dur2,1)
    x = [ swr_dur2(nr_swr,1) swr_dur2(nr_swr,1) swr_dur2(nr_swr,2) swr_dur2(nr_swr,2)]/2500;
    y = [low_patchV2 high_patchV2 high_patchV2 low_patchV2];
    patch(x, y, color_plot2, 'edgecolor', 'none', 'FaceAlpha', .7,'LineWidth',2);
end




% Plot SWR - close SWRs
for nr_swr = 1:size(swr_dur3,1)
    x = [ swr_dur3(nr_swr,1) swr_dur3(nr_swr,1) swr_dur3(nr_swr,2) swr_dur3(nr_swr,2)]/2500;
    y = [low_patchV3 high_patchV3 high_patchV3 low_patchV3];
    patch(x, y, color_plot3, 'edgecolor', 'none', 'FaceAlpha', .7,'LineWidth',2);
end


% Plot SWR - movement & close SWRs
for nr_swr = 1:size(swr_dur4,1)
    x = [ swr_dur4(nr_swr,1) swr_dur4(nr_swr,1) swr_dur4(nr_swr,2) swr_dur4(nr_swr,2)]/2500;
    y = [low_patchV4 high_patchV4 high_patchV4 low_patchV4];
    patch(x, y, color_plot4, 'edgecolor', 'none', 'FaceAlpha', .7,'LineWidth',2);
end

% Adjust Y and X limits
set(gca,'ylim',[ylim_low ylim_high], 'xlim',[0 time_vec(end)]);
labels = {'All SWR', 'locom-SWR', 'clust-SWR', 'lomoc&clust-SWR'};
legend(labels)

ax(2) = subplot(212);
plot(time_vec,sData.daqdata.runSpeed), hold on,
plot(time_vec, loc_vec)
set(gca, 'xlim',[0 time_vec(end)]);
xlabel('Time (sec)');
ylabel('Running speed (cm/sec)')
linkaxes(ax, 'x');
