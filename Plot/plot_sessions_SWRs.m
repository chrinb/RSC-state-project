function plot_sessions_SWRs(sData)

%{
Plot SWR occurrence in session together with state scoring and running
speed.
%}

swr_idx = sData.ephysdata.absRipIdx; 

lfp = sData.ephysdata.lfp;

run_speed = smoothdata(downsample( sData.daqdata.runSpeed, 10), 'gaussian', 250);
% run_speed = downsample(sData.daqdata.runSpeed, 10);

swr_idx_vec          = nan(1, length(lfp));
swr_idx_vec(swr_idx) = 1;

% Binned (1 sec) SWR rate
bin_width_sec         = 1;  
bin_edges             = 0:bin_width_sec:(length(lfp)/2500);
[counts, bin_centers] = histcounts(swr_idx/2500, bin_edges);
x_bins =(bin_centers(1:end-1)+bin_centers(2:end))/2;

%%  Create hypnogram vector
hypnogram_vector = zeros(1, length(lfp));

hypnogram_vector(sData.behavior.quiet_wakefulness == 1)  = 0;
hypnogram_vector(sData.behavior.active_wakefulness == 1) = 1;
hypnogram_vector(sData.behavior.NREM_vector == 1)        = 2;
hypnogram_vector(sData.behavior.REM_vector  == 1)        = 3;
%% Plot
time_vec = (0:length(lfp)-1)./2500;
time_vec_ds = downsample(time_vec, 10);
lfp_ds = downsample(lfp, 10);
% font_size = 16;
y_ax_lim = [-4 max(counts)];


% Downscale running speed if necessary to fit plot
if max(run_speed) > 1
    run_speed = run_speed*0.1;
end
run_speed_scale = [-4 -3.5];

fig = figure; 

hold on
b = bar(x_bins, counts, 1, 'FaceColor',[.5 .5 .5],'EdgeColor',[0 0 0]);
b.ShowBaseLine = "off";
plot([-10 -10], [0 max(counts)], 'k', 'LineWidth',3)

plot(time_vec_ds, run_speed-4, 'k', 'LineWidth',1)
plot([-10 -10], run_speed_scale, 'k', 'LineWidth',3)

plot(time_vec, swr_idx_vec-1.5, 'D', 'Color','b')
plot([-10 -10], [-2.5 -1.5], 'k', 'LineWidth',3)

plot(time_vec,lfp-2, 'k', 'LineWidth',.5 )

% plot(time_vec, hypnogram_vector+3, 'k', 'LineWidth',1 )
% ytickvalues = 3:6;
% x = zeros(size(ytickvalues));
% str = {'QW','AW','NREM', 'REM'};
% 
% for i = 1:4
%     text(-15, ytickvalues(i), str(i), 'HorizontalAlignment', 'right','FontSize',font_size);
% end

% Plot states
if isfield(sData, 'episodes')
    
    REM_episodes = rem_sleep(sData);
    NREM_episodes = nrem_sleep(sData);

    NREM_start_end = NREM_episodes./2500;
    REM_start_end = REM_episodes./2500;

    for i = 1:length(NREM_start_end)
        x = [NREM_start_end(i,1) NREM_start_end(i,1) NREM_start_end(i,2) NREM_start_end(i,2)];
        y = [y_ax_lim flip(y_ax_lim)];
        patch(x, y, 'green', 'edgecolor', 'none', 'FaceAlpha', .2);
    end
    % REM bouts
    if size(REM_start_end(:,1),1) == 1
        a = [REM_start_end(1,1) REM_start_end(1,1) REM_start_end(1,2) REM_start_end(1,2)];
        b = [y_ax_lim flip(y_ax_lim)];
        patch(a, b, 'red', 'edgecolor', 'none', 'FaceAlpha', .2);
    else
        for i = 1:length(REM_start_end)
            a = [REM_start_end(i,1) REM_start_end(i,1) REM_start_end(i,2) REM_start_end(i,2)];
            b = [y_ax_lim flip(y_ax_lim)];
            patch(a, b, 'red', 'edgecolor', 'none', 'FaceAlpha', .2);
        end
    end
end

if ~sum(sData.behavior.active_wakefulness) == 0
    [AW_start, AW_stop]    = findTransitions(sData.behavior.active_wakefulness);
    AW_start_stop = [AW_start', AW_stop']./2500;

    for i = 1:length(AW_start_stop)
        x = [AW_start_stop(i,1) AW_start_stop(i,1) AW_start_stop(i,2) AW_start_stop(i,2)];
        y = [y_ax_lim flip(y_ax_lim)];
        patch(x, y, 'm', 'edgecolor', 'none', 'FaceAlpha', .2);
    end
end

fig.Children.YAxis.Visible = 'off';

set(gca, 'XLim', [-10 time_vec(end)], 'ylim', [-4 max(counts)])
xlabel('Time (s)', FontSize=16)
