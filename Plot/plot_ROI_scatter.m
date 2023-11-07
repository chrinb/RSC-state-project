function plot_ROI_scatter(all_ratios_across_sessions, sData, params )

%{
Scatterplot showing the ratios of single ROI mean activity for different
states vs. quiet wakefulness (which serves as a reference state).
%}

ratios_rel_to_QW = all_ratios_across_sessions;

%%  Plot ratios
jitter_factor = 0.03;
n_rois        = size(ratios_rel_to_QW, 1);

figure, hold on
x_coordinate = 1;
xlabel_tick  = [];

% Loop over AW, NREM and REM
% [sort_ratios, sort_ratios_idx] = deal( zeros( size(ratios_rel_to_QW(:, 1), 1), 3 ));
[sort_ratios, sort_ratios_idx] = deal( cell( 1,3 ));

state_list = [1 3 4];

for state_nr = 1:3
    
    temp_ratios = ratios_rel_to_QW(:, state_nr);

    % Delete zeros, NaNs, and inf
    temp_ratios(temp_ratios == 0)       = [];
    temp_ratios(isnan(temp_ratios))     = [];
    temp_ratios(~isfinite(temp_ratios)) = [];

    idx = state_list(state_nr);
    % Check if state contains ratio values
    if ~sum(temp_ratios) == 0  

        rand_vec_for_plot = (randn(n_rois,1)*jitter_factor) + x_coordinate;
    
        [sort_ratios{:, state_nr}, sort_ratios_idx{:, state_nr} ] = sort(temp_ratios);
    
        below_one_idx = temp_ratios < 1;
        above_one_idx = temp_ratios > 1;
    
        scatter( rand_vec_for_plot(below_one_idx), temp_ratios(below_one_idx), 'filled')
        scatter( rand_vec_for_plot(above_one_idx), temp_ratios(above_one_idx), 'filled')
       
        x_coordinate_2 = x_coordinate + 0.17;

        plot( x_coordinate_2, mean( temp_ratios, 'omitnan'), 'ko', 'MarkerFaceColor','k', 'LineWidth',4)
        errorbar(x_coordinate_2, mean( temp_ratios ), std( temp_ratios), 'LineWidth', 3, 'color', 'k')
    end
    xlabel_tick = [xlabel_tick, x_coordinate];

    x_coordinate = x_coordinate + 0.4;
end
yline(1, 'r--', 'LineWidth',2)
plot_middle = mean( [1, x_coordinate_2]);
plot_edgeL  = plot_middle-0.65;
plot_edgeR  = plot_middle+0.65;
set(gca, 'xlim', [plot_edgeL plot_edgeR])

set(gca, 'xtick', xlabel_tick);
labels = {'AW/QW', 'NREM/QW', 'REM/QW'};
xticklabels(labels)
ylabel('Ratio', FontSize=16)
title(['ROI state preference, ', params.cell_type, ', ', params.signal_type], Interpreter="none")
ax = gca;
ax.XAxis.FontSize = 16;

%% Plot cells with highest & lowest ratios
switch params.beh_state
    case 'AW'        
        ratios = ratios_rel_to_QW(:, 1);
        idx    = sort_ratios_idx{:, 1};
    case 'NREM'
        ratios =  ratios_rel_to_QW(:, 2);
        idx   =  sort_ratios_idx{:, 2};
    case 'REM'
        ratios =  ratios_rel_to_QW(:, 3);
        idx   =  sort_ratios_idx{:, 3};
end

state_vectors_2p = get_state_logicals(sData);

[signal_to_plot, ~,  pc_rois, in_rois] = get_roi_signals_from_sData(sData, params);

if strcmp(params.cell_type, 'pc')
    signal_to_plot = signal_to_plot{1,:};
    roi_nr_in_entire_dataset = pc_rois(idx);

elseif strcmp(params.cell_type, 'in')
    signal_to_plot = signal_to_plot{2,:};
    roi_nr_in_entire_dataset = in_rois(idx);

elseif strcmp(params.cell_type, 'axon')
    signal_to_plot = signal_to_plot{1,1};
    roi_nr_in_entire_dataset = idx;
end


if isfield(sData, 'episodes')

    REM_episodes = rem_sleep(sData);
    NREM_episodes = nrem_sleep(sData);

    NREM_start_end = NREM_episodes./2500;
    REM_start_end = REM_episodes./2500;
end

imaging_sampling_rate = find_imaging_framerate(sData);
time_vector           = linspace(0, size(signal_to_plot,2), size(signal_to_plot,2))/imaging_sampling_rate;

figure;
sgtitle(['Neurons with lowest/highest ', params.beh_state, '/QW ratios'])

list_neuron  = [idx(1:3),  idx(end-2:end)];  
list_neuron2 = [roi_nr_in_entire_dataset(1:3), roi_nr_in_entire_dataset(end-2:end)];
sub_factor  = 0;
running_speed_scale_factor = 0.3;
for rois_to_plot = 1:7

    plot_factor = 60 - sub_factor;
    
    if rois_to_plot < 7
        if rois_to_plot > 3
            col = [0.3010 0.7450 0.9330];
        else
            col = [0.6350 0.0780 0.1840];
        end
        roi_nr            = list_neuron(rois_to_plot);
        roi_nr_in_dataset = list_neuron2(rois_to_plot);
        if strcmp(params.zscore, 'yes')
            signal =  okada( zscore( signal_to_plot(roi_nr,:)), 2 );
        else
            signal =  okada( signal_to_plot(roi_nr,:), 2 );
        end

        plot(time_vector, signal-plot_factor, 'Color',col), hold on
        txt = sprintf(['ROI ', num2str(roi_nr_in_dataset), newline, 'ratio: ',  num2str(ratios(roi_nr))] );
        text(-80, -plot_factor, 0,txt, FontSize=12)

        % Scale bar 10 z-score DF/F
        if rois_to_plot == 6
            plot([0 0],[-plot_factor -plot_factor+10], 'LineWidth', 3, 'Color',[0 0 0])
        end
    elseif rois_to_plot == 7
        run_speed_ds       = downsample(sData.daqdata.runSpeed, 10);
        time_vector_ephys  = linspace(0, size(run_speed_ds,1), size(run_speed_ds,1))/250;
        plot(time_vector_ephys, (smoothdata(run_speed_ds, 'gaussian', 100))*running_speed_scale_factor-plot_factor, 'Color','k');
        txt = sprintf('Running\n speed');
        text(-85, -plot_factor, 0, txt, FontSize=12)

        % Scale bar = 10 cm/s
        scale_bar_ylim = [0 10]*running_speed_scale_factor;
        scale_bar_ylim = scale_bar_ylim-plot_factor;
        plot([0 0], scale_bar_ylim, 'LineWidth', 3, 'Color',[0 0 0])
    end

    set(gca,'ytick',[])
    sub_factor = sub_factor-10;
end
xlabel('Time (s)', FontSize=12)

% Overlay state bouts as colored patches
y_ax_lim = [-50 -120];
if isfield(sData, 'episodes')

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

if ~sum(ratios_rel_to_QW(:, 1)) == 0  
    [AW_start, AW_stop]    = findTransitions(state_vectors_2p{1,4});
    AW_start_stop = [AW_start', AW_stop']./imaging_sampling_rate;

    for i = 1:length(AW_start_stop)
        x = [AW_start_stop(i,1) AW_start_stop(i,1) AW_start_stop(i,2) AW_start_stop(i,2)];
        y = [y_ax_lim flip(y_ax_lim)];
        patch(x, y, 'm', 'edgecolor', 'none', 'FaceAlpha', .2);
    end
end

set(gca, 'xlim', [time_vector_ephys(1), time_vector_ephys(end)])