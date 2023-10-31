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
% switch params.beh_state
%     case 'AW'        
%         ratios = ratios_rel_to_QW(:, 1);
%         idx    = sort_ratios_idx{:, 1};
%     case 'NREM'
%         ratios =  ratios_rel_to_QW(:, 2);
%         idx   =  sort_ratios_idx{:, 2};
%     case 'REM'
%         ratios =  ratios_rel_to_QW(:, 3);
%         idx   =  sort_ratios_idx{:, 3};
% end
% 
% state_vectors_2p = get_state_logicals(sData);
% 
% [signal_to_plot, ~, ~, ~] = get_roi_signals_from_sData(sData, params);
% 
% if strcmp(params.cell_type, 'pc')
%     signal_to_plot = signal_to_plot{1,:};
% elseif strcmp(params.cell_type, 'in')
%     signal_to_plot = signal_to_plot{2,:};
% elseif strcmp(params.cell_type, 'axon')
%     signal_to_plot = signal_to_plot{1,:};
% end
% 
% 
% if isfield(sData, 'episodes')
% 
%     REM_episodes = rem_sleep(sData);
%     NREM_episodes = nrem_sleep(sData);
% 
%     NREM_start_end = NREM_episodes./2500;
%     REM_start_end = REM_episodes./2500;
% end
% 
% imaging_sampling_rate = find_imaging_framerate(sData);
% time_vector  = linspace(0, size(signal_to_plot,2), size(signal_to_plot,2))/imaging_sampling_rate;
% 
% figure,
% sgtitle('3 neurons with lowest/highest state ratios')
% 
% list_neuron = [idx(1:3),  idx(end-2:end)];    
% list_plot  = [1 3 5 2 4 6];
% hAx = zeros(1,6);
% for rois_to_plot = 1:6
% 
%     hAx(rois_to_plot) = subplot(3, 2, list_plot(rois_to_plot));
%     nr                = list_neuron(rois_to_plot);
%     plot(time_vector, zscore( signal_to_plot(nr,:))), hold on
%     title([params.cell_type, ' ROI: ', num2str(nr), ', ratio: ',  num2str(ratios(nr))] )
%     y_ax_lim = [-5 10];
%     
%     % NREM bouts
%     if isfield(sData, 'episodes')
% 
%         for i = 1:length(NREM_start_end)
%             x = [NREM_start_end(i,1) NREM_start_end(i,1) NREM_start_end(i,2) NREM_start_end(i,2)];
%             y = [y_ax_lim flip(y_ax_lim)];
%             patch(x, y, 'green', 'edgecolor', 'none', 'FaceAlpha', .2);
%         end
%         % REM bouts
%         if size(REM_start_end(:,1),1) == 1
%             a = [REM_start_end(1,1) REM_start_end(1,1) REM_start_end(1,2) REM_start_end(1,2)];
%             b = [y_ax_lim flip(y_ax_lim)];
%             patch(a, b, 'red', 'edgecolor', 'none', 'FaceAlpha', .2);
%         else
%             for i = 1:length(REM_start_end)
%                 a = [REM_start_end(i,1) REM_start_end(i,1) REM_start_end(i,2) REM_start_end(i,2)];
%                 b = [y_ax_lim flip(y_ax_lim)];
%                 patch(a, b, 'red', 'edgecolor', 'none', 'FaceAlpha', .2);
%             end
%         end
%         set(gca, 'xlim', [time_vector(1) time_vector(end)], 'ylim', [y_ax_lim])
%     end
%         
% 
%     if ~sum(ratios_rel_to_QW(:, 1)) == 0  
%         [AW_start, AW_stop]    = findTransitions(state_vectors_2p{1,4});
%         AW_start_stop = [AW_start', AW_stop']./imaging_sampling_rate;
% 
%         for i = 1:length(AW_start_stop)
%             x = [AW_start_stop(i,1) AW_start_stop(i,1) AW_start_stop(i,2) AW_start_stop(i,2)];
%             y = [y_ax_lim flip(y_ax_lim)];
%             patch(x, y, 'm', 'edgecolor', 'none', 'FaceAlpha', .2);
%         end
%     end
% end
% 
% linkaxes(hAx, 'x');
