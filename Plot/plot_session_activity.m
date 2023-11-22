function plot_session_activity(sData, params)

% Written by Christoffer Berge | Vervaeke lab

%{
Plot colorplots of DF/F and deconvolved, and their averages for the
session, as well as a state-activity vector (hypnogram)
%}

[signal_to_plot, cmap] = get_roi_signals_from_sData(sData, params);

state_vectors_2p = get_state_logicals(sData);


%% Create hypnogram vector

% Initialize empty vector
hypnogram_vector = zeros(1, size(signal_to_plot{1,:},2) );

% Create hypnogram vector
hypnogram_vector(state_vectors_2p{1,1} == 1) = 2; % NREM sleep
hypnogram_vector(state_vectors_2p{1,2} == 1) = 3; % REM sleep
hypnogram_vector(state_vectors_2p{1,3} == 1) = 0; % Quiet wakefulness
hypnogram_vector(state_vectors_2p{1,4} == 1) = 1; % Active wakefulness

%% Compute cell-by-cell correlation with a particular state

% Determine which state to sort by
switch params.beh_state{:}
    case 'NREM'
        idx = 1;
        txt2 = 'NREM';
    case 'REM'
        idx = 2;
        txt2 = 'REM';
    case 'QW'
        idx = 4;
        txt2 = 'QW';
    case 'AW'
        idx = 3;
        txt2 = 'AW';
end

% Loop over cell-type matrices
for num_cell_types = 1:size(signal_to_plot,1)

%     signal       = zscore(signal_to_plot, 0, 2);
    signal       = signal_to_plot{num_cell_types,:};
    n_rois       = size(signal,1);
    corr_coefs   = zeros(n_rois, 4);
    corr_coefs_p = zeros(n_rois, 4);
    
    for roi_nr = 1:size(signal)
    
        roi_data = signal(roi_nr,:);
    %     roi_data = roi_data./max(roi_data);
    
        for state_nr = 1:4
            current_state_vec = state_vectors_2p{1, state_nr};
            
            if ~isempty(current_state_vec)
                [temp_corr, temp_corr_p]       = corrcoef(current_state_vec, roi_data);
                corr_coefs(roi_nr, state_nr )  = temp_corr(2,1);
                corr_coefs_p(roi_nr, state_nr) = temp_corr_p(2,1);
            end
        end
    end
    [cell_type_sort_coefs{num_cell_types, 1}, cell_type_sort_coefs_idx_all{num_cell_types, 1}] = sort(corr_coefs,1, 'descend');
    cell_type_sort_coefs_idx{num_cell_types}    = cell_type_sort_coefs_idx_all{num_cell_types, 1}(:,idx);
end

% [sorted_coefs, sort_coef_indx_all] = sort(corr_coefs,1, 'descend');
% sorted_coefs_p_all                 = corr_coefs_p(sort_coef_indx_all);

% coefs_to_analyze = sorted_coefs(:,idx);
% sort_coef_indx   = sort_coef_indx_all(:, idx);
% sorted_coefs_p   = sorted_coefs_p_all(:, idx);
%% Plot

% Get 2P frame rate
imaging_sampling_rate = find_imaging_framerate(sData);

time_vector  = linspace(0, length(roi_data), length(roi_data))/imaging_sampling_rate;

if isfield(sData, 'episodes')
    REM_episodes = rem_sleep(sData);
    NREM_episodes = nrem_sleep(sData);

    NREM_start_end = NREM_episodes./2500;
    REM_start_end = REM_episodes./2500;
end

signal_to_plot{1,1} = okada(signal_to_plot{1,1}, 2);
signal_to_plot{2,1} = okada(signal_to_plot{2,1}, 2);
% test1 =signal_to_plot{1,1};
% % test2 =signal_to_plot{2,1};
% 
% signal_to_plot{1,1} = test1 ;
% % % signal_to_plot{2,1} = test2;
% signal_to_plot{1,1} = smoothdata(signal_to_plot{1,1}, 2,'gaussian', 20);
% signal_to_plot{2,1} = smoothdata(signal_to_plot{2,1}, 2, 'gaussian', 20);

%% Plot 

figure, 
% cmap = [0 2]
font_size = 16;

num_signals = size(signal_to_plot,1);

n_plots     =  (num_signals*num_signals)+3;
mean_colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980]};

hAx(1) = subplot(n_plots, 1, 1);
plot(time_vector, hypnogram_vector, 'k','LineWidth',1)
axis off

ytickvalues = 0:3;
str = {'QW', 'AW','NREM', 'REM'};
for i = ytickvalues
    text(-5, ytickvalues(i+1), str(i+1), 'HorizontalAlignment', 'right','FontSize',font_size);
end

for num_cell_types = 1:size(signal_to_plot,1)
    
    signal = signal_to_plot{num_cell_types,:};

    y1 = [1 size(signal,1)];

    hAx(num_cell_types*2) = subplot(n_plots, 1, [2 3]+num_cell_types*num_cell_types-1);
    imagesc(time_vector, y1, signal(cell_type_sort_coefs_idx{num_cell_types},:)),
    ylabel('Neuron #', FontSize=12)
    clim(cmap)
%     colormap(flipud(gray)) 
    c = colorbar;
    ylabel(c, 'z-score DF/F', Fontsize=10)
    c.Position(1) = 0.92;
    c.Position(3) = 0.01;
    hAx(num_cell_types*2).XAxis.Visible = 'off';
    % if num_cell_types == 1 && size(signal,1) > 100
    %     set(gca, 'ylim',[1 100])
    % end

    hAx(num_cell_types*3) = subplot(n_plots, 1, num_cell_types*3+1);
    plot(time_vector, mean(signal,'omitnan'), 'Color',mean_colors{num_cell_types})
    hold on
    y_lims = get(gca, 'ylim');
    hAx(num_cell_types*3).YAxis.Visible = 'off';
    
    if num_signals == 1 || num_cell_types == 2
        hAx(num_cell_types*3).XAxis.Visible = 'on';
        hAx(num_cell_types*3).XAxis.FontSize=10;
    else
         hAx(num_cell_types*3).XAxis.Visible = 'off';
    end

    txt = sprintf('    Mean\nz-score DF/F');
    text(-50, 0, 0, txt)

    % Plot episodes
    if isfield(sData, 'episodes')
        % NREM bouts
        for i = 1:length(NREM_start_end)
            x = [NREM_start_end(i,1) NREM_start_end(i,1) NREM_start_end(i,2) NREM_start_end(i,2)];
            y = [-1 2 2 -1];
            patch(x, y, 'blue', 'edgecolor', 'none', 'FaceAlpha', .2);
        end
        % REM bouts
        if size(REM_start_end(:,1),1) == 1
            a = [REM_start_end(1,1) REM_start_end(1,1) REM_start_end(1,2) REM_start_end(1,2)];
            b = [-1 2 2 -1];
            patch(a, b, 'red', 'edgecolor', 'none', 'FaceAlpha', .2);
        else
            for i = 1:length(REM_start_end)
                a = [REM_start_end(i,1) REM_start_end(i,1) REM_start_end(i,2) REM_start_end(i,2)];
                b = [-1 2 2 -1];
                patch(a, b, 'red', 'edgecolor', 'none', 'FaceAlpha', .2);
            end
        end
    end
    set(gca, 'YLim', y_lims)
    plot([time_vector(end)+5 time_vector(end)+5], [-.5 .5],'k', 'LineWidth',3)
    
end
xlabel('Time (s)', FontSize=16)
linkaxes(hAx, 'x')
set(gca, 'xlim', [time_vector(1) time_vector(end)+5])

%% Histogram 
% figure,
% histogram(coefs_to_analyze, 20)
% set(gca, 'xlim', [-1 1])
% xline(0, 'r--', 'LineWidth',1)
% ylabel('Counts')
% xlabel('Correlation coefficients')
% title([txt, ' - ', txt2])
% axis square


%% Plot cells with highest positive and negative correlations
% if isfield(sData, 'episodes')
% 
%     figure,
%     sgtitle(['3 neurons with highest/lowest correlation: ', state_indicies{idx}])
%     
%     list_neuron = [sort_coef_indx(1), sort_coef_indx(2), sort_coef_indx(3), sort_coef_indx(end), sort_coef_indx(end-1), sort_coef_indx(end-2)];
%     list_corr   = [1, 2, 3, length(sort_coef_indx) length(sort_coef_indx)-1, length(sort_coef_indx)-2];
%     
%     list_plot = [1 3 5 2 4 6];
%     for rois_to_plot = 1:6
%     
%         hAx(rois_to_plot) = subplot(3, 2, list_plot(rois_to_plot));
%         nr      = list_neuron(rois_to_plot);
%         corr_nr = list_corr(rois_to_plot); 
%         plot(time_vector, zscore( signal(nr,:))), hold on
%         title(['ROI: ', num2str(nr), ' Corr. coef. ' num2str(coefs_to_analyze(corr_nr))] )
%         % y_ax_lim = get(gca, 'ylim');
%         y_ax_lim = [-5 6];
%         % NREM bouts
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
%     linkaxes(hAx, 'x');
% end
