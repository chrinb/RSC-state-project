function plot_session_activity(sData, params)

% Written by Christoffer Berge | Vervaeke lab

%{
Plot colorplots of DF/F and deconvolved, and their averages for the
session, as well as a state-activity vector (hypnogram)
%}

%% Get excitatory and inhibitory indices
[pc_rois, in_rois] = remove_cells_longitudinal(sData);

%% Get signal data

switch params.signal_type
    case 'dff'
        txt      = 'newdff';
        axon_txt = 'DffFilt';
        txt2     = '';
    case 'deconv'
        txt = 'ciaDeconvolved';
        axon_txt = 'Dec';
        txt2     = '';
    case 'transients'
        axon_txt = '_sig_transients';
        txt      = '_sig_transients';
        txt2     = params.cell_type;
end

switch params.cell_type
    case {'pc', 'in'}
        signal_to_plot{1,:} = zscore( sData.imdata.roiSignals(2).([txt2, txt])(pc_rois,:), 0, 2);
        signal_to_plot{2,:} = zscore( sData.imdata.roiSignals(2).([txt2, txt])(in_rois,:), 0 ,2);
        cmap                = [-1 2];
    case 'axon'
        signal_to_plot{1,:} = sData.imdata.roiSignals(2).(['mergedAxons',axon_txt]);
        cmap                = [0 .3];
end

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

font_size = 16;

figure, 

num_signals = size(signal_to_plot,1);

hAx(1) = subplot((num_signals*2)+2, 1, 1);
plot(time_vector, hypnogram_vector, 'LineWidth',1)
set(gca, 'xlim', [time_vector(1) time_vector(end)])
axis off
ytickvalues = 0:3;
%     x = zeros(size(ytickvalues));
str = {'QW', 'AW','NREM', 'REM'};
for i = ytickvalues
    text(-5, ytickvalues(i+1), str(i+1), 'HorizontalAlignment', 'right','FontSize',font_size);
end

for num_cell_types = 1:size(signal_to_plot,1)
    
    y1 = [1 size(signal_to_plot{num_cell_types,:},1)];

    hAx(num_cell_types*2) = subplot( (num_signals*2)+2, 1,[num_cell_types, num_cell_types+1]+num_cell_types);
    imagesc(time_vector, y1, signal_to_plot{num_cell_types,:}(cell_type_sort_coefs_idx{num_cell_types},:)),
    ylabel('Neuron #', FontSize=16)
    caxis(cmap)
%     colormap(flipud(gray)) 
    c = colorbar;
    c.Position(1) = 0.92;
    c.Position(3) = 0.01;
    
    hAx(num_cell_types*3) = subplot((num_signals*2)+2, 1, (num_signals*2)+2);
    plot(time_vector, mean(signal_to_plot{num_cell_types,:},'omitnan'))
    hold on
    xlabel('Time (s)', FontSize=16)
    ylabel('Mean DF/F', FontSize=16)
    set(gca, 'xlim', [time_vector(1) time_vector(end)])
    y_lims = get(gca, 'ylim');
end 
    
if isfield(sData, 'episodes')
    % NREM bouts
    for i = 1:length(NREM_start_end)
        x = [NREM_start_end(i,1) NREM_start_end(i,1) NREM_start_end(i,2) NREM_start_end(i,2)];
        y = [-1 2 2 -1];
        patch(x, y, 'green', 'edgecolor', 'none', 'FaceAlpha', .2);
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


set(gca, 'xlim', [time_vector(1) time_vector(end)], 'ylim', (y_lims))

linkaxes(hAx, 'x')

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
