function cell_type_layer_activity = get_mean_activity_piezo(sData, params)

%{
Extract mean activity for different cell types across different layers and
compare the signals.
%}

state_vectors_2p = get_state_logicals(sData);

if isfield(sData, 'episodes')
    % rems = rem_sleep(sData);
    % nrems = nrem_sleep(sData);
 
    try
        [rem_start, rem_stop ] = findTransitions( state_vectors_2p{2});
        REM_start_end  = [rem_start; rem_stop]';   
    catch
        REM_start_end = zeros(0, 2); 
    end

    try 
        [nrem_start, nrem_stop ] = findTransitions( state_vectors_2p{1});
        NREM_start_end = [nrem_start; nrem_stop]';
    catch
        NREM_start_end = zeros(0, 2);
    end
end

% Get 2P frame rate
imaging_sampling_rate = find_imaging_framerate(sData);

% Divide by 4 to get FPS per plane
imaging_sampling_rate = imaging_sampling_rate/4;

%% Find indices of excitatory & inhibitory cells, and axons
[~, ~, ~, all_idx_pc, all_idx_in, ~] = remove_cells_piezo(sData);

% Loop over cell types

cell_type_layer_activity = cell(1,3);

for i = 1:3

    if i == 1
        params.cell_type = 'axon';
    elseif i == 2
        params.cell_type = 'in';
    elseif i == 3
        params.cell_type = 'pc';
    end

    
    total_frames_per_plane = ceil( sData.daqdata.frame_onset_reference_frame(end)/4 );
    
    % It seems like some nr of frames (10-300) are recorded in one plane before
    % piezo start. These "extra" frames are not present in the
    % frame_onset_reference_frame variable. Therefore use that variable to
    % remove extra frames.
    frames_in_rec = size(sData.imdata.roiSignals(2).newdff, 2);
    dff_diff      = frames_in_rec - total_frames_per_plane;
    
    for current_plane = 1:4
    
        % Get DF/F from selected cell type in current plane
        switch params.cell_type
            case 'in'
                % Find the logical indices of the ROIs for plane nr = i
                roi_plane_idx = sData.imdata.plane_indices(all_idx_in) == current_plane;
                % Extract the DF/F for cell type for current plane (and remove
                % piezo artefacts)
                temp_dff      = sData.imdata.roiSignals(2).newdff( all_idx_in(roi_plane_idx), dff_diff+1:end);
                clim = [0 2];
            case 'axon'
                temp_dff = sData.imdata.roiSignals(2).mergedAxonsDffFilt( sData.imdata.plane_indices_merged == current_plane, dff_diff+1:end);
                clim = [0 .5];
            case 'pc'
                roi_plane_idx = sData.imdata.plane_indices(all_idx_pc) == current_plane;
                temp_dff      = sData.imdata.roiSignals(2).newdff( all_idx_pc(roi_plane_idx), dff_diff+1:end);
                clim = [0 2];
        end
            
        % Get nr of ROIs from current plane
        if strcmp(params.zscore, 'yes')
            temp_dff = zscore(temp_dff, 0 ,2);
        end
        
        % Store cell type activity per layer
        layer_activity{current_plane} = temp_dff;
        clear temp_dff
    end
    cell_type_layer_activity{1, i} = layer_activity;
    clear layer_activity
end

all_axon_activity = cell_type_layer_activity{1,1};
all_in_activity   = cell_type_layer_activity{1,2};
all_pc_activity   = cell_type_layer_activity{1,3};

n_axon_rois = cellfun(@(x) size(x,1), all_axon_activity );
n_in_rois = cellfun(@(x) size(x,1), all_in_activity );
n_pc_rois = cellfun(@(x) size(x,1), all_pc_activity );

all_n_rois = horzcat(n_axon_rois, n_in_rois, n_pc_rois);

time_vec = linspace(0, total_frames_per_plane, total_frames_per_plane)/imaging_sampling_rate;

% anon_fun = @(x) okada(x,2);

if strcmp(params.filter, 'yes')
    all_in_activity = cellfun(@(x) okada(x,2), all_in_activity, 'uni', false);
    all_pc_activity = cellfun( @(x) okada(x,2), all_pc_activity, 'uni', false);
end

% Loop over planes and calculate mean
[mean_axon_activity, mean_in_activity, mean_pc_activity] = deal( zeros(4, frames_in_rec-dff_diff));
for i = 1:4
    
    % Axons
    if size(all_axon_activity{1,i},1) > 1
        mean_axon_activity(i,:) = mean( all_axon_activity{1,i} );

    elseif size(all_axon_activity{1,i},1) == 1
         mean_axon_activity(i,:) = all_axon_activity{1,i};

    elseif size(all_axon_activity{1,i},1) == 0
        mean_axon_activity(i,:) = NaN(1, size(time_vec,2));
    end
    
    % Inhibitory cells
    if size(all_in_activity{1,i},1) > 1
        mean_in_activity(i,:) = mean( all_in_activity{1,i} );

    elseif size(all_in_activity{1,i},1) == 1
         mean_in_activity(i,:) = all_in_activity{1,i};

    elseif size(all_in_activity{1,i},1) == 0
        mean_in_activity(i,:) = NaN(1, size(time_vec,2));
    end
    
    % Excitatory cells
    if size(all_pc_activity{1,i},1) > 1
        mean_pc_activity(i,:) = mean( all_pc_activity{1,i} );

    elseif size(all_pc_activity{1,i},1) == 1
         mean_pc_activity(i,:) = all_pc_activity{1,i};

    elseif size(all_pc_activity{1,i},1) == 0
        mean_pc_activity(i,:) = NaN(1, size(time_vec,2));
    end
end

all_signals = vertcat( mean_axon_activity, mean_in_activity, mean_pc_activity );


%% Get states
if isfield(sData, 'episodes')
    REM_episodes = rem_sleep(sData);
    NREM_episodes = nrem_sleep(sData);

    NREM_start_end = round( (NREM_episodes./2500)*imaging_sampling_rate);
    REM_start_end = round( (REM_episodes./2500)*imaging_sampling_rate);
end

% AW_data
% QW_data

nrem_data = arrayfun(@(x) all_signals(:,  NREM_start_end(x,1):NREM_start_end(x,2)), (1:size(NREM_start_end,1)), 'UniformOutput', false);
rem_data = arrayfun(@(x) all_signals(:,  REM_start_end(x,1):REM_start_end(x,2)), (1:size(REM_start_end,1)), 'UniformOutput', false);

nrem_data_cat = horzcat( nrem_data{:});

%% Get states (NREM, REM, QW, AW);
output = get_state_logicals(sData);

nrem_test = all_signals(:, output{1,1});
qw_test   = all_signals(:, output{1,3});

%% All data
axon = vertcat(all_axon_activity{:});
pc = vertcat(all_pc_activity{:});
in = vertcat(all_in_activity{:});

axon_qw_test = axon(:, output{1,3});
axon_sum = sum(axon_qw_test,2);
[axonM, sortID] = sort(axon_sum);


% pc_qw_test = pc(:, output{1,3});
% pc_sum = sum(pc_qw_test,2);
% [~, sortID] = sort(pc_sum);


axon_nrem = arrayfun(@(x) axon(:,  NREM_start_end(x,1):NREM_start_end(x,2)), (1:size(NREM_start_end,1)), 'UniformOutput', false);
axon_nrem_mean = cellfun(@mean, axon_nrem, 'UniformOutput',false);

axon_rem = arrayfun(@(x) axon(:,  REM_start_end(x,1):REM_start_end(x,2)), (1:size(REM_start_end,1)), 'UniformOutput', false);


pc_nrem = arrayfun(@(x) pc(:,  NREM_start_end(x,1):NREM_start_end(x,2)), (1:size(NREM_start_end,1)), 'UniformOutput', false);
pc_rem  = arrayfun(@(x) pc(:,  REM_start_end(x,1):REM_start_end(x,2)), (1:size(REM_start_end,1)), 'UniformOutput', false);

in_nrem = arrayfun(@(x) in(:,  NREM_start_end(x,1):NREM_start_end(x,2)), (1:size(NREM_start_end,1)), 'UniformOutput', false);

%% Optional: plot single ROIs and states
% data = axon;
% 
% scale_factor = 5;
% 
% figure, plot( time_vec, data(sortID,:)+(1:size(data,1))'*scale_factor)
% 
% temp_y_max = max( (1:size(data,1))'*scale_factor)+1 ;
% temp_y_min = -1;
%  if isfield(sData, 'episodes')
% 
%         % NREM bouts
%         for ep_nr = 1:length(NREM_start_end)
%             x = [NREM_start_end(ep_nr,1) NREM_start_end(ep_nr,1) NREM_start_end(ep_nr,2) NREM_start_end(ep_nr,2)]./imaging_sampling_rate;
%             y = [temp_y_min temp_y_max temp_y_max temp_y_min];
%             patch(x, y, 'green', 'edgecolor', 'none', 'FaceAlpha', .2);
%         end
% 
%         % REM bouts
%         if size(REM_start_end(:,1),1) == 1
% 
%             a = [REM_start_end(1,1) REM_start_end(1,1) REM_start_end(1,2) REM_start_end(1,2)]./imaging_sampling_rate;
%             b =[temp_y_min temp_y_max temp_y_max temp_y_min];
%             patch(a, b, 'red', 'edgecolor', 'none', 'FaceAlpha', .2);
%         else
%             for ep_nr = 1:length(REM_start_end)
%                 a = [REM_start_end(ep_nr,1) REM_start_end(ep_nr,1) REM_start_end(ep_nr,2) REM_start_end(ep_nr,2)]./imaging_sampling_rate;
%                 b = [temp_y_min temp_y_max temp_y_max temp_y_min];
%                 patch(a, b, 'red', 'edgecolor', 'none', 'FaceAlpha', .2);
%             end
%         end
%  end
 
%% Plot top and bottom sort
% 
% data = axon;
% scale_factor = 4;
% 
% top_cells = data(sortID(1:10),:);
% bottom_cells = data(sortID(end-10:end),:);
% 
% axon_rem_top    = arrayfun(@(x) top_cells(:,  REM_start_end(x,1):REM_start_end(x,2)), (1:size(REM_start_end,1)), 'UniformOutput', false);
% axon_rem_bottom = arrayfun(@(x) bottom_cells(:,  REM_start_end(x,1):REM_start_end(x,2)), (1:size(REM_start_end,1)), 'UniformOutput', false);
% 
% axon_nrem_top    = arrayfun(@(x) top_cells(:,  NREM_start_end(x,1):NREM_start_end(x,2)), (1:size(NREM_start_end,1)), 'UniformOutput', false);
% axon_nrem_bottom = arrayfun(@(x) bottom_cells(:,  NREM_start_end(x,1):NREM_start_end(x,2)), (1:size(NREM_start_end,1)), 'UniformOutput', false);
% 
% [axon_rem_ep3_mean_top, axon_rem_ep3_SEM_top]  = cellfun(@episode_third_analysis, axon_rem_top, 'uni', 0);
% [axon_rem_ep3_mean_bottom, axon_rem_ep3_SEM_bottom]  = cellfun(@episode_third_analysis, axon_rem_bottom, 'uni', 0);
% 
% [axon_nrem_ep3_mean_top, axon_nrem_ep3_SEM_top]  = cellfun(@episode_third_analysis, axon_nrem_top, 'uni', 0);
% [axon_nrem_ep3_mean_bottom, axon_nrem_ep3_SEM_bottom]  = cellfun(@episode_third_analysis, axon_nrem_bottom, 'uni', 0);
% 
% 
% figure,
% plot(vertcat(axon_rem_ep3_mean_top{:})', 'r'), hold on, 
% plot( mean(vertcat(axon_rem_ep3_mean_top{:})), 'r', 'linew',2);
% 
% plot(vertcat(axon_rem_ep3_mean_bottom{:})', 'b')
% plot( mean(vertcat(axon_rem_ep3_mean_bottom{:})), 'b', 'linew',2);
% title('REM Episode 3 means from top and bottom 10 wake active cells')
% 
% figure, 
% plot(vertcat(axon_nrem_ep3_mean_top{:})', 'r'), hold on, 
% plot( mean(vertcat(axon_nrem_ep3_mean_top{:})), 'r', 'linew',2);
% 
% plot(vertcat(axon_nrem_ep3_mean_bottom{:})', 'b')
% plot( mean(vertcat(axon_nrem_ep3_mean_bottom{:})), 'b', 'linew',2);
% title('NREM Episode 3 means from top and bottom 10 wake active cells')
% 

%% 
% figure, 
% hold on
% plot( time_vec, top_cells+(1:10)'*scale_factor, 'color', 'r')
% plot( time_vec, bottom_cells+(11:21)'*scale_factor, 'color', 'b')
% 
% plot(time_vec, mean( top_cells)-8, 'color', 'r')
% plot(time_vec, mean( bottom_cells)-10, 'Color', 'b')
% 
% temp_y_max = max( 1:20'*scale_factor)+1 ;
% temp_y_min = -20;
%  if isfield(sData, 'episodes')
% 
%         % NREM bouts
%         for ep_nr = 1:length(NREM_start_end)
%             x = [NREM_start_end(ep_nr,1) NREM_start_end(ep_nr,1) NREM_start_end(ep_nr,2) NREM_start_end(ep_nr,2)]./imaging_sampling_rate;
%             y = [temp_y_min temp_y_max temp_y_max temp_y_min];
%             patch(x, y, 'green', 'edgecolor', 'none', 'FaceAlpha', .2);
%         end
% 
%         % REM bouts
%         if size(REM_start_end(:,1),1) == 1
% 
%             a = [REM_start_end(1,1) REM_start_end(1,1) REM_start_end(1,2) REM_start_end(1,2)]./imaging_sampling_rate;
%             b =[temp_y_min temp_y_max temp_y_max temp_y_min];
%             patch(a, b, 'red', 'edgecolor', 'none', 'FaceAlpha', .2);
%         else
%             for ep_nr = 1:length(REM_start_end)
%                 a = [REM_start_end(ep_nr,1) REM_start_end(ep_nr,1) REM_start_end(ep_nr,2) REM_start_end(ep_nr,2)]./imaging_sampling_rate;
%                 b = [temp_y_min temp_y_max temp_y_max temp_y_min];
%                 patch(a, b, 'red', 'edgecolor', 'none', 'FaceAlpha', .2);
%             end
%         end
%  end

 %%


% axon_nrem_cat = horzcat(axon_nrem{:});
% pc_nrem_cat   = 
% in_nrem_cat   = 
% pc_rem = pc(:, REM_start_end(1):REM_start_end(2));
% in_rem = in(:, REM_start_end(1):REM_start_end(2));
%% Plot all
if strcmp(params.plot, 'yes')

figure, 
    plot_order = [1 4 7 10 2 5 8 11 3 6 9 12];
    min_y = 1;
    max_y = -1;
    
    sgtitle(sData.sessionInfo.sessionID, 'interpreter', 'none', fontSize=10)
    for i = 1:12
    
        if i <=4
            clim = 'k';
        elseif i >= 5 && i <=8
            clim = 'r';
        elseif i >=9
            clim = 'b';
        end
            
        h(i) = subplot(4,3,plot_order(i));
    
        signal_to_plot = all_signals(i,:);
        sig_min        = min(signal_to_plot);
        sig_max        = max(signal_to_plot);
    
        hold on
        plot(time_vec, all_signals(i,:), 'Color',clim)
        title(['N rois = ', num2str( all_n_rois(i))])
         if isfield(sData, 'episodes')
    
            % NREM bouts
            for ep_nr = 1:length(NREM_start_end)
                x = [NREM_start_end(ep_nr,1) NREM_start_end(ep_nr,1) NREM_start_end(ep_nr,2) NREM_start_end(ep_nr,2)]./imaging_sampling_rate;
                y = [min_y max_y max_y min_y];
                patch(x, y, 'green', 'edgecolor', 'none', 'FaceAlpha', .2);
            end
    
            % REM bouts
            if size(REM_start_end(:,1),1) == 1
                
                a = [REM_start_end(1,1) REM_start_end(1,1) REM_start_end(1,2) REM_start_end(1,2)]./imaging_sampling_rate;
                b =[min_y max_y max_y min_y];
                patch(a, b, 'red', 'edgecolor', 'none', 'FaceAlpha', .2);
            else
                for ep_nr = 1:length(REM_start_end)
                    a = [REM_start_end(ep_nr,1) REM_start_end(ep_nr,1) REM_start_end(ep_nr,2) REM_start_end(ep_nr,2)]./imaging_sampling_rate;
                    b = [min_y max_y max_y min_y];
                    patch(a, b, 'red', 'edgecolor', 'none', 'FaceAlpha', .2);
                end
            end
         end
         
         if ~isnan(sig_min)
         set(gca, 'ylim', [sig_min sig_max])
         end
    end
    linkaxes(h, 'x')
    
    %% Plot correlation matrix of mean signals
    all_corr_mat = corr(all_signals');
    
    figure,
    sgtitle(sData.sessionInfo.sessionID, 'interpreter', 'none', fontSize=10)
    imagesc(all_corr_mat), 
    yline(4.5,'r', 'linew',2)
    yline(8.5,'r', 'linew',2)
    xline(4.5,'r', 'linew',2)
    xline(8.5,'r', 'linew',2)
    colorbar
    caxis([-1 1])
end
