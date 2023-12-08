function plot_piezo_session_activity(sData, params)

% Written by Christoffer Berge | Vervaeke lab

%{
Plot state-activity vector (hypnogram), roi x time colorplots, and
averages for different cell types in multi-plane recordings (Piezo)
%}


%% Find indices of excitatory & inhibitory cells, and axons
[cell_idx_pc, cell_idx_in, cell_idx_axon, all_idx_pc, all_idx_in, ~] = remove_cells_piezo(sData);

% Check which of the four planes contains ROIs of the cell-type of interest
% and determine how many planes to plot based on that. 
plane_idx = [1,2,3,4];
switch params.cell_type
    case 'in'
        plane_nr = plane_idx( ~cellfun(@isempty, cell_idx_in));
    case 'axon'
        plane_nr = plane_idx( ~cellfun(@isempty, cell_idx_axon));
    case 'pc'
        plane_nr = plane_idx( ~cellfun(@isempty, cell_idx_pc));
end


%% Create hypnogram vector

total_frames_per_plane = ceil( sData.daqdata.frame_onset_reference_frame(end)/4 );

% It seems like some nr of frames (10-300) are recorded in one plane before
% piezo start. These "extra" frames are not present in the
% frame_onset_reference_frame variable. Therefore use that variable to
% remove extra frames.
frames_in_rec = size(sData.imdata.roiSignals(2).newdff, 2);
dff_diff      = frames_in_rec - total_frames_per_plane;

% Initialize empty vector
hypnogram_vector = zeros(1, total_frames_per_plane );

state_vectors_2p = get_state_logicals(sData);

% Create hypnogram vector
hypnogram_vector(state_vectors_2p{1,1} == 1)  = 2; % NREM
hypnogram_vector(state_vectors_2p{1,2} == 1)  = 3; % REM
hypnogram_vector(state_vectors_2p{1,3} == 1)  = 0; % QW
hypnogram_vector(state_vectors_2p{1,4}  == 1) = 1; % AW

% Get 2P frame rate
imaging_sampling_rate = find_imaging_framerate(sData);

% Divide by 4 to get FPS per plane
imaging_sampling_rate = imaging_sampling_rate/4;

% If sleep recording, get start/stop of NREM & REM episodes
if isfield(sData, 'episodes')
    REM_episodes = rem_sleep(sData);
    NREM_episodes = nrem_sleep(sData);

    NREM_start_end = NREM_episodes./2500;
    REM_start_end = REM_episodes./2500;
end
%% Plot 
figure,

time_vector  = linspace(0, total_frames_per_plane, total_frames_per_plane)/imaging_sampling_rate;
font_size    = 14;

h(1) = subplot( numel(plane_nr)*2 + 1, 1, 1);
plot(time_vector, hypnogram_vector, 'LineWidth',1)
set(gca, 'xlim', [time_vector(1) time_vector(end)])
axis off
ytickvalues = 0:3;
str = {'QW', 'AW','NREM', 'REM'};
for current_plane = ytickvalues
    text(-5, ytickvalues(current_plane+1), str(current_plane+1), 'HorizontalAlignment', 'right','FontSize',font_size);
end

subplot_nr = 1;

for current_plane = plane_nr

    % Get DF/F from selected cell type in current plane
    switch params.cell_type
        case 'in'
            % Find the logical indices of the ROIs for plane nr = i
            roi_plane_idx = sData.imdata.plane_indices(all_idx_in) == current_plane;
            % Extract the DF/F for cell type for current plane (and remove
            % piezo artefacts)
            temp_dff      = sData.imdata.roiSignals(2).newdff( all_idx_in(roi_plane_idx), dff_diff+1:end);
            clim = [0 1];
        case 'axon'
            temp_dff = sData.imdata.roiSignals(2).mergedAxonsDffFilt( sData.imdata.plane_indices_merged == current_plane, dff_diff+1:end);
            clim = [0 .5];
        case 'pc'
            roi_plane_idx = sData.imdata.plane_indices(all_idx_pc) == current_plane;
            temp_dff      = sData.imdata.roiSignals(2).newdff( all_idx_pc(roi_plane_idx), dff_diff+1:end);
            clim = [0 1];
    end
        
    % Get nr of ROIs from current plane
    if strcmp(params.zscore, 'yes')
        temp_dff = zscore(temp_dff, 0 ,2);
    end
    % temp_dff = okada(temp_dff,2);

    y1       = [1, size(temp_dff,1) ];

    h(subplot_nr*2) = subplot( numel(plane_nr)*2 + 1, 1, subplot_nr*2);
    imagesc(time_vector, y1, temp_dff)
    ylabel(['Plane ',num2str(current_plane)], FontSize=14)
    caxis(clim)
    c = colorbar;
    c.Position(1) = 0.92;
    c.Position(3) = 0.01;

    signal_to_plot = mean(temp_dff, 'omitnan');

    h(subplot_nr*2+1) = subplot( numel(plane_nr)*2 + 1, 1, subplot_nr*2+1) ;
    plot(time_vector, signal_to_plot)
    set(gca, 'xlim', [time_vector(1) time_vector(end)])
    
    % Plot states on top of mean activity traces 
    max_y = max( signal_to_plot, [], 2);
    min_y = min ( signal_to_plot, [], 2);

    if isfield(sData, 'episodes')
        % NREM bouts
        for ep_nr = 1:length(NREM_start_end)
            x = [NREM_start_end(ep_nr,1) NREM_start_end(ep_nr,1) NREM_start_end(ep_nr,2) NREM_start_end(ep_nr,2)];
            y = [min_y max_y max_y min_y];
            patch(x, y, 'green', 'edgecolor', 'none', 'FaceAlpha', .2);
        end
        % REM bouts
        if size(REM_start_end(:,1),1) == 1
            a = [REM_start_end(1,1) REM_start_end(1,1) REM_start_end(1,2) REM_start_end(1,2)];
            b =[min_y max_y max_y min_y];
            patch(a, b, 'red', 'edgecolor', 'none', 'FaceAlpha', .2);
        else
            for ep_nr = 1:length(REM_start_end)
                a = [REM_start_end(ep_nr,1) REM_start_end(ep_nr,1) REM_start_end(ep_nr,2) REM_start_end(ep_nr,2)];
                b = [min_y max_y max_y min_y];
                patch(a, b, 'red', 'edgecolor', 'none', 'FaceAlpha', .2);
            end
        end
    end

    subplot_nr = subplot_nr + 1;

end

linkaxes(h, 'x')