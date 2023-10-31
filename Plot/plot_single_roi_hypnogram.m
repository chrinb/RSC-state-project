function plot_single_roi_hypnogram(sData, roi_nr)

%{
Plot ROI activity trace and overlay different states in colored shades
%}

% Choose parameters
params.signal_type = 'Dff';
params.cell_type   = 'pc';
params.zscore      = 'no';

% Get signal to plot
[signal_to_plot, cmap, pc_rois, in_rois] = get_roi_signals_from_sData(sData, params);

if strcmp(params.cell_type, 'pc')
    signal_to_plot = signal_to_plot{1,:};
    txt = 'Excitatory cells';
elseif strcmp(params.cell_type, 'in')
    signal_to_plot = signal_to_plot{2,:};
    txt = 'Inhibitory cells';
elseif strcmp(params.cell_type, 'axon')
    txt = 'Axons';
end

% Get state vectors in 2P time
state_vectors_2p = get_state_logicals(sData);


imaging_sampling_rate = find_imaging_framerate(sData);
roi_data = signal_to_plot(roi_nr,:) ;
time_vector  = linspace(0, length(roi_data), length(roi_data))/imaging_sampling_rate;

if isfield(sData, 'episodes')
    REM_episodes = rem_sleep(sData);
    NREM_episodes = nrem_sleep(sData);

    NREM_start_end = NREM_episodes./2500;
    REM_start_end = REM_episodes./2500;
end

font_size = 16;

figure, 
plot(time_vector, roi_data ), hold on
patch_edges = [-.2 1 1 -1];

if isfield(sData, 'episodes')
    % NREM bouts
    for i = 1:length(NREM_start_end)
        x = [NREM_start_end(i,1) NREM_start_end(i,1) NREM_start_end(i,2) NREM_start_end(i,2)];
        y = patch_edges;
        patch(x, y, 'green', 'edgecolor', 'none', 'FaceAlpha', .2);
    end
    % REM bouts
    if size(REM_start_end(:,1),1) == 1
        a = [REM_start_end(1,1) REM_start_end(1,1) REM_start_end(1,2) REM_start_end(1,2)];
        b = patch_edges;
        patch(a, b, 'red', 'edgecolor', 'none', 'FaceAlpha', .2);
    else
        for i = 1:length(REM_start_end)
            a = [REM_start_end(i,1) REM_start_end(i,1) REM_start_end(i,2) REM_start_end(i,2)];
            b = patch_edges;
            patch(a, b, 'red', 'edgecolor', 'none', 'FaceAlpha', .2);
        end
    end
end