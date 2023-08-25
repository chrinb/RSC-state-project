function plot_piezo_session_activity(varargin)

% Written by Christoffer Berge | Vervaeke lab

% Plot state-activity vector (hypnogram), roi x time colorplots, and
% averages for different cell types in multi-plane recordings (Piezo)

sData  = varargin{1,1};
params = varargin{1,2};

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
for i = ytickvalues
    text(-5, ytickvalues(i+1), str(i+1), 'HorizontalAlignment', 'right','FontSize',font_size);
end

subplot_nr = 1;

for i = plane_nr

    % Get DF/F from selected cell type in current plane
    switch params.cell_type
        case 'in'
            % Find the logical indices of the ROIs for plane nr = i
            roi_plane_idx = sData.imdata.plane_indices(all_idx_in) == i;
            % Extract the DF/F for cell type for current plane (and remove
            % piezo artefacts)
            temp_dff      = sData.imdata.roiSignals(2).newdff( all_idx_in(roi_plane_idx), dff_diff+1:end);
        case 'axon'
            temp_dff = sData.imdata.roiSignals(2).mergedAxonsDffFilt( sData.imdata.plane_indices_merged == i, dff_diff+1:end);
        case 'pc'
            roi_plane_idx = sData.imdata.plane_indices(all_idx_pc) == i;
            temp_dff      = sData.imdata.roiSignals(2).newdff( all_idx_pc(roi_plane_idx), dff_diff+1:end);
    end
        
    % Get nr of ROIs from current plane
    y1       = [1, size(temp_dff,1) ];

    h(subplot_nr*2) = subplot( numel(plane_nr)*2 + 1, 1, subplot_nr*2);
    imagesc(time_vector, y1, temp_dff)
    ylabel('ROI #', FontSize=14)
    caxis([0 .5])
    c = colorbar;
    c.Position(1) = 0.92;
    c.Position(3) = 0.01;

    h(subplot_nr*2+1) = subplot( numel(plane_nr)*2 + 1, 1, subplot_nr*2+1) ;
    plot(time_vector, mean(temp_dff, 'omitnan'))
    set(gca, 'xlim', [time_vector(1) time_vector(end)])

    subplot_nr = subplot_nr + 1;

end

linkaxes(h, 'x')