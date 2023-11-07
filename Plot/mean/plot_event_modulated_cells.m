function plot_event_modulated_cells(sData)

%% Variables
params.cell_type          = 'pc';
params.signal_type        = 'Dff';
params.zscore             =  'yes';
params.event_type         = 'SWR';
params.beh_state          = 'awake';
params.select_swr         = 'all';
params.remSpin            = 'all';
params.baseSub            = 'No subtract'; 
params.swr_for_analysis   = 4; 
params.swr_cluster_thresh = 2;
% Define variables
nr_of_seconds  = 3;
nr_of_frames   = (nr_of_seconds*31*2)+1;
frames         = sData.daqdata.frame_onset_reference_frame;
time           = (-(31*nr_of_seconds):(31*nr_of_seconds))./31;
sessionID      = sData.sessionInfo.sessionID;


%% Get all signals
[signal_to_plot, ~, pc_rois, in_rois] = get_roi_signals_from_sData(sData, params);

switch params.cell_type
    case 'pc'
        neg_roi_idx = sData.imdata.SWR_pc_neg_cells;
        pos_roi_idx = sData.imdata.SWR_pc_pos_cells;
        signal  = signal_to_plot{1,:};
    case 'in'
        neg_roi_idx = sData.imdata.SWR_in_neg_cells;
        pos_roi_idx = sData.imdata.SWR_in_pos_cells;        
        signal  = signal_to_plot{2,:};
    case 'axon'
        pos_roi_idx = sData.imdata.SWR_in_neg_cells;
        signal  = signal_to_plot{1,1};
end

%% Get event idx
switch params.event_type
    case 'SWR'

        if strcmp(params.beh_state, 'awake')
            select_swr_idx = sData.ephysdata.absRipIdx; 
        elseif strcmp(params.beh_state, 'sleep')
            select_swr_idx = sort([sData.ephysdata.NREM_spindle_uncoupled_swr, sData.ephysdata.spindle_coupled_swr]);
        end

        % Get the indicies of user specified SWR types
        event_idx = get_swr_idx(params.swr_for_analysis ,sData,select_swr_idx, params);
        
        % Convert SWR time stamps from e-phys to 2P time
        event_idx = frames(round(event_idx));
        
        % Preallocate 
        temp_mean_event_activity = zeros(size(pos_roi_idx,1), nr_of_frames);
        peri_event_activity = zeros(length(event_idx), nr_of_frames) ;
        temp_all_data            = zeros(size(pos_roi_idx,1), size(event_idx,2), size(time,2));                

        time_win = (62:108); % -1s to +500ms after SWR peak
    case 'Spindle'
        spindle_select     = [];
        spin_start_end_str = strcat('NREMspindleStartEnd', spindle_select);
        
        % convert sleep spindle time stamps from e-phys to 2P time
        event_idx = sData.ephysdata2.(spin_start_end_str)(:,1);
    
        % Get time points for sleep spindle onsets
        event_idx = get_spindle_idx(params , event_idx);
        
        time_win             = (93:124); % 0s to +1s after spindle onset
    case 'SWA'
%         output = so_sig_an(sData, params);
end

%% Extract signals
roi_idx = {pos_roi_idx; neg_roi_idx};


for i = 1:2
    t = 1;
    temp_roi_idx = roi_idx{i};

    for roi_nr = 1:size(temp_roi_idx,1)

        roi = temp_roi_idx(roi_nr);
    
        % Extract signal from currrent ROI
        roi_signal        = okada( signal(roi, :),2); 
        % roi_signal        = signal(roi, :); 

        % Get all peri-SWR activity windows from current ROI
        peri_event_activity = extract_avg_activity(event_idx, roi_signal, peri_event_activity, nr_of_seconds, params);
    
        % store all roi x time x SWR data
        temp_all_data(t,:,:)   = peri_event_activity;
        
        % store mean SWR-activity
        temp_mean_event_activity(roi_nr,:)  = mean(peri_event_activity, 'omitnan');
        t = t+1;

    end

    all_data{i}            = temp_all_data;
    mean_event_activity{i} = temp_mean_event_activity;
    clear temp_all_data temp_mean_event_activity
end

%% Loop over positively-modulated cells to find good examples
figure,

switch params.cell_type
    case 'pc'
        all_rois_idx = pc_rois;
    case 'in'
        all_rois_idx = in_rois;
end
store_roi = zeros(1, size(neg_roi_idx,1));
for i = 1:size(neg_roi_idx,1)
    subplot(211)
    imagesc( squeeze (all_data{1,2}(i,:,:)))
    axis square
    c = colorbar;
    c.Position(1) = 0.8;
    clim([0 2])
    subplot(212)
    plot(mean( squeeze(all_data{1,2}(i,:,:))))
    axis square
    title(['ROI nr ', num2str(all_rois_idx(neg_roi_idx(i))) ])
    prompt = sprintf('Store ROI %d', i);
    x = input(prompt,'s');
    
    if strcmp(x,'y')
        store_roi(i) = i;       
        clear figure
    end
end
         


%% Plot

figure,
h(1) = subplot(3,2, [1,3]);
y_lim1 = [1 size(mean_event_activity{1,1},1)];
imagesc(time, y_lim1, mean_event_activity{1,1});

h(2) = subplot(3,2, [2,4]);
y_lim2 = [1 size(mean_event_activity{1,2},1)];
imagesc(time, y_lim2, mean_event_activity{1,2});

h(3) = subplot(3,2,5);
plot(time,  mean(mean_event_activity{1,1}))

h(4) = subplot(3,2,6);
plot(time,  mean( mean_event_activity{1,2}))

linkaxes(h, 'x')

set(gca, 'xlim', [time(1) time(end)])
