function plot_subset_mod_neurons(sData, roi_list)

% Written by Christoffer Berge | Vervaeke Lab

% User inputs sData and list of modulated neurons. Function loops over
% neurons and plot their peri-SWR-activity for visualization.

nr_of_seconds = 3;

% Add roisignals folder to path to get roi array
try
    folder_name = dir;
    addpath( [folder_name(1).folder '\roisignals'])
catch
end


dF_F   = sData.imdata.roiSignals(2).newdff;
deconv = sData.imdata.roiSignals(2).ciaDeconvolved;


%% Select SWRs for analysis
% prompt     = sprintf('Select SWRs to analyze: 1 = All | 2 = Awake | 3 Spindle-Uncoupled | 4 Spindle-Coupled ');
swr_select = 1;
swr_idx    = get_sleep_swr_idx(sData, swr_select);

% prompt = sprintf('All SWRs (1) | Rem. locom. SWRs (2) | Rem. clust. SWRs (3) | Rem- locom. & clust. SWRs (4) ');
swr_for_analysis      = 4;

% Get the indicies of user specified SWR types
swr_idx = get_swr_idx(swr_for_analysis, sData, swr_idx);

%% 
win_length = (nr_of_seconds*31)*2+1;
time       = linspace(-nr_of_seconds,nr_of_seconds,win_length);

% find SWR indicies in 2P imaging frame time
frames      = sData.daqdata.frame_onset_reference_frame;
swr_idx = frames(round(swr_idx));

% Initialize variables
[npil_swr_activity,dFF_swr_activity,deconv_swr_activity, dFF_swr_activity_zscore,deconv_swr_activity_zscore] ...
    = deal(zeros(length(swr_idx), win_length));

% loop over ROIs
for roinr = roi_list'
    
    % Select the activity for a given ROI
    % ROI_npil_signal   = npil(roinr, :);
    roi_dff           = dF_F(roinr,:);
    roi_dff_zscore    = zscore(dF_F(roinr,:));
    roi_deconv        = deconv(roinr,:);
    roi_deconv_zscore = zscore(deconv(roinr,:));

    % Loop over SWRs and extract 2P segments
    for swr_nr = 1:length(swr_idx)
        
        % Find idx of first and last frame in SWR-window
        swr_window_start = swr_idx(swr_nr) - (nr_of_seconds*31); 
        swr_window_end   = swr_idx(swr_nr) + (nr_of_seconds*31); 
    
        % skip SWRs at the beginning or end with a window shorter than
        % length specified in seconds by user above
        if swr_window_start > 1 && swr_window_end < length(dF_F)
    
    %             if length(swr_window_start:swr_window_end) == win_length
    %         npil_swr_activity(swr_nr, :) = ...
    %             ROI_npil_signal(swr_window_start:swr_window_end); 
            dFF_swr_activity(swr_nr, :) = ...
                roi_dff(swr_window_start:swr_window_end); 
            dFF_swr_activity_zscore(swr_nr,:) = ...
                roi_dff_zscore(swr_window_start:swr_window_end); 
    
            deconv_swr_activity(swr_nr, :) = ...
                roi_deconv(swr_window_start:swr_window_end); 
            deconv_swr_activity_zscore(swr_nr, :) = ...
                roi_deconv_zscore(swr_window_start:swr_window_end); 
        end
    end

    deconv_swr_activity        = smoothdata(deconv_swr_activity,2, 'gaussian', 5);
    deconv_swr_activity_zscore = smoothdata(deconv_swr_activity_zscore,2, 'gaussian', 5);
    
    c_lim_dff = [-3 3];
    c_lim_dec = [0 2];
    sessionID = sData.sessionInfo.sessionID;
    
    
    figure,
    x1 = [time(1), time(end)];
    y1 = [1 size(dFF_swr_activity,1)];
    
    sgtitle([ sessionID, ', ROI # ' num2str(roinr)], 'Interpreter', 'none')
    
    subplot(241);
    imagesc(x1, y1, dFF_swr_activity)
    ylabel('SWR #')
    title('DF/F')
    c(1) = colorbar;
    c(1).Position(1) = 0.289;
    c(1).Position(2) = 0.552;
    c(1).Position(3) = 0.006;
    
    subplot(242);
    imagesc(x1, y1,dFF_swr_activity_zscore)
    % ylabel('SWR #')
    title('DF/F (z-score)')
    c(2) = colorbar;
    c(2).Position(1) = 0.496;
    c(2).Position(2) = 0.552;
    c(2).Position(3) = 0.006;
    
    caxis(c_lim_dff)
    
    subplot(243)
    imagesc(x1, y1, deconv_swr_activity)
    % ylabel('SWR #')
    title('Deconv. DF/F')
    c(3) = colorbar;
    c(3).Position(1) = 0.703;
    c(3).Position(2) = 0.552;
    c(3).Position(3) = 0.006;
    
    subplot(244)
    imagesc(x1, y1, deconv_swr_activity_zscore)
    % ylabel('SWR #')
    title('Deconv. DF/F (z-score)')
    c(4) = colorbar;
    c(4).Position(1) = 0.909;
    c(4).Position(2) = 0.552;
    c(4).Position(3) = 0.006;
    caxis(c_lim_dec)
    
    subplot(245),
    shadedErrorBar(time, mean(dFF_swr_activity,'omitnan'), ...
        std(dFF_swr_activity)/sqrt(numel(dFF_swr_activity(:, 1))) ,'lineprops', 'b');
    xlabel('Time from SWR peak (s)')
    ylabel('Mean DF/F')
    set(gca, 'xlim',[time(1) time(end)])
    
    subplot(246),
    shadedErrorBar(time, mean(dFF_swr_activity_zscore, 'omitnan'), ...
        std(dFF_swr_activity_zscore)/sqrt(numel(dFF_swr_activity_zscore(:, 1))) ,'lineprops', 'b');
    xlabel('Time from SWR peak (s)')
    ylabel('Mean DF/F (z-score)')
    set(gca, 'xlim',[time(1) time(end)])
    
    subplot(247),
    shadedErrorBar(time, mean(deconv_swr_activity,'omitnan'), ...
        std(deconv_swr_activity)/sqrt(numel(deconv_swr_activity(:, 1))) ,'lineprops', 'b');
    xlabel('Time from SWR peak (s)')
    ylabel('Mean deconv. DF/F')
    set(gca, 'xlim',[time(1) time(end)])
    
    subplot(248),
    shadedErrorBar(time, mean(deconv_swr_activity_zscore,'omitnan'), ...
        std(deconv_swr_activity_zscore)/sqrt(numel(deconv_swr_activity_zscore(:, 1))) ,'lineprops', 'b');
    xlabel('Time from SWR peak (s)')
    ylabel('Mean deconv. DF/F (z-score)')
    set(gca, 'xlim',[time(1) time(end)], 'ylim',[-.1 1])


    prompt = sprintf('Next ROI? ');
    x       = input(prompt,'s');
    if isempty(x)
        close 
        clc
    end
end



