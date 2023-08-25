function plot_cluster_axon(sData)

% Written by Christoffer Berge | Vervaeke Lab

% Loop over merged axons in sData, plot the individual ROI signals, their
% average, and the Okada-filtered signals. (Note that only the average and
% filtered time series are scaled to match the DF/F scale bar. The single
% ROIs are scaled separately).

% Load data
[pc_rois, ~] = remove_cells(sData.imdata.roi_arr);
dff          = sData.imdata.roiSignals(2).newdff(pc_rois,:);
roiClustIDs  = sData.analysis.roiClustIDs;
clusterID    = sData.analysis.clusterID;
% Create time vector for plotting
time_imaging = (0:length(dff)-1)/31;
Color        = [0 0 0];
label_coord  = 0;
fontsize     = 10;
dff_scale_bar = 0.7;
% Loop over nr of clusters
for i = 1:length(roiClustIDs)
    
    % Locate clustered axon data in sData
    temp                  = ismember(clusterID(2,:),i);
    temp2                 = clusterID(1,:);
    roi_idx               = temp2(temp);

    merged_axon_data      = sData.imdata.roiSignals(2).mergedAxonsDff(roi_idx, :);
    merged_axon_data_filt = sData.imdata.roiSignals(2).mergedAxonsDffFilt(roi_idx, :);

    max_lim       = max(  [max( max(merged_axon_data)) max(max(merged_axon_data_filt))] );
    min_lim       = min(  [min( min(merged_axon_data)) min(min(merged_axon_data_filt))] );

    figure(1), clf
    n_rois_to_plot = length(roiClustIDs(i).rois);
    % loop over nr of ROIs in current cluster
    for ii = 1:length(roiClustIDs(i).rois)
        
        h(ii) = subplot( n_rois_to_plot+3 ,1,ii);
        plot(time_imaging, dff(roiClustIDs(i).rois(ii),:), 'Color', Color)
        set(gca,'visible','off', 'xlim', [time_imaging(1) time_imaging(end)])
        roi_string = string(roiClustIDs(i).rois(ii));
        text(1, label_coord, roi_string, 'HorizontalAlignment', 'right','FontSize',fontsize, 'Position',[-5 label_coord 0]);
    end
    
    % Plot average of clustered ROIs
    h(ii+1) = subplot( n_rois_to_plot+3 ,1,ii+1);
    plot(time_imaging, merged_axon_data,'r')
    set(gca,'visible','off', 'xlim', [time_imaging(1) time_imaging(end)],'ylim', [min_lim dff_scale_bar])
    str1 = 'Average';
    text(1, label_coord, str1, 'HorizontalAlignment', 'right','FontSize',fontsize, 'Position',[-5 label_coord 0]);
    
    % Plot okada-filtered average of clustered ROIs
    h(ii+2) = subplot( n_rois_to_plot+3 ,1,ii+2);
    plot(time_imaging, merged_axon_data_filt,'m')
    set(gca,'visible','off', 'xlim', [time_imaging(1) time_imaging(end)], 'ylim', [min_lim dff_scale_bar])
    str2 = 'Filtered'; 
    text(1, label_coord, str2, 'HorizontalAlignment', 'right','FontSize',fontsize, 'Position',[-5 label_coord 0]);
    
    h(ii+3) = subplot( n_rois_to_plot+3 ,1,ii+3);
    plot(time_imaging, ones(1, size(dff,2)), 'Color', [1 1 1] ), hold on
    plot([0 0], [dff_scale_bar 0], 'LineWidth',3, 'Color',[0 0 0])
    plot([0 30], [0 0], 'LineWidth', 3, 'Color',[0 0 0])
    set(gca,'visible','off', 'xlim', [time_imaging(1) time_imaging(end)], 'ylim', [min_lim dff_scale_bar])
    text(1, 2, [num2str(dff_scale_bar*100), '% DF/F'],'HorizontalAlignment', 'right','FontSize',fontsize, 'Position',[-5 .2 0]);
    text(1, 2, '30 s','HorizontalAlignment', 'right','FontSize',fontsize, 'Position',[25 min_lim-.01 0]);

%     xlabel('Time (sec)');
    linkaxes(h, 'x')
    
    prompt = sprintf('\nCluster # %d, press "enter" for next cluster: ', i);
    x       = input(prompt,'s');
    if isempty(x)
        close 
        clc
    end
    clear data
end
