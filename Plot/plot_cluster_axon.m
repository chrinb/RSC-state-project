function plot_cluster_axon(sData)

% Written by Christoffer Berge | Vervaeke Lab

% Loop over merged axons in sData, plot the individual ROI signals, their
% average, and the Okada-filtered signals. (Note that only the average and
% filtered time series are scaled to match the DF/F scale bar. The single
% ROIs are scaled separately).
params.use_roi_classification = 'grid';

% Load data
[pc_rois, ~] = remove_cells(sData, params);
dff          = sData.imdata.roiSignals(2).newdff(pc_rois,:);

% Find "true" ROI idx 
roi_arr                  = sData.imdata.roi_arr;
rois_after_grid_analysis = roi_arr(pc_rois);
fov_dim                  = roi_arr(1, 1).imagesize;

ID_of_ROIs_in_cluster  = sData.imdata.roiClustIDs;
clusterID     = sData.imdata.clusterID;

% Create time vector for plotting
time_imaging  = (0:length(dff)-1)/31;
x_label_pos   = -20;
fontsize      = 14;
sub_factor    = 10;

% Loop over nr of clusters
for i = 1:length(ID_of_ROIs_in_cluster)
    
    % Locate clustered axon data in sData
    temp                  = ismember(clusterID(2,:),i);
    temp2                 = clusterID(1,:);
    roi_idx               = temp2(temp);

    merged_roi_idx = ID_of_ROIs_in_cluster(i).rois;
    true_id_loc    = zeros(numel(merged_roi_idx),1);
    for j = 1:numel(merged_roi_idx)
        temp_id = rois_after_grid_analysis(merged_roi_idx(j)).uid;
        
        for k = 1:numel(roi_arr)
            if strcmp(roi_arr(k).uid, temp_id)
            true_id_loc(j) = k;
            end
        end
    end

    merged_axon_data      = sData.imdata.roiSignals(2).mergedAxonsDff(roi_idx, :);
    merged_axon_data_filt = sData.imdata.roiSignals(2).mergedAxonsDffFilt(roi_idx, :);

    figure(1), clf
    
    h(1) = subplot(121);
    hold on
    % loop over nr of ROIs in current cluster
    tmp_roi_array = zeros(fov_dim);
    for ii = 1:length(ID_of_ROIs_in_cluster(i).rois)
     
        plot_factor = 60 - sub_factor;

        tmp_roi_idx = ID_of_ROIs_in_cluster(i).rois(ii);
        
        

        plot(time_imaging, dff(tmp_roi_idx, :)-plot_factor, 'Color', 'k')
       
        plot( [610 610], [-plot_factor -plot_factor+0.5], 'LineWidth',2, 'Color', 'k');

        roi_string = ['ROI # ', num2str(true_id_loc(ii))];

        text(x_label_pos, -plot_factor, 0, roi_string, 'HorizontalAlignment', 'right','FontSize',fontsize);   
        
        sub_factor = sub_factor-1;
        
        tmp_roi_array        = tmp_roi_array + roi_arr(true_id_loc(ii)).mask;
        tmp_roi_position{ii} = roi_arr(true_id_loc(ii)).center;
    end
    
    plot_factor = 60 - sub_factor;
    plot(time_imaging, merged_axon_data-plot_factor,'r')
    plot( [610 610],  [-plot_factor -plot_factor+0.5], 'LineWidth',2, 'Color', 'k');
    str1 = 'Average';
    text(x_label_pos,  -plot_factor, 0, str1, 'HorizontalAlignment', 'right','FontSize',fontsize);
     
    sub_factor = sub_factor-1;
    plot_factor = 60 - sub_factor;
    plot(time_imaging, merged_axon_data_filt-plot_factor,'m')
    plot( [610 610],  [-plot_factor -plot_factor+0.5], 'LineWidth',2, 'Color', 'k');
    str2 = 'Filtered'; 
    text(x_label_pos,  -plot_factor, 0, str2, 'HorizontalAlignment', 'right','FontSize',fontsize);
    
    % Make scaling look better
    y_min      = min(merged_axon_data_filt-plot_factor)-0.1;
    tmp_lim    = h(1).YLim;
    tmp_lim(1) = y_min;
    % h(1).XLim = [time_imaging(1) 620];
    % h(1).YLim = tmp_lim;
    set(gca, 'xlim', [time_imaging(1) 620], 'ylim', tmp_lim)
    h(1).YAxis.Visible = "off";
    hold off 
    
    % Plot ROI masks
    h(2) = subplot(122);
    hold on
    h(2).XLim = [1, size(tmp_roi_array,2)];
    h(2).YLim = [1, size(tmp_roi_array,1)];
    imagesc(tmp_roi_array)
    colormap gray
    h(2).XAxis.Visible = "off";
    h(2).YAxis.Visible = 'off';

    % for txt_idx = 1:length(ID_of_ROIs_in_cluster(i).rois)
    % 
    %     if rem(txt_idx,2) == 0 || rem(txt_idx,2) == 1
    %         txt_position = 'right';
    %         displacement_f = -10;
    %     elseif rem(txt_idx,2) == 2
    %         txt_position = 'left';
    %         displacement_f = -10;
    %     end
    %     text( tmp_roi_position{1, txt_idx}(1)+displacement_f, tmp_roi_position{1, txt_idx}(2), 0, num2str(true_id_loc(txt_idx)), {txt_position}, 'FontSize',10, 'Color', 'w')
    % end

    % Plot average of clustered ROIs
    % h(ii+1) = subplot( n_rois_to_plot+3 ,1,ii+1); hold on
    % plot(time_imaging, merged_axon_data,'r')
    % plot( [610 610], [0 0.5], 'LineWidth',2, 'Color', 'k');
    % set(gca,'visible','off', 'xlim', [time_imaging(1) 620],'ylim', [min_lim dff_scale_bar])
    % str1 = 'Average';
    % text(x_label_pos, y_label_pos, 0, str1, 'HorizontalAlignment', 'right','FontSize',fontsize);
    % 
    % % Plot okada-filtered average of clustered ROIs
    % h(ii+2) = subplot( n_rois_to_plot+3 ,1,ii+2); hold on
    % plot(time_imaging, merged_axon_data_filt,'m')
    % plot( [610 610], [0 0.5], 'LineWidth',2, 'Color', 'k');
    % % set(gca,'visible','off', 'xlim', [time_imaging(1) 620], 'ylim', [min_lim dff_scale_bar])
    % set(gca, 'xlim', [time_imaging(1) 620], 'ylim', [min_lim dff_scale_bar])
    % h(ii+2).YAxis.Visible = 'off';
    % str2 = 'Filtered'; 
    % text(x_label_pos, y_label_pos, 0, str2, 'HorizontalAlignment', 'right','FontSize',fontsize);
    % 
    % plot([0 5], [-10 -10 ], 'LineWidth',2, 'Color', 'k');
    % h(ii+3) = subplot( n_rois_to_plot+3 ,1,ii+3);
    % plot(time_imaging, ones(1, size(dff,2)), 'Color', [1 1 1] ), hold on
    % plot([0 0], [dff_scale_bar 0], 'LineWidth',3, 'Color',[0 0 0])
    % plot([0 30], [0 0], 'LineWidth', 3, 'Color',[0 0 0])
    % set(gca,'visible','off', 'xlim', [time_imaging(1) time_imaging(end)], 'ylim', [min_lim dff_scale_bar])
    % text(1, 2, [num2str(dff_scale_bar*100), '% DF/F'],'HorizontalAlignment', 'right','FontSize',fontsize, 'Position',[-5 .2 0]);
    % text(1, 2, '30 s','HorizontalAlignment', 'right','FontSize',fontsize, 'Position',[25 min_lim-.01 0]);

%     xlabel('Time (sec)');
    % linkaxes(h, 'x')
    
    prompt = sprintf('\nCluster # %d, press "enter" for next cluster: ', i);
    x       = input(prompt,'s');
    if isempty(x)
        clear fig
        clc
    end
    clear data
end
