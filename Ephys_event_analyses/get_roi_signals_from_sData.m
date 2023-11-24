function [signal_to_plot, cmap, pc_rois, in_rois] = get_roi_signals_from_sData(sData, params)

% [signal, params, label3, rois_for_an, roiClustIDs, cells_to_exclude] = get_roi_signals_from_sData(sData, params)

% Written by Christoffer Berge | Vervaeke Lab

% Function that (1) selects which ROI signal to analyze, (2) experiment
% type (bulk, axon, all ROIs), (3) whether to split roi array (e.g., into
% excitatory vs. inhibitory cells) or access a subset of ROIs (NOT
% FINISHED)

[pc_rois, in_rois] = remove_cells_longitudinal(sData, params);


switch params.signal_type
    case 'Dff'
        txt      = 'newdff';
        axon_txt = 'DffFilt';
        txt2     = '';
    case 'deconv'
        txt = 'ciaDeconvolved';
        axon_txt = 'Dec';
        txt2     = '';
    case 'transients'
        axon_txt = '_sig_transients';
        txt      = 'all_sig_transients';
%         txt2     = params.cell_type;
end

switch params.cell_type
    case {'pc', 'in'}
        signal_to_plot{1,:} = sData.imdata.roiSignals(2).(txt)(pc_rois,:);
        signal_to_plot{2,:} = sData.imdata.roiSignals(2).(txt)(in_rois,:);
        cmap                = [0 2];
    case 'axon'
        signal_to_plot{1,:} = sData.imdata.roiSignals(2).(['mergedAxons',axon_txt]);
        cmap                = [0 2 ];
end

% Z-score normalize
if strcmp(params.zscore, 'yes') 
    signal_to_plot{1,:} = zscore(signal_to_plot{1,:} , 0 ,2);
    try signal_to_plot{2,:} = zscore(signal_to_plot{2,:} , 0 ,2); catch, end
% elseif strcmp(params.zscore, 'yes') && strcmp(params.cell_type, 'axon')
%     signal_to_plot{1,:} = zscore(signal_to_plot{1,:} , 0 ,2);
end

% 
% %% Select which ROI signal to use
% checkParameter = @(param, n, str) (isnumeric(param) && param==n) || strcmp(param, str);
% 
% if checkParameter(params.signal_type, 1, 'deconv')
%     signal = sData.imdata.roiSignals(2).ciaDeconvolved;
%     text   = 'deconv.';
% elseif checkParameter(params.signal_type, 2, 'dff')
%     signal = sData.imdata.roiSignals(2).newdff;
%     text   = '';
% elseif checkParameter(params.signal_type, 3, '')
%     prompt = sprintf('Type structname: '); 
%     signal = input(prompt);
%     text = [];
% end

%% Select which experiment type 

% Bulk dendritic analysis: average over ROIs (4x4 grid, 16 ROIs) to get average
% activity in FOV

% Axonal analysis: merge correlated ROIs

% Otherwise run analysis on all ROIs


%% Check if ROI array is longitudinal and if yes use ROI classification index to select ROIs
% if isfield(sData.imdata, 'roi_classification')
% 
%     % Get indicies of cells present in current session
%     roi_idx = sData.imdata.roi_classification == 1;
% 
%     % Get signals from present cells
%     signal = signal(roi_idx,:);
% 
%     % Find ROIs to exclude (because not present in current session)
%     temp  = 1:length(sData.imdata.roi_arr);
%     temp2 = ~roi_idx;
%     cells_to_exclude = temp(temp2);
% 
%     % Update ROI array to reflect subset of ROIs present in current session
%     roi_arr = sData.imdata.roi_arr(roi_idx);
% 
% if isfield(sData.imdata, 'roi_arr')
%     roi_arr = sData.imdata.roi_arr;
% else
%     roi_arr = [];
% end
% 
% roiClustIDs = [];
% if checkParameter(params.exp_type, 1, 'bulk') % bulk analysis
% 
%     % In the case of axon imaging sessions, check if there is a bulk axon
%     % signal and use that.  
%     try 
%         signal = mean(sData.imdata.roiSignals(2).axonBulkDff);
%         label3 = '# SWR';
%         cells_to_exclude = [];
%     catch
%         signal = mean(signal);
%         label3 = '# SWR';
%         cells_to_exclude = [];
%     end
% 
% elseif checkParameter(params.exp_type, 3, 'axon')   % axons
%     
%     % Check if there is a field in sData called 'mergedAxons', if not,
%     % merge.
%     cells_to_exclude = [];
%     if ~isfield(sData.imdata.roiSignals(2), 'mergedAxonsDffFilt')
%        
%         % Select principal/excitatory cells
%         [pc_rois, ~]       = remove_cells(sData);
% 
%         % Select DF/F for axon clustering
%         ROImtx = sData.imdata.roiSignals(2).newdff(pc_rois,:);
%         % Cluster axons
%         [merged_ROI_dff,merged_ROI_deconv, roiClustIDs] = hierClust_axons(ROImtx,1, sData, pc_rois);
%         % Select either the merged DF/F or deconvolved signal
%         if checkParameter(params.signal_type, 1, 'deconv')
%             signal = merged_ROI_deconv;
%         elseif checkParameter(params.signal_type, 2, 'dff')
%             signal = merged_ROI_dff;
%         end
% 
%     else
%         if checkParameter(params.signal_type, 1, 'deconv')
%             signal = sData.imdata.roiSignals(2).mergedAxonsDec;
%         elseif checkParameter(params.signal_type, 2, 'dff')
%             signal = sData.imdata.roiSignals(2).mergedAxonsDff;
% %             signal = sData.analysis.transients.sig_transients;
%         end
% 
%     end
%     label3 = '# ROI';
% else
%     label3 = '# ROI';
% end
% 
% %% Select neuron subtype for analysis
% 
% % Get indicies of excitatory and inhibitory cells
% [pc_rois, in_rois] = remove_cells_longitudinal(sData);
% 
% % Population analysis + select PCs
% if  checkParameter(params.exp_type, 2, 'default') && checkParameter(params.cell_type, 1 , 'pc')
% %     [pc_rois, in_rois] = remove_cells(sData); % Select principal cells for analysis
%     rois_for_an        = pc_rois;
%     cells_to_exclude   = in_rois;
% % Population analysis + select INs
% elseif checkParameter(params.exp_type, 2, 'default') && checkParameter(params.cell_type, 2 , 'in')
% %     [pc_rois, in_rois] = remove_cells(sData); % Select inhibitory cells for analysis
%     rois_for_an        = in_rois;
%     cells_to_exclude   = pc_rois;
% 
% % Axon population analysis + get IN indicies
% elseif checkParameter(params.exp_type, 3, 'axon')
% %     [~, in_rois]     = remove_cells(roi_arr); % Get IN indicies
% %     cells_to_exclude = in_rois;
%     rois_for_an      = 1:size(signal,1); 
% 
% elseif checkParameter(params.exp_type, 1, 'bulk')
%     rois_for_an = 1:size(signal,1); % Select all ROIs (1 in the case of bulk analysis)
% end
% 
% 
% 
% 
