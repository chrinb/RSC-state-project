function plot_transients(varargin)

% Written by Christoffer Berge | Vervaeke lab

% Plot a random subset of ROI activity patterns.
% User inputs sData, number of ROIs to plot, and signal type
% (e.g., 'pc' for principal cells, 'in' for inhibitory cells, 'axon' for
% axonal data).

% Load data
sData              = varargin{1,1};
% dff                = sData.imdata.roiSignals(2).newdff;

% Find excitatory/inhibitory indicies 
[pc_rois, in_rois] = remove_cells(sData);

% If session is part of multi-day recordings, remove ROIs not present in
% current session
if isfield(sData.imdata, 'roi_classification')

    roi_classification = sData.imdata.roi_classification;

    pc_roi_idx = roi_classification(pc_rois); % of all ROIs, index out cell type
    in_roi_idx = roi_classification(in_rois); % of all ROIs, index out cell type
    log_idx_pc = pc_roi_idx == 1;
    log_idx_in = in_roi_idx == 1;
    pc_rois    = pc_rois(log_idx_pc);
    in_rois    = in_rois(log_idx_in);
end
% try
%     [pc_rois, in_rois] = remove_cells(sData);
%     roi_idx = sData.imdata.roi_classification(pc_rois);
%     roi_idx = roi_idx == 1;
%     pc_rois = pc_rois(roi_idx);
% catch
%     [pc_rois, in_rois] = remove_cells(sData);
% end
n_rois_to_plot = varargin{1,2};
n_frames       = size(sData.imdata.roiSignals(2).newdff,2);
% If this is a running session, also plot wheel position
if isfield(sData, 'behavior.wheelPosDs')
    n_rois_to_plot = n_rois_to_plot + 1;

    if n_frames > size(sData.behavior.wheelPosDs,1)
         temp                                      = NaN(1,n_frames);
         temp(1:length(sData.behavior.wheelPosDs)) = sData.behavior.wheelPosDs;
         sData.behavior.wheelPosDs                 = temp;
    elseif n_frames < size(sData.behavior.wheelPosDs,1)
        sData.behavior.wheelPosDs = sData.behavior.wheelPosDs(1:n_frames);
    end
end

switch varargin{1,3}
    case 'pc'
    rois_to_plot = randsample(size(pc_rois,2), n_rois_to_plot);
    dff          = sData.imdata.roiSignals(2).newdff(pc_rois,:);
    deconv       = sData.imdata.roiSignals(2).ciaDeconvolved(pc_rois,:);
%     dff          = dff(roi_idx,:);
    case 'in'
    rois_to_plot = randsample(size(in_rois,2), n_rois_to_plot);
    dff          = sData.imdata.roiSignals(2).newdff(in_rois,:);
    deconv       = sData.imdata.roiSignals(2).ciaDeconvolved(in_rois,:);

%     dff          = dff(roi_idx,:);
    case 'axon'
    dff          = sData.imdata.roiSignals(2).mergedAxonsDffFilt;
    rois_to_plot = randsample(size(dff,1), n_rois_to_plot);

end

figure,

time_imaging  = linspace(1, size(dff,2), size(dff,2) )/31;
max_lim       = max( max(dff(rois_to_plot,:)));
min_lim       = min( min(dff(rois_to_plot,:)));
dff_scale_bar = .5; % 1 / 100% DF/F 
grayColor     = [.2 .2 .2];
% timescale = dsearchn(time_imaging', [0 30]');
for i = 1:n_rois_to_plot

    roi = rois_to_plot(i);
    % Plot scale bars
    if i == n_rois_to_plot 
        h(i) = subplot( n_rois_to_plot ,1,i);
        plot(time_imaging, ones(1, size(dff,2)), 'Color', [1 1 1] ), hold on
        plot([0 0], [dff_scale_bar 0], 'LineWidth',3, 'Color',[0 0 0]), hold on
        plot([0 30], [0 0], 'LineWidth', 3, 'Color',[0 0 0])
        set(gca,'visible','off', 'ylim', [min_lim max_lim], 'xlim', [time_imaging(1) time_imaging(end)])
    elseif isfield(sData, 'behavior.wheelPosDs') && i == n_rois_to_plot-1        
        h(i) = subplot( n_rois_to_plot ,1,i);
        plot(time_imaging, sData.behavior.wheelPosDs, 'Color', grayColor)
        set(gca,'visible','off', 'xlim', [time_imaging(1) time_imaging(end)])
        text(-85, 80, 0, ['Position on wheel'])
    else
    h(i) = subplot( n_rois_to_plot ,1,i);
    plot(time_imaging, dff(roi,:), 'Color', grayColor)
    text(-55, 0.2, 0, ['ROI ', num2str(roi)])
    
    hold on
    plot(time_imaging, (deconv(roi,:)*.5)-.2,'r')

    set(gca,'visible','off', 'xlim', [time_imaging(1) time_imaging(end)],...
        'ylim', [min_lim max_lim])
    end

end
linkaxes(h, 'x')