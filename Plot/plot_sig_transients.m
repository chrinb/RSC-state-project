function plot_sig_transients(varargin)

% Written by Christoffer Berge | Vervaeke lab

% Plot a random subset of ROI activity patterns and their significant
% transients. User inputs sData, number of ROIs to plot, and signal type
% (e.g., 'pc' for principal cells, 'in' for inhibitory cells, 'axon' for
% axonal data).

% Load data
sData              = varargin{1,1};
% dff                = sData.imdata.roiSignals(2).newdff;

% Find excitatory/inhibitory indicies 
[pc_rois, in_rois] = remove_cells_longitudinal(sData);

% try
%     [pc_rois, in_rois] = remove_cells(sData);
%     roi_idx = sData.imdata.roi_classification(pc_rois);
%     roi_idx = roi_idx == 1;
%     pc_rois = pc_rois(roi_idx);
% catch
%     [pc_rois, in_rois] = remove_cells(sData);
% end
n_rois_to_plot     = varargin{1,2};
% time_imaging       = linspace(1, size(dff,2), size(dff,2) )/31;

switch varargin{1,3}
    case 'pc'
    rois_to_plot = randsample(size(pc_rois,2), n_rois_to_plot);
    dff          = sData.imdata.roiSignals(2).newdff(pc_rois,:);
    transients   = sData.imdata.roiSignals(2).all_sig_transients(pc_rois,:);
%     dff          = dff(roi_idx,:);
    case 'in'
    rois_to_plot = randsample(size(in_rois,2), n_rois_to_plot);
    dff          = sData.imdata.roiSignals(2).newdff(in_rois,:);
    transients   = sData.imdata.roiSignals(2).all_sig_transients(in_rois,:);

%     dff          = dff(roi_idx,:);
    case 'axon'
    dff          = sData.imdata.roiSignals(2).mergedAxonsDffFilt;
    rois_to_plot = randsample(size(dff,1), n_rois_to_plot);
%     transients   = sData.imdata.roiSignals(2).all_sig_transients(pc_rois,:);

end


figure,

time_imaging  = linspace(1, size(dff,2), size(dff,2) )/31;
max_lim       = max( max(dff(rois_to_plot,:)));
min_lim       = min( min(dff(rois_to_plot,:)));
dff_scale_bar = 1; % 1 / 100% DF/F 
grayColor     = [.6 .6 .6];
% timescale = dsearchn(time_imaging', [0 30]');
for i = 1:n_rois_to_plot       

    roi = rois_to_plot(i);
    % Plot scale bars
    if i == n_rois_to_plot
        h(i) = subplot( n_rois_to_plot ,1,i);
        plot(time_imaging, ones(1, size(dff,2)), 'Color', [1 1 1] ), hold on
        plot([0 0], [dff_scale_bar 0], 'LineWidth',3, 'Color',[0 0 0]), 
        plot([0 30], [0 0], 'LineWidth', 3, 'Color',[0 0 0])
        set(gca,'visible','off', 'ylim', [min_lim max_lim])
    else
    h(i) = subplot( n_rois_to_plot ,1,i);
    plot(time_imaging, dff(roi,:), 'Color', grayColor)
    hold on
    transient = transients(roi,:);
    transient(transient == 0) = NaN;
    plot(time_imaging, transient,'r')
    text(-55, 0.2, 0, ['ROI ', num2str(roi)])
    set(gca,'visible','off', 'xlim', [time_imaging(1) time_imaging(end)],...
        'ylim', [min_lim max_lim])
    end

end
linkaxes(h, 'x')