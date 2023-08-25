function plot_sample_signals(varargin)

% Written by Christoffer Berge | Vervaeke lab

% Plot a random subset of ROI activity patterns (DF/F and deconvolved DFF).
% User inputs sData, number of ROIs to plot, and signal type
% (e.g., 'pc' for principal cells, 'in' for inhibitory cells, 'axon' for
% axonal data).

% Load data
sData              = varargin{1,1};
% dff                = sData.imdata.roiSignals(2).newdff;
try
    [pc_rois, in_rois] = remove_cells(sData.imdata.roi_arr);
catch
    [pc_rois, in_rois] = remove_cells(sData.imdata.roiArr);
end
n_rois_to_plot     = varargin{1,2};
% time_imaging       = linspace(1, size(dff,2), size(dff,2) )/31;

switch varargin{1,3}
    case 'pc'
    rois_to_plot = randsample(pc_rois, n_rois_to_plot);
    dff          = sData.imdata.roiSignals(2).newdff;
    deconv       = sData.imdata.roiSignals(2).ciaDeconvolved;
    case 'in'
    rois_to_plot = randsample(in_rois, n_rois_to_plot);
    dff          = sData.imdata.roiSignals(2).newdff;
    deconv       = sData.imdata.roiSignals(2).ciaDeconvolved;
    case 'axon'
    dff          = sData.imdata.roiSignals(2).mergedAxonsDffFilt;
    deconv       = sData.imdata.roiSignals(2).mergedAxonsDec;
    rois_to_plot = randsample(size(dff,1), n_rois_to_plot);

end


figure,

time_imaging  = linspace(1, size(dff,2), size(dff,2) )/62;
max_lim       = max( max(dff(rois_to_plot,:)));
deconv_multiplication_factor = 1;
min_lim       = min ((deconv(1,:)-.4)*deconv_multiplication_factor);
dff_scale_bar = 1; % 1 / 100% DF/F 
blackColor     = [0 0 0];
fontsize       = 12;
% timescale = dsearchn(time_imaging', [0 30]');
for i = 1:n_rois_to_plot

    roi = rois_to_plot(i);
    % Plot scale bars
    if i == n_rois_to_plot
        h(i) = subplot( n_rois_to_plot ,1,i);
        plot(time_imaging, ones(1, size(dff,2)), 'Color', [1 1 1] ), hold on
        plot([0 0], [dff_scale_bar 0], 'LineWidth',3, 'Color',[0 0 0]), hold on
        plot([0 30], [0 0], 'LineWidth', 3, 'Color',[0 0 0])
        set(gca,'visible','off', 'ylim', [min_lim max_lim])
        text(1, 2, [num2str(dff_scale_bar*100), '% DF/F'],'HorizontalAlignment', 'right','FontSize',fontsize, 'Position',[-5 .2 0]);
        text(1, 2, '30 s','HorizontalAlignment', 'right','FontSize',fontsize, 'Position',[25 min_lim-.01 0]);
    else
        h(i) = subplot( n_rois_to_plot ,1,i);
        plot(time_imaging, dff(roi,:), 'Color', blackColor)
        hold on
        plot(time_imaging, (deconv(roi,:)-.4)*deconv_multiplication_factor,'r', 'LineWidth',1)
    
        set(gca,'visible','of', 'xlim', [time_imaging(1) time_imaging(end)],...
            'ylim', [min_lim max_lim])
    end

end
linkaxes(h, 'x')