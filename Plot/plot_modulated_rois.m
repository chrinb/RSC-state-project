function plot_modulated_rois(sorted_rois)

% Written by Christoffer Berge | Vervaeke Lab

% Function that plots modulated ROIs from "event_mod_an" function by
% checking whether input contains activated and/or suppressed ROIs and
% plotting the data + mean & standard error.

if size(sorted_rois) == [1,1]
    sorted_rois = sorted_rois{1,1};
end


%% Assign variables
ROIs_activated_dff         = sorted_rois{1,1};
ROIs_activated_dff_zscore  = sorted_rois{1,2};
ROIs_activated_dec         = sorted_rois{1,3};
ROIs_activated_dec_zscore  = sorted_rois{1,4};
ROIs_activated_idx         = sorted_rois{1,5};
ROIs_suppressed_dff        = sorted_rois{1,6};
ROIs_suppressed_dff_zscore = sorted_rois{1,7};
ROIs_suppressed_dec        = sorted_rois{1,8};
ROIs_suppressed_dec_zscore = sorted_rois{1,9};
ROIs_suppressed_idx        = sorted_rois{1,10};

time_vec = (-(31*3):(31*3))./31;
y_label      = 'ROI #';
y_label_dff  = 'Mean DF/F';
y_label_dffZ = 'Mean DF/F z-score';
y_label_dec  = 'Mean deconv. DF/F';
y_label_decZ = 'Mean deconv. DF/F z-score';
x_label      = 'Time from SWR peak (sec)';
c_lim_dff    = [-.5 .5];
c_lim_dec    = [-.05 .5];
%% Calulcate mean and SE of modulated ROIs
if ~isempty(ROIs_activated_idx) && size(ROIs_activated_idx,1) > 1
    ROIs_activated_dff_mean            = mean(ROIs_activated_dff,'omitnan'); 
    ROIs_activated_dff_zscore_mean     = mean(ROIs_activated_dff_zscore,'omitnan'); 
    ROIs_activated_dec_mean        = mean(ROIs_activated_dec,'omitnan'); 
    ROIs_activated_dec_zscore_mean = mean(ROIs_activated_dec_zscore,'omitnan'); 
    
    ROIs_aDFF_mean_SE  = std(ROIs_activated_dff,[], 1)./ sqrt(size(ROIs_activated_dff,1));
    ROIs_aDFFz_mean_SE = std(ROIs_activated_dff_zscore,[], 1)./ sqrt(size(ROIs_activated_dff_zscore,1));
    ROIs_aDec_mean     = std(ROIs_activated_dec,[], 1)./ sqrt(size(ROIs_activated_dec_zscore,1));
    ROIs_aDecz_mean    = std(ROIs_activated_dec_zscore,[], 1)./ sqrt(size(ROIs_activated_dec_zscore,1));
end

if ~isempty(ROIs_suppressed_idx) && size(ROIs_suppressed_idx,1) > 1
    ROIs_suppressed_dff_mean            = mean(ROIs_suppressed_dff,'omitnan'); 
    ROIs_suppressed_dff_zscore_mean     = mean(ROIs_suppressed_dff_zscore,'omitnan'); 
    ROIs_suppressed_dec_mean        = mean(ROIs_suppressed_dec,'omitnan'); 
    ROIs_suppressed_dec_zscore_mean = mean(ROIs_suppressed_dec_zscore,'omitnan'); 
    
    ROIs_sDFF_mean_SE  = std(ROIs_suppressed_dff,[], 1)./ sqrt(size(ROIs_suppressed_dff,1));
    ROIs_sDFFz_mean_SE = std(ROIs_suppressed_dff_zscore,[], 1)./ sqrt(size(ROIs_suppressed_dff_zscore,1));
    ROIs_sDec_mean     = std(ROIs_suppressed_dec,[], 1)./ sqrt(size(ROIs_suppressed_dec_zscore,1));
    ROIs_sDecz_mean    = std(ROIs_suppressed_dec_zscore,[], 1)./ sqrt(size(ROIs_suppressed_dec_zscore,1));
end

%% Plot results
if ~isempty(ROIs_activated_idx) && size(ROIs_activated_idx,1) > 1
    x1 = [time_vec(1) time_vec(end)];
    y2 = 1:size(ROIs_activated_idx,1);

    figure, 
    sgtitle('Activated ROIs')

    subplot(3,4, [1 5])
    imagesc(x1, y2, ROIs_activated_dff)
    ylabel(y_label)
    xlabel(x_label)
    colorbar
    title('DF/F')


    subplot(3,4, [2 6])
    imagesc(x1, y2, ROIs_activated_dff_zscore)
    ylabel(y_label)
    xlabel(x_label)
    caxis(c_lim_dff)
    colorbar
    title('DF/F z-score')

    subplot(3,4, [3 7])
    imagesc(x1, y2, ROIs_activated_dec)
    ylabel(y_label)
    xlabel(x_label)
%     colormap(flipud(gray))
    colorbar
    title('Deconvolved DF/F')

    subplot(3,4, [4 8])
    imagesc(x1, y2, ROIs_activated_dec_zscore)
    ylabel(y_label)
    xlabel(x_label)
    caxis(c_lim_dec)
    colorbar
    title('Deconvolved DF/F z-score')

    H(1) = subplot(3,4, 9);
    shadedErrorBar(time_vec, ROIs_activated_dff_mean, ROIs_aDFF_mean_SE)
    ylabel(y_label_dff)
    xlabel(x_label)

    H(2) = subplot(3,4, 10);
    shadedErrorBar(time_vec, ROIs_activated_dff_zscore_mean, ROIs_aDFFz_mean_SE)
    ylabel(y_label_dffZ)
    xlabel(x_label)

    H(3) = subplot(3,4, 11);
    shadedErrorBar(time_vec, ROIs_activated_dec_mean, ROIs_aDec_mean)
    ylabel(y_label_dec)
    xlabel(x_label)

    H(4) = subplot(3,4, 12);
    shadedErrorBar(time_vec, ROIs_activated_dec_zscore_mean, ROIs_aDecz_mean)
    ylabel(y_label_decZ)
    xlabel(x_label)
    linkaxes(H,'x')
    set(gca, 'xlim',[time_vec(1) time_vec(end)])
end

if ~isempty(ROIs_suppressed_idx) && size(ROIs_suppressed_idx,1) > 1
    x1 = [time_vec(1) time_vec(end)];
    y2 = 1:size(ROIs_suppressed_idx,1);

    figure, 
    sgtitle('Suppressed ROIs')

    subplot(3,4, [1 5])
    imagesc(x1, y2, ROIs_suppressed_dff)
    ylabel(y_label)
    xlabel(x_label)
    colorbar
    title('DF/F')

    subplot(3,4, [2 6])
    imagesc(x1, y2, ROIs_suppressed_dff_zscore)
    ylabel(y_label)
    xlabel(x_label)
    caxis(c_lim_dec)
    colorbar
    title('DF/F z-score')

    subplot(3,4, [3 7])
    imagesc(x1, y2, ROIs_suppressed_dec)
    ylabel(y_label)
    xlabel(x_label)
%     colormap(flipud(gray))
    colorbar
    title('Deconvolved DF/F')

    subplot(3,4, [4 8])
    imagesc(x1, y2, ROIs_suppressed_dec_zscore)
    ylabel(y_label)
    xlabel(x_label)
    caxis(c_lim_dff)
    colorbar
    title('Deconvolved DF/F z-score')

    H(1) = subplot(3,4, 9);
    shadedErrorBar(time_vec, ROIs_suppressed_dff_mean, ROIs_sDFF_mean_SE)
    ylabel(y_label_dff)
    xlabel(x_label)

    H(2) = subplot(3,4, 10);
    shadedErrorBar(time_vec, ROIs_suppressed_dff_zscore_mean, ROIs_sDFFz_mean_SE)
    ylabel(y_label_dffZ)
    xlabel(x_label)

    H(3) = subplot(3,4, 11);
    shadedErrorBar(time_vec, ROIs_suppressed_dec_mean, ROIs_sDec_mean)
    ylabel(y_label_dec)
    xlabel(x_label)

    H(4) = subplot(3,4, 12);
    shadedErrorBar(time_vec, ROIs_suppressed_dec_zscore_mean, ROIs_sDecz_mean)
    ylabel(y_label_decZ)
    xlabel(x_label)

    linkaxes(H,'x')
    set(gca, 'xlim',[time_vec(1) time_vec(end)])
end

