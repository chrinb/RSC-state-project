function rois_to_keep = find_static_grid_rois(RoiSignals_MeanF, threshold)

% Written by Christoffer Berge | Vervaeke lab

%{
Take mean grid ROI fluorescence from channel 1 and determine which ROI has
fluorescence changes exceeding threshold
%}

F = RoiSignals_MeanF';

zf = zscore(F, 0 ,2);

zftrend = detrend(zf);

% frames = sData.daqdata.frame_onset_reference_frame;
% rem_times = rem_sleep(sData);
% rem_times = frames(rem_times);

% data_rem         = zf(:, rem_times(1):rem_times(2));
% data_outside_rem = zf(:, [1:rem_times(1), rem_times(2):end]) ;
%% Find grid ROIS to keep
n_rois = size(F,1);

% threshold = 1.5;
upper_threshold = threshold;
lower_threshold = -threshold;

potential_rois_to_keep = zeros( 1, n_rois);
for roi_nr = 1:n_rois

    data = zftrend(roi_nr,:);
    min_val = min(data);
    max_val = max(data);

    if min_val > lower_threshold && max_val < upper_threshold

        potential_rois_to_keep(roi_nr) = roi_nr;
    end
end

rois_to_keep = potential_rois_to_keep(potential_rois_to_keep > 0);
tpm = 1:n_rois;
worst_rois   = tpm(potential_rois_to_keep == 0);

figure, 

sgtitle(['Threshold = ' num2str(threshold), ' std'])
subplot(231)
imagesc(F(rois_to_keep,:))
ylabel('ROI #', FontSize=16)
title('Mean fluorescence')

subplot(232)
imagesc(zf(rois_to_keep,:))
caxis([-1 1 ])
title('z-score mean F')

subplot(233)
imagesc(zftrend(rois_to_keep,:))
caxis([-1 1 ])
title('z-score mean F detrended')


subplot(234)
plot(mean(F(rois_to_keep,:)))

subplot(235)
plot(mean(zf(rois_to_keep,:)))

subplot(236)
plot(mean(zftrend(rois_to_keep,:)))

% figure, 
% histogram(zftrend(10,:))
% xline(threshold, 'r--', 'LineWidth',2)
% xline(-threshold, 'r--', 'LineWidth',2)
% title(['ROI 10, threhold ', num2str(threshold), ' STD'],FontSize=16)
% rois_to_keep