function [rois_to_keep_peaks, zf, zftrend, roi_pdf] = find_static_grid_rois(sessionObjects, threshold)

% Written by Christoffer Berge | Vervaeke lab

%{
Take mean grid ROI fluorescence from channel 1 and determine which ROI has
fluorescence changes exceeding threshold
%}
F   = sessionObjects.loadData('RoiSignals_MeanF');

F = table2array(F)';

zf = zscore(F, 0 ,2);

zftrend = detrend(zf);

% frames = sData.daqdata.frame_onset_reference_frame;
% rem_times = rem_sleep(sData);
% rem_times = frames(rem_times);

% data_rem         = zf(:, rem_times(1):rem_times(2));
% data_outside_rem = zf(:, [1:rem_times(1), rem_times(2):end]) ;

%% Strategy 1: Apply threshold to the distribution of fluorescence values per ROI

% Some ROIs have noticable changes in mean fluorescence during REM. This
% somtimes result in activity distributions that are skewed. 
n_rois = size(F,1);

% threshold = 1.5;
upper_threshold = threshold;
lower_threshold = -threshold;

potential_rois_to_keep = zeros( 1, n_rois);
for roi_nr = 1:n_rois

    data = zf(roi_nr,:);
    min_val = min(data);
    max_val = max(data);

    if min_val > lower_threshold && max_val < upper_threshold

        potential_rois_to_keep(roi_nr) = roi_nr;
    end
end

rois_to_keep_meanF = potential_rois_to_keep(potential_rois_to_keep > 0);
% tpm = 1:n_rois;
% worst_rois   = tpm(potential_rois_to_keep == 0);



%% Strategy 2: Keep ROIs with one peak in distribution

% Create kernel PDF for each ROI
[roi_pdf, roi_pdf_x, peaks] = deal( cell( n_rois, 1));
for roi_nr = 1:n_rois

    data_min = min(zf(roi_nr,:));
    data_max = max(zf(roi_nr,:));

    x = linspace(data_min, data_max);
    tmp_kernel = fitdist(zf(roi_nr,:)', 'kernel');
    roi_pdf_x{roi_nr,1} = x;
    roi_pdf{roi_nr,1} = pdf(tmp_kernel, x);

%     peaks{roi_nr} = findpeaks(roi_pdf{roi_nr,1}, 'MinPeakHeight', 0.15, 'MinPeakProminence', .02);
    peaks{roi_nr} = findpeaks(roi_pdf{roi_nr,1});

end

%% Plot example roi
roi = 40;
figure, 
subplot(131)
plot(zf(roi,:))
axis square
subplot(132)
histogram(zf(roi,:))
axis square
title(['ROI ', num2str(roi)], FontSize=16)
subplot(133)
findpeaks(roi_pdf{roi,1})
axis square
%% Calc
is_smaller_than = @(z) numel(z)<2;
keep_1peak_rois_idx = cellfun(is_smaller_than, peaks);
n_rois_list = 1:n_rois;
rois_to_keep_peaks = n_rois_list(keep_1peak_rois_idx);
%% Plot
figure, 

sgtitle(['Threshold = ' num2str(threshold), ' std'])
subplot(231)
imagesc(F(rois_to_keep_meanF,:))
ylabel('ROI #', FontSize=16)
title('Mean fluorescence')

subplot(232)
imagesc(zf(rois_to_keep_meanF,:))
caxis([-1 1 ])
title('z-score mean F')

subplot(233)
imagesc(zftrend(rois_to_keep_meanF,:))
caxis([-1 1 ])
title('z-score mean F detrended')


subplot(234)
plot(mean(F(rois_to_keep_meanF,:)))

subplot(235)
plot(mean(zf(rois_to_keep_meanF,:)))

subplot(236)
plot(mean(zftrend(rois_to_keep_meanF,:)))

% figure, 
% histogram(zftrend(10,:))f
% xline(threshold, 'r--', 'LineWidth',2)
% xline(-threshold, 'r--', 'LineWidth',2)
% title(['ROI 10, threhold ', num2str(threshold), ' STD'],FontSize=16)
% rois_to_keep