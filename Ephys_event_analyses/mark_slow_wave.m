function [nrem_SOs, nrem_delta_waves] = mark_slow_wave(varargin)

% Written by Christoffer Berge | Vervaeke Lab

% Find slow oscillations & delta waves. Code is mainly based on method
% outlined in Kim et al. (2019) "Competing Roles of Slow Oscillations and 
% Delta Waves in Memory Consolidation versus Forgetting". 

sData = varargin{1,1};
%% Delta waves
ecog               = sData.ephysdata2.lfp';
ecog_delta         = sData.ephysdata2.deltaband';
ecog_sigma         = sData.ephysdata2.sigmaband;
ecog_filt_zscore   = zscore(ecog_delta);
time               = (0:length(ecog)-1)/2500;
signal_prime       = diff(ecog_filt_zscore);

nrem_times = nrem_sleep(sData);
raw_nrem_bouts = [];

%% Find threshold
filt_nrem_bouts = [];

if size(nrem_times,1) > 1
    for i = 1:length(nrem_times)
        filt_nrem_bouts = horzcat(filt_nrem_bouts, ecog_delta(nrem_times(i, 1):...
            nrem_times(i,2)));
        raw_nrem_bouts = horzcat(raw_nrem_bouts, ecog_delta(nrem_times(i, 1):...
            nrem_times(i,2)));
    end
else
    i               = 1;
    filt_nrem_bouts = ecog_delta(nrem_times(i, 1):nrem_times(i,2));
end


% Find all peaks and troughs in filtered signal
[neg_peaks, ~] = islocalmin(filt_nrem_bouts); 
[val_pos_peaks, pos_peaks_locs] = findpeaks(filt_nrem_bouts);
% test1 = zeros(1,length(filt_nrem_bouts));
% test1(pos_peaks_locs) = 0.1;

val_neg_peaks = filt_nrem_bouts(neg_peaks);

% Define the positive threshold as the top 15 percentile of
% the peaks and the negative threshold as the bottom 40 percentile of the
% troughs
thresh_up  = prctile(val_pos_peaks, 85);
thresh_low = prctile(val_neg_peaks, 40);

%% 
nrem_time_vec = (0:length(ecog_delta)-1)/2500;

% Find all peaks and troughs in filtered signal
[neg_peaks, ~] = islocalmin(ecog_delta); 
[val_pos_peaks, pos_peaks_locs] = findpeaks(ecog_delta);
% test1 = zeros(1,length(signal_filt));
% test1(pos_peaks_locs) = 0.1;

val_neg_peaks = ecog_delta(neg_peaks);

peak_idx_pos = ismember(ecog_delta, val_pos_peaks );
peak_idx_neg = ismember(ecog_delta, val_neg_peaks);

peak_pos_times = find(peak_idx_pos);
peak_neg_times = find(peak_idx_neg);

%% Find zero-crossings

% Find positive-to-negative zero crossings
DownZCi                      = @(signal_filt) find(signal_filt(1:end-1) >= 0 & signal_filt(2:end) < 0);    % Returns Down Zero-Crossing Indices
dzc                          = DownZCi(ecog_delta);
dzc(ecog_delta(dzc) > 0.02) = [];
% % Find negative-to-positive zero crossings
UpZCi                         = @(signal_filt) find(signal_filt(1:end-1) <= 0 & signal_filt(2:end) > 0);	
yzc                           = UpZCi(ecog_delta);
yzc(ecog_delta(yzc) < -0.02) = [];

%%
% start with first trough after first peak
% dzc = dzc(dzc > peak_pos_times(1)  );
slow_wave_mat = zeros(length(dzc),4);
for i = 1:length(dzc)
    
    % Remove peaks after, and troughs before, current positive-to-negative zero crossing  
    preceding_peaks = peak_pos_times(peak_pos_times < dzc(i));
    succeeding_troughs = peak_neg_times(peak_neg_times > dzc(i));
    % Remove negative-to-positive zero crossings before and after current
    % zero-crossing in two different vectors- 
    preceding_neg_to_pos_zc  = yzc(yzc < dzc(i));
    succeeding_neg_to_pos_zc = yzc(yzc > dzc(i));
    
    % Find index of peak preceding current trough
    [~, preceding_zc_idx]      = min(abs(dzc(i) - preceding_neg_to_pos_zc));
    [~, preceding_peak_idx]    = min(abs( dzc(i) - preceding_peaks));
    [~, succeeding_trough_idx] = min(abs(succeeding_troughs - dzc(i) ));
    [~, succeeding_zc_idx]     = min(abs(succeeding_neg_to_pos_zc - dzc(i) ));
    
    try
        slow_wave_mat(i, 1) = preceding_neg_to_pos_zc(preceding_zc_idx);
        slow_wave_mat(i, 2) = preceding_peaks(preceding_peak_idx);
        slow_wave_mat(i, 3) = succeeding_troughs(succeeding_trough_idx);
        slow_wave_mat(i, 4) = succeeding_neg_to_pos_zc(succeeding_zc_idx);
    end
end

% remove rows with zeros (typically first or last negative-to-positive zero
% crossing where one or more events are missing)
slow_wave_mat = slow_wave_mat(all(slow_wave_mat,2),:);


% remove peaks and troughs below and above respective thresholds
peaks_above_up_thresh_idx  = ecog_delta( slow_wave_mat(:,2)) > thresh_up;
peaks_below_low_thresh_idx = ecog_delta( slow_wave_mat(:,3)) < thresh_low;

% Store remaining peaks and add/subtract to find SO/delta waves indicies
peaks_mat = [peaks_above_up_thresh_idx', peaks_below_low_thresh_idx'];
temp_mat1 = peaks_mat(:,1) + peaks_mat(:,2);
temp_mat2 = peaks_mat(:,1) - peaks_mat(:,2);

SO_idx         = temp_mat1 > 1;
delta_wave_idx = temp_mat2 == -1;

% Min/max duration of SO & delta waves = 150ms - 500ms
min_duration   = 375;
max_duration   = 1250;
putative_SOs   = slow_wave_mat(SO_idx,:);
putative_delta = slow_wave_mat(delta_wave_idx,:);
% check duration between SO peak and trough
SO_duration    = putative_SOs(:,3)-putative_SOs(:,2);
Delta_duration = putative_delta(:,3)-putative_delta(:,2);

SO_keep_idx    = SO_duration    > min_duration & SO_duration    < max_duration;
Delta_keep_idx = Delta_duration > min_duration & Delta_duration < max_duration;
SOs            = putative_SOs(SO_keep_idx,:);
delta_waves    = putative_delta(Delta_keep_idx,:);

%% Discard slow osillations/delta waves outside of NREM bouts
% nrem_SOs = zeros(size(SOs,1),1);
nrem_SOs_idx   = zeros( size(SOs,1), length(nrem_times));
nrem_delta_idx = zeros( size(delta_waves,1), length(nrem_times));

if size(nrem_times, 1) > 1
    for i = 1:length(nrem_times)
        nrem_bout           = nrem_times(i,1):nrem_times(i,2);
        nrem_SOs_idx(:,i)   = ismember(SOs(:,3), nrem_bout);
        nrem_delta_idx(:,i) = ismember(delta_waves(:,3), nrem_bout); 
    end
else
    i                   = 1;
    nrem_bout           = nrem_times(i,1):nrem_times(i,2);
    nrem_SOs_idx(:,i)   = ismember(SOs(:,3), nrem_bout);
    nrem_delta_idx(:,i) = ismember(delta_waves(:,3), nrem_bout); 
end

% Sum rows to obtain all NREM SO indicies and convert to logicals
nrem_SOs_idx   = logical( sum(nrem_SOs_idx,2));
nrem_delta_idx = logical( sum(nrem_delta_idx,2)); 
% nrem_SOs_idx1 = repmat(nrem_SOs_idx,1, 4);
nrem_SOs         = SOs(nrem_SOs_idx,:);
nrem_delta_waves = delta_waves(nrem_delta_idx,:);

%% Plot results
NREM_start_end = nrem_times./2500;


% Remove slow waves inside spindles for better visualization
% Remove detected slow waves that are inside sleep spindles
ECoG_spindle          = sData.ephysdata2.NREMspindleStartEnd;
spindle_center        = round(sData.ephysdata2.NREMAbsSpindleIdx);
delta_wave_start_stop = [nrem_delta_waves(:,2), nrem_delta_waves(:,4)];
so_start_stop         = [nrem_SOs(:,2), nrem_SOs(:,4)];

log_idx_delta = [];
log_idx_so    = [];
for i = 1:length(spindle_center)
   log_idx_delta = horzcat(log_idx_delta, delta_wave_start_stop(:,1) > ECoG_spindle(i,1) & ...
       delta_wave_start_stop(:,2) < ECoG_spindle(i,2) );
    log_idx_so = horzcat(log_idx_so, so_start_stop(:,1) > ECoG_spindle(i,1) & ...
       so_start_stop(:,2) < ECoG_spindle(i,2) );
end
delta_wave_removed       = ~logical(sum(log_idx_delta,2));
so_removed               = ~logical(sum(log_idx_so,2));
nrem_delta_waves_removed = nrem_delta_waves(delta_wave_removed,: );
nrem_so_removed          = nrem_SOs(so_removed,: );

% slow_wave = slow_wave(sortId,:); 


% Downsample ephys signals for less laggy visualization
duration_sec    = time(end) - time(1);
downsampled_n   = 250*duration_sec;
downsampled_idx = round(linspace(1,length(ecog),downsampled_n));
ecog_ds         = ecog(downsampled_idx);
ecog_delta_ds   = ecog_delta(downsampled_idx);
ecog_sigma_ds   = ecog_sigma(downsampled_idx);
ecog_time_ds    = time(downsampled_idx);

if length(varargin) > 1
    
    displacement_factor = max(ecog) + 0.1;
    ymax                = max(ecog);
    ymin                = min(ecog-3);  
    patch_lim           = [ymin, ymax, ymax, ymin];

    figure,
    plot(ecog_time_ds, ecog_ds)
%     yline(thresh_low)
%     yline(thresh_up)
    hold on
    plot(ecog_time_ds, ecog_delta_ds-displacement_factor)
    plot(ecog_time_ds, ecog_sigma_ds-displacement_factor-.5)
    % Loop over slow oscillations and plot them as individual patches
    for i = 1:length(nrem_so_removed)
        x = [nrem_so_removed(i,2) nrem_so_removed(i,2) nrem_so_removed(i,4) nrem_so_removed(i,4)]/2500;
        y = patch_lim;
        h1 = patch(x, y, 'magenta', 'edgecolor', 'none', 'FaceAlpha', .2);
    end
    
    % Loop over delta waves and plot them as individual patches
    for i = 1:length(nrem_delta_waves_removed)
        x = [nrem_delta_waves_removed(i,2) nrem_delta_waves_removed(i,2) nrem_delta_waves_removed(i,4) nrem_delta_waves_removed(i,4)]/2500;
        y = patch_lim;
        h2 = patch(x, y, 'cyan', 'edgecolor', 'none', 'FaceAlpha', .2);
    end
    
    % Loop over sleep spindles and plot them as individual patches
    for k = 1:length(sData.ephysdata2.NREMspindleStartEnd)
        z = [ECoG_spindle(k,1) ECoG_spindle(k,1) ECoG_spindle(k,2) ECoG_spindle(k,2) ]/2500;
        v = patch_lim;
        h3 = patch(z,v, 'red', 'edgecolor', 'none', 'FaceAlpha', .2);
    end

    % Loop over NREM bouts and plot them as individual patches
    for i = 1:length(NREM_start_end)
        x = [NREM_start_end(i,1) NREM_start_end(i,1) NREM_start_end(i,2) NREM_start_end(i,2)];
        y = patch_lim;
        h4 = patch(x, y, 'black', 'edgecolor', 'none', 'FaceAlpha', .1);
    end
     set(gca, 'ylim', [min(ecog_sigma_ds-displacement_factor-.5) ymax])
     legend( [h1, h2, h3, h4], 'SO','Delta', 'Spindle','NREM');
     set(gca,'xlim', [0 nrem_time_vec(end)])
     legend( [h1, h2, h3, h4], 'SO','Delta', 'Spindle','NREM');

end
% legend( [h1, h2, h3, h4], 'SO','Delta', 'Spindle','NREM');
% set(gca,'xlim', [0 nrem_time_vec(end)])

