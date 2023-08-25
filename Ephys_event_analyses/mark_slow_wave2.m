function nrem_delta_waves = mark_slow_wave2(varargin)

% Written by Christoffer Berge | Vervaeke Lab

% Alternative code to "slow_wave_an" to find delta waves based on method
% outlined in Maingret et al. (2016) "Hippocampo-cortical coupling mediates 
% memory consolidation during sleep".

sData = varargin{1,1};
%% Delta waves
signal             = sData.ephysdata2.lfp';
signal_filt        = sData.ephysdata2.deltaband';
signal_filt_zscore = zscore(signal_filt);
time               = (0:length(signal)-1)/2500;
signal_prime       = diff(signal_filt_zscore);
time_prime         = time(1:end-1);

nrem_times = nrem_sleep(sData);
raw_nrem_bouts = [];

% Zero-crossings
DownZCi = @(signal_prime) find(signal_prime(1:end-1) >= 0 & signal_prime(2:end) < 0);    % Returns Down Zero-Crossing Indices
dzc = DownZCi(signal_prime);
% dzc_alt(signal_prime(dzc_alt) > 0.02) = [];

% % Find negative-to-positive zero crossings
UpZCi = @(signal_prime) find(signal_prime(1:end-1) <= 0 & signal_prime(2:end) > 0);	
yzc = UpZCi(signal_prime);
% yzc_alt(signal_filt(yzc_alt) < -0.02) = [];


% start with first trough after first peak
% dzc = dzc(dzc > peak_pos_times(1)  );
slow_wave_mat = zeros(length(yzc),3);
downward_zc = yzc;
upward_zc   = dzc;

for i = 1:length(yzc)
    try
        slow_wave_mat(i,1) = upward_zc(i);
        slow_wave_mat(i,3) = upward_zc(i+1);
        temp_id = downward_zc > slow_wave_mat(i,1) & downward_zc < slow_wave_mat(i,3);
        slow_wave_mat(i,2) = downward_zc(temp_id);
    end
end

% Min/max duration of SO = 150ms - 500ms
min_duration = 375;
max_duration = 1250;

delta_duration = slow_wave_mat(:,3)-slow_wave_mat(:,1);
delta_keep_idx = delta_duration > min_duration & delta_duration < max_duration;
putative_delta_waves = slow_wave_mat(delta_keep_idx,:);

delta_peaks = putative_delta_waves(:,2);
delta_end   = putative_delta_waves(:,3);
peaks_above_threshold = signal_filt_zscore(delta_peaks) < -2;
peaks_above_threshold2 = signal_filt_zscore(delta_peaks) < -1 & signal_filt_zscore(delta_end) < 1.5;

test = [peaks_above_threshold; peaks_above_threshold2];
sum_test = logical( sum(test,1));

delta_waves_to_keep = putative_delta_waves(sum_test, :);

%% Discard slow osillations/delta waves outside of NREM bouts
% nrem_SOs = zeros(size(SOs,1),1);
nrem_delta_idx = zeros( size(delta_waves_to_keep,1), length(nrem_times));

for i = 1:length(nrem_times)
    nrem_bout           = nrem_times(i,1):nrem_times(i,2);
    nrem_delta_idx(:,i) = ismember(delta_waves_to_keep(:,3), nrem_bout); 
end

% Sum rows to obtain all NREM SO indicies and convert to logicals
nrem_delta_idx = logical( sum(nrem_delta_idx,2)); 

nrem_delta_waves = delta_waves_to_keep(nrem_delta_idx,:);

%% Plot results
NREM_start_end = nrem_times./2500;
ECoG_spindle      = sData.ephysdata2.NREMspindleStartEnd/2500;

if length(varargin) > 1
    
    figure,
    plot(time ,signal_filt)
    hold on
    plot(time,signal*.4-.2)
    

    for i = 1:length(nrem_delta_waves)
        x = [nrem_delta_waves(i,1) nrem_delta_waves(i,1) nrem_delta_waves(i,3) nrem_delta_waves(i,3)]/2500;
        y = [-1 1 1 -1];
        patch(x, y, 'magenta', 'edgecolor', 'none', 'FaceAlpha', .3);
    end

    for k = 1:length(sData.ephysdata2.NREMspindleStartEnd)
        z = [ECoG_spindle(k,1) ECoG_spindle(k,1) ECoG_spindle(k,2) ECoG_spindle(k,2) ];
        v = [-1 1 1 -1];
        patch(z,v, 'red', 'edgecolor', 'none', 'FaceAlpha', .3);
    end

    %NREM bouts
    for i = 1:length(NREM_start_end)
        x = [NREM_start_end(i,1) NREM_start_end(i,1) NREM_start_end(i,2) NREM_start_end(i,2)];
        y = [-1 1 1 -1];
        patch(x, y, 'black', 'edgecolor', 'none', 'FaceAlpha', .2);
    end
     set(gca, 'ylim', [-.6 .4])
end

% 
% %% Plot
% figure, plot(time_prime, signal_prime)
% hold on
% plot(time_prime(dzc), signal_prime(dzc), 'bp')
% plot(time_prime(yzc), signal_prime(yzc), 'ro')
% 
% 
% figure, plot(signal_filt)
% hold on
% plot(signal_prime*6)