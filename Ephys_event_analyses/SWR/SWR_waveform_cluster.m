%% Filter signal in 10-1000 Hz band
LFP = sData.ephysdata.lfp;
signalLength = length(LFP);
srate = 2500;
fs = 2500;
nyquist = srate/2;
transw = 0.1
shape = [0 0 1 1 0 0];
delta_t = 1/srate;
% tt = delta_t:delta_t:signalLeng0th*delta_t;
ripple_freq = [10 1000];
ripple_freq_shape = [0 ripple_freq(1)-ripple_freq(1)*transw, ripple_freq(1) ripple_freq(2), ripple_freq(2)+ripple_freq(2)*transw nyquist]/nyquist;
order_ripple = round( 50*srate/ripple_freq(1));
filtkern_ripple = fir1(order_ripple, ripple_freq/nyquist);
% filtkern_ripple = firls(order_ripple, ripple_freq_shape, [0 0 1 1 0 0]);
filt_ripple_pow = abs(fft(filtkern_ripple).^2);
hz_ripple = linspace(0, nyquist, floor(length(filtkern_ripple)/2)+1);
filtpow_ripple = filt_ripple_pow(1:length(hz_ripple));

figure,
plot(hz_ripple, filtpow_ripple, 'b', 'linew', 2), hold on
plot([0 ripple_freq(1) ripple_freq ripple_freq(2) nyquist], shape, 'r','linew', 2)
set(gca, 'xlim', [0, ripple_freq(2)+50])
xlabel('Frequency (Hz)'), ylabel('Filtergain')
legend({'Filterkernel'; 'Ideal filter'})
title('Ripple frequency')

ripple_filtsig_10_1000 = filtfilt(filtkern_ripple, 1, LFP);


for i = 1:length(sData.ephysdata.absRipIdx)
lfpPeakIdx = sData.ephysdata.absRipIdx(i);
lfpStartIdx = lfpPeakIdx(1) - (0.5*fs);
%if the ripple timepoint is near the beginning of the trace
if lfpStartIdx < 0; lfpStartIdx = 1; end
lfpEndIdx = lfpPeakIdx(1) + (0.5*fs);
%if the ripple timepoint is near the end of the trace
if lfpEndIdx > length(LFP); lfpEndIdx = length(LFP); end
rippleSnips(i).lfp = ripple_filtsig_10_1000(lfpStartIdx:lfpEndIdx);
rippleIdx(i) = lfpPeakIdx(1);
%     if runSignal(rippleIdx(i)) ~= 0; rippleIdx(i) = NaN; end %take out timepoints when animal is walking
% %     if runSignal(moving_mean_move(i)) ~= 0; rippleIdx(i) = NaN; end %take out timepoints when animal is walking
%
end

%% Stuff
n_swrs = length(sData.ephysdata.absRipIdx);

first_peak = zeros(1,n_swrs);
[dataZ_short, dataZ_short_raw] = deal( zeros( n_swrs, 401));
fs = 2500;

% 10-1000Hz filtered LFP data
for i = 1:n_swrs
    dataZ_short(i,:) = zscore( rippleSnips(i).lfp(1050:1450));
    [~,locs,~,~] = findpeaks( dataZ_short(i,:), fs, 'MinPeakHeight',0.5, 'MinPeakDistance',0.004 );
%     first_peak(i) = locs(1);
end

% Raw LFP data
for i = 1:n_swrs
    dataZ_short_raw(i,:) = zscore( sData.ephysdata.rippleSnips(i).lfp(1050:1450));
    [~,locs,~,~] = findpeaks( dataZ_short_raw(i,:), fs, 'MinPeakHeight',0.5, 'MinPeakDistance',0.004 );
%     first_peak(i) = locs(1);
end

swr_to_plot = 113;
figure, hold on
plot(dataZ_short(swr_to_plot,:))
plot(dataZ_short_raw(swr_to_plot,:))
legend('10-1000Hz', 'raw')
%% plot spectrogram
window  = 120;              % Window size for computing the spectrogram (FFT) [# samples]
overlap = 119;              % Overlap of the windows for computing the spectrogram [# samples]
nFFT    = 100:1:250;          % Vector defining the frequencies for computing the FFT
fs = 2500;
x = rippleSnips(6).lfp;
[~,F,T,P] = spectrogram(x,window,overlap,nFFT,fs);

time = linspace(0,1000,2501);

figure, 
subplot(411)
plot(time, x),

subplot(412),
s = surf(T,F, (abs(P)),'edgecolor','none');
          colormap(jet)
          view([0 0])
          
center = length(s.ZData)/2; 

subplot(413),
plot(mean(s.ZData,1)), xline(center+200), xline(center-200), xline(center)

subplot(414),
findpeaks( mean(s.ZData,1), fs, 'MinPeakHeight',6e-06), yline(6e-06)

[~,locscenter,~,~] = findpeaks( mean(s.ZData,1), fs, 'MinPeakHeight',6e-06);


centers = repmat( (center/2500), 1, length(locscenter));
compare = diff( [centers; locscenter]);
[~, id] = min(abs(compare));
centerpeak = locscenter(id);


%% plot peaks in filtered signal

filtered_SWRsnips      = zeros( n_swrs, 401);
filtered_SWRsnips_full = zeros( n_swrs, 2501);
time2                  = linspace(0, 0.16, 401);

for i = 1:n_swrs
    start = sData.ephysdata.absRipIdx(i) - 200;
    stop = sData.ephysdata.absRipIdx(i) + 200; 
    filtered_SWRsnips(i,:) = sData.ephysdata.ripplefreq(start:stop);
    filtered_SWRsnips_full(i,:) = sData.ephysdata.ripplefreq( (start-1050):(stop+1050));
end

swr = 139;
figure, 
subplot(311),
plot(time2, sData.ephysdata.rippleSnips(swr).lfp(1050:1450))

subplot(312),
plot(time2,filtered_SWRsnips(swr,:))

subplot(313),
findpeaks(filtered_SWRsnips(swr,:),  'MinPeakHeight',0.05, 'MinPeakDistance',0.004 );

[~,locs,~,~] = findpeaks( filtered_SWRsnips(swr,:), fs, 'MinPeakHeight',0.05, 'MinPeakDistance',0.004 );

%% find first peak

[first_peak, last_peak, locsum] = deal( zeros(1,n_swrs) );

for i = 1:length(sData.ephysdata.absRipIdx)
    [~,locs,~,~] = findpeaks( filtered_SWRsnips(i,:), fs, 'MinPeakHeight',0.05, 'MinPeakDistance',0.004 );
    locsum(i) = sum(diff(locs) > 0.02);
    first_peak(i) = locs(1);
    last_peak(i) = locs(end);
end

[min_peak, idmin] = min(first_peak); 
[max_peak, idmax] = max(last_peak);


%% find peaks 

[first_peak, last_peak] = deal( zeros(1,n_swrs) );
peaks_to_keep = cell(1, n_swrs);
% CHANGE TO 0.08 IF USING 'fs' IN findpeaks!!!!
window_center = 200;
% CHANGE TO 0.02 IF USING 'fs' IN findpeaks!!!!
peak_distance = 50;

for i = 1:n_swrs
    
    [~,locs,~,~] = findpeaks( filtered_SWRsnips(i,:), 'MinPeakHeight',0.05, 'MinPeakDistance',0.004 );
    locszero = zeros(1, length(locs));  
    [~, centerIdx] = min( abs(locs-window_center));
    centerIdx_copy = centerIdx; 
    locszero(centerIdx) = 1;
    
    while centerIdx > 1

        % find peak closest to center

        t = 1;
        if locs(centerIdx)-locs(centerIdx-t) < peak_distance 
            locszero(centerIdx-t) = 1;
            centerIdx = centerIdx - 1;

        else
            centerIdx = 0;

        end
    end

    while centerIdx_copy < length(locs)
        
        if locs(centerIdx_copy+ t) - locs(centerIdx_copy) < peak_distance
            locszero(centerIdx_copy+t) = 1;
            centerIdx_copy = centerIdx_copy + 1;
    %         if centerIdx_copy == length(locs)

        else
            centerIdx_copy = length(locs);

        end
    end

    peaks_to_keep{i} = locs(locszero > 0); 
    
    % timing of first peak in seconds in 400ms window
    first_peak(i) = min(peaks_to_keep{i});
    last_peak(i) = max(peaks_to_keep{i});
end

[align_to_peak, id]         = min(first_peak);
[last_peak_to_keep, idl]    = max(last_peak);

%% align SWR-snippets to earliest peak

% determine where the latest peak of a SWR occurs 
shift_for_latest_peak = first_peak(idl) - align_to_peak;
last_peak_shifted = last_peak(idl) - shift_for_latest_peak;

% initialize variables
shifts = zeros(1, n_swrs);

fig = figure; hold on
for i = 1:n_swrs
    shifts(i) = first_peak(i) - align_to_peak;
    plot(time2-shifts(i), filtered_SWRsnips(i,:))
    
end
xline(last_peak_shifted)

window = [align_to_peak, last_peak_shifted];

% fig = figure; hold on
% for i = 1:n_swrs
%     shifts(i) = first_peak(i) - align_to_peak;
%     plot(time2-shifts(i), sData.ephysdata.rippleSnips(i).lfp(1050:1450))
% end
% xline(last_peak_shifted)

%% find aligned time vectors (if using 'fs' in findpeaks)
cutoff = [0 0.1];

aligned_time_vec = zeros( n_swrs, 251); 
for i = 1:n_swrs
    time3 = time2 - shifts(i);
    [~, cutoffId1] = min( abs( time3 - cutoff(1) )); 
    [~, cutoffId2] = min( abs( time3 - cutoff(2) ));
    temp = time3(cutoffId1:cutoffId2);
    if length(temp) < 251
        zero_pad = zeros(1, 251);
       temp(numel(zero_pad)) = 0;
    end
    aligned_time_vec(i,:) = temp;
end

%% do some plotting

% 100-250 Hz filtered data
aligned_SWRs_filt = zeros( n_swrs, length(filtered_SWRsnips));

fig = figure; 
subplot(311), hold on
for i = 1:n_swrs
    aligned_SWRs_filt(i,:) = circshift( filtered_SWRsnips(i,:), -shifts(i), 2);
    plot( aligned_SWRs_filt(i,:))
end
xline(last_peak_shifted)
plot( mean(aligned_SWRs_filt,1), 'linew',2)

% SWR with latest peak:
% plot(aligned_SWRs_filt(99,:), 'linew',2)

% trim aligned SWRs
aligned_SWRs_trim1 = aligned_SWRs_filt(:, 1:237);

% 10-1000 Hz filtered data
aligned_SWRs = zeros( n_swrs, length(filtered_SWRsnips));

subplot(312); hold on
for i = 1:n_swrs
    aligned_SWRs(i,:) = circshift( dataZ_short(i,:), -shifts(i), 2);
    plot( aligned_SWRs(i,:))
end
xline(last_peak_shifted)
plot( mean(aligned_SWRs,1), 'linew',2)
% SWR with latest peak:
% plot(aligned_SWRs(99,:), 'linew',2)

% trim aligned SWRs
aligned_SWRs_trim2 = aligned_SWRs(:, 1:237);

% Raw data
aligned_SWRs_raw = zeros( n_swrs, length(filtered_SWRsnips));

subplot(313); hold on
for i = 1:n_swrs
    aligned_SWRs_raw(i,:) = circshift( dataZ_short_raw(i,:), -shifts(i), 2);
    plot( aligned_SWRs_raw(i,:))
end
xline(last_peak_shifted)
plot( mean(aligned_SWRs_raw,1), 'linew', 2, 'MarkerFaceColor', 'r')

% SWR with latest peak:
% plot(aligned_SWRs(99,:), 'linew',2)

% trim aligned SWRs
aligned_SWRs_trim3 = aligned_SWRs_raw(:, 1:237);

%% Plot cluster results

% Cluster1 = labels_freq == 0;
% Cluster2 = labels_freq == 1;
% Cluster3 = labels_freq == 2;
% 
% nrOfSeconds = 3;
% nrOfFrames = (nrOfSeconds*31*2)+1;
% time = (-(31*nrOfSeconds):(31*nrOfSeconds))./31;
% 
% % Find average DF/F signal across all ROIs
% signal = sData.imdata.roiSignals(2).newdff;
% meanDFF = mean(signal);
% 
% prompt = sprintf('Do baseline subtraction? (y = yes | everything else = no) ');
% baseSub = input(prompt,'s');
% 
% % loop over awake SWRs
% figure(1), 
% for i = 1:3
%     
%     if i == 1
%         RippleIdx = sData.ephysdata.frameRipIdx(Cluster1);
%     elseif i == 2
%         RippleIdx = sData.ephysdata.frameRipIdx(Cluster2);
%     elseif i == 3
%         RippleIdx = sData.ephysdata.frameRipIdx(Cluster3);
%     end
%     SWRactivity = zeros(length(RippleIdx),nrOfFrames);
%     
%     for awake_SWRnr = 1:length(RippleIdx) % creates vector of nr of ripples during recording
%         
%         %gives frame nr X  before ripple peak 
%         x = RippleIdx(awake_SWRnr) - (nrOfSeconds*31); 
%         
%         % check if first time point in SWR-aligned time series begins
%         % before first imaging frame, or ends after last imaging frame in
%         % in session. If so, set time point to 1 or last imaging frame,
%         % respectively. 
%         if x <= 0
%             x = 1;
%         end
%         
%         %gives frame nr X after ripple peak
%         y = RippleIdx(awake_SWRnr) + (nrOfSeconds*31); 
%         if y > length(signal)
%             y = length(signal);
%         end
% 
%         if length(x:y) == (nrOfSeconds*31)*2+1 %if the number of frames from 2s prior to 2s post ripple occurence equals 125 frames, do this
%             SWRactivity(awake_SWRnr, :) = meanDFF(x:y); %index every individual ripple frame segment for the mean DF/F )
%         elseif length(x:y) < (nrOfSeconds*31)*2+1 && length(x:y) > (nrOfSeconds*31)*2+1 %skip ripples that occurred so early/late in recording that there isn't 2 full seconds before/after ripple peak 
%             SWRactivity(awake_SWRnr, :) = NaN;
%         end
% 
%         % Baseline subtraction
%         if strcmp(baseSub, 'y')
%             baselineDFF = nanmean(SWRactivity(awake_SWRnr,1:31));
%             SWRactivity(awake_SWRnr,:) = SWRactivity(awake_SWRnr,:)-baselineDFF;
%             clear baselineDFF
%         end
%     end
%     
%     zscore_SWRactivity = zscore(SWRactivity,0,2);
%     
%     subplot(4,3,i),
%     imagesc(SWRactivity)
%     
%     subplot(4,3,i+3),
%     plot(mean(SWRactivity,1))
% %     set(gca, 'ylim', [min(mean(SWRactivity,1),[],2), max(mean(SWRactivity,1),[],2)])
%     
%     subplot(4,3,i+6),
%     imagesc(zscore_SWRactivity)
%     
%     subplot(4,3,i+9),
%     plot(mean(zscore_SWRactivity))
% 
% end
