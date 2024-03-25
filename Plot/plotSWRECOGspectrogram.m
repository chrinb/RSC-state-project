function [time, F,T, P] = plotSWRECOGspectrogram(sData, params)


%% Run analysis

ripple_idx = sData.ephysdata.absRipIdx;
n_ripples  = length(sData.ephysdata.absRipIdx);

% ecog_z = zscore(sData.ephysdata2.lfp');
% 
% wide_pass_freq = [10 500];
% nyquistFs = 2500/2;
% filter_kernel_bandpass = fir1(600, [wide_pass_freq(1) wide_pass_freq(2)]./nyquistFs);
% ecog_bp_freq = filtfilt(filter_kernel_bandpass, 1, ecog_z); % Fil

if strcmp(params.signal_type, 'ECoG')
    % signal = zscore(sData.ephysdata2.lfp');
    signal = sData.ephysdata2.lfp';

    % signal = ecog_bp_freq;
    % signal = sData.ephysdata2.lfp-mean(sData.ephysdata2.lfp);

elseif strcmp(params.signal_type,'LFP')
    signal = sData.ephysdata.lfp';
end

% find RSC ECoG snippets
for k = 1:n_ripples
    snippet_start = ripple_idx(k)-750;
    snippet_end   = ripple_idx(k)+750;

    if snippet_start > 0 && snippet_end < length(signal)
        snippet(k,:)  = signal(snippet_start:snippet_end);
    end
end

% window  = 128;              % Window size for computing the spectrogram (FFT) [# samples]
% overlap = 127;

% Window, overlap, and nFFT parameters chosen bc they produced a "smooth" spectrogram
% without increases at certain regular frequencies (visible as bands in the spectrogram.
% This was particularly prominent in ECoG data).

window  = 120;              % Window size for computing the spectrogram (FFT) [# samples]
overlap = 119; % Overlap of the windows for computing the spectrogram [# samples]
nFFT    = 50:20:300;        % Vector defining the frequencies for computing the FFT
Fs      = 2500;             % Signal sampling frequency.

clear P
for z = 1:size(snippet,1)
    fprintf('Calculating spectrogram...ripple %d of %d\n',z,size(snippet,1));
    
    x = snippet(z,:);
    % x = zscore(x,0,2);
    x = x - mean(x);

    [~,~,~,cP] = spectrogram(x,window,overlap,nFFT,Fs);
    P(:,:,z) = cP;
end


% P = mean(P,3);
% zscoreP = zscore(P,0,2);
[~,F,T,~] = spectrogram(x,window,overlap,nFFT,Fs);
% Plot spectrogram
time = (-(length(T)/2):(length(T)/2)-1)./2500;

% figure(2), clf
% sgtitle(sessionID),
% subplot(121)
% contourf(time,F,10*log10(abs(P)),200,'edgecolor','none');
% colormap(viridis)
% colorbar
% title('Avg HPC-SWR aligned ECoG spectrogram')
% ylabel('Frequency (Hz)')
% 
% subplot(212)
% contourf(time,F,10*log10(abs(zscoreP)),200,'edgecolor','none');
% colormap(jet)
% caxis([-5 5])
% colorbar
% title('Avg HPC-SWR aligned ECoG spectrogram (z-score)')
% ylabel('Frequency (Hz)')
% xlabel('Time from HPC SWR peak (s)')