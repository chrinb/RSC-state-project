function sData = lfp_deltaband(sData)

% Written by Christoffer Berge | Vervaeke lab

% Filter LFP and ECoG in slow oscillation, delta & sigma band (part of this
% was missing in the original pipeline...)


lfp           = sData.ephysdata.lfp;
ecog          = sData.ephysdata2.lfp;
signal_length = length(sData.ephysdata.lfp);
srate         = 2500;
nyquist       = srate/2;

%% Define frequency-domain shape and cut-offs
shape  = [0 0 1 1 0 0];
transw = 0.1;

so_freq       = [0.1 1];
sofreq_shape = [0 so_freq(1)-so_freq(1)*transw, so_freq(1) so_freq(2), so_freq(2)+so_freq(2)*transw nyquist]/nyquist;

delta_freq       = [0.5 4];
delta_freq_shape = [0 delta_freq(1)-delta_freq(1)*transw, delta_freq(1) delta_freq(2), delta_freq(2)+delta_freq(2)*transw nyquist]/nyquist;

sigma_freq        = [8 18];
sigma_freq_shape  = [0 sigma_freq(1)-sigma_freq(1)*transw, sigma_freq(1) sigma_freq(2), sigma_freq(2)+sigma_freq(2)*transw nyquist]/nyquist;

%% Generate filter kernels

% Specify filter order
order_so     = round( 3*srate/so_freq(1) );
filtkern_so  = fir1(order_so, so_freq/nyquist);

order_delta     = round( 3*srate/delta_freq(1) );
filtkern_delta  = fir1(order_delta, delta_freq/nyquist);

order_sigma     = round( 20*srate/sigma_freq(1) );
filtkern_sigma  = fir1(order_sigma, sigma_freq/nyquist);

%% Evaluate kernels and their power spectrum

% Compute the power spectrum of the filter kernels
filt_so_pow     = abs(fft(filtkern_so).^2);
filt_delta_pow  = abs(fft(filtkern_delta).^2);
filt_sigma_pow  = abs(fft(filtkern_sigma).^2);

% Compute the frequencies vector and remove negative frequencies
hz_so      = linspace(0, nyquist, floor(length(filtkern_so)/2)+1);
filtpow_so = filt_so_pow(1:length(hz_so));

hz_delta      = linspace(0, nyquist, floor(length(filtkern_delta)/2)+1);
filtpow_delta = filt_delta_pow(1:length(hz_delta));

hz_sigma      = linspace(0, nyquist, floor(length(filtkern_sigma)/2)+1);
filtpow_sigma = filt_sigma_pow(1:length(hz_sigma));

%% Plot kernels
% figure, 
% 
% subplot(131)
% plot(hz_so, filtpow_so, 'b', 'linew', 2), hold on
% plot([0 so_freq(1) so_freq so_freq(2) nyquist], shape, 'r','linew', 2)
% set(gca, 'xlim', [0, so_freq(2)*4])
% xlabel('Frequency (Hz)'), ylabel('Filtergain')
% legend({'Filterkernel'; 'Ideal filter'})
% title('SO (0.1-1 HZ)')
% 
% subplot(132)
% plot(hz_delta, filtpow_delta, 'b', 'linew', 2), hold on
% plot([0 delta_freq(1) delta_freq delta_freq(2) nyquist], shape, 'r','linew', 2)
% set(gca, 'xlim', [0, delta_freq(2)*4])
% xlabel('Frequency (Hz)'), ylabel('Filtergain')
% legend({'Filterkernel'; 'Ideal filter'})
% title('Delta (0.1-4 HZ)')
% 
% subplot(133)
% plot(hz_sigma, filtpow_sigma, 'b', 'linew', 2), hold on
% plot([0 sigma_freq(1) sigma_freq sigma_freq(2) nyquist], shape, 'r','linew', 2)
% set(gca, 'xlim', [0, sigma_freq(2)*4])
% xlabel('Frequency (Hz)'), ylabel('Filtergain')
% legend({'Filterkernel'; 'Ideal filter'})
% title('Sigma (8-18 Hz)')
%% Apply filter kernels to data

% Filter LFP and ECoG in slow oscillation frequency band
tic;
so_lfp_filtsig  = filtfilt(filtkern_so, 1, lfp);
so_ecog_filtsig = filtfilt(filtkern_so, 1, ecog);

% Filter LFP in delta & sigma band
delta_LFP_filtsig = filtfilt(filtkern_delta, 1, lfp);
sigma_LFP_filtsig = filtfilt(filtkern_sigma, 1, lfp);

toc;
% Save filtered signals in sData
sData.ephysdata.soband    = so_lfp_filtsig;
sData.ephysdata.deltaband = delta_LFP_filtsig;
sData.ephysdata.sigmaband = sigma_LFP_filtsig;

sData.ephysdata2.soband  = so_ecog_filtsig;
