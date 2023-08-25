function [sData] = bandpassfilter(sData)

% Function used for 4 channel ephys experiments ? (PFC, S1, RSC + EMG)

EMG = sData.ephysdata.lfp;
LFP2 = sData.ephysdata2.lfp;
LFP3 = sData.ephysdata3.lfp;
LFP4 = sData.ephysdata4.lfp;
signalLength = length(sData.ephysdata.lfp);
srate = 2500;
nyquist = srate/2;
delta_t = 1/srate; 
tt = delta_t:delta_t:signalLength*delta_t;

%% Define frequency-domain shape and cut-offs
shape = [0 0 1 1 0 0];
transw = 0.1;

LFP_freq = [0.5 30];

EMG_freq = [100 1000];
EMG_freq_shape = [0 EMG_freq(1)-EMG_freq(1)*transw, EMG_freq(1) EMG_freq(2), EMG_freq(2)+EMG_freq(2)*transw nyquist]/nyquist;

delta_freq = [0.5 4];
delta_freq_shape = [0 delta_freq(1)-delta_freq(1)*transw, delta_freq(1) delta_freq(2), delta_freq(2)+delta_freq(2)*transw nyquist]/nyquist;

theta_freq = [5 9];
theta_freq_shape = [0 theta_freq(1)-theta_freq(1)*transw, theta_freq(1) theta_freq(2), theta_freq(2)+theta_freq(2)*transw nyquist]/nyquist;

sigma_freq = [8 18];
sigma_freq_shape = [0 sigma_freq(1)-sigma_freq(1)*transw, sigma_freq(1) sigma_freq(2), sigma_freq(2)+sigma_freq(2)*transw nyquist]/nyquist;

ripple_freq = [100 250];

%% Generate filter kernels
% Generate filter order
order_LFP = round(3*srate/LFP_freq(1));
order_EMG = round( 10*srate/EMG_freq(1));
order_delta = round( 3*srate/delta_freq(1));
order_theta = round( 25*srate/theta_freq(1));
order_sigma = round( 20*srate/sigma_freq(1));
order_ripple = round( 50*srate/ripple_freq(1)); 

filtkern_LFP = fir1(order_LFP, LFP_freq/nyquist);
filtkern_EMG = fir1(order_EMG, EMG_freq/nyquist);
filtkern_delta = fir1(order_delta, delta_freq/nyquist);
filtkern_theta = fir1(order_theta, theta_freq/nyquist);
filtkern_sigma = fir1(order_sigma, sigma_freq/nyquist);
%filtkern_ripple = fir1(order_ripple, ripple_freq/nyquist);

%% Evaluate kernels and their power spectrum
%compute the power spectrum of the filter kernels
filt_LFP_pow = abs(fft(filtkern_LFP).^2);
filt_EMG_pow = abs(fft(filtkern_EMG).^2);
filt_delta_pow = abs(fft(filtkern_delta).^2);
filt_theta_pow = abs(fft(filtkern_theta).^2);
filt_sigma_pow = abs(fft(filtkern_sigma).^2);
%filt_ripple_pow = abs(fft(filtkern_ripple).^2);

%compute the frequencies vector and remove negative frequencies
hz_LFP = linspace(0, nyquist, floor(length(filtkern_LFP)/2)+1);
filtpow_LFP = filt_LFP_pow(1:length(hz_LFP));

hz_EMG = linspace(0, nyquist, floor(length(filtkern_EMG)/2)+1);
filtpow_EMG = filt_EMG_pow(1:length(hz_EMG));

hz_delta = linspace(0, nyquist, floor(length(filtkern_delta)/2)+1);
filtpow_delta = filt_delta_pow(1:length(hz_delta));

hz_theta = linspace(0, nyquist, floor(length(filtkern_theta)/2)+1);
filtpow_theta = filt_theta_pow(1:length(hz_theta));

hz_sigma = linspace(0, nyquist, floor(length(filtkern_sigma)/2)+1);
filtpow_sigma = filt_sigma_pow(1:length(hz_sigma));

%hz_ripple = linspace(0, nyquist, floor(length(filtkern_ripple)/2)+1);
%filtpow_ripple = filt_ripple_pow(1:length(hz_ripple));

figure,
subplot(321)
plot(hz_LFP, filtpow_LFP, 'b', 'linew', 2), hold on
plot([0 LFP_freq(1), LFP_freq, LFP_freq(2), nyquist],shape, 'r', 'linew', 2), hold off
set(gca, 'xlim', [0 LFP_freq(2)+10])
xlabel('Frequency (Hz)'), ylabel('Filtergain')
legend({'Filterkernel'; 'Ideal filter'})
title('0.5-30 Hz filtered LFP')

subplot(322)
plot(hz_EMG, filtpow_EMG, 'b', 'linew', 2), hold on
plot([0 EMG_freq(1) EMG_freq EMG_freq(2) nyquist], shape, 'r','linew', 2)
set(gca, 'xlim', [0, EMG_freq(2)+100])
xlabel('Frequency (Hz)'), ylabel('Filtergain')
legend({'Filterkernel'; 'Ideal filter'})
title('EMG (100-1000 Hz)')

subplot(323)
plot(hz_delta, filtpow_delta, 'b', 'linew', 2), hold on
plot([0 delta_freq(1) delta_freq delta_freq(2) nyquist], shape, 'r','linew', 2)
set(gca, 'xlim', [0, delta_freq(2)*4])
xlabel('Frequency (Hz)'), ylabel('Filtergain')
legend({'Filterkernel'; 'Ideal filter'})
title('Delta (0.1-4 HZ)')

subplot(324)
plot(hz_theta, filtpow_theta, 'b', 'linew', 2), hold on
plot([0 theta_freq(1) theta_freq theta_freq(2) nyquist], shape, 'r','linew', 2)
set(gca, 'xlim', [0, theta_freq(2)*4])
xlabel('Frequency (Hz)'), ylabel('Filtergain')
legend({'Filterkernel'; 'Ideal filter'})
title('Theta (5-9 Hz)')

subplot(325)
plot(hz_sigma, filtpow_sigma, 'b', 'linew', 2), hold on
plot([0 sigma_freq(1) sigma_freq sigma_freq(2) nyquist], shape, 'r','linew', 2)
set(gca, 'xlim', [0, sigma_freq(2)*4])
xlabel('Frequency (Hz)'), ylabel('Filtergain')
legend({'Filterkernel'; 'Ideal filter'})
title('Sigma (8-18 Hz)')

%subplot(326)
%plot(hz_ripple, filtpow_ripple, 'b', 'linew', 2), hold on
%plot([0 ripple_freq(1) ripple_freq ripple_freq(2) nyquist], shape, 'r','linew', 2)
%set(gca, 'xlim', [0, ripple_freq(2)+50])
%xlabel('Frequency (Hz)'), ylabel('Filtergain')
%legend({'Filterkernel'; 'Ideal filter'})
%title('Ripple frequency')


prompt = sprintf('Keep filtern kernels? (yes = 1 no = 2): ');
keepKernel = input(prompt);
if keepKernel == 1
    % do nothing
else
    error('Try changing order of the kernel(s)');
end

%% Apply filter kernels to data
% create 0.5-30 Hz bandpass-filtered signal
LFP2_filtsig = filtfilt(filtkern_LFP, 1, LFP2);
LFP3_filtsig = filtfilt(filtkern_LFP, 1, LFP3);
LFP4_filtsig = filtfilt(filtkern_LFP, 1, LFP4);

% create 100-1000 Hz bandpass-filtered EMG signal
EMG_filtsig = filtfilt(filtkern_EMG, 1, EMG);

% create delta-frequency filtered signal
delta_LFP2_filtsig = filtfilt(filtkern_delta, 1, LFP2);
delta_LFP3_filtsig = filtfilt(filtkern_delta, 1, LFP3);
delta_LFP4_filtsig = filtfilt(filtkern_delta, 1, LFP4);

% create theta-frequency filtered signal
theta_LFP2_filtsig = filtfilt(filtkern_theta, 1, LFP2);
theta_LFP3_filtsig = filtfilt(filtkern_theta, 1, LFP3);
theta_LFP4_filtsig = filtfilt(filtkern_theta, 1, LFP4);

% create sigma-frequency filtered signal
sigma_LFP2_filtsig = filtfilt(filtkern_sigma, 1, LFP2);
sigma_LFP3_filtsig = filtfilt(filtkern_sigma, 1, LFP3);
sigma_LFP4_filtsig = filtfilt(filtkern_sigma, 1, LFP4);

% create ripple-frequency filtered signal
%ripple_filtsig = filtfilt(filtkern_ripple, 1, LFP3);
%ripple_2_filtsig = filtfilt(filtkern_ripple, 1 , LFP2);

sData.ephysdata2.lfpFilt = LFP2_filtsig;
sData.ephysdata3.lfpFilt = LFP3_filtsig;
sData.ephysdata4.lfpFilt = LFP4_filtsig;

sData.ephysdata.EMGfilt = EMG_filtsig;

sData.ephysdata2.deltaband = delta_LFP2_filtsig;
sData.ephysdata3.deltaband = delta_LFP3_filtsig;
sData.ephysdata4.deltaband = delta_LFP4_filtsig;

sData.ephysdata2.thetaband = theta_LFP2_filtsig;
sData.ephysdata3.thetaband = theta_LFP3_filtsig;
sData.ephysdata4.thetaband = theta_LFP4_filtsig;

sData.ephysdata2.sigmaband = sigma_LFP2_filtsig;
sData.ephysdata3.sigmaband = sigma_LFP3_filtsig;
sData.ephysdata4.sigmaband = sigma_LFP4_filtsig;
%sData.ephysdata.ripplefreq = ripple_filtsig;
%sData.ephysdata.ripplefreq_2 = ripple_2_filtsig;