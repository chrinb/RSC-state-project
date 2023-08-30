function sData = filter_ephys_sig(varargin)

% Written by Christoffer Berge | Vervaeke Lab

% Filter LFP/ECOG/EMG e-phys signals in various channels into specific 
% frequency bands. 

sData = varargin{1,1};

% Load signals and define variables 
lfp           = sData.ephysdata.lfp;
ecog          = sData.ephysdata2.lfp;
emg           = sData.ephysdata3.lfp;
signal_length = length(sData.ephysdata.lfp);
srate         = 2500;
nyquist       = srate/2;


%% Define frequency-domain shape and cut-offs
shape  = [0 0 1 1 0 0];
transw = 0.1;

ECOG_freq = [0.5 30];

EMG_freq       = [100 1000];
EMG_freq_shape = [0 EMG_freq(1)-EMG_freq(1)*transw, EMG_freq(1) EMG_freq(2), EMG_freq(2)+EMG_freq(2)*transw nyquist]/nyquist;

delta_freq       = [0.5 4];
delta_freq_shape = [0 delta_freq(1)-delta_freq(1)*transw, delta_freq(1) delta_freq(2), delta_freq(2)+delta_freq(2)*transw nyquist]/nyquist;

theta_freq       = [5 9];
theta_freq_shape = [0 theta_freq(1)-theta_freq(1)*transw, theta_freq(1) theta_freq(2), theta_freq(2)+theta_freq(2)*transw nyquist]/nyquist;

sigma_freq        = [8 18];
sigma_freq2       = [10 16];
sigma_freq_shape  = [0 sigma_freq(1)-sigma_freq(1)*transw, sigma_freq(1) sigma_freq(2), sigma_freq(2)+sigma_freq(2)*transw nyquist]/nyquist;
sigma_freq_shape2 = [0 sigma_freq2(1)-sigma_freq2(1)*transw, sigma_freq2(1) sigma_freq2(2), sigma_freq2(2)+sigma_freq2(2)*transw nyquist]/nyquist;

ripple_freq = [100 250];

%% Generate filter kernels
% Generate filter order
order_ECOG   = round( 3*srate/ECOG_freq(1) );
order_EMG    = round( 10*srate/EMG_freq(1) );
order_delta  = round( 3*srate/delta_freq(1) );
order_theta  = round( 25*srate/theta_freq(1) );
order_sigma  = round( 20*srate/sigma_freq(1) );
order_sigma2 = round( 50*srate/sigma_freq2(1) );
order_ripple = round( 50*srate/ripple_freq(1) ); 

filtkern_ECOG   = fir1(order_ECOG, ECOG_freq/nyquist);
filtkern_EMG    = fir1(order_EMG, EMG_freq/nyquist);
filtkern_delta  = fir1(order_delta, delta_freq/nyquist);
filtkern_theta  = fir1(order_theta, theta_freq/nyquist);
filtkern_sigma  = fir1(order_sigma, sigma_freq/nyquist);
filtkern_sigma2 = fir1(order_sigma2, sigma_freq2/nyquist);
filtkern_ripple = fir1(order_ripple, ripple_freq/nyquist);

%% Evaluate kernels and their power spectrum
%compute the power spectrum of the filter kernels
filt_ECOG_pow   = abs(fft(filtkern_ECOG).^2);
filt_EMG_pow    = abs(fft(filtkern_EMG).^2);
filt_delta_pow  = abs(fft(filtkern_delta).^2);
filt_theta_pow  = abs(fft(filtkern_theta).^2);
filt_sigma_pow  = abs(fft(filtkern_sigma).^2);
filt_sigma_pow2 = abs(fft(filtkern_sigma2).^2);
filt_ripple_pow = abs(fft(filtkern_ripple).^2);

%compute the frequencies vector and remove negative frequencies
hz_LFP      = linspace(0, nyquist, floor(length(filtkern_ECOG)/2)+1);
filtpow_LFP = filt_ECOG_pow(1:length(hz_LFP));

hz_EMG      = linspace(0, nyquist, floor(length(filtkern_EMG)/2)+1);
filtpow_EMG = filt_EMG_pow(1:length(hz_EMG));

hz_delta      = linspace(0, nyquist, floor(length(filtkern_delta)/2)+1);
filtpow_delta = filt_delta_pow(1:length(hz_delta));

hz_theta      = linspace(0, nyquist, floor(length(filtkern_theta)/2)+1);
filtpow_theta = filt_theta_pow(1:length(hz_theta));

hz_sigma      = linspace(0, nyquist, floor(length(filtkern_sigma)/2)+1);
filtpow_sigma = filt_sigma_pow(1:length(hz_sigma));

hz_sigma2      = linspace(0, nyquist, floor(length(filtkern_sigma2)/2)+1);
filtpow_sigma2 = filt_sigma_pow2(1:length(hz_sigma2));

hz_ripple      = linspace(0, nyquist, floor(length(filtkern_ripple)/2)+1);
filtpow_ripple = filt_ripple_pow(1:length(hz_ripple));

if length(varargin) > 1
    figure,
    subplot(321)
    plot(hz_LFP, filtpow_LFP, 'b', 'linew', 2), hold on
    plot([0 ECOG_freq(1), ECOG_freq, ECOG_freq(2), nyquist],shape, 'r', 'linew', 2), hold off
    set(gca, 'xlim', [0 ECOG_freq(2)+10])
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
    plot(hz_sigma2, filtpow_sigma2, 'b', 'linew', 2), hold on
    plot([0 sigma_freq2(1) sigma_freq2 sigma_freq2(2) nyquist], shape, 'r','linew', 2)
    set(gca, 'xlim', [0, sigma_freq2(2)*4])
    xlabel('Frequency (Hz)'), ylabel('Filtergain')
    legend({'Filterkernel'; 'Ideal filter'})
    title('Sigma (10-16 Hz)')

    subplot(326)
    plot(hz_ripple, filtpow_ripple, 'b', 'linew', 2), hold on
    plot([0 ripple_freq(1) ripple_freq ripple_freq(2) nyquist], shape, 'r','linew', 2)
    set(gca, 'xlim', [0, ripple_freq(2)+50])
    xlabel('Frequency (Hz)'), ylabel('Filtergain')
    legend({'Filterkernel'; 'Ideal filter'})
    title('Ripple frequency')


    prompt = sprintf('Keep filtern kernels? (yes = 1 no = 2): ');
    keepKernel = input(prompt);
    if keepKernel == 1
        % do nothing
    else
        error('Try changing order of the kernel(s)');
    end
end
%% Apply filter kernels to data
% create 0.5-30 Hz bandpass-filtered signal
ECOG_filtsig = filtfilt(filtkern_ECOG, 1, ecog);
%LFP_filtsig = filtfilt(filtkern_LFP, 1, LFP);
%LFP4_filtsig = filtfilt(filtkern_LFP, 1, LFP4);

% create 100-1000 Hz bandpass-filtered EMG signal
EMG_filtsig = filtfilt(filtkern_EMG, 1, emg);

% create delta-frequency filtered signal
delta_ecog_filtsig = filtfilt(filtkern_delta, 1, ecog);
delta_lfp_filtsig  = filtfilt(filtkern_delta, 1, lfp);
%delta_LFP3_filtsig = filtfilt(filtkern_delta, 1, LFP);
%delta_LFP4_filtsig = filtfilt(filtkern_delta, 1, LFP4);

% create theta-frequency filtered signal
theta_ecog_filtsig = filtfilt(filtkern_theta, 1, ecog);
theta_lfp_filtsig = filtfilt(filtkern_theta, 1, lfp);
%theta_LFP4_filtsig = filtfilt(filtkern_theta, 1, LFP4);

% create sigma-frequency filtered signal
sigma_ecog_filtsig = filtfilt(filtkern_sigma, 1, ecog);
sigma2_ecog_filtsig = filtfilt(filtkern_sigma2, 1, ecog);
%sigma_LFP_filtsig = filtfilt(filtkern_sigma, 1, LFP);
%sigma_LFP4_filtsig = filtfilt(filtkern_sigma, 1, LFP4);

% create ripple-frequency filtered signal
ripple_filtsig = filtfilt(filtkern_ripple, 1, lfp);
ripple_2_filtsig = filtfilt(filtkern_ripple, 1 , ecog);

% 0.5-30 Hz ECoG
sData.ephysdata2.lfpFilt = ECOG_filtsig;

% EMG 
sData.ephysdata3.EMGfilt = EMG_filtsig;

% Delta
sData.ephysdata2.deltaband = delta_ecog_filtsig;
sData.ephysdata.deltaband  = delta_lfp_filtsig;
% Theta
sData.ephysdata.thetaband = theta_lfp_filtsig;
sData.ephysdata2.thetaband = theta_ecog_filtsig;

% Sigma
sData.ephysdata2.sigmaband = sigma_ecog_filtsig;
sData.ephysdata2.sigmaband2 = sigma2_ecog_filtsig;

% Ripple 
sData.ephysdata.ripplefreq = ripple_filtsig;
sData.ephysdata2.ripplefreq = ripple_2_filtsig;