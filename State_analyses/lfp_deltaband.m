function sData = lfp_deltaband(sData)

% Written by Christoffer Berge | Vervaeke lab

% Quick fix to filter out delta band activity from LFP signal (this was not
% originally implemented...


lfp           = sData.ephysdata.lfp;
signal_length = length(sData.ephysdata.lfp);
srate         = 2500;
nyquist       = srate/2;

%% Define frequency-domain shape and cut-offs
shape  = [0 0 1 1 0 0];
transw = 0.1;

delta_freq       = [0.5 4];
delta_freq_shape = [0 delta_freq(1)-delta_freq(1)*transw, delta_freq(1) delta_freq(2), delta_freq(2)+delta_freq(2)*transw nyquist]/nyquist;

%% Generate filter kernels
% Generate filter order

order_delta     = round( 3*srate/delta_freq(1) );
filtkern_delta  = fir1(order_delta, delta_freq/nyquist);

%% Evaluate kernels and their power spectrum
%compute the power spectrum of the filter kernels

filt_delta_pow  = abs(fft(filtkern_delta).^2);


% Compute the frequencies vector and remove negative frequencies
hz_delta      = linspace(0, nyquist, floor(length(filtkern_delta)/2)+1);
filtpow_delta = filt_delta_pow(1:length(hz_delta));

%% Apply filter kernels to data

delta_LFP_filtsig  = filtfilt(filtkern_delta, 1, lfp);

% Delta
sData.ephysdata.deltaband  = delta_LFP_filtsig;
