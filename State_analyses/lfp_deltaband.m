function sData = lfp_deltaband(sData, params)

% Written by Christoffer Berge | Vervaeke lab

% Filter LFP and ECoG in slow oscillation, delta & sigma band (part of this
% was missing in the original pipeline...)


signal           = sData.ephysdata.(params.signal);
% ecog          = sData.ephysdata2.lfp;
signal_length = length(sData.ephysdata.lfp);
srate         = 2500;
nyquist       = srate/2;

%% Define frequency-domain shape and cut-offs
shape  = [0 0 1 1 0 0];
transw = 0.1;

switch params.freq_band
    case 'slow'
        freq       = [0.1 1];
        label = 'slowband';
        order_nr = round( 3*srate/freq(1) );
    case 'delta'
        freq       = [0.5 4];
        label = 'deltaband';
        order_nr = round( 3*srate/freq(1) );

    case 'sigma'
        freq        = [8 18];
        label = 'sigmaband';
        order_nr     = round( 20*srate/freq(1) );
    case 'theta'
        freq       = [5 12];
        label = 'thetaband';
        order_nr  = round( 25*srate/freq(1) );
end
filter_req_shape = [0 freq(1)-freq(1)*transw, freq(1) freq(2), freq(2)+freq(2)*transw nyquist]/nyquist;


%% Generate filter kernels

% Specify filter order

filtkern  = fir1(order_nr, freq/nyquist);

% order_delta     = round( 3*srate/delta_freq(1) );
% filtkern_delta  = fir1(order_delta, delta_freq/nyquist);
% 
% order_sigma     = round( 20*srate/sigma_freq(1) );
% filtkern_sigma  = fir1(order_sigma, sigma_freq/nyquist);

%% Evaluate kernels and their power spectrum

% Compute the power spectrum of the filter kernels
filt_pow     = abs(fft(filtkern).^2);



% Compute the frequencies vector and remove negative frequencies
hz_vec      = linspace(0, nyquist, floor(length(filtkern)/2)+1);
filtpow_trim = filt_pow(1:length(hz_vec));

% Plot kernels
% figure, 
% plot(hz_vec, filtpow_trim, 'b', 'linew', 2), hold on
% plot([0 freq(1) freq freq(2) nyquist], shape, 'r','linew', 2)
% set(gca, 'xlim', [0, freq(2)*4])
% xlabel('Frequency (Hz)'), ylabel('Filtergain')
% legend({'Filterkernel'; 'Ideal filter'})
% title([params.freq_band, ' ', num2str(freq)])


%% Apply filter kernels to data

% Filter LFP and ECoG in slow oscillation frequency band
tic;
filtered_signal  = filtfilt(filtkern, 1, signal);
% so_ecog_filtsig = filtfilt(filtkern, 1, ecog);

% Filter LFP in delta & sigma band
% delta_LFP_filtsig = filtfilt(filtkern_delta, 1, lfp);
% sigma_LFP_filtsig = filtfilt(filtkern_sigma, 1, lfp);

toc;
% Save filtered signals in sData
sData.ephysdata.(label)    = filtered_signal;
% sData.ephysdata.deltaband = delta_LFP_filtsig;
% sData.ephysdata.sigmaband = sigma_LFP_filtsig;
% 
% sData.ephysdata2.soband  = so_ecog_filtsig;
