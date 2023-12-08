function [r1, c1, lags1] = mean_inh_vs_rem_thetaDeltaRatio(sData)

params.cell_type              = 'in';
params.signal_type            = 'Dff';
params.use_roi_classification = 'grid';
params.zscore                 = 'yes';
[signal_to_plot, ~]           = get_roi_signals_from_sData(sData, params);

%% Theta signals
theta             = sData.ephysdata.thetaband;
delta             = sData.ephysdata.deltaband;
theta_ampl        = abs(hilbert(theta));
delta_ampl        = abs(hilbert(delta));
theta_ampl_smooth = smoothdata(theta_ampl, 'gaussian', 5000);
delta_ampl_smooth = smoothdata(delta_ampl, 'gaussian', 5000);
theta_delta_R     = theta_ampl_smooth./delta_ampl_smooth;

%% Stuff
% time_vector  = linspace(0, length(roi_data), length(roi_data))/imaging_sampling_rate;
% time_ephys   = (0:length(theta)-1)./2500;

REM_episodes = rem_sleep(sData);

frames = sData.daqdata.frame_onset_reference_frame;

signal = mean( okada(signal_to_plot{2,1}, 2));

imag_srate = find_imaging_framerate(sData);
ephys_srate = 2500;
ds_factor = round(ephys_srate/imag_srate);

rem_start_stop = frames(REM_episodes);

% Downsample to match imaging data
thetaR_ds     = downsample(theta_delta_R, ds_factor);
theta_ampl_ds = downsample(theta_ampl_smooth, ds_factor);

snip_theta_ampl = theta_ampl_ds( rem_start_stop(1):rem_start_stop(2));
snip_thetaR     = thetaR_ds( rem_start_stop(1):rem_start_stop(2));
snip_signal     = signal( rem_start_stop(1):rem_start_stop(2));

%% Find mean inhibitory activity frequency via FFT
% imag_srate = find_imaging_framerate(sData);
% nyquist = imag_srate/2;
% 
% mean_inh_pow = abs( fft(snip_signal).^2);
% 
% hz_vec = linspace(0, nyquist, floor(length(snip_signal)/2)+1);
% min_inh_pow_trim = mean_inh_pow(1:length(hz_vec));
% 
% figure, 
% plot(hz_vec, min_inh_pow_trim)
%% Correlation % cross-correlation
r1         = corrcoef(snip_signal, snip_thetaR');
r2        = corrcoef(snip_signal, snip_theta_ampl);
[c1, lags1] = xcorr(  snip_thetaR, snip_signal, 'normalized');
[c2, lags2] = xcorr(  snip_theta_ampl, snip_signal, 'normalized');

%% Plot
figure, 
sgtitle(sData.sessionInfo.sessionID, 'interpreter', 'none')
subplot(311)
hold on
plot(snip_thetaR)
plot(snip_signal)
legend('Theta/delta ratio', 'Mean DF/F inhibitory cells')
xlabel('Time (frames)')
title( ['Correlation coef ' , num2str( r1(2))] )

subplot(312)
hold on
plot(snip_theta_ampl*20)
plot(snip_signal)
legend('Theta ampl.', 'Mean DF/F inhibitory cells')
xlabel('Time (frames)')
title( ['Correlation coef ' , num2str( r2(2))] )

subplot(313)
hold on
plot(lags1, c1, 'k')
plot(lags2, c2, 'b')
title('Cross correlation')
legend('theta/delta ratio', 'theta ampl.')
ylabel('Corr coef.')
xlabel('Lag (frames)')
