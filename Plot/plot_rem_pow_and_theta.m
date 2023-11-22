function plot_rem_pow_and_theta(sData)

params.ephys_signal = 'lfp';
params.freq_band    = 'theta';
params.zscore       = 'yes';

output = freq_band_analysis(sData, params);

time_vec = (0:length(output.state_ephys{1, 1})-1)/2500;

%% Plot
figure, 
% h(1) = subplot(311);
plot(time_vec, output.state_ephys_filt_ampl{1,1}, 'k'  )
ylabel('Z-score theta power', FontSize=16)
title('CA1 LFP theta (5-9Hz) power during REM episode')

% h(2) = subplot(312);
% plot(time_vec, output.state_ephys_filt{1,1}, 'k'  )
% ylabel('Filtered CA1 LFP theta (5-9Hz)', FontSize=16)
% 
% h(3) = subplot(313);
% plot(time_vec, output.state_ephys{1,1}, 'k'  )
% ylabel('Raw CA1 LFP theta ', FontSize=16)

xlabel('Time (s)', FontSize=16)
% linkaxes(h,'x');
