function [mean_pow_episode, mean_power_across_session, nrem_start_times, nrem_start_times_uncorrected] = get_swa_power_in_session(sData)

% Written by Christoffer Berge | Vervaeke lab

%{
Calculate mean slow-wave power per NREM episode in session, the mean power
of episodes, and when they occurred.
%}

params.freq_band    = 'delta';
params.ephys_signal = 'ecog';

freq = [.5, 4];

% Get state/frequency band-specific data
state_data = freq_band_analysis(sData, params);

% Mean power of this frequency band during session
% freq_band_mean_pow = bandpower(state_data.f, 2500, freq );

% Mean power per episode
mean_pow_episode = cellfun(@mean, state_data.state_ephys_filt_ampl, 'UniformOutput',false);

% Mean power over session
mean_power_across_session = mean(horzcat( mean_pow_episode{:}), 'omitnan');
%% Get NREM episode start times relative to experiment start

% Get NREM start stop times
nrem_start_stop_ephys_time  = nrem_sleep(sData);

% Experiment start time
exp_date               = datetime(sData.sessionInfo.date, InputFormat='yyyyMMdd');
exp_start_time         = datetime(sData.sessionInfo.sessionStartTime, InputFormat='HH:mm:ss');
time_portion           = timeofday(exp_start_time);
correct_exp_start_time = exp_date + time_portion;

% Convert NREM start/stop times from ephys time to duration object in seconds
nrem_start_stop_time_seconds = seconds(nrem_start_stop_ephys_time./2500);

% Add the start times in seconds to experiment start time
nrem_start_times             = correct_exp_start_time + nrem_start_stop_time_seconds(:,1);
nrem_start_times_uncorrected = exp_start_time + nrem_start_stop_time_seconds(:,1);

% Convert back to characters
% event_time_string = datestr(nrem_start_times, 'HH:MM:SS');
