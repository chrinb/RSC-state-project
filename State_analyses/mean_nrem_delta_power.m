function output = mean_nrem_delta_power(sData)

% Written by Christoffer Berge | Vervaeke lab

%{
Compute mean power in the delta frequency band (0.5-4Hz) for each NREM
sleep episode in session, and then average to get session mean.
%}

% Get NREM start stop times
nrem_start_stop_ephys_time  = nrem_sleep(sData);

% Get state/frequency band-specific data
params.freq_band    = 'delta';
params.ephys_signal = 'ecog';
state_data          = freq_band_analysis(sData, params);

for episode_nr = 1:size(nrem_start_stop_ephys_time,1)
    mean_pow_episode


