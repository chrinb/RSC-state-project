function output = freq_band_analysis(sData, params)

% Written by Christoffer Berge | Vervaeke lab

%{
Extract raw and filtered NREM/REM snippet(s) from CA1 LFP and cortical 
ECoG channel in sleep sessions and store snippets in cell array. Exclude 
NREM/REM episodes that are (1) ongoing at recording start, (2) continues 
past recording end, or (3) are < certain duration threshold. 
%} 

% Frequency band parameter input specifices which frequency and states to 
% analyze. 
switch params.freq_band
    case 'SO'
        freq_band    = 'soband';
        state_times  = nrem_sleep(sData);
        duration_lim = 10; % min 10 s episode for NREM
    case 'delta'
        freq_band    = 'deltaband';
        state_times  = nrem_sleep(sData);
        duration_lim = 10; % min 10 s episode for NREM

    case 'theta'
        freq_band = 'thetaband';
        state_times = rem_sleep(sData);
        duration_lim     = 30; % min 30 s episode for REM
    case 'sigma'
        freq_band    = 'sigmaband';
        state_times  = nrem_sleep(sData);
        duration_lim = 10; % min 10 s episode for NREM

end

% Select ephys signal
switch params.ephys_signal
    case 'lfp'
        signal           = sData.ephysdata.lfp;
        signal_filt      = sData.ephysdata.(freq_band);
    case 'ecog'
        signal           = sData.ephysdata2.lfp;
        signal_filt      = sData.ephysdata2.(freq_band);
end

if strcmp(params.zscore, 'yes')
    signal      = zscore(signal);
    signal_filt = zscore(signal_filt);
end
% Compute the amplitude envelope (power) of filtered signal
signal_filt_ampl = abs( hilbert(signal_filt));

% Ephys sample rate
srate = 2500;

% Set a 3s threshold to skip REM episodes that extends before/after
% recording start/end. 
threshold = srate*3;

% Set minimum duration for episodes to be analyzed
ep_duration_min = duration_lim*srate;

% Preallocate
[state_ephys, state_ephys_filt, state_ephys_filt_ampl, state_snippet] ...
    = deal( cell( size(state_times,1), 1));

% Loop over episodes
for state_ep_nr = 1:size(state_times,1)
    
    % Get start/stop times of current episode
    tmp_state_times = [state_times(state_ep_nr, 1) state_times(state_ep_nr, 2)];
    
    % NEED BETTER WAY TO SORT EPISODES

    % Check if episode fits criteria
    if tmp_state_times(1) > threshold && tmp_state_times(2) < size(signal,1)-threshold && size(tmp_state_times(1):tmp_state_times(2),2) > ep_duration_min
        
        % Get time points of episode
        state_snippet{state_ep_nr} = tmp_state_times(1):tmp_state_times(2);
        
        % Extract raw, frequency band, and frequency amplitude signal from
        % signal
        state_ephys{state_ep_nr}           = signal(state_snippet{state_ep_nr});
        state_ephys_filt{state_ep_nr}      = signal_filt(state_snippet{state_ep_nr});
        state_ephys_filt_ampl{state_ep_nr} = signal_filt_ampl(state_snippet{state_ep_nr});

    end
end

% Store results in struct
output = struct();

output.state_ephys           = state_ephys;
output.state_ephys_filt      = state_ephys_filt;
output.state_ephys_filt_ampl = state_ephys_filt_ampl;

output.ep_duration         = state_snippet;
output.state_times         = state_times;
output.signal              = signal;


    