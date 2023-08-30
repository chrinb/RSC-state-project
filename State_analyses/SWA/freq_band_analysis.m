function output = freq_band_analysis(sData, params)

% Written by Christoffer Berge | Vervaeke lab

% Extract raw and filtered NREM/REM snippet(s) from CA1 LFP and cortical 
% ECoG channel in sleep sessions. Exclude NREM/REM episodes that are (1) 
% ongoing at recording start, (2) continues past recording end, or (3) 
% are < 30s in duration. 

% Frequency band parameter input specifices which frequency and states to 
% analyze. 
switch params.freq_band
    case 'SO'
        freq_band = 'soband';
        state_times = nrem_sleep(sData, 1);
    case 'delta'
        freq_band = 'deltaband';
        state_times = nrem_sleep(sData, 1);
    case 'theta'
        freq_band = 'thetaband';
        state_times = rem_sleep(sData);
    case 'sigma'
        freq_band = 'sigmaband';
        state_times = nrem_sleep(sData, 1);
end


% Load variables
lfp         = sData.ephysdata.lfp;
ecog        = sData.ephysdata2.lfp;
lfp_filt    = sData.ephysdata.(freq_band);
ecog_filt   = sData.ephysdata2.(freq_band);
srate       = 2500;

% Set a 5s threshold to skip REM episodes that extends before/after
% recording start/end. 
threshold = srate*5;

% Preallocate
[state_lfp, state_lfp_filt, state_lfp_filt_ampl, state_ecog, state_ecog_filt, ...
    state_ecog_filt_ampl, state_snippet] = deal( cell( size(state_times,1), 1));

% Loop over episodes
for state_ep_nr = 1:size(state_times,1)
    
    % Get start/stop times of current episode
    tmp_state_times = [state_times(state_ep_nr, 1) state_times(state_ep_nr, 2)];
    
    % Check if episode fits criteria
    if tmp_state_times(1) > threshold && tmp_state_times(2) < size(lfp,1)-threshold && size(tmp_state_times(1):tmp_state_times(2),2) > 30*srate
        
        % Get time points of episode
        state_snippet{state_ep_nr} = tmp_state_times(1):tmp_state_times(2);
        
        % Extract raw, theta band, and theta amplitude signal from LFP
        state_lfp{state_ep_nr}           = lfp(state_snippet{state_ep_nr});
        state_lfp_filt{state_ep_nr}      = lfp_filt(state_snippet{state_ep_nr});
        state_lfp_filt_ampl{state_ep_nr} = abs( hilbert(lfp_filt(state_snippet{state_ep_nr})));

        % Extract raw, theta band, and theta amplitude signal from ECoG
        state_ecog{state_ep_nr}           = ecog(state_snippet{state_ep_nr});
        state_ecog_filt{state_ep_nr}      = ecog_filt(state_snippet{state_ep_nr});
        state_ecog_filt_ampl{state_ep_nr} = abs( hilbert(ecog_filt(state_snippet{state_ep_nr})));

    end
end

% Store results in struct
output = struct();

output.state_lfp           = state_lfp;
output.state_lfp_filt      = state_lfp_filt;
output.state_lfp_filt_ampl = state_lfp_filt_ampl;

output.state_ecog           = state_ecog;
output.state_ecog_filt      = state_ecog_filt;
output.state_ecog_filt_ampl = state_ecog_filt_ampl;

output.ep_duration         = state_snippet;
output.state_times         = state_times;


    