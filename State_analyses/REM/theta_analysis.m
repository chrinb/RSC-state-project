function output = theta_analysis(sData)

% Written by Christoffer Berge | Vervaeke lab

% Extract REM snippet from CA1 LFP and cortical ECoG channel in sleep 
% sessions. Exclude REM episodes that are (1) ongoing at recording start, 
% (2) continues past recording end, or (3) are < 30s in duration. 

% Load variables
lfp           = sData.ephysdata.lfp;
ecog          = sData.ephysdata2.lfp;
lfp_theta     = sData.ephysdata.thetaband;
ecog_theta    = sData.ephysdata2.thetaband;
srate         = 2500;

% Get REM sleep start/stop times
rem_theta_times = rem_sleep(sData);

% Set a 5s threshold to skip REM episodes that extends before/after
% recording start/end. 
threshold = srate*5;

% Preallocate
[rem_lfp, rem_lfp_theta, rem_lfp_theta_ampl, rem_ecog, rem_ecog_theta, ...
    rem_ecog_theta_ampl, rem_snippet] = deal( cell( size(rem_theta_times,1), 1));

% Loop over REM episodes
for rem_ep_nr = 1:size(rem_theta_times,1)
    
    % Get start/stop times of current REM episode
    tmp_rem_times = [rem_theta_times(rem_ep_nr, 1) rem_theta_times(rem_ep_nr, 2)];
    
    % Check if REM episode fits criteria
    if tmp_rem_times(1) > threshold && tmp_rem_times(2) < size(lfp,1)-threshold && size(tmp_rem_times(1):tmp_rem_times(2),2) > 30*srate
        
        rem_snippet{rem_ep_nr} = tmp_rem_times(1):tmp_rem_times(2);
        
        % Extract raw, theta band, and theta amplitude signal from LFP
        rem_lfp{rem_ep_nr}            = lfp(rem_snippet{rem_ep_nr});
        rem_lfp_theta{rem_ep_nr}      = lfp_theta(rem_snippet{rem_ep_nr});
        rem_lfp_theta_ampl{rem_ep_nr} = abs( hilbert(lfp_theta(rem_snippet{rem_ep_nr})));

        % Extract raw, theta band, and theta amplitude signal from ECoG
        rem_ecog{rem_ep_nr}            = ecog(rem_snippet{rem_ep_nr});
        rem_ecog_theta{rem_ep_nr}      = ecog_theta(rem_snippet{rem_ep_nr});
        rem_ecog_theta_ampl{rem_ep_nr} = abs( hilbert(lfp_theta(rem_snippet{rem_ep_nr})));


    end
end

output = struct();

output.rem_lfp            = rem_lfp;
output.rem_lfp_theta      = rem_lfp_theta;
output.rem_lfp_theta_ampl = rem_lfp_theta_ampl;

output.rem_ecog            = rem_ecog;
output.rem_ecog_theta      = rem_ecog_theta;
output.rem_ecog_theta_ampl = rem_ecog_theta_ampl;

output.ep_duration         = rem_snippet;


    