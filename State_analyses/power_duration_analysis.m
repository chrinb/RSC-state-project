function [initial_mean_power_vec, episode_duration_vec] = power_duration_analysis(sData, params)

% Written by Christoffer Berge | Vervaeke lab

% For each episode in session, compute mean power in first 1s/4s and extraxt
% corresponding episode duration to examine relationship. 

%% Compute mean power during window

% Window size for power average
win_size = params.win_size*2500;

% Get state/frequency band-specific data
state_data = freq_band_analysis(sData, params);

switch params.ephys_signal
    case 'lfp'
        txt = 'state_lfp_filt_ampl';
    case 'ecog'
        txt = 'state_ecog_filt_ampl';
end

% Preallocate
[initial_mean_power, episode_duration] = deal( cell( size(state_data.(txt), 1), 1));

% Loop over episodes
for ep_nr = 1:size(state_data.(txt), 1)
    
    % Get data
    tmp_data_ephys  = state_data.(txt){ep_nr, 1};
    tmp_ep_duration = length([state_data.state_times(ep_nr,1):state_data.state_times(ep_nr,2)] )/ 2500;
    % Check if state data is empty (indicating that episode didn't fit
    % criteria). 
    if ~isempty(tmp_data_ephys) 

        % Mean power in beginning of episode
        initial_mean_power{ep_nr, 1} = mean( tmp_data_ephys(1:win_size));
        
        % Get the duration of corresponding episode (in seconds)
        episode_duration{ep_nr, 1} = tmp_ep_duration;
    end

end

% Concatenate mean binned power and DF/F values into vector shape
initial_mean_power_vec = horzcat( initial_mean_power{:});
episode_duration_vec   = horzcat( episode_duration{:});


%% Calculate linear regression
% x = initial_mean_power_vec';
% y = episode_duration_vec';
% 
% X = [ ones(length(x),1 ), x];
% 
% b = X\y;
% 
% yCalc = X*b;
% 
% figure, 
% scatter(x, y)
% hold on
% plot(x, yCalc, '-')
% xlabel('Binned theta power', FontSize=16)
% ylabel('Binned mean DF/F', FontSize=16)