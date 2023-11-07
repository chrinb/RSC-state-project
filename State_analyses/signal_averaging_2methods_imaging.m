function  output = signal_averaging_2methods_imaging(all_snippets_current_state, all_times_current_state)

% Written by Christoffer Berge | Vervaeke lab

%{
Function that computes average of *imaging* data (ROI x time array) of various lengths
(for example due to sleep episodes of different lengths) using two
strategies: First find the longest data vector (i.e. longest episode) and
then either (2) pad the shorter ones with NaNs, or (2) interpolate the
shorter ones so all have the same length (normalized).
%}

n_episodes = size(all_times_current_state,1);

state_durations =  cell2mat(cellfun(@(c) size(c,2), all_times_current_state, 'uni', 0));

% Find maximum state duration
[max_state_duration, max_dur_idx] = max(state_durations);

%% Two different strategies: (1) Pad shorter state episodes with NaNs, and (2) interpolate shorter state episodes

% Initialize matrices
[padded_state_mat, interp_state_mat] = deal( []);

tic;
% Loop over episodes
for ep_nr = 1:n_episodes
    
    % Get data from one session
    temp_ep_data = all_snippets_current_state{ep_nr, 1};

    [tmp_state_pad, temp_state_interp] = deal( temp_ep_data);
    if ~(ep_nr == max_dur_idx) && ~isempty(temp_ep_data)

        % Pad with NaNs
        tmp_state_pad(:,end:max_state_duration) = NaN;

          % Interpolate short state
        x  = 1:size(temp_ep_data,2);
        v  = temp_state_interp';
        xq = (linspace(1, size(temp_state_interp,2), max_state_duration));
        temp_state_interp = interp1(x, v, xq, 'spline' );

    else
        temp_state_interp = temp_ep_data'; 

    end
   
    % Store padded and interpolated state vector separate matrices
    padded_state_mat = [padded_state_mat; tmp_state_pad];
    interp_state_mat = [interp_state_mat; temp_state_interp'];
end
t = toc;
fprintf('\n Averaged signals for state %.1f seconds',t)

% Compute mean over all state episodes
mean_padded_state           = mean(padded_state_mat, 1 , 'omitnan');
mean_interp_state           = mean(interp_state_mat, 1, 'omitnan');
mean_interp_state_SE        = std(interp_state_mat, 1, 'omitnan')/sqrt( size(interp_state_mat,1));
mean_interp_state_DF        = std(interp_state_mat, 1, 'omitnan');

% Save results in struct
output = struct();
output.padded_state_mat     = padded_state_mat;
output.interp_state_mat     = interp_state_mat;
output.mean_padded_state    = mean_padded_state;
output.mean_interp_state    = mean_interp_state;
output.mean_interp_state_SE = mean_interp_state_SE;
output.mean_interp_state_DF = mean_interp_state_DF;
