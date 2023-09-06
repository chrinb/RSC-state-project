function  output = signal_averaging_2methods(num_sessions, data)

% Written by Christoffer Berge | Vervaeke lab

%{
Function that computes average of data of various lengths
(for example due to sleep episodes of different lengths) using two
strategies: First find the longest data vector (i.e. longest episode) and
then either (2) pad the shorter ones with NaNs, or (2) interpolate the
shorter ones so all have the same length (normalized).
%}

% Preallocate
state_durations = zeros( num_sessions, 1 );

% Loop over sessions
for session_nr = 1:num_sessions
    
    % Get data from one session
    tmp = data{session_nr, 1};
    
    % Store the length of longest episode vector per session
    state_durations(session_nr, 1) = max( cellfun(@(c) size(c,2), tmp.ep_duration), [], 'all');
end

% Find maximum state duration
max_state_duration = max(state_durations);

%% Two different strategies: (1) Pad shorter state episodes with NaNs, and (2) interpolate shorter state episodes

% Initialize matrices
[padded_state_mat, padded_state_mat_smooth, interp_state_mat, interp_state_mat_smooth] = deal( []);

% Loop over sessions
for session_nr = 1:num_sessions
    
    % Get data from one session
    tmp = data{session_nr, 1};

    % Find indices of cell arrays shorter than max state duration (or
    % matching max duration)
    log_idx1 = cellfun(@(c) size(c,2), tmp.ep_duration) <= max_state_duration;

    % Find indices of non-zero cell arrays
    log_idx2 = any(cellfun(@(c) size(c,2), tmp.ep_duration), 2);
    
    % Find overlap between these two indices
    final_idx = (log_idx1 == log_idx2);
    
    % Use this index/indices to find the state vector 
    tmp_state = tmp.state_ephys_filt_ampl(final_idx);
    
    % Loop over state episodes and pad with NaNs
    for ep_nr = 1:size(tmp_state,1)

        % Create two state vectors, one for padding, one for
        % interpolation
        [tmp_state_pad, tmp_state_interp] = deal( tmp_state{ep_nr});
        
        % Check if state vector is shorter than max, if so, pad it
        if size(tmp_state_pad,1) < max_state_duration
            
            % Pad short state
            tmp_state_pad(end:max_state_duration) = NaN;

            % Interpolate short state
            x  = 1:size(tmp_state_interp,1);
            v  = tmp_state_interp';
            xq = linspace(1, size(tmp_state_interp,1), max_state_duration);

            tmp_state_interp = interp1(x, v, xq, 'spline' );
        else
            tmp_state_interp = tmp_state_interp'; 

        end
   
        % Store padded and interpolated state vector separate matrices
        padded_state_mat = [padded_state_mat; tmp_state_pad'];
        interp_state_mat = [interp_state_mat; tmp_state_interp];

        padded_state_mat_smooth = [padded_state_mat_smooth; smoothdata(tmp_state_pad', 'Gaussian', 2500)];
        interp_state_mat_smooth = [interp_state_mat_smooth; smoothdata(tmp_state_interp, 'Gaussian', 2500)];

    end
end

% Compute mean over all state episodes
mean_padded_state           = mean(padded_state_mat, 1 , 'omitnan');
mean_padded_state_smooth    = mean(padded_state_mat_smooth, 1, 'omitnan');
mean_interp_state           = mean(interp_state_mat, 1, 'omitnan');
mean_interp_state_SE        = std(interp_state_mat, 1, 'omitnan')/sqrt( size(interp_state_mat,1));
mean_interp_state_smooth    = mean(interp_state_mat_smooth, 1, 'omitnan');
mean_interp_state_smooth_SE = std(interp_state_mat_smooth, 1, 'omitnan')/sqrt( size(interp_state_mat_smooth,1));

% Save results in struct

output = struct();

output.padded_state_mat            = padded_state_mat;
output.interp_state_mat            = interp_state_mat;
output.padded_state_mat_smooth     = padded_state_mat_smooth;
output.interp_state_mat_smooth     = interp_state_mat_smooth;
output.mean_padded_state           = mean_padded_state;
output.mean_padded_state_smooth    = mean_padded_state_smooth;
output.mean_interp_state           = mean_interp_state;
output.mean_interp_state_SE        = mean_interp_state_SE;
output.mean_interp_state_smooth    = mean_interp_state_smooth;
output.mean_interp_state_smooth_SE = mean_interp_state_smooth_SE;