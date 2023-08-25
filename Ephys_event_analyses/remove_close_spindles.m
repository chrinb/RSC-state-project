function adjusted_spindle_idx = remove_close_spindles(spindle_idx, sec_between_spindles)

% Written by Christoffer Berge | Vervaeke Lab

% Function that identifies clustered spindles (i.e. inter-spindle interval
% smaller than a user-defined threshold) and removes all but the first
% spindle in that cluster. This code therefore does not capture "true" single spindles.

%% Adjust interval between sleep spindles

% Set inter-spindle-interval = 2 seconds unless otherwise specified
if isempty(sec_between_spindles)
    sec_between_spindles = 1.5;
end

samples_between_spindles = sec_between_spindles*2500;
adjusted_spindle_idx     = spindle_idx;

% Find nr of e-phys samples between each SWR
pre_spindle_interval = diff(spindle_idx);
spindle_locs          = pre_spindle_interval < samples_between_spindles;
% test             = pre_swr_interval > samplesBetweenRipples;
% testA            = [true(1,1) test];
adjusted_spindle_locs  = [false(1,1); spindle_locs]; %adds a logical 0 to make the vector equal in length to adjustedrippleIdx

% Remove SWRs that occurs to close to the previous SWR. 
adjusted_spindle_idx(adjusted_spindle_locs) = [];
