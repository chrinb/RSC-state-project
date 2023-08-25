function spindle_idx = get_spindle_idx(opts, select_spindle_idx)

% Written by Christoffer Berge | Vervaeke Lab

% Function that selects which spindles to analyze, either (1) all spindles, 
% or (2) removing clustered SWRs (but keeping the first in each cluster)

% Anon function for evaluating options
checkParameter = @(param, n, str) (isnumeric(param) && param==n) || strcmp(param, str);

% Keep all spindles
if checkParameter(opts.remSpin, 1 , 'all')
    spindle_idx = select_spindle_idx;

% Remove clustered spindles
elseif checkParameter(opts.remSpin, 2 , 'remove close')
    [spindle_idx] = remove_close_spindles(select_spindle_idx, []);
end
