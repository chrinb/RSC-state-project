function opts = get_multi_session_opts(keyword)

% Written by Christoffer Berge || Vervaeke lab

% Function that gets options for multi session analysis code using Nansen

% Find right data variable
% if strcmp(keyword.exp_type, 'bulk')
%     txt = 'bulk'
% elseif 

strcmp(keyword.beh_state, 'sleep')
% Select appropriate text string for experiment type
if strcmp(keyword.exp_type, 'default')
    opts.signal_text  = 'signal_swr_activity_mean';
    opts.signal_text2 = 'signal_swr_activity_meanZ';
elseif strcmp(keyword.exp_type, 'bulk')
    opts.signal_text  = 'signal_swr';
    opts.signal_text2 = 'signal_swr';
end

% Check behavioral state
if strcmp(keyword.beh_state, 'awake')
    lol = 1;
elseif strcmp(keyword.beh_state, 'sleep')
    lol = 2;
end

