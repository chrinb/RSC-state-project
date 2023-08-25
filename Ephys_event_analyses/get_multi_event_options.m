function opts = get_multi_event_options(keyword)

% Written by Christoffer Berge | Vervaeke Lab

% Get various options for multi-event analysis

opts = getDefaultOptions();

%% Specify experiment type to use for analysis

% Check if user specified experiment type as second input
select_opt = contains(keyword, ["default", "axon", "bulk"]);
if sum(select_opt) > 0
    try
        switch keyword(select_opt)
            case 'default'
                opts = getDefaultOptions();
            case 'axon'
                opts.exp_type = 3;
            case 'bulk'
                opts.exp_type = 1;
                opts.baseSub  = 1; 
        end
    catch
    end
end

%% Specify signal type to use for analysis

% Check if user specified signal type as third input
select_opt = contains(keyword, ["dff", "deconv", "bulk"]);
try
    switch keyword(select_opt)
        case 'dff' 
            opts.signal_type = 2;
        case 'deconv'
            opts.signal_type = 1;
        case 'bulk'
            opts.signal_type = 3;   
    end
catch
end

%% Specify ROIs to use for analsis

% Check if user specified which ROI's to use as fourth input
select_opt = contains(keyword, ["pc", "in", "select"]);
try
    switch keyword(select_opt)
        case 'pc'
            opts.split_rois = 1;
        case 'in'
            opts.split_rois = 2;
        case 'select'
            opts.split_rois = 3;
    end
catch
end

%% Specify ROIs to use for analsis

% Check if user specified behavioral state
select_opt = contains(keyword, ["awake", "sleep"]);
try
    switch keyword(select_opt)
        case 'awake'
            opts.beh_state = 1;
        case 'sleep'
            opts.beh_state = 2;
    end
catch
end

end

% Default settings to use if user doesn't input any of the settings above
function opts = getDefaultOptions()
    opts = struct;
    opts.exp_type    = 2; % 1 = bulk analysis, 2 population analysis, 3 = axons
    opts.signal_type = 2; % Analyze DF/F
    opts.plotting    = 0; % Plot
    opts.split_rois  = 1; % Split ROI array: 1 = principal cells, 2 = inhibitory cells
    opts.baseSub     = 0; % No baseline subtraction
    opts.beh_state   = 1; % awake recording
    opts.predefined  = 1;
end

