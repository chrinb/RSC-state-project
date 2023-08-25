function opts = get_event_analysis_options(keyword)

% Written by Christoffer Berge | Vervaeke Lab

% Get various options for analysis

opts = getDefaultOptions();

%% Specify experiment type to use for analysis

% Check if user specified experiment type as second input
try
    switch keyword{1,2}
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

%% Specify signal type to use for analysis

% Check if user specified signal type as third input
try
    switch keyword{1,3}
        case 'dff' 
            opts.signal_type = 2;
        case 'deconv'
            opts.signal_type = 1;
%         case 'bulk'
%             opts.signal_type = 3;   
    end
catch
end

%% Specify cell type to use for analsis

% Check if user specified which ROI's to use as fourth input
try
    switch keyword{1,4}
        case 'pc'
            opts.split_rois = 1;
        case 'in'
            opts.split_rois = 2;
%         case 'select'
%             opts.split_rois = 3;
    end
catch
end

%% Specify criteria to filter SWRs 
try
    opts.swr_for_analysis = keyword{1,5};
catch
end

%% Specify threshold for SWRs occurring in close temporal proximity
try
    opts.swr_cluster_thresh = keyword{1,6};
catch
end

end


% Default settings to use if user doesn't input any of the settings above
function opts = getDefaultOptions()
    opts = struct;
    opts.exp_type           = 2; % 1 = bulk analysis, 2 population analysis, 3 = axons
    opts.signal_type        = 2; % Analyze DF/F
    opts.split_rois         = 1; % Split ROI array: 1 = principal cells, 2 = inhibitory cells
    opts.swr_for_analysis   = 4; % remove clustered and close-to-movement SWRs 
    opts.plotting           = 1; % Plot
    opts.baseSub            = 0; % No baseline subtraction
    opts.swr_for_analysis   = 4; % remove clustered and close-to-movement SWRs 
    opts.swr_cluster_thresh = 2; % Clustered SWRs occurring within 2 sec of each other (keep only first)
    opts.remSpin            = 1; % Use all spindles
end

