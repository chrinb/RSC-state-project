function analysis_output =  spindle_mod_an(varargin)

% Written by Christoffer Berge | Vervaeke Lab

% Function that finds spindle-modulated ROIs 


sData = varargin{1,1};

checkParameter = @(param, n, str) (isnumeric(param) && param==n) || strcmp(param, str);

% Check if user specified particular options for analysis. If not use
% default options. 
if nargin > 1
    if isstruct( varargin{1,2} )
        opts = varargin{1,2};
    else
        try
            opts.exp_type           = varargin{1,2};
            opts.signal_type        = varargin{1,3};
            opts.split_rois         = varargin{1,4};
%             opts.swr_for_analysis   = varargin{1,5};
%             opts.swr_cluster_thresh = varargin{1,6};
%             opts.swr_awake_or_sleep = varargin{1,7};
        catch
        end
    end
end

% Get ROI signals
if nargin > 1 && isfield(opts, 'preDef')
    [signal, text, opts, label3, rois_for_an, ~, cells_to_exclude] = ...
    get_roi_signals_from_sData(sData, opts );
else
    [signal, text, opts, label3, rois_for_an, ~, cells_to_exclude] = ...
    get_roi_signals_from_sData(sData, get_event_analysis_options(varargin) );
end

% Z-score data. 
zscore_signal = bsxfun(@minus, signal, mean(signal,2) ); 
zscore_signal = bsxfun(@rdivide, zscore_signal, std(signal,0,2));

%% Define variables for analysis
nr_of_seconds  = 3; % Before/after sindle onset
nr_of_shuffles = 1000; % nr of shuffle iterations
min_nr_events  = 5; 

frames         = sData.daqdata.frame_onset_reference_frame;
% time           = (-(31*nr_of_seconds):(31*nr_of_seconds))./31;
% sessionID      = sData.sessionInfo.sessionID;
win_length           = (nr_of_seconds*31)*2+1;
All_ROIs             = zeros(nr_of_shuffles,1);

% time2                 = linspace(-nr_of_seconds,nr_of_seconds,win_length);
% swr_idx              = frames(round(swr_idx));
time_win             = (62:108); % -1s to +500ms after SWR peak
% time_win = 1:win_length;

%% Select sleep spindles for analysis

spindle_select     = [];
spin_start_end_str = strcat('NREMspindleStartEnd', spindle_select);

% convert sleep spindle time stamps from e-phys to 2P time
spindle_idx_start = sData.ephysdata2.(spin_start_end_str)(:,1);

% Get time points for sleep spindle onsets
spindle_idx_start = get_spindle_idx(opts , spindle_idx_start);

%% Preallocate 
[signal_spindle_activity, signal_spindle_activity_zscore] ...
    = deal( zeros(length(spindle_idx_start), win_length) );

[signal_spindle_activity_mean,  signal_spindle_activity_meanZ,...
signal_spindle_median, signal_spindle_medianZ] = deal( zeros(size(signal,1), win_length) );

% convert sleep spindle onset time stamps from e-phys to 2P time
spindle_idx_start = frames(round(spindle_idx_start));

%% Run main analysis

% Loop over rois
for roinr = rois_for_an

    roi_signal        = signal(roinr, :); 
    roi_signal_zscore = zscore_signal(roinr,:); 

    % Get all peri-spindle activity windows from current ROI
    [signal_spindle_activity,signal_spindle_activity_zscore] = ...
        extract_avg_activity(spindle_idx_start,roi_signal, roi_signal_zscore,signal_spindle_activity,...
        signal_spindle_activity_zscore, nr_of_seconds, signal, opts.baseSub, roinr);

    % Find mean/median spindle-aligned activity for current ROI
    signal_spindle_activity_mean(roinr,:)  = mean(signal_spindle_activity, 'omitnan');
    signal_spindle_activity_meanZ(roinr,:) = mean(signal_spindle_activity_zscore, 'omitnan');

%     signal_swr_activity_median(roinr,:)  = median(signal_swr_activity, 'omitnan');
%     signal_swr_activity_medianZ(roinr,:) = median(signal_swr_activity_zscore, 'omitnan');

    % Run shuffle analysis
    All_ROIs = shuffle_analysis(win_length,...
        nr_of_shuffles, signal_spindle_activity, roinr, All_ROIs, signal_spindle_activity_mean,...
        time_win, min_nr_events);
end

%% Find activated and suppressed ROIs

% Get activated and suppressed ROIs
sorted_rois         = split_modulated_rois(All_ROIs, win_length, signal_spindle_activity_mean, ...
                    signal_spindle_activity_meanZ, opts, sData, spindle_idx_start, nr_of_seconds, 'spindle');

% Store variables used in the analysis
analysis_varables   = {All_ROIs, win_length, signal_spindle_activity_mean, ...
                    signal_spindle_activity_meanZ, opts, spindle_idx_start, nr_of_seconds};

% Find percentage of activated or suppressed ROIs
prc_mod = prc_mod_cells(rois_for_an, sorted_rois);

%% Store output in struct
analysis_output = struct;

analysis_output.sorted_rois       = sorted_rois;
analysis_output.analysis_varables = analysis_varables;
analysis_output.opts              = opts;
analysis_output.prc_mod           = prc_mod;
analysis_output.rois              = rois_for_an;