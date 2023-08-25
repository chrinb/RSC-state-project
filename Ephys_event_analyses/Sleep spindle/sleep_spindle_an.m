function [analysis_output] = sleep_spindle_an(varargin)

%Written by Christoffer Berge | Vervaeke Lab

% Find mean/median sleep spindle onset-aligned DF/F for bulk or population 
% recordings in sleep sessions. Takes as input a sData struct in the form of
% "sData.imdata.roiSignals(2).matrix" where "matrix" corresponds to a field
% containing a ROI x time matrix for a 2P recording session. User specifies
% whether to analyze bulk signals (averages all ROIs to a single mean
% trace) or population signals (loops over ROIs and calculates mean
% sleep spindle onset-aligned activity). Results can be plotted in "plot_sleep_spindle_an"
% function. 

sData = varargin{1,1};

% Anon function for evaluating options
checkParameter = @(param, n, str) (isnumeric(param) && param==n) || strcmp(param, str);

% Check if user specified particular options for analysis. If not use
% default options.(check "getSwrAnalysisOptions" to see what they are). 
if nargin > 1
    if isstruct( varargin{1,2} )
        opts = varargin{1,2};
    else
        try
        opts.exp_type     = varargin{1,2};
        opts.signal_type  = varargin{1,3};
        opts.split_rois   = varargin{1,4};
        opts.spindle_type = varargin{1,5};
%         opts.preDef       = 1;
        opts.baseSub      = 1;
        opts.plotting     = 1;
        catch
        end
    end
end

% Add roisignals folder to path to get roi array
try
    folder_name = dir;
    addpath( [folder_name(1).folder '\roisignals'])
catch
end

% Get ROI signals
if nargin > 1 && isfield(opts, 'preDef')
    [signal, text, opts, label3, rois_for_an, ~, cells_to_exclude] = ...
    get_roi_signals_from_sData(sData, opts );
else
    [signal, text, opts, label3, rois_for_an, ~, cells_to_exclude] = ...
    get_roi_signals_from_sData(sData, get_event_analysis_options(varargin) );
end

% Z-score data
zscore_signal = bsxfun(@minus, signal, mean(signal,2) ); 
zscore_signal = bsxfun(@rdivide, zscore_signal, std(signal,0,2));

% Create different variables
nr_of_seconds = 3; % nr of seconds before/after spindle onset
nr_of_frames  = (nr_of_seconds*31*2)+1; % nr of recording frames in spindle window
time          = linspace(-nr_of_seconds,nr_of_seconds,nr_of_frames); % time vector in seconds
frames        = sData.daqdata.frame_onset_reference_frame;
sessionID     = sData.sessionInfo.sessionID;

%% Select spindle freq
% prompt = sprintf('8-18Hz (1) or 10-16hz (2)? ');
% spindle_band_select = input(prompt);
% 
% if spindle_band_select == 1
%     spindle_select = [];
% elseif spindle_band_select == 2
%     spindle_select = '1016';
% end
spindle_select     = [];
% spin_idx_str       = strcat('NREMAbsSpindleIdx', spindle_select);
spin_start_end_str = strcat('NREMspindleStartEnd', spindle_select);

% spindleIdx = frames(round(sData.ephysdata2.(spin_idx_str)));

% convert sleep spindle time stamps from e-phys to 2P time
spindle_idx_start = sData.ephysdata2.(spin_start_end_str)(:,1);

% Get time points for sleep spindle onsets
spindle_idx_start = get_spindle_idx(opts , spindle_idx_start);
%% Preallocate 
[signal_spindle_activity, signal_spindle_activity_zscore] = deal( zeros(length(spindle_idx_start), nr_of_frames) );

[signal_spindle_activity_mean,  signal_spindle_activity_meanZ,...
signal_spindle_median, signal_spindle_medianZ] = deal( zeros(size(signal,1), nr_of_frames) );

% convert sleep spindle onset time stamps from e-phys to 2P time
spindle_idx_start = frames(round(spindle_idx_start));
%% Run main analysis

% loop over ROIs
for roinr = rois_for_an(1:size(rois_for_an,2)) 

    roi_signal        = signal(roinr, :); % extract time series data for current ROI
    roi_signal_zscore = zscore_signal(roinr,:); % extract z-scored time series data for current ROI

    % Get all peri-spindle activity windows from current ROI
    [signal_spindle_activity,signal_spindle_activity_zscore] = ...
        extract_avg_activity(spindle_idx_start, roi_signal, roi_signal_zscore, signal_spindle_activity,...
        signal_spindle_activity_zscore, nr_of_seconds, opts);

    % average all awake spindle-aligned snippets for a given ROI
    signal_spindle_activity_mean(roinr, :)  = mean(signal_spindle_activity,'omitnan'); 
    signal_spindle_activity_meanZ(roinr, :) = mean(signal_spindle_activity_zscore, 'omitnan'); 

    signal_spindle_median(roinr, :)  = median(signal_spindle_activity, 'omitnan'); 
    signal_spindle_medianZ(roinr, :) = median(signal_spindle_activity_zscore, 'omitnan'); 
end

%% Plot results

% Label cell type analyzed
if opts.split_rois == 1
    celltype = 'PC';
elseif opts.split_rois == 2
    celltype = 'IN';
else
    celltype = '';
end


if ~isempty(opts.plotting)
sort_idx = plot_sleep_spindle_an(sessionID,                      opts,...
                                 time,                           signal_spindle_activity,...
                                 signal_spindle_activity_zscore, signal_spindle_activity_mean,...
                                 signal_spindle_activity_meanZ,  celltype,...
                                 spindle_idx_start,              text,...
                                 label3,                         cells_to_exclude);
end


%% Save output as cell array

if checkParameter(opts.exp_type, 1 , 'bulk')
    analysis_output.signal_spindle_activity        = signal_spindle_activity;
    analysis_output.signal_spindle_activity_zscore = signal_spindle_activity_zscore;
    analysis_output.opts                           = opts;
    analysis_output.sort_idx                       = sort_idx;

elseif checkParameter(opts.exp_type, 2 , 'default') || checkParameter(opts.exp_type, 3 , 'axon')
    analysis_output.signal_spindle_activity_mean = signal_spindle_activity_mean;
    analysis_output.signal_spindle_activity_meanZ = signal_spindle_activity_meanZ;
    analysis_output.signal_spindle_median         = signal_spindle_median;
    analysis_output.signal_spindle_medianZ        = signal_spindle_medianZ;
    analysis_output.opts                          = opts;
    analysis_output.sort_idx                      = sort_idx;

end
