function analysis_output = swr_awake_an(varargin)

% Written by Christoffer Berge | Vervaeke Lab

% Code to analyse mean SWR-aligned activity during awake recording sessions

% Input(s):    First and required input argument is sData struct. If no other
%              inputs are provided the analysis will use default options.
%              User can, however, specify session type (input #2), signal type
%              (input #3), and whether to split ROI array or not (input
%              #4). See "getSwrAnalysisOptions" for options. 

sData = varargin{1,1};

% Anon function for evaluating options
checkParameter = @(param, n, str) (isnumeric(param) && param==n) || strcmp(param, str);

% Check if user specified particular options for analysis. If not use
% default options.(check "get_event_analysis_options" to see what they are). 
if nargin > 1
    if isstruct( varargin{1,2} )
        opts = varargin{1,2};
    else
        try
            opts.exp_type           = varargin{1,2};
            opts.signal_type        = varargin{1,3};
            opts.split_rois         = varargin{1,4};
            opts.swr_for_analysis   = varargin{1,5};
            opts.swr_cluster_thresh = varargin{1,6};
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
    [signal, text, opts, label3, rois_for_an, roiClustIDs, cells_to_exclude] = ...
    get_roi_signals_from_sData(sData, opts );
else
    [signal, text, opts, label3, rois_for_an, roiClustIDs, cells_to_exclude] = ...
    get_roi_signals_from_sData(sData, get_event_analysis_options(varargin) );
end

% Z-score data. 
zscore_signal = bsxfun(@minus, signal, mean(signal,2) ); 
zscore_signal = bsxfun(@rdivide, zscore_signal, std(signal,0,2));
% zscore_signal = okada(zscore_signal, 2);

% Binarize signal
% if checkParameter(opts.signal_type, 1 , 'deconv')
%     temp_sig      = zscore_signal > 0;
%     zscore_signal = double(temp_sig);
% end

% Various variables
nr_of_seconds  = 3;
frames         = sData.daqdata.frame_onset_reference_frame;
nr_of_frames   = (nr_of_seconds*31*2)+1;
time           = (-(31*nr_of_seconds):(31*nr_of_seconds))./31;
sessionID      = sData.sessionInfo.sessionID;
select_swr_idx = sData.ephysdata.absRipIdx; 

% Get the indicies of user specified SWR types
RippleIdx = get_swr_idx(opts.swr_for_analysis ,sData,select_swr_idx, opts);

% Convert SWR time stamps from e-phys to 2P time
RippleIdx = frames(round(RippleIdx));

% Preallocate 
[signal_swr_activity, signal_swr_activity_zscore] = deal( zeros(length(RippleIdx), nr_of_frames) );

[signal_swr_activity_mean, signal_swr_activity_meanZ,...
signal_swr_activity_median, signal_swr_activity_medianZ] = deal( zeros(size(signal,1), nr_of_frames) );

[all_data, all_data_z] = deal( zeros(size(rois_for_an,2), size(RippleIdx,2), size(time,2)));

t = 1;
%% Run main analysis

% Loop over nr of ROIs
for roinr = rois_for_an 
    
    % Extract signal from currrent ROI
    roi_signal        = signal(roinr, :); 
    roi_signal_zscore = zscore_signal(roinr,:); 

    % Get all peri-SWR activity windows from current ROI
    [signal_swr_activity,signal_swr_activity_zscore] = ...
        extract_avg_activity(RippleIdx, roi_signal, roi_signal_zscore, signal_swr_activity,...
        signal_swr_activity_zscore, nr_of_seconds, opts);

    % store all roi x time x SWR data
    all_data(t,:,:)   = signal_swr_activity;
    all_data_z(t,:,:) = signal_swr_activity_zscore;
    
    % store mean SWR-activity
    signal_swr_activity_mean(roinr,:)  = mean(signal_swr_activity, 'omitnan');
    signal_swr_activity_meanZ(roinr,:) = mean(signal_swr_activity_zscore, 'omitnan');
    % store median SWR-activity
    signal_swr_activity_median(roinr,:)  = median(signal_swr_activity, 'omitnan');
    signal_swr_activity_medianZ(roinr,:) = median(signal_swr_activity_zscore, 'omitnan');
    t = t+1;
%     fprintf('ROI %d\', roinr)
end

%% Plot results

% Label cell type analyzed
if checkParameter(opts.split_rois, 1 , 'principal cells')
    celltype = 'PC';
elseif checkParameter(opts.split_rois, 2 , 'inhibitory cells')
    celltype = 'IN';
else
    celltype = '';
end

sort_idx = [];
if checkParameter(opts.plotting, 1 , 'plot')
sort_idx = plot_swr_awake_an(sessionID,        opts, ...
                  time,                        signal_swr_activity,...
                  signal_swr_activity_zscore,  signal_swr_activity_mean, ...
                  signal_swr_activity_meanZ,   celltype,...
                  RippleIdx,                   text,...
                  label3,                      cells_to_exclude);
end

%% Save output as struct
analysis_output = struct;

if checkParameter(opts.exp_type, 1 , 'bulk')
    analysis_output.signal_swr_activity   = signal_swr_activity;
    analysis_output.signal_swr_activityZ  = signal_swr_activity_zscore;
    analysis_output.opts                  = opts;
    analysis_output.sort_idx              = sort_idx;
    analysis_output.roiClustIDs           = roiClustIDs;

elseif checkParameter(opts.exp_type, 2 , 'default') || checkParameter(opts.exp_type, 3 , 'axon')
    analysis_output.signal_swr_activity_mean    =  signal_swr_activity_mean(rois_for_an,:);
    analysis_output.signal_swr_activity_meanZ   = signal_swr_activity_meanZ(rois_for_an,:);
    analysis_output.signal_swr_activity_median  = signal_swr_activity_median(rois_for_an,:);
    analysis_output.signal_swr_activity_medianZ = signal_swr_activity_medianZ(rois_for_an,:);
    analysis_output.opts                        = opts;
    analysis_output.sort_idx                    = sort_idx;
    analysis_output.roiClustIDs                 = roiClustIDs;
    analysis_output.all_data_z                  = all_data_z;
    analysis_output.all_data                    = all_data;
    analysis_output.roi_idx                     = rois_for_an;
end