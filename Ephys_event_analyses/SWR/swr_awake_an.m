function analysis_output = swr_awake_an(varargin)

% Written by Christoffer Berge | Vervaeke Lab

%{
Code to analyse mean SWR-aligned activity during awake recording sessions
%} 

sData  = varargin{1,1};
params = varargin{1,2};

% Get ROI signals
[signal, params, label3, rois_for_an, roiClustIDs, cells_to_exclude] = get_roi_signals_from_sData(sData, params );

% Z-score data. 
zscore_signal = bsxfun(@minus, signal, mean(signal,2) ); 
zscore_signal = bsxfun(@rdivide, zscore_signal, std(signal,0,2));

% Define variables
nr_of_seconds  = 3;
frames         = sData.daqdata.frame_onset_reference_frame;
nr_of_frames   = (nr_of_seconds*31*2)+1;
time           = (-(31*nr_of_seconds):(31*nr_of_seconds))./31;
sessionID      = sData.sessionInfo.sessionID;
select_swr_idx = sData.ephysdata.absRipIdx; 

% Get the indicies of user specified SWR types
RippleIdx = get_swr_idx(params.swr_for_analysis ,sData,select_swr_idx, params);

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
        signal_swr_activity_zscore, nr_of_seconds, params);

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

sort_idx = [];
if strcmp(params.plotting,'plot')
sort_idx = plot_swr_awake_an(sessionID,        params, ...
                  time,                        signal_swr_activity,...
                  signal_swr_activity_zscore,  signal_swr_activity_mean, ...
                  signal_swr_activity_meanZ,   ...
                  RippleIdx,                   ...
                  label3,                      cells_to_exclude);
end

%% Save output as struct
analysis_output = struct;

if strcmp(params.exp_type, 'bulk')
    analysis_output.signal_swr_activity   = signal_swr_activity;
    analysis_output.signal_swr_activityZ  = signal_swr_activity_zscore;
    analysis_output.params                  = params;
    analysis_output.sort_idx              = sort_idx;
    analysis_output.roiClustIDs           = roiClustIDs;

elseif strcmp(params.exp_type, 'default') || strcmp(params.exp_type, 'axon')
    analysis_output.signal_swr_activity_mean    =  signal_swr_activity_mean(rois_for_an,:);
    analysis_output.signal_swr_activity_meanZ   = signal_swr_activity_meanZ(rois_for_an,:);
    analysis_output.signal_swr_activity_median  = signal_swr_activity_median(rois_for_an,:);
    analysis_output.signal_swr_activity_medianZ = signal_swr_activity_medianZ(rois_for_an,:);
    analysis_output.params                      = params;
    analysis_output.sort_idx                    = sort_idx;
    analysis_output.roiClustIDs                 = roiClustIDs;
    analysis_output.all_data_z                  = all_data_z;
    analysis_output.all_data                    = all_data;
    analysis_output.roi_idx                     = rois_for_an;
end