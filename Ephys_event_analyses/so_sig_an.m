function analysis_output = so_sig_an(varargin)

% Written by Christoffer Berge | Vervaeke Lab

% Code that extracts mean signal (e.g. DF/F or deconvolved DF/F) activity
% during slow oscillations and delta waves, as defined in function
% "mark_slow_wave". Can be used for bulk imaging sessions and population
% sessions. 

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
%             opts.swr_for_analysis   = varargin{1,5};
%             opts.swr_cluster_thresh = varargin{1,6};
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

% Various variables
nr_of_seconds  = 1;
frames         = sData.daqdata.frame_onset_reference_frame;
nr_of_frames   = (nr_of_seconds*31*2)+1;
time           = (-(31*nr_of_seconds):(31*nr_of_seconds))./31;
sessionID      = sData.sessionInfo.sessionID;

% Get slow oscillation & delta wave time points
[SOs, delta_waves] = mark_slow_wave(sData);
% delta_waves = mark_slow_wave2(sData);

% Convert to imaging time
frames       = sData.daqdata.frame_onset_reference_frame;

% Extract time stamps for SO/delta events in 2P time. Index 3 corresponds
% to the 
SO_times     = frames( SOs(:,3)); 
delta_times  = frames( delta_waves(:,3));

% Merge SOs and delta waves?
% prompt = sprintf('Merge SOs and delta waves (1)? ');
% merge_select = input(prompt);
% if merge_select == 1
%     % Designate the variable SO_times to contain both SO and delta troughs
%     % and use the output of SO analysis for all slow wave activity
%     SO_times = horzcat(SO_times, delta_times);
% end

delta_idx = delta_times;
SO_idx    = SO_times;
% preallocate 
[signal_delta_activity,signal_delta_activity_zscore]      = deal( zeros(length(delta_idx), nr_of_frames));
[signal_delta_activity_mean, signal_delta_activity_meanZ] = deal( zeros( size(signal,1), nr_of_frames));

[signal_SO_activity, signal_SO_activity_zscore]     = deal( zeros(length(SO_idx), nr_of_frames));
[signal_SO_activity_mean, signal_SO_activity_meanZ] = deal( zeros( size(signal,1), nr_of_frames));

[all_data_SO, all_data_SOZ]       = deal( zeros(size(rois_for_an,2), size(SO_idx,2), size(time,2)));
[all_data_delta, all_data_deltaZ] = deal( zeros(size(rois_for_an,2), size(delta_idx,2), size(time,2)));
%% Run main analysis
t = 1;
%loop over ROIs
for roinr = rois_for_an
 
    % Extract signal from currrent ROI
    roi_signal        = signal(roinr, :); 
    roi_signal_zscore = zscore_signal(roinr,:); 


    % Get all peri-delta activity windows from current ROI
    [signal_delta_activity,signal_delta_activity_zscore] = ...
    extract_avg_activity(delta_idx, roi_signal, roi_signal_zscore,signal_delta_activity,...
    signal_delta_activity_zscore, nr_of_seconds, opts);

    % Get all peri-SO activity windows from current ROI
    [signal_SO_activity,signal_SO_activity_zscore] = ...
    extract_avg_activity(SO_idx, roi_signal, roi_signal_zscore,signal_SO_activity,...
    signal_SO_activity_zscore, nr_of_seconds, opts);

    % store all roi x time x SWR data
    all_data_SO(t,:,:)     = signal_SO_activity;
    all_data_SOZ(t,:,:)    = signal_SO_activity_zscore;
    all_data_delta(t,:,:)  = signal_delta_activity;
    all_data_deltaZ(t,:,:) = signal_delta_activity_zscore;

    signal_SO_activity_mean(roinr,:)  = mean(signal_SO_activity,'omitnan');
    signal_SO_activity_meanZ(roinr,:) = mean(signal_SO_activity_zscore,'omitnan');

    signal_delta_activity_mean(roinr,:)        = mean(signal_delta_activity,'omitnan');
    signal_delta_activity_meanZ(roinr,:) = mean(signal_delta_activity_zscore,'omitnan');
    t = t+1;
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
    sort_idx = plot_so_sig_an(sessionID,                    opts, ...
                              time,                         signal_SO_activity,...
                              signal_SO_activity_zscore,    signal_SO_activity_mean, ...
                              signal_SO_activity_meanZ,     signal_delta_activity,...
                              signal_delta_activity_zscore, signal_delta_activity_mean,...
                              signal_delta_activity_meanZ,  celltype,...
                              label3,                       text,...
                              cells_to_exclude,             delta_idx,...
                            SO_idx,                         nr_of_frames);
end


%% Save output as struct
analysis_output = struct;

if checkParameter(opts.exp_type, 1 , 'bulk')
    analysis_output.signal_delta_activity        = signal_delta_activity;
    analysis_output.signal_delta_activity_zscore = signal_delta_activity_zscore;
    analysis_output.signal_SO_activity           = signal_SO_activity;
    analysis_output.signal_SO_activity_zscore    = signal_SO_activity_zscore;
    analysis_output.opts                         = opts;
    analysis_output.sort_idx                     = sort_idx;
%     analysis_output.roiClustIDs                  = roiClustIDs;

elseif checkParameter(opts.exp_type, 2 , 'default') || checkParameter(opts.exp_type, 3 , 'axon')
    analysis_output.signal_delta_activity_mean  = signal_delta_activity_mean(rois_for_an,:);
    analysis_output.signal_delta_activity_meanZ = signal_delta_activity_meanZ(rois_for_an,:);
    analysis_output.signal_SO_activity_mean     = signal_SO_activity_mean(rois_for_an,:);
    analysis_output.signal_SO_activity_meanZ    = signal_SO_activity_meanZ(rois_for_an,:);
    analysis_output.all_data_SO                 = all_data_SO;
    analysis_output.all_data_SOZ                = all_data_SOZ;
    analysis_output.all_data_delta              = all_data_delta;
    analysis_output.all_data_deltaZ             = all_data_deltaZ;
    analysis_output.opts                        = opts;
    analysis_output.sort_idx                    = sort_idx;
    analysis_output.roi_idx                     = rois_for_an;

end

