function [analysis_output] = swr_sleep_an(varargin)

% Written by Christoffer Berge | Vervaeke Lab

% Find mean/median SWR-aligned DF/F for bulk or population recordings in sleep 
% sessions. Takes as input a sData struct in the form of
% "sData.imdata.roiSignals(2).matrix" where "matrix" corresponds to a field
% containing a ROI x time matrix for a 2P recording session. User specifies
% whether to analyze bulk signals (averages all ROIs to a single mean
% trace) or population signals (loops over ROIs and calculates mean
% SWR-aligned activity). Results can be plotted in "plot_swr_sleep_an"
% function. 

sData = varargin{1,1};

% Anon function for evaluating options
checkParameter = @(param, n, str) (isnumeric(param) && param==n) || strcmp(param, str);

% Check if user specified particular options for analysis. If not use
% default options (check "getSwrAnalysisOptions" to see what they are). 
if nargin > 1
    if isstruct( varargin{1,2} )
        opts = varargin{1,2};
    else
        try
            opts.exp_type           = varargin{1,2};
            opts.signal_type        = varargin{1,3};
            opts.split_rois         = varargin{1,4};
            opts.swr_to_include     = varargin{1,5};
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

%% Select signal to be analyzed

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

% Various variables
nr_of_seconds  = 3;
frames         = sData.daqdata.frame_onset_reference_frame;
nr_of_frames   = (nr_of_seconds*31*2)+1;
time           = (-(31*nr_of_seconds):(31*nr_of_seconds))./31;
sessionID      = sData.sessionInfo.sessionID;
select_swr_idx = sData.ephysdata.absRipIdx; 

%% Select spindle freq
% prompt = sprintf('8-18Hz (1) or 10-16hz (2)? ');
% spindle_band_select = input(prompt);
% 
% if spindle_band_select == 1
%     spindle_select = [];
% elseif spindle_band_select == 2
%     spindle_select = '1016';
% end
spindle_select = [];
unclassified_swr_str           = strcat('unclassified_swr', spindle_select);
NREM_spindle_uncoupled_swr_str = strcat('NREM_spindle_uncoupled_swr', spindle_select);
spindle_coupled_swr_str        = strcat('spindle_coupled_swr', spindle_select);
spin_idx_str                   = strcat('NREMAbsSpindleIdx', spindle_select);
nr_of_spindles_in_analysis     = length(sData.ephysdata2.(spin_idx_str) );

%% Select SWRs for analysis

% Get specified SWR indicies
RippleIdx = get_swr_idx(opts.swr_for_analysis ,sData, select_swr_idx, opts);

unclassified_SWR_idx = ismember(RippleIdx, sData.ephysdata.(unclassified_swr_str) );
single_swr_idx       = ismember(RippleIdx, sData.ephysdata.(NREM_spindle_uncoupled_swr_str));
spindle_coupled_idx  = ismember(RippleIdx, sData.ephysdata.(spindle_coupled_swr_str));

% convert SWR time stamps from e-phys to 2P time
RippleIdx = frames(round(RippleIdx));

%% Preallocate variables

[signal_swr_activity,signal_swr_activity_zscore]  = deal(zeros( length(RippleIdx), nr_of_frames));

[signal_swr_activity_mean_unclassified,      signal_swr_activity_mean_unclassified_zscore ,...
    signal_swr_activity_median_unclassified, signal_swr_activity_mean_single,...
    signal_swr_activity_mean_single_zscore,  signal_swr_activity_median_single, ...
    signal_swr_activity_mean_coupled,        signal_swr_activity_mean_coupled_zscore,...
    signal_swr_activity_median_coupled]      = deal(zeros( size(signal,1), nr_of_frames));

all_data_zscored_unclassified = zeros(size(rois_for_an,2),  sum(unclassified_SWR_idx), size(time,2));
all_data_zscored_single       = zeros(size(rois_for_an,2),  sum(single_swr_idx),       size(time,2));
all_data_zscored_coupled      = zeros(size(rois_for_an,2),  sum(spindle_coupled_idx),  size(time,2));
%% Run main analysis
% signal = zscore(signal,0,2);
% loop over ROIs 
t = 1;
for roinr = rois_for_an
    roi_signal        = signal(roinr, :); %creates vector of the frame signals for a particualr roi
    roi_signal_zscore = zscore_signal(roinr,:);

    % Get all peri-SWR activity windows from current ROI
    [signal_swr_activity,signal_swr_activity_zscore] = ...
        extract_avg_activity(RippleIdx,roi_signal, roi_signal_zscore,signal_swr_activity,...
        signal_swr_activity_zscore,nr_of_seconds, signal, opts);

    % Take all SWR x time snippets for a ROI and sort them into different
    % matrices depending on their classification
    signal_swr_unclassified        = signal_swr_activity(unclassified_SWR_idx,:);
    signal_swr_unclassified_zscore = signal_swr_activity_zscore(unclassified_SWR_idx,:);

    signal_swr_single              = signal_swr_activity(single_swr_idx,:);
    signal_swr_single_zscore       = signal_swr_activity_zscore(single_swr_idx,:);

    signal_swr_coupled             = signal_swr_activity(spindle_coupled_idx,:);
    signal_swr_coupled_zscore      = signal_swr_activity_zscore(spindle_coupled_idx,:);
    
    % Store all PETHs for each ROI
    all_data_zscored_unclassified(t,:,:) = signal_swr_unclassified_zscore;
    all_data_zscored_single(t,:,:)       = signal_swr_single_zscore;
    all_data_zscored_coupled(t,:,:)      = signal_swr_coupled_zscore;
    t = t+1;
    % Calculate mean & median for each ROI for each SWR category
    signal_swr_activity_mean_unclassified(roinr,:)          = mean(signal_swr_unclassified, 'omitnan');
    signal_swr_activity_median_unclassified(roinr,:)        = median(signal_swr_unclassified, 'omitnan');

    signal_swr_activity_mean_unclassified_zscore(roinr,:)   = mean(signal_swr_unclassified_zscore, 'omitnan');
    signal_swr_activity_median_unclassified_zscore(roinr,:) = median(signal_swr_unclassified_zscore, 'omitnan');    

    signal_swr_activity_mean_single(roinr,:)                = mean(signal_swr_single, 'omitnan');
    signal_swr_activity_median_single(roinr,:)              = median(signal_swr_single, 'omitnan');

    signal_swr_activity_mean_single_zscore(roinr,:)         = mean(signal_swr_single_zscore, 'omitnan');
    signal_swr_activity_median_single_zscore(roinr,:)       = median(signal_swr_single_zscore, 'omitnan');

    signal_swr_activity_mean_coupled(roinr,:)               = mean(signal_swr_coupled, 'omitnan');
    signal_swr_activity_median_coupled(roinr,:)             = median(signal_swr_coupled, 'omitnan');

    signal_swr_activity_mean_coupled_zscore(roinr,:)        = mean(signal_swr_coupled_zscore, 'omitnan');
    signal_swr_activity_median_coupled_zscore(roinr,:)      = median(signal_swr_coupled_zscore, 'omitnan');
end

%% Plot results
plot_z = 1; % Plot z-scored data by default. Set to 0 for non-z-scored plots.

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
sort_idx = plot_swr_sleep_an(sessionID,                     opts, ...
                  time,                                     signal_swr_unclassified, ...
                  signal_swr_single,                        signal_swr_coupled, ...
                  signal_swr_unclassified_zscore,           signal_swr_single_zscore, ...
                  signal_swr_coupled_zscore,                signal_swr_activity_mean_unclassified,...
                  signal_swr_activity_mean_single,               signal_swr_activity_mean_coupled, ...
                  signal_swr_activity_mean_unclassified_zscore,  signal_swr_activity_mean_single_zscore,...
                  signal_swr_activity_mean_coupled_zscore,       opts.signal_type,...
                  plot_z,                                   celltype,...
                  cells_to_exclude,                         label3,...
                  text);
end

%% Store output as struct

analysis_output = struct;

if checkParameter(opts.exp_type, 1 , 'bulk')
    analysis_output.signal_swr_activity_unclassified        = signal_swr_unclassified;
    analysis_output.signal_swr_activity_unclassified_zscore = signal_swr_unclassified_zscore;
    analysis_output.signal_swr_activity_single              = signal_swr_single;
    analysis_output.signal_swr_activity_single_zscore       = signal_swr_single_zscore;
    analysis_output.signal_swr_activity_coupled             = signal_swr_coupled;
    analysis_output.signal_swr_activity_coupled_zscore      = signal_swr_coupled_zscore;
    analysis_output.sort_idx                       = sort_idx;
    analysis_output.roiClustIDs                    = roiClustIDs;
    analysis_output.opts                           = opts;

elseif checkParameter(opts.exp_type, 2 , 'default') || checkParameter(opts.exp_type, 3 , 'axon')
    analysis_output.signal_swr_activity_mean_unclassified        = signal_swr_activity_mean_unclassified(rois_for_an,:);
    analysis_output.signal_swr_activity_mean_unclassified_zscore = signal_swr_activity_mean_unclassified_zscore(rois_for_an,:);
    analysis_output.signal_swr_activity_mean_single              = signal_swr_activity_mean_single(rois_for_an,:);
    analysis_output.signal_swr_activity_mean_single_zscore       = signal_swr_activity_mean_single_zscore(rois_for_an,:);
    analysis_output.signal_swr_activity_mean_coupled             = signal_swr_activity_mean_coupled(rois_for_an,:);
    analysis_output.signal_swr_activity_mean_coupled_zscore      = signal_swr_activity_mean_coupled_zscore(rois_for_an,:);
    analysis_output.sort_idx                                = sort_idx;
    analysis_output.roiClustIDs                             = roiClustIDs;
    analysis_output.all_data_zscored_unclassified           = all_data_zscored_unclassified;
    analysis_output.all_data_zscored_single                 = all_data_zscored_single;
    analysis_output.all_data_zscored_coupled                = all_data_zscored_coupled;
    analysis_output.opts                                    = opts;
    analysis_output.roi_idx                                 = rois_for_an;
     
end
