function rip_clus_an(varargin)

% Calculate mean activity (e.g., DF/F or deconvolved DF/F) for single vs.
% clustered SWRs.

sData                 = varargin{1,1};
swr_idx               = sData.ephysdata.absRipIdx;
nr_nrem_swr           = sort([sData.ephysdata.NREM_spindle_uncoupled_swr, ...
                        sData.ephysdata.spindle_coupled_swr]);
cluster_idx           = find_swr_clusters(sData, swr_idx);
idx_single_swr        = ~cluster_idx;
idx_clustered_swr     = swr_idx(cluster_idx);
nr_nrem_clust_swr_idx = ismember(nr_nrem_swr, idx_clustered_swr);

%% Specify imaging data
prompt = sprintf('Bulk (1) or single ROI (2) data? ');
select = input(prompt);

prompt = sprintf('All (1) or NREM (2) SWRs? ');
swr_select = input(prompt);

prompt = sprintf('Which signal? (1 = deconvolved | 2 = DF/F | 3 = other) ');
Signal_select = input(prompt);      
if Signal_select == 1
    signal = sData.imdata.roiSignals(2).ciaDeconvolved;
    text = 'deconvolved';
elseif Signal_select == 2
    signal = sData.imdata.roiSignals(2).newdff;
    text = [];
elseif Signal_select == 3
    prompt = sprintf('Type structname: '); 
    signal = input(prompt);
    text = [];
end

% If bulk imaging average across all ROIs before analysis
if select == 1
    signal = mean(signal);
    label1 = ['Mean ' text, 'DF/F'];
    label2 = '# SWR';
else
    label1 = ['Mean ' text, 'DF/F'];
    label2 = '# ROI';
end

nr_of_seconds = 3;
frames        = sData.daqdata.frame_onset_reference_frame;
swr_idx       = sData.ephysdata.absRipIdx;
nr_of_frames  = (nr_of_seconds*31*2)+1;
time          = (-(31*nr_of_seconds):(31*nr_of_seconds))./31;

% get indicies of first nrem SWR in cluster
idx_first_swr_in_clust = diff(nr_nrem_clust_swr_idx);
idx_first_swr_in_clust = [false(1,1), idx_first_swr_in_clust];
idx_first_swr_in_clust(idx_first_swr_in_clust < 0) = 0;

% get indicies of first SWR in cluster (all SWRs, awake & NREM)
idx_first_swr_in_clust_all = diff(cluster_idx);
idx_first_swr_in_clust_all = [false(1,1), idx_first_swr_in_clust_all];
idx_first_swr_in_clust_all(idx_first_swr_in_clust_all < 0) = 0;


if swr_select == 1
    RippleTimes_clust = swr_idx(logical(idx_first_swr_in_clust_all));
    RippleIdx_clust   = idx_first_swr_in_clust_all;
else
    RippleTimes_clust = nr_nrem_swr(logical(idx_first_swr_in_clust));
    RippleIdx_clust   = idx_first_swr_in_clust;

end

% convert to 
RippleTimes_clust = frames(round(RippleTimes_clust));
% RippleIdx_sing  = frames(round(idx_single_swr));
baseSub = [];
% Specify whether to threshold deconvolved dF/F or not
if Signal_select == 1
    prompt = sprintf('Threshold deconvolved dF/F? (1) ');
    dothreshold = input(prompt);
    if dothreshold == 1
        threshold = true(1,1);
    else
        threshold = false(1,1);
    end
elseif select == 1
    prompt = sprintf('Do baseline subtraction? (1) | everything else = no) ');
    baseSub = input(prompt);
end
RippleIdx = remove_movement_swr(sData, swr_idx(idx_single_swr));

signal_clust_swr_activity     = zeros(length(RippleTimes_clust), nr_of_frames);
signal_clust_swr_activity_pop = zeros(size(signal,1), nr_of_frames);

%% Run main analysis
for roinr = 1:size(signal,1) %runs through the nr of coloumns (corresponding to nr of rois)
    roi_signal = signal(roinr, :); %creates vector of the frame signals for a particualr roi

    for swr_nr = 1:length(RippleTimes_clust)
        % find SWR window start end
        swr_window_start = RippleTimes_clust(swr_nr) - (nr_of_seconds*31); 
        swr_window_end   = RippleTimes_clust(swr_nr) + (nr_of_seconds*31); 

        % skip SWRs at the beginning or end with a window shorter than
        % length specified in seconds by user above
        if swr_window_start > 1 && swr_window_end < size(signal,2)
            signal_clust_swr_activity(swr_nr, :) = ...
                roi_signal(swr_window_start:swr_window_end); 
        end
        
        % Baseline subtraction
        if baseSub == 1
            % Baseline = mean activity in -3 to -2 sec before ripple peak
            baselineDFF = nanmean(signal_clust_swr_activity(swr_nr,1:31));
            signal_clust_swr_activity(swr_nr,:) = signal_clust_swr_activity(swr_nr,:)-baselineDFF;
        end
    end
    % Threshold signal using K-means clustering to separate "noise" vs
    % true deconvolved events.
    if Signal_select == 1 && threshold == 1 
        [threshold_percentile, ~] = threshold_deconvolved(signal, roi_signal, roinr);

        % Remove deconvolved values below threshold
        ROI_idx = reshape(signal_clust_swr_activity, 1,[]);
        ROI_idx(ROI_idx < threshold_percentile(roinr)) = 0;
        signal_clust_swr_activity = reshape(ROI_idx, size(signal_clust_swr_activity));
    end
        
    signal_clust_swr_activity_pop(roinr,:) = nanmean(signal_clust_swr_activity);
end



% for plotting bulk analysis
signal_swr_activity_zscore = zscore(signal_clust_swr_activity,0,2);
signal_bulk_SE             = std(signal_clust_swr_activity, 'omitnan') ./ ...
    sqrt(size(signal_clust_swr_activity,1));
signal_bulk_SE_zscore      = std(signal_swr_activity_zscore, 'omitnan') ./ ...
    sqrt(size(signal_swr_activity_zscore,1));

% for plotting population analysis
signal_swr_activity_awakeZ = zscore( signal_clust_swr_activity_pop, 0,2);
signal_pop_SE             = std(signal_clust_swr_activity_pop, 'omitnan') ./ ...
    sqrt(size(signal_clust_swr_activity_pop,1));
signal_pop_SE_zscore      = std(signal_swr_activity_awakeZ, 'omitnan') ./ ...
    sqrt(size(signal_swr_activity_awakeZ,1));

if select == 1
    mean_signal = {signal_clust_swr_activity, signal_swr_activity_zscore};
elseif select == 2
    mean_signal = {signal_clust_swr_activity_pop, signal_swr_activity_awakeZ};
end


if length(varargin) > 1
    figure,
    if select == 1

        subplot(221)
        x1 = [time(1), time(end)];
        y1 = [1 size(signal_clust_swr_activity,1)];
        imagesc(x1, y1, signal_clust_swr_activity) 
        ylabel(label2)
        xlabel('Time from SWR peak (sec)')

        subplot(222)
        imagesc(x1, y1, signal_swr_activity_zscore) 
        ylabel(label2)
        xlabel('Time from SWR peak (sec)')

        subplot(223)
        shadedErrorBar(time,nanmean(signal_clust_swr_activity),  signal_bulk_SE);
        xlabel('Time from SWR peak (sec)')
        ylabel(label1)
        set(gca, 'xlim',[time(1) time(end)])

        subplot(224)
        shadedErrorBar(time,nanmean(signal_swr_activity_zscore), signal_bulk_SE_zscore);
        xlabel('Time from SWR peak (sec)')
        ylabel(label1)
        set(gca, 'xlim',[time(1) time(end)])

    elseif select == 2

        subplot(221)
        x1 = [time(1), time(end)];
        y1 = [1 size(signal_clust_swr_activity_pop,1)];
        imagesc(x1, y1, signal_clust_swr_activity_pop) 
        ylabel(label2)
        xlabel('Time from SWR peak (sec)')

        subplot(222)
        imagesc(x1, y1, signal_swr_activity_awakeZ) 
        ylabel(label2)
        xlabel('Time from SWR peak (sec)')

        subplot(223)
        shadedErrorBar(time,nanmean(signal_clust_swr_activity_pop),  signal_pop_SE);
        xlabel('Time from SWR peak (sec)')
        ylabel(label1)
        set(gca, 'xlim',[time(1) time(end)])

        subplot(224)
        shadedErrorBar(time,nanmean(signal_swr_activity_awakeZ), signal_pop_SE_zscore);
        xlabel('Time from SWR peak (sec)')
        ylabel(label1)
        set(gca, 'xlim',[time(1) time(end)])
    end
end

