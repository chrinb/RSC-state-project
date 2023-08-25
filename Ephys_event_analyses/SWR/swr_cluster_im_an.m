function swr_cluster_im_an(sData)

% Written by Christoffer Berge | Vervaeke Lab

% Function that plots the mean imaging data from SWR-clusters.

% Get SWR-cluster indicies
swr_cluster_idx = swr_cluster_an(sData,[]);

%% Select signal to be analyzed
prompt = sprintf('Bulk (1) or single ROI (2) data? ');
bulk_or_pop = input(prompt);

prompt = sprintf('Which signal? (1 = deconvolved | 2 = DF/F | 3 = other) ');
DFF_or_deconv = input(prompt); 

if DFF_or_deconv == 1
    signal = sData.imdata.roiSignals(2).ciaDeconvolved;
elseif DFF_or_deconv == 2
    signal = sData.imdata.roiSignals(2).newdff;
elseif DFF_or_deconv == 3
    prompt = sprintf('Type structname: '); 
    signal = input(prompt);
end

% If bulk signal average all ROIs and set nr of ROIs = 1
if bulk_or_pop == 1
    signal  = mean(signal); % If bulk analysis, average all ROIs
    correct_nr_rois = 1; % In bulk sessions all 16 ROIs are averaged together
else
    % if population analysis remove ROIs tagged as overexpressing cells/INs
    % NOTE: "remove_cells" function require "roi_arr" to be added to path
    correct_nr_rois = remove_cells; 
end

% Create variables for analysis
nr_of_seconds = 3; % nr of seconds before/after SWR peak
nr_of_frames  = (nr_of_seconds*31*2)+1; % nr of recording frames in SWR window
time          = linspace(-nr_of_seconds,nr_of_seconds,nr_of_frames); % time vector in seconds
frames        = sData.daqdata.frame_onset_reference_frame;
sessionID     = sData.sessionInfo.sessionID;

%% Select spindle freq
prompt = sprintf('8-18Hz (1) or 10-16hz (2)? ');
spindle_band_select = input(prompt);

if spindle_band_select == 1
    spindle_select = [];
elseif spindle_band_select == 2
    spindle_select = '1016';
end

unclassified_swr_str           = strcat('unclassified_swr', spindle_select);
NREM_spindle_uncoupled_swr_str = strcat('NREM_spindle_uncoupled_swr', spindle_select);
spindle_coupled_swr_str        = strcat('spindle_coupled_swr', spindle_select);
spin_idx_str                   = strcat('NREMAbsSpindleIdx', spindle_select);
nr_of_spindles_in_analysis     = length(sData.ephysdata2.(spin_idx_str) );

% If signal to be analyzed is deconvolved, specify whether to threshold or not
threshold = false(1,1);
if DFF_or_deconv == 1
    prompt = sprintf('Threshold deconvolved dF/F? (1) ');
    dothreshold = input(prompt);
    if dothreshold == 1
        threshold = true(1,1);
    end
end
    
% If signal to be analyzed is bulk signal, specify whether to do baseline 
% subtraction or not
baseSub = [];
if bulk_or_pop == 1
    prompt = sprintf('Do baseline subtraction? (1) | everything else = no) ');
    baseSub = input(prompt);
end

idx_vec = swr_cluster_idx{3, 1}-1;

% Loop over nr of SWR-clusters
for swr_cluster = idx_vec(~isnan(idx_vec))
    RippleIdx = swr_cluster_idx{8,swr_cluster};
    
    % convert SWR time stamps from e-phys to 2P time
    RippleIdx = frames(round(RippleIdx));

    % loop over ROIs 
    for roinr = correct_nr_rois
        roi_signal = signal(roinr, :); %creates vector of the frame signals for a particualr roi
        
        % loop over SWRs
        for swr_nr = 1:length(RippleIdx)
            swr_window_start = RippleIdx(swr_nr) - (nr_of_seconds*31); 
            swr_window_end   = RippleIdx(swr_nr) + (nr_of_seconds*31); 
    
            % skip SWRs occurring too early or late (i.e. the SWR-window
            % extends before or beyond recording length) 
            if swr_window_start > 1 && swr_window_end < size(signal,2)
                signal_swr_activity(swr_nr, :) = ...
                    roi_signal(swr_window_start:swr_window_end); 
            end
    
            % Baseline subtraction
            if baseSub == 1
                % Baseline = mean activity in -3 to -2 sec before ripple peak
                baselineDFF                   = nanmean(signal_swr_activity(swr_nr,1:31));
                signal_swr_activity(swr_nr,:) = signal_swr_activity(swr_nr,:)-baselineDFF;
            end
        
            % Threshold signal using K-means clustering to separate "noise" vs
            % true deconvolved events.
            if DFF_or_deconv == 1 && threshold == 1 
                [threshold_percentile, ~] = threshold_deconvolved(signal, roi_signal, roinr);
    
                % Remove deconvolved values below threshold
                ROI_idx = reshape(signal_swr_activity, 1,[]);
                ROI_idx(ROI_idx < threshold_percentile(roinr)) = 0;
                signal_swr_activity = reshape(ROI_idx, size(signal_swr_activity));
            end
        end
            % Calculate mean & median for each ROI for each SWR category
        signal_swr_pop_mean(roinr,:)   = nanmean(signal_swr_activity);
        signal_swr_pop_median(roinr,:) = nanmedian(signal_swr_activity);
    end
    swr_cluster_bulk{swr_cluster} = signal_swr_activity;
    swr_cluster_pop{swr_cluster} = signal_swr_pop_mean;
end

swr_cluster_bulk
swr_cluster_pop