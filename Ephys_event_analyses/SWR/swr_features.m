function swr_feat = swr_nrem_features(varargin)

% Get various SWR features from data: nr of SWRs in awake vs. sleep,
% duration, amplitude, cycles, clusters. 


sData = varargin{1,1};

swr_idx      = sData.ephysdata.absRipIdx;
nr_nrem_swr  = sort([sData.ephysdata.NREM_spindle_uncoupled_swr, ...
               sData.ephysdata.spindle_coupled_swr]);
nrem_swr_idx = ismember(swr_idx, nr_nrem_swr   );

% NREM duration in recording
mean_sleep_features = sleepfeatures(sData);
total_nrem = mean_sleep_features{7};
nrem_swr_rate = length(nr_nrem_swr)/total_nrem;

% Duration
[~, swr_length]   = mark_ripple_onset_offset(sData);
nrem_swr_duration = swr_length(nrem_swr_idx)/2500;

% Amplitude
rippleAmpResp = rippleAmplitudeAnalysis([],sData,1,1,nrem_swr_idx);

% Clusters
cluster_idx = find_swr_clusters(sData);
clustered_swr = swr_idx(cluster_idx);

nr_nrem_clust_swr_idx = ismember(nr_nrem_swr, clustered_swr);
nr_nrem_clust_swr = sum(nr_nrem_clust_swr_idx);
proportion_clust_nrem_swr = nr_nrem_clust_swr/length(clustered_swr);
proportion_clust_swr = length(clustered_swr)/length(swr_idx);
% get indicies of first SWR in cluster
temp_var = diff(nr_nrem_clust_swr_idx);
temp_var(temp_var < 0) = 0;
idx_first_swr_in_clust = circshift(temp_var,1);
idx_first_swr_in_clust(1) = 0;

% test2 = diff(clustered_swr);
% test3 = test2 > 500;
% adTest3 = [false(1,1), test3];
% swr_feat = {nr_nrem_swr,
