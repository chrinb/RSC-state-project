function swr_awake_feat = swr_features_awake(varargin)

% Written by Christoffer Berge | Vervaeke Lab

% Get various awake SWR features
sData = varargin{1,1};

% Find time points and indicies of SWRs
swr_idx                         = sData.ephysdata.absRipIdx;
abs_swr_nr_in_rec               = length(sData.ephysdata.absRipIdx);
awake_swr_rate_per_min = abs_swr_nr_in_rec/ (length(sData.ephysdata.lfp)/2500/60);

% Duration
[~, swr_length] = mark_ripple_onset_offset(sData);
swr_dur_in_sec  = swr_length/2500;

% Amplitude
rippleAmpResp = rippleAmplitudeAnalysis([],sData,1,1,[]);

% Clusters
cluster_idx = find_swr_clusters(sData);
clustered_swr = swr_idx(cluster_idx);

% Determine proportion of clustered SWRs  
proportion_clust_swr      = length(clustered_swr)/length(swr_idx);

%% Now do a more detailed SWR cluster analysis
swr_clust                     = swr_cluster_an(sData);
[true_sing_idx, true_sing_nr] = true_singlet(sData);


% SWR cluster types (e.g., 2 = doublets, 3 = triplets etc.)
swr_cluster_types = swr_clust{3,1};

% Proportion of different SWR clusters in the session (nr of SWRs in each
% cluster divided by total nr of SWRs)
prop_swr_in_each_cluster = swr_clust{4, 1};

% Nr of SWRs cluster types (e.g, 14 doublets etc.)
nr_of_swr_cluster_types = swr_clust{5,1};

% Inter-SWR-intervals
ISI = swr_clust{7, 1};
%% Output
swr_awake_feat = {'Total nr of SWR in rec.',        'SWR per min', ...
                  'Proportion clustered SWRs',      'SWR duration (sec)', ...
                  'SWR amplitude (z-score)',        'SWR cluster types',...
                  'Nr of SWR cluster types',        'Proportion SWR in each cluster',... 
                  'ISI',                            'Nr of true singlets';
                  abs_swr_nr_in_rec,                 awake_swr_rate_per_min, ...
                  proportion_clust_swr,              swr_dur_in_sec, ...
                  rippleAmpResp.zScoreAmp,           swr_cluster_types, ...
                  nr_of_swr_cluster_types,           prop_swr_in_each_cluster, ...
                  ISI,                               true_sing_nr};
