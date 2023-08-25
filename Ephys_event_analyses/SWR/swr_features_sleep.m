function swr_nrem_feat = swr_features_sleep(varargin)

% Written by Christoffer Berge | Vervaeke Lab

% Function that computes various features of SWRs in sleep recordings:
% (1) Total nr of SWRs in recording
% (2) Proportion of clustered SWRs during NREM sleep
% (3) Overall proportion of clustered SWRs in recording
% (4) Length of NREM SWRs in seconds
% (5) NREM SWR per min
% (6) SWR amplitude (z-score)
% (7) Nr of spindle-coupled SWRs (SWR nested inside spindle)
% (8) Proportion if spindle-coupled SWR
% (9) Different SWR cluster types (e.g., doublets, triplets, etc.)
% (10) Nr of SWRs in each cluster type
% (11) Proportion of all SWRs in each cluster type 
% (12) Nr of clustered SWRs in NREM sleep
% (13) Proportion SWR clusters during NREM sleep
% (14) Nr of clustered SWRs during sleep spindles
% (15) Proportion SWR clusters during sleep spindles
% (16) Inter-SWR-intervals of all SWRs

% Load data 
sData = varargin{1,1};

% Find time points and indicies of SWRs
swr_idx                         = sData.ephysdata.absRipIdx;
nrem_swr_times                  = sort([sData.ephysdata.NREM_spindle_uncoupled_swr, ...
                                  sData.ephysdata.spindle_coupled_swr]);
nrem_swr_idx                    = ismember(swr_idx, nrem_swr_times   );
abs_swr_nr_in_rec               = length(sData.ephysdata.absRipIdx);
nr_spindle_coupled_swrs         = length(sData.ephysdata.spindle_coupled_swr);
proportion_spindle_coupled_swrs = nr_spindle_coupled_swrs/length(swr_idx);

% NREM duration in recording
mean_sleep_features   = sleepfeatures(sData);
total_nrem            = mean_sleep_features{7};
nrem_swr_rate_per_min = length(nrem_swr_times)/total_nrem;

% Duration
[~, swr_length]     = mark_ripple_onset_offset(sData);
nrem_swr_dur_in_sec = swr_length(nrem_swr_idx)/2500;

% Amplitude
rippleAmpResp = rippleAmplitudeAnalysis([],sData,1,1,nrem_swr_idx);

% Clusters
cluster_idx   = find_swr_clusters(sData);
clustered_swr = swr_idx(cluster_idx);

nrem_clust_swr_idx = ismember(nrem_swr_times, clustered_swr);
nr_nrem_clust_swr  = sum(nrem_clust_swr_idx);

% Determine proportion of clustered SWRs occurring during NREM and how
% many of the overall SWRs are clustered. 
proportion_clust_nrem_swr = nr_nrem_clust_swr/length(clustered_swr);
proportion_clust_swr      = length(clustered_swr)/length(swr_idx);

%% Now do a more detailed SWR cluster analysis
swr_clust                     = swr_cluster_an(sData);
[true_sing_idx, true_sing_nr] = true_singlet(sData);

nrem_clust                       = {};
nr_cluster_swr_in_nrem     = NaN(1,9);
nr_cluster_swr_in_spindles = NaN(1,9);

% Because the vector containing nr of SWR cluster types likely contain NaNs
% (for missing SWR cluster types), find the indicies of the SWR cluster
% types that are in the data
nr_of_cluster_types_idx  = ~isnan(swr_clust{3,1});
correct_cluster_type_idx = (swr_clust{3,1}(nr_of_cluster_types_idx))-1;
% loop over different SWR cluster types
for i = correct_cluster_type_idx

    % Check if any NREM SWR/spindle-coupled SWR also belongs to a cluster
    nrem_clust_log_idx = ismember(nrem_swr_times, swr_clust{1,i});
    spin_clust_log_idx = ismember(sData.ephysdata.spindle_coupled_swr, swr_clust{1,i} );
    
    % Store logical indicies as is in one cell array and in the other sum
    % to get the overall nr
    nrem_clust{i,1} = nrem_clust_log_idx;
    nrem_clust{i,2} = sum(nrem_clust_log_idx); 

    spin_clust{i,1} = spin_clust_log_idx;
    spin_clust{i,2} = sum(spin_clust_log_idx);

    % Find nr of SWR in each cluster during NREM sleep
    nr_cluster_swr_in_nrem(1,i)     = sum(nrem_clust_log_idx);
    nr_cluster_swr_in_spindles(1,i) = sum(spin_clust_log_idx);
end


% Proportion of different SWR clusters in NREM sleep (nr of SWR in each
% cluster during NREM sleep divided by total nr of NREM SWRs)
prop_swr_cluster_in_nrem  = nr_cluster_swr_in_nrem ./ length(nrem_swr_times);

% Proportion of different SWR clusters in sleep spindles (nr of SWR in each
% cluster during spindles divided by total spindle-coupled SWRs
prop_swr_cluster_in_spindles = nr_cluster_swr_in_spindles./ length(sData.ephysdata.spindle_coupled_swr);

% SWR cluster types (e.g., 2 = doublets, 3 = triplets etc.)
swr_cluster_types = swr_clust{3,1};

% Proportion of different SWR clusters in the session (nr of SWRs in each
% cluster divided by total nr of SWRs)
prop_swr_in_each_cluster = swr_clust{4, 1};

% Nr of SWRs cluster types (e.g, 14 doublets etc.)
nr_of_swr_cluster_types = swr_clust{5,1};

% Inter-SWR-intervals
ISI = swr_clust{7, 1};

%% Store output in cell array
swr_nrem_feat = {'Total nr of SWR in rec.',             'Proportion clustered NREM SWR', ...
                'Proportion clustered SWR',             'NREM SWR duration (sec)', ...
                'NREM SWR per min',                     'SWR amplitude (z-score)', ...
                'Nr of spindle-coupled SWRs',           'Proportion spindle-coupled SWR', ...
                'SWR cluster types',                    'Nr of SWR cluster types',...
                'Proportion SWR in each cluster',       'Nr of clustered SWRs in NREM',...
                'Proportion SWR clusters in NREM',      'Nr of clustered SWRs in spindles',...
                'Proportion SWR clusters in spindles',  'ISI',...
                'Nr of true singlets';
                 abs_swr_nr_in_rec,                  proportion_clust_nrem_swr, ...
                 proportion_clust_swr,               nrem_swr_dur_in_sec, ...
                 nrem_swr_rate_per_min,              rippleAmpResp.zScoreAmp,...
                 nr_spindle_coupled_swrs,            proportion_spindle_coupled_swrs, ...
                 swr_cluster_types,                  nr_of_swr_cluster_types,...
                 prop_swr_in_each_cluster,           nr_cluster_swr_in_nrem, ...
                 prop_swr_cluster_in_nrem,           nr_cluster_swr_in_spindles, ...
                 prop_swr_cluster_in_spindles,       ISI,...
                 true_sing_nr};