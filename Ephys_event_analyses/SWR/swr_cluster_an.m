function swr_clust = swr_cluster_an(sData, plot_data)

% Written by Christoffer Berge | Vervaeke Lab

% Find different SWR-cluster types (e.g., singles, doubles etc.) in 
% recording 

%% Select SWRs for analysis
% prompt = sprintf('All SWRs (1) | Rem. locom. SWRs (2) | Rem. clust. SWRs (3) | Rem- locom. & clust. SWRs (4) ');
% swr_for_analysis = input(prompt);
% 
% RippleIdx     = sData.ephysdata.absRipIdx;
% 
% if swr_for_analysis == 2
%     [RippleIdx,~] = riprun2(sData, RippleIdx);
% elseif swr_for_analysis == 3
%     [RippleIdx]   = RemoveRip(RippleIdx,[],[]);
% elseif swr_for_analysis == 4
%     [RippleIdx,~] = riprun2(sData, RippleIdx);
%     [RippleIdx]   = RemoveRip(RippleIdx,[],[]);
% end

RippleIdx = get_swr_idx(swr_remove, sData, select_swr_idx, opts)

% Find logical indices of clustered SWRs (occurring within 200ms of each other)
[all_swr_cluster_idx, inter_swr_intervals] = find_swr_clusters(sData,RippleIdx);

% Find the inter-SWR-intervals of the clustered SWRs and their time stamps 
swr_cluster_idx          = RippleIdx(all_swr_cluster_idx);
cluster_swr_inter_interv = [1, diff(swr_cluster_idx)];

tempV          = [];
i              = 1;
cluster_id_vec = [];
while i <= length(swr_cluster_idx)

    % Check if the clustered inter-SWR-interval is smaller than 500 samples
    % (200 ms). If true, add the index of that interval to vector tempV and
    % increment i with 1 to make while loop progress.
    if cluster_swr_inter_interv(i) < 500 
        tempV = [tempV, i];
        i = i+1;
    % If not true (indicating end of cluster): 
    % (1) find nr of SWR in cluster
    % (2) create a vector with elements having the value
    % same as length of cluster, 
    % (3) add this vector to a vector containing the
    % identity of all SWR clusters
    % (4) set the value of the clustered
    % inter-SWR-interval index (i) = 1 such that while loop now continues
    % from that index and looks for inter-SWR-intervals > 500 from that
    % index
    % (5) go back to the start of the current cluster and set the first
    % element to be equal to the nr of SWR in cluster
    else
        nr_swr_in_clus              = length(tempV);
        cluster_identity            = ones(1 ,nr_swr_in_clus)*(nr_swr_in_clus);
        cluster_id_vec              = [cluster_id_vec cluster_identity];
        cluster_swr_inter_interv(i) = 1;
        cluster_swr_inter_interv(i - nr_swr_in_clus) = nr_swr_in_clus;
        tempV = [];
    end
    
    % In the case that the last SWR is clustered do the following to add
    % their indicies to the "cluster_id_vec" vector
    if i == (length(swr_cluster_idx) + 1) && cluster_swr_inter_interv(i-1) < 500 
        nr_swr_in_clus              = length(tempV);
        cluster_identity            = ones(1 ,nr_swr_in_clus)*(nr_swr_in_clus);
        cluster_id_vec              = [cluster_id_vec cluster_identity];
        cluster_swr_inter_interv(i - nr_swr_in_clus) = nr_swr_in_clus;
        tempV = [];
    end
    
end


% Find the indicies of the first SWR in a cluster (locs) and the width of the peak (w)
% corresponding to how many SWRs occuring in that cluster (i.e. 4 clustered
% SWRs = width = 4)
% 
% [~, locs, w,~] = findpeaks(double(swr_cluster_idx), 'minPeakHeight',0.5);


% Find unique elements in cluster_id_vec (corresponding to cluster types), 
% count them using "accumarray" (see Matlab page on "unique" for for info).
[swr_cluster_types_in_rec,~,cluster_idx] = unique(cluster_id_vec');
n_counts                                 = accumarray(cluster_idx,1);

% Create two zero-vectors that can contain data on up to 9 SWR cluster
% types (e.g., from doublets to decets). Hard limit set to 9 because more 
% than 10 clustered SWRs are extremely rare/unlikely.
[swr_cluster_types, nr_swr_per_cluster] = deal( zeros(1,9) );

% Types of SWRs in recording (e.g. 2 = doublets, 3 = triplets, etc.)
swr_cluster_types( swr_cluster_types_in_rec-1 ) = swr_cluster_types_in_rec;
swr_cluster_types( swr_cluster_types == 0 )     = NaN; 

% Nr of SWRs per cluster type, starting with doublets
nr_swr_per_cluster( swr_cluster_types_in_rec-1 )     = n_counts;
nr_swr_per_cluster( nr_swr_per_cluster ==0 ) = NaN; 

% Nr of different SWR cluster types
nr_of_swr_cluster_types = nr_swr_per_cluster./swr_cluster_types;

% Proportion of SWRs in each cluster type (all clustered SWRs / total nr of
% SWRs)
prop_swr_in_each_cluster = nr_swr_per_cluster ./ length( sData.ephysdata.absRipIdx);

%% Now find the indicies of all SWRs in a cluster

% First, determine nr of SWRs to add to the clusters (1 less than widt of
% peak w)
% nr_of_swrs_to_add  = w-1;

% Create new variable to contain the indicies of all clustered SWRs (not just
% the first in a cluster)
% locs_with_all_swrs = locs;
% 
% % loop over nr of 1st SWRs in a cluster 
% swrs_to_add = [];
% IC          = [];
% for i = 1:length(locs)
% 
%     % determine nr of SWRs to add (scalar/vector)
%     nr_swrs_to_add = 1:nr_of_swrs_to_add(i);
%     test = ones(1, nr_of_swrs_to_add(i))*nr_of_swrs_to_add(i);
%     % create new variable containing the indicies of SWRs to add
%     swrs_to_add = horzcat(swrs_to_add, locs_with_all_swrs(i) + [nr_swrs_to_add]);
%     IC          = [IC, horzcat(ic(i), test)];
% end
% 
% % Sort the new variable containing the indicies of clustered SWRs
% locs_with_all_swrs = sort([locs_with_all_swrs, swrs_to_add]);

swr_cluster_enumerate = 1:length(all_swr_cluster_idx);
locs_with_all_swrs    = swr_cluster_enumerate(all_swr_cluster_idx);

%% Find SWR onset/offset
[swr_start_stop,~, ~] = mark_ripple_onset_offset(sData);

% Find time stamps and duration of SWR clusters
clust_swr_times = sData.ephysdata.absRipIdx(locs_with_all_swrs);
clust_swr_dur   = swr_start_stop(locs_with_all_swrs, :);

% Find the duration of non-cluster SWRs: (1) check which SWR-durations are
% not part of the cluster-SWR durations, then (2) index using those
% non-clustered SWR indicies. Check if there are any clustered SWRs. If
% not, set single SWR idx = RippleIdx
if ~isempty(clust_swr_dur)
    single_swr_dur_indx = ~ismember(swr_start_stop,clust_swr_dur);
    single_swr_dur      = swr_start_stop(single_swr_dur_indx(:,1),:);

else 
    single_swr_dur_indx = ismember(sData.ephysdata.absRipIdx, RippleIdx);
    single_swr_dur      = swr_start_stop(single_swr_dur_indx',:);
end


cluster_swr_idx = cluster_id_vec-1;
swr_clust = [];

% Loop over nr of cluster types in session
for i = (swr_cluster_types_in_rec-1)'

    % Select a particular cluster tyoe
    clust_swr = cluster_swr_idx' == i;
    % store time stamps of different SWR clusters
    swr_clust{1,i} = clust_swr_times(clust_swr);
    % store the onset/offset of different SWR clusters
    swr_clust{2,i} = clust_swr_dur(clust_swr', :);
    % Store the time stamp of the first SWR in the current cluster type
    swr_clust{8,i} = swr_cluster_idx(cluster_swr_inter_interv == (i+1));
end

% Store the nr of different clusters, proportion of SWRs that belonged to 
% specific SWR clusters, the nr of SWRs in each cluster, indicies of 
% non-clustered SWRs (NOTE: this is not "true singlets"), and the overall 
% inter-SWR-interval (ISI) in recording 
swr_clust{3,1} = swr_cluster_types;
swr_clust{4,1} = prop_swr_in_each_cluster;
swr_clust{5,1} = nr_of_swr_cluster_types;
swr_clust{6,1} = sData.ephysdata.absRipIdx(single_swr_dur_indx(:,1));
swr_clust{7,1} = inter_swr_intervals';
%% Plot results

% Find start/end of NREM bouts and convert to seconds.
% [~, NREMstartend] = swrspindle(sData);
% NREMstartend      = NREMstartend./2500;


if ~isempty(plot_data)
    % Plotting variables
    time_vec   = (0:length(sData.ephysdata.lfp)-1)./2500;
    low_patchV = abs(mean(sData.ephysdata.lfp));
    ylim_low   = mean(sData.ephysdata.lfp)-1;
    ylim_high  = mean(sData.ephysdata.lfp)+1; 
    patchC     = ["r", "g", "b", "c", "m", "y", "k", ""];
    legendName = ["singlet", "doublet", "triplet", "quartet", "quintet", "sextet", "septet", "octet", "nonet", "decet" ]; 
    h          = [];
    
    figure,
    plot(time_vec, sData.ephysdata.lfp), hold on
    % plot(upper_thresh)
    
%     % Loop over NREM bouts and plot them as individual patches
%     for i = 1:length(NREMstartend)
%         x = [NREMstartend(i,1) NREMstartend(i,1) NREMstartend(i,2) NREMstartend(i,2)];
%         y = [-.6 .6 .6 -.6];
%         patch(x, y, 'green', 'edgecolor', 'none', 'FaceAlpha', .3);
%     end

    % Plot SWR singlets
    for nr_swr_per_cluster = 1:size(single_swr_dur,1)
        x = [ single_swr_dur(nr_swr_per_cluster,1) single_swr_dur(nr_swr_per_cluster,1) single_swr_dur(nr_swr_per_cluster,2) single_swr_dur(nr_swr_per_cluster,2)]/2500;
        y = [low_patchV 1 1 low_patchV];
        h{1} = patch(x, y, [0 0.4470 0.7410], 'edgecolor', 'none', 'FaceAlpha', .7,'LineWidth',2);
    end
    
    % Plot SWR clusters
    correct_idx = ~isnan(swr_cluster_types);
    for i = (swr_cluster_types(correct_idx))-1
        
        for swr_cluster_start_end = 1:size(swr_clust{2,i},1)
        x = [ swr_clust{2,i}(swr_cluster_start_end,1) swr_clust{2,i}(swr_cluster_start_end,1) ...
              swr_clust{2,i}(swr_cluster_start_end,2) swr_clust{2,i}(swr_cluster_start_end,2)]/2500;
        y = [low_patchV 1 1 low_patchV];
        h{i+1} = patch(x, y, patchC(i), 'edgecolor', 'none', 'FaceAlpha', .7,'LineWidth',2);
        end
    end
    
    % Adjust Y and X limits
    set(gca,'ylim',[ylim_low ylim_high], 'xlim',[0 time_vec(end)]);
    
    % Add legends
    i = [1, (swr_cluster_types(correct_idx))];
    legend( [h{i}], legendName(i) ), hold on

end

