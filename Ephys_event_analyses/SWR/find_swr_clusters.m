function [swr_cluster_idx, inter_swr_intervals] = find_swr_clusters(sData, RippleIdx)

% Written by Christoffer Berge | Vervaeke Lab

% Function to find indicies of SWR clusters.

if nargin == 1
    RippleIdx = sData.ephysdata.absRipIdx;
end

% Find inter-SWR-interval (in e-phys time samples)
inter_swr_intervals = diff(RippleIdx); 

% Find indicies of intervals > 500 samples (200ms). (The result is a vector
% that is 1 element shorter than "nr of SWR in rec vector" and it misses
% the indicies of the first SWR in a cluster. Some adjustments therefore
% has to be made.)
cluster_swr_idx = inter_swr_intervals < 500;

% Add one element to make vector of indicies equal in length to the nr of
% SWRs. This element is set = false, but that can be changed in the next
% step.
adjusted_cluster_swr_idx1 = [false(1,1), cluster_swr_idx];

% Now add one (false) element to the end of the cluster-SWR index vector in
% a new vector. This vector is now equal length as nr of SWRs in recording,
% and it captures the index of the first SWR in each cluster
adjusted_cluster_swr_idx2 = [cluster_swr_idx, false(1,1)];

% Concatenate the two cluster-SWR index vectors
temp_var = vertcat(adjusted_cluster_swr_idx1, adjusted_cluster_swr_idx2);

% Sum across columns and convert results to logical. The result is a vector
% of cluster-SWR indicies 
swr_cluster_idx = logical(sum(temp_var));
