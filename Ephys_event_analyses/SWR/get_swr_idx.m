function RippleIdx = get_swr_idx(swr_remove, sData, select_swr_idx, opts)

% Written by Christoffer Berge | Vervaeke Lab

% Function that selects SWRs for analysis: either (1) all SWRs, (2)
% removing SWRs during/close to movement, (3) removing clustered SWRs (but
% keeping the first in each cluster), or (4) removing movement & clustered
% SWRs. 

% Input:    Numerical (1-4) 
% Output:   SWR indicies in ephys time
proximity_thresh = opts.swr_cluster_thresh;
% Keep all SWRs
if swr_remove == 1
    RippleIdx = select_swr_idx;

% Remove SWRs occuring in/close to movement
elseif swr_remove == 2
    [RippleIdx,~] = remove_movement_swr(sData, select_swr_idx);

% Remove clustered SWRs
elseif swr_remove == 3
    [RippleIdx] = remove_clustered_swr(select_swr_idx, proximity_thresh);

% First remove movement-associated SWRs, then the remaining clustered SWRs
elseif swr_remove == 4
    [RippleIdx]   = remove_clustered_swr(select_swr_idx, proximity_thresh);
    [RippleIdx,~] = remove_movement_swr(sData, RippleIdx);
end