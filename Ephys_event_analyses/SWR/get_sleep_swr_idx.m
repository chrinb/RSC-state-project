function swr_idx = get_sleep_swr_idx(sData, swr_select)

% Written by Christoffer Berge | Vervaeke Lab

% Function that that selects either all or a subset of SWRs in sleep
% recordings: awake SWRs, spindle-uncoupled SWRs or spindle-coupled SWRs.

if swr_select == 1 
    swr_idx = sData.ephysdata.absRipIdx;
elseif swr_select == 2
    swr_idx = sData.ephysdata.awake_swr;
elseif swr_select == 3 
    swr_idx = sData.ephysdata.NREM_spindle_uncoupled_swr;
elseif swr_select == 4 
    swr_idx = sData.ephysdata.spindle_coupled_swr;
end