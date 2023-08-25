function select_swr_idx = select_swr_for_mod_an(sData, opts)

% Written by Christoffe Berge | Vervaeke lab

% Select SWR indicies for significant modulation analysis. User inputs
% sData and opts struct with a subfield called "select_swr" which contains
% the string specifying which SWR indicies to use. 

if  strcmp(opts.select_swr, 'all')
    select_swr_idx = sData.ephysdata.absRipIdx;

elseif strcmp(opts.select_swr, 'awake only')
    select_swr_idx = sData.ephysdata.awake_swr; 

elseif  strcmp(opts.select_swr, 'all NREM')
    select_swr_idx = sort([sData.ephysdata.NREM_spindle_uncoupled_swr, sData.ephysdata.spindle_coupled_swr]); 

elseif  strcmp(opts.select_swr, 'spindle-uncoupled')
        select_swr_idx = sData.ephysdata.NREM_spindle_uncoupled_swr; 

elseif  strcmp(opts.select_swr, 'spindle-coupled')
    select_swr_idx = sData.ephysdata.spindle_coupled_swr; 

elseif strcmp(opts.select_swr, 'randomize')
    select_swr_idx = sort( randsample(size(sData.ephysdata.lfp,1), size(sData.ephysdata.absRipIdx,2) ))';
end