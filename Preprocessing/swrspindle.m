function [sData, NREMepisodes] = swrspindle(varargin)

% Written by Christoffer Berge || Vervaeke lab

% Function that classifies SWRs as (1) undetermined (if occurring outside
% of defined NREM episode, i.e. awake, drowsiness, awake-NREM transitions,
% REM sleep), (2) spindle-coupled, or (3) spindle-uncoupled.

sData = varargin{1,1};

if length(varargin) < 2
    prompt = sprintf('8-18Hz (1) or 10-16hz (2)? ');
    spindle_band_select = input(prompt);

    if spindle_band_select == 1
        spindle_select = [];
    elseif spindle_band_select == 2
        spindle_select = '1016';
    end
else
    spindle_select = varargin{1,2};
end

% find all NREM episodes
NREMep       = sData.episodeMerge.state == 'NREM' & sData.episodeMerge.state_duration > 30;
NREMarr      = table2array(sData.episodeMerge(NREMep,2:3));
NREMepisodes = NREMarr*2500;


spin_idx_str = strcat('NREMAbsSpindleIdx', spindle_select);
spin_start_end_str = strcat('NREMspindleStartEnd', spindle_select);

swrIdx = sData.ephysdata.absRipIdx;
spindle_coupled_swrIdx = [];
for k = 1:length(swrIdx)
    % loop over NREM episodes
    for l = 1:size(NREMepisodes,1)
        
        % check if SWR occurs within a NREM episode
        if swrIdx(k) > NREMepisodes(l,1) && swrIdx(k) < NREMepisodes(l,2)
            NREM_swrIdx(k) = swrIdx(k);
        end

        % loop over spindles
        for i = 1:length(sData.ephysdata2.(spin_idx_str) )
            
            % check if SWR occurs between NREM spindle onset/offset
           if swrIdx(k) > sData.ephysdata2.(spin_start_end_str)(i,1) && ...
           swrIdx(k) < sData.ephysdata2.(spin_start_end_str)(i,2) 
       
           spindle_coupled_swrIdx(k) = swrIdx(k);
           end
        end
    end
end
    
% find awake SWR
unclassified_swrIdx  = ~ismember(swrIdx, NREM_swrIdx);
unclassified_swrIdx1 = sData.ephysdata.absRipIdx(unclassified_swrIdx);

% separate spindle-coupled SWRs from remaining NREM SWRs
NREM_swrIdx(NREM_swrIdx == 0)                       = [];
spindle_coupled_swrIdx(spindle_coupled_swrIdx == 0) = [];
NREM_swrIdxLog                                      = ~ismember(NREM_swrIdx,spindle_coupled_swrIdx);
NREM_swrIdx                                         = NREM_swrIdx(NREM_swrIdxLog);

% separate spindle-coupled SWRs from remaining awake SWRs
unclassified_swrIdx2 = ~ismember(unclassified_swrIdx1, spindle_coupled_swrIdx);

unclassified_swr_str           = strcat('unclassified_swr', spindle_select);
NREM_spindle_uncoupled_swr_str = strcat('NREM_spindle_uncoupled_swr', spindle_select);
spindle_coupled_swr_str        = strcat('spindle_coupled_swr', spindle_select);

sData.ephysdata.(unclassified_swr_str)           = unclassified_swrIdx1(unclassified_swrIdx2);
sData.ephysdata.(NREM_spindle_uncoupled_swr_str) = NREM_swrIdx;
sData.ephysdata.(spindle_coupled_swr_str)        = spindle_coupled_swrIdx;

        

           
            
