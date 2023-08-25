function sData = delete_duplicate_swr(sData)


swr_abs_ripIdx   = sData.ephysdata.absRipIdx;
swr_snips        = sData.ephysdata.rippleSnips;
swr_frame_ripIdx = sData.ephysdata.frameRipIdx;


[~, doubleLocs] = unique(swr_abs_ripIdx);
x = 1:length(swr_abs_ripIdx);
x(doubleLocs) = [];

if length(x) > 1
    for i = 1:length(x)
        i = x(i);
        swr_abs_ripIdx(i) = 0;
        swr_snips(i).lfp = NaN(2501,1);
        swr_frame_ripIdx(i) = 0;
    end
    % delete duplicate spindles
    swr_abs_ripIdx( swr_abs_ripIdx == 0) =[];
    swr_frame_ripIdx( swr_frame_ripIdx == 0) = [];

    sData.ephysdata.absRipIdx   = swr_abs_ripIdx;
    sData.ephysdata.rippleSnips = swr_snips;  
    sData.ephysdata.frameRipIdx = swr_frame_ripIdx;

elseif length(x) == 1
    swr_abs_ripIdx(x) = [];
    swr_snips(x).lfp = NaN(2501,1);
    swr_frame_ripIdx(x) = [];

    
    sData.ephysdata.absRipIdx   = swr_abs_ripIdx;
    sData.ephysdata.rippleSnips = swr_snips;  
    sData.ephysdata.frameRipIdx = swr_frame_ripIdx;
end
