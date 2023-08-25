function sData = merge_ripples(sData)

% Written by Christoffer Berge | Vervaeke Lab

% Function that merges close SWRs. Close SWRs are here defined as SWRs with
% peaks located so close together such that, when setting all values in the
% smoothed, filtered LFP envelop below the upper threshold equal to zero,
% those peaks occur on the same "bump" or increase in the filtered LFP
% signal.

% Find SWR onset/offset times
[swr_start_stop, ~, ~] = mark_ripple_onset_offset(sData);

RipIdx = sData.ephysdata.absRipIdx;
test1 = sData.ephysdata.absRipIdx;
rippleSnips = sData.ephysdata.rippleSnips;

% Loop over SWRs starting from the second SWR
for swr_nr = 2:length(RipIdx)

    % check if SWR occurs within onset/ofset of previous SWR
    if RipIdx(swr_nr) > swr_start_stop(swr_nr-1,1) && RipIdx(swr_nr) < swr_start_stop(swr_nr-1,2)
        RipIdx(swr_nr) = 0;
        rippleSnips(swr_nr).lfp = NaN(2501,1);
    end
end
% test2 = test1;
frames = sData.daqdata.frame_onset_reference_frame;
% Delete overlapping SWRs
RipIdx(RipIdx ==0) = [];
% x = 3
sData.ephysdata.absRipIdx = RipIdx;
sData.ephysdata.frameRipIdx = frames(RipIdx);
sData.ephysdata.rippleSnips = rippleSnips;

