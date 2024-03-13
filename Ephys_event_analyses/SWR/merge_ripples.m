function sData = merge_ripples(sData)

% Written by Christoffer Berge | Vervaeke Lab

% Function that merges close SWRs. Close SWRs are here defined as SWRs with
% peaks located so close together such that, when setting all values in the
% smoothed, filtered LFP envelop below the upper threshold equal to zero,
% those peaks occur on the same "bump" or increase in the filtered LFP
% signal.

% Find SWR onset/offset times
[swr_start_stop, ~, ~] = mark_ripple_onset_offset(sData);

swr_idx      = sData.ephysdata.absRipIdx;
rippleSnips = sData.ephysdata.rippleSnips;

test = zeros(1, length(sData.ephysdata.lfp));

% for swr_nr = 1:length(swr_idx)
% 
%     tmp_start = swr_start_stop(swr_nr ,1);
%     tmp_stop  = swr_start_stop(swr_nr, 2);
% 
%     test(tmp_start:tmp_stop) = ones;
% end
% 
% figure, hold on
% plot(sData.ephysdata.lfp)
% plot(test)
% 
% for swr_nr = 2:length(swr_idx)
% 
% 
%     if swr_start_stop(swr_nr, 1) < swr_start_stop(swr_nr-1, 2) && (swr_start_stop(swr_nr, 1)-swr_start_stop(swr_nr-1, 2)) > 2500*.2
% 
%         swr_start_stop(swr_nr-1, 1) =  swr_start_stop(swr_nr-1, 1);
%         swr_start_stop(swr_nr-1, 2) =  swr_start_stop(swr_nr, 1);
%     end
% end


% Loop over SWRs starting from the second SWR
for swr_nr = 2:length(swr_idx)
    % check if SWR occurs within onset/ofset of previous SWR
    if swr_idx(swr_nr) > swr_start_stop(swr_nr-1,1) && swr_idx(swr_nr) < swr_start_stop(swr_nr-1,2)
        swr_idx(swr_nr) = 0;
        rippleSnips(swr_nr).lfp = NaN(2501,1);
    end
end
% test2 = test1;
frames = sData.daqdata.frame_onset_reference_frame;
% Delete overlapping SWRs
swr_idx(swr_idx ==0) = [];
% x = 3
sData.ephysdata.absRipIdx = swr_idx;
sData.ephysdata.frameRipIdx = frames(swr_idx);
sData.ephysdata.rippleSnips = rippleSnips;

