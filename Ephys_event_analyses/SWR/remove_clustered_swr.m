function [adjustedRippleIdx] = remove_clustered_swr(swr_idx, secBetweenRipples)

% Written by Christoffer Berge | Vervaeke Lab

% Function that identifies clustered SWRs (i.e. inter-ripple interval
% smaller than a user-defined threshold) and removes all but the first SWR
% in that cluster. This code therefore does not capture "true" single SWRs.

%% Adjust interval between ripples

% Set inter-SWR-interval = 2 seconds unless otherwise specified
if isempty(secBetweenRipples)
    secBetweenRipples = 2;
end

samplesBetweenRipples = secBetweenRipples*2500;
adjustedRippleIdx     = swr_idx;

% Find nr of e-phys samples between each SWR
pre_swr_interval = diff(swr_idx);
riplocs          = pre_swr_interval < samplesBetweenRipples;
% test             = pre_swr_interval > samplesBetweenRipples;
% testA            = [true(1,1) test];
adjustedRiplocs  = [false(1,1), riplocs]; %adds a logical 0 to make the vector equal in length to adjustedrippleIdx

% Remove SWRs that occurs to close to the previous SWR. 
adjustedRippleIdx(adjustedRiplocs) = [];

%%  Plotting check
% [swr_start_stop, ~, ~] = mark_ripple_onset_offset(sData);
% 
% 
% swr_to_include = swr_start_stop(testA',:);
% swr_to_exclude = swr_start_stop(adjustedRiplocs',:);
% 
% time_vec = (0:length(sData.ephysdata.lfp)-1)/2500;
% 
% figure, plot(time_vec, sData.ephysdata.lfp), hold on
% 
% for l = 1:length(swr_to_include)
%     x = [ swr_to_include(l,1) swr_to_include(l,1) swr_to_include(l,2) swr_to_include(l,2)]/2500;
%     y = [-1 1 1 -1];
%     patch(x, y, [0 0.4470 0.7410], 'edgecolor', 'none', 'FaceAlpha', .7,'LineWidth',2);
% end
% 
% for l = 1:length(swr_to_exclude)
%     x = [ swr_to_exclude(l,1) swr_to_exclude(l,1) swr_to_exclude(l,2) swr_to_exclude(l,2)]/2500;
%     y = [-1 1 1 -1];
%     patch(x, y, 'green', 'edgecolor', 'none', 'FaceAlpha', .7,'LineWidth',2);
% end