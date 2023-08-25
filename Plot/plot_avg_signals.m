function plot_avg_signals(varargin)

% Function that plots e-phys event (e.g., SWRs, sleep spindles, SOs)
% aligned data as (1) an event-by-event (trial-by-trial) heatmap and
% associated mean + SE. This could for example be imaging or ECoG/LFP data
% aligned to SWRs. 

% Different scenarios:
% awake: 2x2 plo
% sleep spindles
% sleep either 2x3 (3x3 for size adjustments) or 2x6 depending on whether
% to show both "raw" and normalized data side by side for the three SWR
% types

