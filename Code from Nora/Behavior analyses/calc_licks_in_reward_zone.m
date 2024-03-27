function fraction_licks_reward_zone = calc_licks_in_reward_zone(sData)

% Written by Christoffer Berge | Vervaeke lab

% Calculate the fraction of licks in reward zone vs. outside of reward
% zone.

% Load lick data
lick_data = sData.behavior.binning.lickBinned;
n_bins    = sData.behavior.meta.nBins;

% Define reward zone: Mice are given 2 drops of water. The first is always
% given in bin nr 1, the second may vary a little. 
bin_size = sData.behavior.meta.binSize;
bin_vec  = 1:n_bins;

% First bin were water was given
min_bin = 1;

% Last bin (across entire session) were water was given
mean_water_given        = mean(sData.behavior.binning.waterGivenBinned);
max_bin                 = find(mean_water_given, 1, 'last');

% Determine reward center as the median of the first and last bin, rounding
% downwards because most rewards are given in the first bin. 
reward_center_bin = floor( median( [min_bin, max_bin] ));
reward_center_cm  = reward_center_bin*bin_size;

% Nr of bins before and after reward zone center that together consistutes
% the entire reward zone.
reward_zone_halfsize    = 10;
reward_center_extension = floor( reward_zone_halfsize/bin_size);

bin_vec_shifted = circshift(bin_vec, -(reward_center_bin-1) );

reward_zone_start_bin = circshift(bin_vec_shifted, reward_center_extension+5 );
reward_zone_start_bin = reward_zone_start_bin(1) ;
reward_zone_start_cm  = reward_zone_start_bin*bin_size;

reward_zone_end_bin   = circshift(bin_vec_shifted, -reward_center_extension );
reward_zone_end_bin   = reward_zone_end_bin(1);
% reward_zone_end_bin   = 1;

reward_zone_end_cm    = reward_zone_end_bin*bin_size;

% Calculate nr of binned licks in reward zone, outside of reward zone, and
% total nr of binned licks
reward_zone_binned_licks = [ lick_data(:, reward_zone_start_bin:n_bins), lick_data(:, 1:reward_zone_end_bin) ];
remaining_binned_licks   = lick_data(:, reward_zone_end_bin+1:reward_zone_start_bin-1);

total_binned_licks           = sum(lick_data, 'all');
sum_reward_zone_binned_licks = sum(reward_zone_binned_licks, 'all');
sum_remaining_binned_licks   = sum(remaining_binned_licks, 'all');


fraction_licks_reward_zone = sum_reward_zone_binned_licks/total_binned_licks;