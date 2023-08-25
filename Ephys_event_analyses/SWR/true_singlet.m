function [true_sing_idx, true_sing_nr] = true_singlet(sData, RippleIdx, plot1)

% Written by Christoffer Berge | Vervaeke Lab

% Function that identifies true singlet SWRs, i.e. SWRs that are not
% flanked by other SWRs either before or after it, in a time window
% specified by user. 
nr_of_seconds = 2;
fs            = 2500;
threshold     = fs*nr_of_seconds;

if isempty(RippleIdx)
    swr_loc       = sData.ephysdata.absRipIdx;
else
    swr_loc = RippleIdx;
end

nr_swr        = length(swr_loc);
true_sing_idx = false(1, nr_swr );
for i = 1:nr_swr
    
    % first SWR
    if i == 1
        inter_swr_interval = swr_loc(i+1)-swr_loc(i);
        if inter_swr_interval > threshold
            true_sing_idx(i) = true;
        end
    end
    % last SWR
    if  i == nr_swr
        inter_swr_interval  = swr_loc(i)-swr_loc(i-1);
        if inter_swr_interval < threshold
            true_sing_idx(i) = true;
        end
    end
    
    if i > 1 && i < nr_swr
    
        pre_swr_interval  = swr_loc(i)   - swr_loc(i-1);
        post_swr_interval = swr_loc(i+1) - swr_loc(i); 

        if pre_swr_interval > threshold && post_swr_interval > threshold
            true_sing_idx(i) = true;
        end
    end

end

% Nr of true singlets in rec
true_sing_nr = sum(true_sing_idx);

% Find SWR onset/offset
[swr_start_stop, ~, ~] = mark_ripple_onset_offset(sData);

%% Plot results


% Plot variables
time_vec              = (0:length(sData.ephysdata.lfp)-1)/2500;
true_singlets_to_plot = swr_start_stop(true_sing_idx,:);
other_swrs_to_plot    = swr_start_stop(~true_sing_idx,:);

if ~isempty(plot1)
    
    figure, 
    plot(time_vec, sData.ephysdata.lfp), hold on
    
    for l = 1:length(true_singlets_to_plot)
        x = [ true_singlets_to_plot(l,1) true_singlets_to_plot(l,1) true_singlets_to_plot(l,2) true_singlets_to_plot(l,2)]/2500;
        y = [-1 1 1 -1];
        h1 = patch(x, y, [0 0.4470 0.7410], 'edgecolor', 'none', 'FaceAlpha', .7,'LineWidth',2);
    end
    
    for l = 1:length(other_swrs_to_plot)
        x = [ other_swrs_to_plot(l,1) other_swrs_to_plot(l,1) other_swrs_to_plot(l,2) other_swrs_to_plot(l,2)]/2500;
        y = [-1 1 1 -1];
        h2 = patch(x, y, 'magenta', 'edgecolor', 'none', 'FaceAlpha', .7,'LineWidth',2);
    end
    
    set(gca,'xlim',[0 time_vec(end)]);
    legend( [h1 h2], 'True singlets', 'Non-singlets');
end
     