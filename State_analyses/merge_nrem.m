function sData = merge_nrem(sData)

% Written by Christoffer Berge || Vervaeke lab

% Merge NREM bouts separated by micro awakenings/microarousals lasting less
% than 5 seconds. 

srate     = 2500;
threshold = 5;
% Get NREM sleep episodes
nrem_episodes = nrem_sleep(sData);

% Number of NREM episodes
n_nrem_ep     = length(nrem_episodes);

% Create two vectors, one containing all NREM episode onsets (except
% first), one containing all NREM episode offsets (except last). Subtract
% the offsets from the onsets and check if difference is smaller than
% microarousal threshold
nrem_start    = nrem_episodes(2:n_nrem_ep, 1);
nrem_end      = nrem_episodes(1:n_nrem_ep-1, 2);
micro_arousal = [nrem_start-nrem_end < srate*threshold; 0];

merged_nrem = [0, 0];
nrem_ep_nr = 1;
counter    = 1;

while nrem_ep_nr <= n_nrem_ep 

    log_idx = micro_arousal(nrem_ep_nr);
    
    % Check if no merge and episode does not end inside previous merge
    if log_idx == 0 && nrem_episodes(nrem_ep_nr, 1) >  merged_nrem(end)
        merged_nrem(counter, 1) = nrem_episodes(nrem_ep_nr, 1);
        merged_nrem(counter, 2) = nrem_episodes(nrem_ep_nr, 2);

        nrem_ep_nr = nrem_ep_nr + 1;
        counter    = counter + 1;
    elseif log_idx == 1 

        if nrem_episodes(nrem_ep_nr, 1) > merged_nrem(end)
            merged_nrem(counter, 1) = nrem_episodes(nrem_ep_nr, 1);
            merged_nrem(counter, 2) = nrem_episodes(nrem_ep_nr+1, 2);
            counter = counter + 1;
        else
            merged_nrem(end, 2) = nrem_episodes(nrem_ep_nr+1,2);
            counter = counter + 1;
        end
        nrem_ep_nr = nrem_ep_nr + 1;
    else
        nrem_ep_nr = nrem_ep_nr + 1;
    end
end

% Remove any zeros in array
merged_nrem = merged_nrem(~(merged_nrem == 0));
merged_nrem = reshape(merged_nrem, size(merged_nrem,1)/2, []);

% Loop over nr of NREM episodes (except the last)
% for nrem_ep_nr = 1:length(nrem_episodes)-1
%     
%     % Check if the duration between NREM episode end and beginning of next
%     % NREM episode is < 5s. If true, create a merged NREM episode
%     if length( nrem_episodes(nrem_ep_nr,2):nrem_episodes(nrem_ep_nr+1,1) ) < 2500*5 && ...
%             ~isnan(nrem_episodes(nrem_ep_nr,2))
%         merged_nrem(nrem_ep_nr,1) = nrem_episodes(nrem_ep_nr,1);
%         merged_nrem(nrem_ep_nr,2) = nrem_episodes(nrem_ep_nr+1,2);
%         nrem_episodes(nrem_ep_nr+1,:) = NaN;
%     else
%         merged_nrem(nrem_ep_nr,:) = nrem_episodes(nrem_ep_nr,:);
%     end
%     
%     % check if second-to-last and last NREM episode have been merged. If
%     % not, store last NREM episode in merged NREM episode array
%     if nrem_ep_nr == length(nrem_episodes)-1 && ... 
%             length( nrem_episodes(nrem_ep_nr,2):nrem_episodes(nrem_ep_nr+1,1) ) > 2500*5 && ...
%             ~isnan(nrem_episodes(nrem_ep_nr+1,2))
%         merged_nrem(nrem_ep_nr+1,:) = nrem_episodes(nrem_ep_nr+1,:);
%     end
% end
% 
% % Remove NaNs
% merged_nrem = rmmissing(merged_nrem);

%% Create new sData state table containing merged NREM episodes

% First find all non-NREM episodes (e.g., IS, REM episodes)
non_nrem_idx  = sData.episodes.state ~= 'NREM';
non_nrem_data = sData.episodes(non_nrem_idx, :);

% Create new table containing the merged NREM episodes
state          = categorical("NREM");
state          = repmat(state, length(merged_nrem),1);
state_start    = merged_nrem(:,1)/2500;
state_end      = merged_nrem(:,2)/2500;
state_duration = state_end-state_start;
new_tab        = table(state, state_start,state_end, state_duration);
new_tab        = [non_nrem_data; new_tab];

% Sort table to get the states in the correct temporal order
[~, sort_idx]      = sortrows(new_tab(:,2));

% Store new table in struct
sData.episodeMerge = new_tab(sort_idx,:);