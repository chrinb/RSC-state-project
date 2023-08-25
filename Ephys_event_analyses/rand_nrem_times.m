function [rand_times, rand_window] = rand_nrem_times(sData)

% Written by Christoffer Berge

% Function that extracts snippets of mean DF/F from bulk recordings during random, 
% non-spindle epochs during NREM. The nr of snippets are equal 
% to nr of sleep spindles in recordings. The windoes are specified to
% fall inside NREM bouts but exclude sleep spindle epochs to allow for a 
% comparison of the mean DF/F signal during sleep spindles and random NREM time points. 

% random time window length
nr_of_seconds = 3;

% NREM episodes > 30s in recording
[~, NREMepisodes] = swrspindle(sData);

% Sleep spindle onset offset times
spindle_start_end = sData.ephysdata2.NREMspindleStartEnd;

% Create vector of zeros
nrem_nonspindle_times = false(1,length(sData.ephysdata2.lfp));

% Loop over NREM bouts and set them = 1
for i = 1:length(NREMepisodes)
    nrem_nonspindle_times(1,NREMepisodes(i,1):NREMepisodes(i,2)) = true;
end

% time = (0:length(sData.ephysdata2.lfp)-1) / 2500;

% loop over spindle start/end at set those epochs = 0
for i = 1:length(spindle_start_end)
    nrem_nonspindle_times(1, spindle_start_end(i,1):spindle_start_end(i,2)) = false;
end

% The resulting vector now contains 1s in non-spindle NREM bouts and 0s
% everywhere else. 

t = 0;
rand_times    = [];
rand_window   = [];
while t < 20

    % select random time point
    rand_time = randi(length(nrem_nonspindle_times));
    % extract snippet 
    rand_win =  [rand_time - (2500*nr_of_seconds), rand_time + (2500*nr_of_seconds)];
    
    % check that the window doesn't exceed start/end of recording
    if rand_win(1) > 0 && rand_win(2) < length(nrem_nonspindle_times)

        % extract the random time window
        rand_snippets = nrem_nonspindle_times( rand_win(1):rand_win(2) );
        
        % Check if vector contains 0 (indicating that it occurs outside
        % NREM bout or inside sleep spindle)
        if ismember(0, rand_snippets)
        
        else % if not, store the start/end of random window 
        rand_times  = vertcat(rand_times, rand_time);
        rand_window = vertcat(rand_window, rand_win);

        nrem_nonspindle_times( rand_win(1):rand_win(2) ) = 0;

        t = t + 1;
        end
    end
end
