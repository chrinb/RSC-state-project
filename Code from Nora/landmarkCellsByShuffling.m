function lcs = landmarkCellsByShuffling(rmaps,varargin)
%% landmarkCellsByShuffling
%
% INPUT
%   -- Required: rmaps of place cells
%
% OUTPUT
%      IDs of obcject location cells
%

% set seed to make reproducible
rng(5)

%% Default parameters
if ~isempty(varargin)
    params = varargin{1};
else
    params.first_landmark_bins = 23:43;
    params.second_landmark_bins = 63:83;
end

params.percentile = 95;


% Update parameters
% Init empty var
lcs = [];

%% Shuffle data
for c = 1:size(rmaps,3)
    rmap = rmaps(:,:,c);
    response = nanmean(rmap);

    % Shuffle 1000 times
    shuffled_avg = zeros(1000,size(response,2));
    for i = 1:1000

        % Init zeroed matrix
        shuffled_rmap = zeros(size(rmap,1),size(rmap,2));

        % Shuffle for each lap
        for l = 1:size(rmap,1)
           shuffled_rmap(l,:) = circshift(rmap(l,:),randi(size(response,2))); 
        end

        shuffled_avg(i,:) = nanmean(shuffled_rmap);            

    end

    baseline_threshold = prctile(shuffled_avg,params.percentile);
    
    % First landmark
    response_first_position = [];
    for b = params.first_landmark_bins
        
        if response(b)>baseline_threshold(b)
            response_first_position = [response_first_position, b];
        end
    end
    
    % Second landmark
    response_second_position = [];
    for b = params.second_landmark_bins
        if response(b)>baseline_threshold(b)
            response_second_position = [response_second_position, b];
        end
    end
    
    response_outside = [];
    
    for b = setdiff(1:105,unique([params.first_landmark_bins,params.second_landmark_bins]))
        if response(b)>baseline_threshold(b)
            response_outside = [response_outside, b];
        end
    end
    
    
    if length(response_first_position)>1
        if length(response_second_position)>1
                lcs = [lcs,c];
        end
    end



end



