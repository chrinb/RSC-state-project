function idx_modulated_rois = shuffle_analysis(varargin)

% Written by Christoffer Berge | Vervaeke Lab

% Shuffle analysis to find ROIs that are modulated by some event (e.g. SWRs
% or sleep spindles). Uses an adapted version of the algorithm
% outlined in Jadhav et al.(2012), Rothschild et al.(2016), 
% Alexander et al.(2018). 

win_length                     = varargin{1,1};
nr_of_shuffles                 = varargin{1,2};
peri_event_activity_single_ROI = varargin{1,3};
roinr                          = varargin{1,4};
idx_modulated_rois             = varargin{1,5};
% mean_activity                  = varargin{1,6};
time_win                       = varargin{1,7};
min_nr_events                  = varargin{1,8};

% initialize variables
mean_activity_all_shuffle_iterations = zeros(nr_of_shuffles, win_length);

% Exclude ROIs with fewer than 5 deconvolved dF/F events in the PETH
check_nr_of_events = peri_event_activity_single_ROI > 0.05;
sum_events         = sum(check_nr_of_events,2);
event_check        = sum(sum_events >=1); 

if event_check >= min_nr_events

    % Smooth convolved PETH: 1 = no smoothing. 
    smoothing_factor                    = 10;
    peri_event_activity_single_ROI      = smoothdata( peri_event_activity_single_ROI, 2, 'gaussian',smoothing_factor);
    mean_peri_event_activity_single_ROI = mean(peri_event_activity_single_ROI);

    % Create shuffled matrix and its mean
    for shuffle_iteration = 1:nr_of_shuffles

        % vector of random integers specifying the shifts to be applied
        % to each row of the PETH
        id = randi( win_length, 1, size(peri_event_activity_single_ROI,1) );

        % shuffle the deconvolved activity vector in each row in the PETH
        % by the vector of random integers (preserves the 
        % temporal structure of the activity
        shuffled_event_activity = cell2mat( arrayfun(@(x) circshift(peri_event_activity_single_ROI(x,:),[1 id(x)]), (1:numel(id))', 'UniformOutput', false));

        % Calculate and store the mean of the each shuffled PETH activity 
        mean_activity_all_shuffle_iterations(shuffle_iteration,:) = mean(shuffled_event_activity, 'omitnan');
    end

    % calculate the mean of all shuffled mean PETHs
    mean_shuffled_data = mean(mean_activity_all_shuffle_iterations); 

    % FIRST: calculate the MSE between activation between the average of the real 
    % PETH and the average of all shuffle iterations in user-specified time
    % window
    MSE_true = immse(mean_peri_event_activity_single_ROI(time_win), mean_shuffled_data(time_win)); 
    
    % Preallocate
    event_modulated_ROIs = zeros(nr_of_shuffles,1);

    % THEN: loop over nr of shuffle iterations and calculate the MSE between 
    % each shuffled PETH and the mean of all shuffle iterations. For each 
    % iteration check whether the MSE of the shuffled PETH to mean shuffled PETH is 
    % smaller than the MSE of the real PETH to the mean shuffled PETH 
    % (lower MSE value indicates greater similarity)
    for shuffle_nr = 1:nr_of_shuffles
            event_modulated_ROIs(shuffle_nr) ...
                = (immse(mean_activity_all_shuffle_iterations(shuffle_nr, time_win), ...
                mean_shuffled_data(time_win)) < MSE_true); 
    end
    
    % store all ROIs whose MSE is larger than 95% of the MSE's obtain from the
    % shuffle procedure
    if sum(event_modulated_ROIs)/numel(event_modulated_ROIs) > 0.95 
        idx_modulated_rois(roinr) = roinr;
    end
end
