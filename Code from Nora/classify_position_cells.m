function sData = classify_position_cells(sData)

%{
Classify position-tuned cells. Based on shuffle method outlined in Lande et al 2023.
%}


bin_size    = sData.behavior.meta.binSize; 
binned_data = sData.imdata.binned.RoidFF;
n_bins      = sData.behavior.meta.nBins;
bin_vector  = 1:n_bins;

%% Shuffle analysis

n_rois = sData.imdata.nROIs;

position_cells_mean         = zeros(n_rois, 1);
position_cells_median       = zeros(n_rois, 1);

for roi_n = 1:n_rois 
        
        % Activity map (trials x bins) and its mean for current ROI
        activity_map = binned_data{1,roi_n};
        activity_map_mean   = mean( activity_map, 'omitnan');
        activity_map_median = median( activity_map, 'omitnan');

        n_shuffles = 1000;
        mean_shuffled_mat = zeros(1000,80);
        for shuffle_n = 1:n_shuffles
            
            % Random integers between 1 and 80, one for every trial
            rand_int = randi( 80, 1, size(activity_map,1) );
    
            % Shuffle the activity map
            shuffled_mat = cell2mat( arrayfun(@(x) circshift(activity_map(x,:),[1 rand_int(x)]), (1:numel(rand_int))', 'UniformOutput', false));
    
            % Store mean shuffled activity map 
            mean_shuffled_mat(shuffle_n,:) = mean(shuffled_mat, 'omitnan');
        end
    
        % Check whether position bins for original tuning curve is larger
        % than 95th percentile of the shuffled responses
        threshold1 = prctile(mean_shuffled_mat, 99);
        threshold2 = prctile(mean_shuffled_mat, 99.9);

        % sig_bins1            = bin_vector(activity_map_mean > threshold1);

        logical_vector1 = activity_map_mean > threshold1;
        logical_vector2 = activity_map_median > threshold2;

        % If there are neighboring significant bins, it's a
        % significant cue rsponse to that cue
        % if any(sig_bins(1:end-1) & sig_bins(2:end))
        % if any(test(1:end-1) & test(2:end))
        if any( conv(double(logical_vector1), ones(1, 3), 'valid') == 3) 
            position_cells_mean(roi_n) = roi_n;
            
            figure(1), clf
            hold on
            plot(threshold1)
            plot(activity_map_mean)
            fprintf('\n Significant position response for ROI # %d', roi_n)

        end

         if any( conv(double(logical_vector2), ones(1, 3), 'valid') == 3) 
            position_cells_median(roi_n) = roi_n;
            
            % figure(1), clf
            % hold on
            % plot(threshold2)
            % plot(activity_map_median)
            % fprintf('\n Significant position response for ROI # %d', roi_n)

        end
        fprintf('\n Calculating for ROI # %d', roi_n)
end

%%

% Remove non-cue cells
all_pos_cells = position_cells_mean;
all_pos_cells(all_pos_cells==0) = [];

all_pos_cells_med = position_cells_median;
all_pos_cells_med(all_pos_cells_med==0) = [];

activity_map_all        = sData.imdata.binned.MeanRoiAct;
activity_map_all_zscore = zscore( activity_map_all, 0, 2);
activity_map_pos_cells  = activity_map_all_zscore(all_pos_cells,:);

%% Sort cue cells
pc_nr           = size(all_pos_cells,1);
MaxData         = max(activity_map_pos_cells, [], 2); % search the max in each row
MaxDataBin      = NaN(pc_nr,2);
MaxDataBin(:,2) = 1:pc_nr; % search Max position in bins
for i = 1:pc_nr

    if find( activity_map_pos_cells(i,:)== MaxData(i) ) % needed if the max number is twice in the dataset
        MaxDataBin(i,1) = find(activity_map_pos_cells(i,:)== MaxData(i));   % search Max position in bins
        continue
    end
end
SortingOrder  = sortrows(MaxDataBin, 1); % second column is the sorted ROI order. plot in these order
sorting       = SortingOrder(:,2);

% Store in sData
sData.imdata.pos_activity.pos_cells_mean         = all_pos_cells; % Cue cell idx (out of all cells)
sData.imdata.pos_activity.pos_cell_mean_sorting  = sorting; % Interal sorting of cue cells based on peak activity
sData.imdata.pos_activity.pos_cells_median       = all_pos_cells_med; % Cue cell idx (out of all cells)
%%
distance = bin_size*n_bins;

time_vec = linspace(0, distance, n_bins);

% Set the y-axis location where you want to place the symbol (outside of the plot)
y_location = -10;  % You can adjust this value as needed to position the symbol

% Set the symbol to be placed outside the plot
symbol = '∇';  % You can change this to any other symbol

% Soft cues
cue_text_size = 12;

x_soft1 = median( [sData.imdata.cues.C1A,  sData.imdata.cues.C1B]) ;
x_soft2 = median( [sData.imdata.cues.C4A,  sData.imdata.cues.C4B]) ;
x_hard1 = median( [sData.imdata.cues.C2A,  sData.imdata.cues.C2B]) ;
x_hard2 = median( [sData.imdata.cues.C3A,  sData.imdata.cues.C3B]) ;


figure, 
imagesc(time_vec, [], smoothdata( activity_map_pos_cells(sorting,:), 'gaussian', 5))
colormap viridis

text(x_soft1, y_location, symbol, 'FontSize', cue_text_size, 'FontWeight', 'bold', 'Color', 'red', ...
 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Interpreter', 'none');    
text(x_soft2, y_location, symbol, 'FontSize', cue_text_size, 'FontWeight', 'bold', 'Color', 'red', ...
 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Interpreter', 'none');

% Hard cues
text(x_hard1, y_location, symbol, 'FontSize', cue_text_size, 'FontWeight', 'bold', 'Color', 'black', ...
 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Interpreter', 'none');
text(x_hard2, y_location, symbol, 'FontSize', cue_text_size, 'FontWeight', 'bold', 'Color', 'black', ...
 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Interpreter', 'none');