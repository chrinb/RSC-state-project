function sData = classify_cue_cells(sData)

%{
Classify cue cells. Based on shuffle method outlined in Lande et al 2023.
%}


bin_size    = sData.behavior.meta.binSize; 
binned_data = sData.imdata.binned.RoidFF;
n_bins      = sData.behavior.meta.nBins;
bin_vector  = 1:n_bins;




% Area +- cue position to search for activity (in cm)
cue_zone = 7;

% Cue positions (in cm)
x_soft1 = median( [sData.imdata.cues.C1A,  sData.imdata.cues.C1B]) ;
x_hard1 = median( [sData.imdata.cues.C2A,  sData.imdata.cues.C2B]) ;
x_hard2 = median( [sData.imdata.cues.C3A,  sData.imdata.cues.C3B]) ;
x_soft2 = median( [sData.imdata.cues.C4A,  sData.imdata.cues.C4B]) ;

% Area around each cue to search for activity
% cue_soft1_bins = round( [x_soft1-cue_zone, x_soft1+cue_zone] /bin_size);
% cue_soft2_bins = round( [x_soft2-cue_zone, x_soft2+cue_zone]/bin_size);
% cue_hard1_bins = round( [x_hard1-cue_zone, x_hard1+cue_zone]/bin_size);
% cue_hard2_bins = round ( [x_hard2-cue_zone, x_hard2+cue_zone]/bin_size);

cue_shift = 3;
cue_soft1_bins    = round( [x_soft1 + cue_shift - cue_zone, x_soft1 + cue_shift + cue_zone]/bin_size); 
cue_soft1_bins(2) = cue_soft1_bins(2) + 1; % Add + 1 to get even nr of bins for all cues
cue_hard1_bins = round( [x_hard1 + cue_shift - cue_zone, x_hard1 + cue_shift + cue_zone]/bin_size)+1; 
cue_hard1_bins(2) = cue_hard1_bins(2) + 1; % Add + 1 to get even nr of bins for all cues
cue_hard2_bins = round ([x_hard2 + cue_shift - cue_zone, x_hard2 + cue_shift + cue_zone]/bin_size);
cue_soft2_bins = round( [x_soft2 + cue_shift - cue_zone, x_soft2 + cue_shift + cue_zone]/bin_size);

all_cue_zones = [cue_soft1_bins; cue_hard1_bins; cue_hard2_bins; cue_soft2_bins ];

%% Shuffle analysis

n_rois = sData.imdata.nROIs;

cue_cells         = zeros(n_rois, 1);
cue_cell_activity = false(n_rois, 4);

sData.imdata.cue_activity.cue_cells         = [];
sData.imdata.cue_activity.cue_cell_sorting  = [];
sData.imdata.cue_activity.cue_cell_activity = []; %  

InOutRatio = 1.4;

%% Loop
cue_template = zeros(1, n_bins);
for i = 1:size(all_cue_zones, 1)
    % Set the elements between start and stop indices to 1 in the vector
    cue_template(all_cue_zones(i, 1):all_cue_zones(i, 2)) = 1;
end


figure(5),clf
colormap viridis
for roi_n = 1:n_rois 
        %% 
        
        fprintf('\n Calculating for ROI # %d', roi_n)

        if sData.imdata.roi_classification(roi_n) == 1

            % Activity map (trials x bins) and its mean for current ROI
            activity_map        = binned_data{1,roi_n};

            activity_map_mean   = mean( activity_map, 'omitnan');
            % activity_map_median = median( activity_map, 'omitnan');

            n_shuffles = 1000;
            mean_shuffled_mat = zeros(1000,80);
            % max_shuffle_xcorr = zeros( 1000, 1);
            for shuffle_n = 1:n_shuffles
                
                % Random integers between 1 and 80, one for every trial
                rand_int = randi( 80, 1, size(activity_map,1) );
        
                % Shuffle the activity map
                shuffled_mat = cell2mat( arrayfun(@(x) circshift(activity_map(x,:),[1 rand_int(x)]), (1:numel(rand_int))', 'UniformOutput', false));
        
                % Store mean shuffled activity map 
                mean_shuffled_mat(shuffle_n,:) = mean(shuffled_mat, 'omitnan');

                % max_shuffle_xcorr(shuffle_n) = max( xcorr( mean( shuffled_mat, 'omitnan'), cue_template ));
            end
        
            % Corr-correlation val
            % max_true_xcorr = max( xcorr( activity_map_mean, cue_template) );

            % Check whether position bins for original tuning curve is larger
            % than 95th percentile of the shuffled responses
            threshold_top    = prctile(mean_shuffled_mat, 99, 1);
            threshold_bottom = prctile(mean_shuffled_mat, 1, 1);

            % shuffle_mean = mean(mean_shuffled_mat);
            
            % test = activity_map_mean;
            % test(activity_map_mean < shuffle_mean) =NaN;
           
            % activity_map_mean_smooth = smoothdata(activity_map_mean, 'Gaussian', 1);

            sig_bin_idx     = activity_map_mean > threshold_top;
            % sig_bins = bin_vector(activity_map_mean_smooth > threshold_top);
    
            % Reliability check (peak activity per trial must happen inside
            % cue zones 25% of trials

            % [max_val, max_id] = max( activity_map,[], 2);
            % test = prctile(max_val, 50);
            % threshold_all_activity_map = prctile(activity_map(:), 95);
            % test = 2;

            % is_within_range = zeros(1,4);
            n_trials = size(activity_map,1);

            % reliability_check = zeros(n_trials, 1);

            peac_activity_above_threshold= zeros(n_trials, 1);
            peak_activity_in_cue = zeros(n_trials, 4);
            reliability_threshold = 0;
    
            % Concatenate all trials into a vector to determine baseline
            % level of activity for thresholding each trial

            % activity_map_cat    = reshape(activity_map, 1, []);
            % threshold_per_trial = 2*std(activity_map_cat, 'omitnan');
            % % Loop over trials
            % for trial_nr = 1:n_trials
            % 
            %     % Check if trial contains peak activity above threshold AND
            %     % that activity falls within cue zone
            % 
            %     trial_data       = smoothdata( activity_map(trial_nr,:), 'gaussian', 5);
            % 
            %     % Test different thresholds: compute threshold based on
            %     % activity in trial, or compute threshold over all trials
            %     % (either entire activity map or concatenated activity map)
            % 
            %     % trial_threshold  = 3*std(activity_map(trial_nr,:));
            %     trial_threshold  = threshold_all_activity_map;
            %     trial_threshold  = threshold_per_trial;
            % 
            %     [trial_max, trial_max_id]  = max(trial_data);
            % 
            %     peac_activity_above_threshold(trial_nr) = trial_max > trial_threshold ;
            % 
            %     if trial_max > trial_threshold
            %         if sum( trial_max_id >= all_cue_zones(:, 1) & trial_max_id <= all_cue_zones(:, 2) ) > 0
            %             peak_activity_in_cue(trial_nr, :) = trial_max_id >= all_cue_zones(:, 1) & trial_max_id <= all_cue_zones(:, 2);
            %         end
            %     end
            % 
            %     % if sum(is_within_range) > 0
            %     %     reliability_check(trial_nr) = 1;
            %     % end
            % end
            % % trials_to_include = peak_activity_in_cue(peac_activity_above_threshold==1,:);
            % 
            % if any( sum(peak_activity_in_cue) / n_trials > 0.1 )
            % 
            %     reliability_threshold = 1;
            % end

            % Loop over place fields and check whether there are at least 3
            % consecutive bins inside cue zone


            % Max place field width in cm
            pf_width_max = 35;
            pf_width_min = 1;

            % Place field start/stop bin
            [pf_start_bin, pf_stop_bin] = findTransitions(sig_bin_idx);
            
            if ~isempty(pf_start_bin)
                
                % Nr of potential place fields
                potential_pfs = size(pf_start_bin,2);

                % Place field size
                pf_width_bins = pf_stop_bin - pf_start_bin;
                pf_width_cm   = pf_width_bins*bin_size;
                
                crition_width = pf_width_cm > pf_width_min & pf_width_cm < pf_width_max;
                
                % Calc activity in-out ratio
                ratio_in_vs_out = [];
                for pf_n = 1:potential_pfs
                
                    pf_bins      =  pf_start_bin(pf_n):pf_stop_bin(pf_n);
                    pf_bin_idx   = ismember(bin_vector, pf_bins);
                    activity_in  = mean( activity_map_mean( pf_bin_idx ), 'omitnan');
                    % activity_out = mean( activity_map_mean( ~pf_log_idx ), 'omitnan');
                
                    activity_out = mean( activity_map_mean( ~sig_bin_idx ), 'omitnan');
                
                    ratio_in_vs_out(pf_n) = activity_in/activity_out;
                end
                
                criterion_activity = ratio_in_vs_out > InOutRatio;
                
                all_criteria = sum( [crition_width; criterion_activity], 1);
                all_criteria = all_criteria == 2;
    
                sig_peaks_in_cue_zone = false(pf_n, 4);
    
                % Loop over place/cue fields and check if criteria are met
                for pf_n = 1:potential_pfs
    
                        if all_criteria(pf_n) > 0
    
                            % If true, check which of the 4 cue zones the
                            % peak occurs in. If more than one, assign
                            % field to cue zone containing most significant
                            % position bins
                            pf_bins      =  pf_start_bin(pf_n):pf_stop_bin(pf_n);
                            % pf_bin_idx   = ismember(bin_vector, pf_bins);
        
                            sig_peaks_b_bins = zeros(1, 4);
                            for cue_nr = 1:4 
                                
                                tmp_cue_zone = all_cue_zones(cue_nr,:);
                               
                                sig_bin_in_cue_idx = pf_bins > tmp_cue_zone(1) & pf_bins < tmp_cue_zone(2);
        
                                sig_bins_in_cue     = pf_bins(sig_bin_in_cue_idx);
                                sig_bins_in_cue_idx = ismember(bin_vector, sig_bins_in_cue);
                                
                                n_successive_samples = 2;

                                if any( conv(double(sig_bins_in_cue_idx), ones(1,n_successive_samples), 'valid') == n_successive_samples)
                                    sig_peaks_in_cue_zone(pf_n, cue_nr) = true;
                                    sig_peaks_b_bins(cue_nr) = numel(sig_bins_in_cue);
                                end
                            end
                            
                            % Check if place fields have significant bins
                            % in more than one cue zone, if so, set
                            % significant response to be the cue with most
                            % bins
                            if sum( sig_peaks_in_cue_zone(pf_n, :)) > 1
                                [~, max_id ]                                               = max(sig_peaks_b_bins);
                                sig_peaks_in_cue_zone(pf_n, :)                             = false(1,4);
                                sig_peaks_in_cue_zone(pf_n, max_id) = true;
                            end
                        end
                end
  
            % Classify as cue cell if there are minimum two significant cue
            % responses AND reliability criterion met
           
            % if sum(sig_peaks_in_cue_zone) >= 2 && reliability_threshold == 1
            sig_peaks_in_cue_zone = any(sig_peaks_in_cue_zone, 1);  

            if sum( sig_peaks_in_cue_zone) >= 2

            % if max_true_xcorr > max(max_shuffle_xcorr)
                cue_cells(roi_n) = roi_n;                  
                cue_cell_activity(roi_n,:) = sig_peaks_in_cue_zone;
    
                fprintf('\n Significant cue response for ROI # %d', roi_n)
                
                figure(5), clf
                h(1) = subplot(211);
                imagesc(activity_map)
                colormap viridis
                axis square
                h(2) = subplot(212);
                hold on
                plot(activity_map_mean, 'k', 'LineWidth',1)
                ax = gca;   
                y_lims = ax.YLim;

                % Thresholds shaded
                x_vec   = 1:length(activity_map_mean);
                x_patch = [x_vec, fliplr(x_vec)];
                y_patch = [threshold_top, fliplr(threshold_bottom)]; 
                patch(x_patch, y_patch, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.2);  


                % Plot cues
                cue_y_lims = [-2, -2, 2, 2];

                cue1_x_patch = [all_cue_zones(1,:), fliplr(all_cue_zones(1,:))];
                cue2_x_patch = [all_cue_zones(2,:), fliplr(all_cue_zones(2,:))];
                cue3_x_patch = [all_cue_zones(3,:), fliplr(all_cue_zones(3,:))];
                cue4_x_patch = [all_cue_zones(4,:), fliplr(all_cue_zones(4,:))];


                patch( cue1_x_patch, cue_y_lims,'r', 'EdgeColor', 'none', 'FaceAlpha', 0.2);  
                patch( cue2_x_patch, cue_y_lims,'b', 'EdgeColor', 'none', 'FaceAlpha', 0.2);  
                patch( cue3_x_patch, cue_y_lims,'b', 'EdgeColor', 'none', 'FaceAlpha', 0.2);  
                patch( cue4_x_patch, cue_y_lims,'r', 'EdgeColor', 'none', 'FaceAlpha', 0.2);  

                set(gca,'YLim', y_lims)


                
                axis square
                title([num2str(roi_n), ' !!'])
                linkaxes(h, 'x')
                pause(1.5)

                % prompt = sprintf('Continue?' );
                % x = input(prompt,'s');
                % if strcmp(x, '')
            else
            
            % end

                figure(5), clf
                h(1) = subplot(211);
                imagesc(activity_map)
                colormap viridis
                axis square
                h(2) = subplot(212);
                hold on
                plot(activity_map_mean, 'k', 'LineWidth',1)
    
                ax = gca;
                y_lims = ax.YLim;
                % Plot upper/lower threshold
                x_vec = 1:length(activity_map_mean);
    
                x_patch = [x_vec, fliplr(x_vec)];
                y_patch = [threshold_top, fliplr(threshold_bottom)]; 
                patch(x_patch, y_patch, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.2);  
    
                % plot(threshold_top)
                % plot(threshold_bottom)
                cue_y_lims = [-2, -2, 2, 2];
    
                cue1_x_patch = [all_cue_zones(1,:), fliplr(all_cue_zones(1,:))];
                cue2_x_patch = [all_cue_zones(2,:), fliplr(all_cue_zones(2,:))];
                cue3_x_patch = [all_cue_zones(3,:), fliplr(all_cue_zones(3,:))];
                cue4_x_patch = [all_cue_zones(4,:), fliplr(all_cue_zones(4,:))];
    
    
                patch( cue1_x_patch, cue_y_lims,'r', 'EdgeColor', 'none', 'FaceAlpha', 0.2);  
                patch( cue2_x_patch, cue_y_lims,'b', 'EdgeColor', 'none', 'FaceAlpha', 0.2);  
                patch( cue3_x_patch, cue_y_lims,'b', 'EdgeColor', 'none', 'FaceAlpha', 0.2);  
                patch( cue4_x_patch, cue_y_lims,'r', 'EdgeColor', 'none', 'FaceAlpha', 0.2);  
    
                set(gca,'YLim', y_lims)
                axis square
                title([num2str(roi_n) ' not cue cell'])
                linkaxes(h, 'x')
             
                % prompt = sprintf('Continue?' );
                % x = input(prompt,'s');
                % pause(.3)
            end
            end

        end
        %%
end

%%

% Remove non-cue cells
all_cue_cells = cue_cells;
all_cue_cells(all_cue_cells==0) = [];

activity_map_all        = sData.imdata.binned.MeanRoiAct;
activity_map_all_zscore = zscore( activity_map_all, 0, 2);
activity_map_cue_cells  = activity_map_all_zscore(all_cue_cells,:);

%% Sort cue cells
pc_nr           = size(all_cue_cells,1);
MaxData         = max(activity_map_cue_cells, [], 2); % search the max in each row
MaxDataBin      = NaN(pc_nr,2);
MaxDataBin(:,2) = 1:pc_nr; % search Max position in bins
for i = 1:pc_nr

    if find( activity_map_cue_cells(i,:)== MaxData(i) ) % needed if the max number is twice in the dataset
        MaxDataBin(i,1) = find(activity_map_cue_cells(i,:)== MaxData(i));   % search Max position in bins
        continue
    end
end
SortingOrder  = sortrows(MaxDataBin, 1); % second column is the sorted ROI order. plot in these order
sorting       = SortingOrder(:,2);

% Store in sData
sData.imdata.cue_activity.cue_cells         = all_cue_cells; % Cue cell idx (out of all cells)
sData.imdata.cue_activity.cue_cell_sorting  = sorting; % Interal sorting of cue cells based on peak activity
sData.imdata.cue_activity.cue_cell_activity = cue_cell_activity; % Which of the cues contained significant response
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
imagesc(time_vec, [], smoothdata( activity_map_cue_cells(sorting,:), 'gaussian', 5))
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