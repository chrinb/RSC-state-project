function [binned_means, binned_sd]  = episode_third_analysis(varargin)
% [mean_vec, SD_vec, SEM_vec, plot_coordinates, e, means]
% Written by Christoffer Berge | Vervaeke lab

%{
Episode third analysis inspired by Grosmark et al. (2012). Function
splits a matrix into three equal-sized parts and computes mean and SD/SEM
for each part. 
%}

data = varargin{1,1};

% Flip if column vector
if size(data,2) == 1
    data = data';
end

if nargin <= 1
    bins = 3;
else
    bins = varargin{1,2};
end

if ~isempty(data)
    % Split the time-normalized episodes into thirds (or specified nr of
    % bins
    num_elements = size(data,2);
    
    samples_per_part = floor( num_elements / bins);
    
    % Split dataset into bins, compute mean and standard deviation per bin
    n_cells = size(data,1);
    binned_data = cell(1, bins);
    bin_boundaries = round( linspace(1, size(data, 2), bins+1) );

    % Split the matrix into sextiles
    for i = 1:bins
        rowIndices = bin_boundaries(i):bin_boundaries(i + 1) - 1;
        binned_data{i} = data(:, rowIndices);
    end
    
    binned_means = cellfun(@(x) mean(x, 'all'), binned_data, 'UniformOutput',true);
    binned_sd    = cellfun(@(x) std(x, 0, 'all'), binned_data, 'UniformOutput',true);

    % part_1 = data(:, 1:samples_per_part);
    % part_2 = data(:, samples_per_part:samples_per_part*2);
    % part_3 = data(:, samples_per_part*2:samples_per_part*3);
    % 
    % % Compute mean and SEM for each part
    % mean_pt1 = mean(part_1, 'all');
    % mean_pt2 = mean(part_2, 'all');
    % mean_pt3 = mean(part_3, 'all');
    % means = [mean_pt1, mean_pt2, mean_pt3];
%     means = struct();
%     means.mean_pt1 = mean_pt1;
%     means.mean_pt2 = mean_pt2;
%     means.mean_pt3 = mean_pt3;
% 
%     SEM_pt1 = std( mean(part_1,2))/ sqrt( numel( mean(part_1,2)));
%     SEM_pt2 = std( mean(part_2,2))/ sqrt(  numel( mean(part_3,2)));
%     SEM_pt3 = std( mean(part_3,2))/ sqrt(  numel( mean(part_3,2)));
% % 
%     SEM = struct();
%     SEM.pt1 = SEM_pt1;
%     SEM.pt2 = SEM_pt2;
% %     SEM.pt3 = SEM_pt3;
%      SEM = [SEM_pt1, SEM_pt2, SEM_pt3];
%     % Collect mean and SD in separate vectors
%     mean_vec = [mean_pt1, mean_pt2, mean_pt3];
%     SEM_vec  = [SEM_pt1, SEM_pt2, SEM_pt3];
%     SD_vec   = [std( mean(part_1,2)), std( mean(part_2,2)), std( mean(part_3,2))];
% 
%     % Calculate plotting variables
%     plot_factor      = round(samples_per_part/2);
%     plot_coordinates = round([plot_factor, (samples_per_part*2)-plot_factor,  (samples_per_part*3)-plot_factor]);
% 
    %% Compute test-statistic
    
    % If there are sufficient nr of observations, compute a test-statistic to
    % examine whether there are significant differences between the three parts
    % e = [];
    
%     if size(data,1) > 20
%      
%         mean_pt1_dim2 = mean( part_1, 2);
%         mean_pt2_dim2 = mean( part_2, 2);
%         mean_pt3_dim2 = mean( part_3, 2);
%     
%     %     % Kruskal-Wallis
%     %     [p, tbl, stats] = kruskalwallis( [mean_pt1_dim2, mean_pt2_dim2, mean_pt3_dim2]);
%     %     
%     %     c = multcompare(stats);
%     % 
%     %     % One-way ANOVA
%     %     [pA, tblA, statsA] = anova1( [mean_pt1_dim2, mean_pt2_dim2, mean_pt3_dim2]);
%     %     d = multcompare(statsA);
%     
%         % As scores within each third not necessarily are normally distributed, and 
%         % the scores in each third are not independent (i.e. they are dependent, repeated
%         % measures), use non-paramatric Friedmans ANOVA to assess group
%         % differences
%         [pB, tblB, statsB] = friedman( [mean_pt1_dim2, mean_pt2_dim2, mean_pt3_dim2])
%             e = multcompare(statsB)
%     end
else
    binned_means = [];
    binned_sd   = [];

end
