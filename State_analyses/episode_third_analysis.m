function [means, SEM]  = episode_third_analysis(data)
% [mean_vec, SD_vec, SEM_vec, plot_coordinates, e, means]
% Written by Christoffer Berge | Vervaeke lab

%{
Episode third analysis inspired by Grosmark et al. (2012). Function
splits a matrix into three equal-sized parts and computes mean and SD/SEM
for each part. 
%}

if ~isempty(data)
    % Split the time-normalized REM episodes into thirds
    num_elements = size(data,2);
    
    elements_per_part = floor( num_elements / 3);
    
    % Parts for the interpolated data
    part_1 = data(:, 1:elements_per_part);
    part_2 = data(:, elements_per_part:elements_per_part*2);
    part_3 = data(:, elements_per_part*2:elements_per_part*3);
    
    % Compute mean and SEM for each part
    mean_pt1 = mean(part_1, 'all');
    mean_pt2 = mean(part_2, 'all');
    mean_pt3 = mean(part_3, 'all');
    means = [mean_pt1, mean_pt2, mean_pt3];
%     means = struct();
%     means.mean_pt1 = mean_pt1;
%     means.mean_pt2 = mean_pt2;
%     means.mean_pt3 = mean_pt3;
    
    SEM_pt1 = std( mean(part_1,2))/ sqrt( numel( mean(part_1,2)));
    SEM_pt2 = std( mean(part_2,2))/ sqrt(  numel( mean(part_3,2)));
    SEM_pt3 = std( mean(part_3,2))/ sqrt(  numel( mean(part_3,2)));
% 
%     SEM = struct();
%     SEM.pt1 = SEM_pt1;
%     SEM.pt2 = SEM_pt2;
%     SEM.pt3 = SEM_pt3;
     SEM = [SEM_pt1, SEM_pt2, SEM_pt3];
    % Collect mean and SD in separate vectors
    mean_vec = [mean_pt1, mean_pt2, mean_pt3];
    SEM_vec  = [SEM_pt1, SEM_pt2, SEM_pt3];
    SD_vec   = [std( mean(part_1,2)), std( mean(part_2,2)), std( mean(part_3,2))];
    
    % Calculate plotting variables
    plot_factor      = round(elements_per_part/2);
    plot_coordinates = round([plot_factor, (elements_per_part*2)-plot_factor,  (elements_per_part*3)-plot_factor]);
    
    %% Compute test-statistic
    
    % If there are sufficient nr of observations, compute a test-statistic to
    % examine whether there are significant differences between the three parts
    e = [];
    
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
    means = [];
    SEM   = [];

end
