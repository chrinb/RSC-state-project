function plot_across_day_swr_changes(varargin)

%{
Plot differences in SWR feature in pre-vs post-RUN sessions over days8
%}

data_to_plot     = varargin{1,1};
swr_type = varargin{1,2};
RF_idx   = varargin{1,3};
GOL_idx  = varargin{1,4};
%% Select datatype
if strcmp(swr_type, 'absolute')
    figure(5); clf
    sgtitle('Absolute difference in SWR nr pre- vs- post RUN')
    ylabel_text = 'Absolute nr of SWRs';

elseif strcmp(swr_type, 'relative')
    figure(6); clf
    sgtitle('SWR occurrence (Hz) pre- vs- post RUN')
    ylabel_text = 'SWR occurrence (Hz)';

elseif strcmp(swr_type, 'amplitude')
    figure(7); clf
    sgtitle('SWR z-scored amplitude pre- vs- post RUN')
    ylabel_text = 'SWR amplitude (z-scored)';

elseif strcmp(swr_type, 'duration')
    figure(8); clf
    sgtitle('SWR durations pre- vs- post RUN')
    ylabel_text = 'SWR duration';
elseif strcmp(swr_type, 'cluster')
    figure(9); clf
    sgtitle('SWR clusters pre- vs- post RUN')
    ylabel_text = 'SWR cluster n';
elseif strcmp(swr_type, 'qw')
    figure(9); clf
    sgtitle('QW bout length pre- vs- post RUN')
    ylabel_text = 'Seconds';
end

data_pre  = data_to_plot(1,:);
data_post = data_to_plot(2,:); 

data_pre = reshape(data_pre, 7,[]);
data_post = reshape(data_post, 7,[]);
%% Plot

% Plot absolute nr of SWRs over days
h(1) = subplot(1,1,1);
hold on

cmap = viridis(7);

for pre_post = 1:2

    if pre_post == 1
        data = data_pre;
    elseif pre_post == 2
        data = data_post;
    end

    for i = 1:7
        
        if iscell(data_to_plot(1,:))
            data_plot = vertcat( data{i,:});
            data_plot(data_plot == 0) = [];
            dim = size(data_plot, 1);
        else
           
            % dim = dim;
           data_plot = data(i,:);
           dim = size(data_plot, 2);
        end

        if pre_post == 1
            plot_x = ((1:7).*2)-1;
            x_vec = ones(1, dim ) *plot_x(i);
        elseif pre_post == 2
            plot_x = ((1:7).*2);
            x_vec = ones(1, dim ) *plot_x(i);
        end

        
        % Color by day
        % if ~(rem(i, 2) == 0)
            tmp_cmap =  cmap(i,:);
        % end

        b = boxchart(x_vec, data_plot, 'BoxFaceColor', tmp_cmap );

        b.MarkerColor = 'white';
        
        jitter_factor1 = (randn( size(x_vec))*0.06) + plot_x(i);
        
        s1 = scatter(jitter_factor1, data_plot, 'filled', 'MarkerFaceColor',[.5 .5 .5]);
        s1.MarkerFaceAlpha = .3;
        ylabel(ylabel_text)
        h(1).XTick = (1.5:2:14);
        h(1).XTickLabel = {'Day 1', 'Day 2', 'Day 3', 'Day 4', 'Day 5', 'Day 6', 'Day 7'};
        set(gca, 'FontSize', 12)    

    end
end

%% Compare means between days
merged_data =  [data_pre, data_post];

idx1 = [1 3 5 7 9 11];
idx2 = [2 4 6 8 10 12];

rf_idx = 1:4;
gol_idx = 5:7;


% Mean across days
if iscell(data_to_plot(1,:))
    merged_data_sorted = cell(size(merged_data));
    merged_data_sorted(:, idx1) = data_pre;
    merged_data_sorted(:, idx2) = data_post;

    merged_data_sorted = merged_data_sorted';

    for i = 1:7
    merged_data_sorted_2{i} =  vertcat(merged_data_sorted{i,:});
    end
else
    merged_data_sorted = zeros( size(merged_data));
    merged_data_sorted(:, idx1) = data_pre;
    merged_data_sorted(:, idx2) = data_post;


    data_mean = mean(merged_data_sorted,2);
    merged_data_sorted = merged_data_sorted';
end

%% RF vs GOL days 
rf_values = merged_data_sorted(:, rf_idx );
gol_values = merged_data_sorted(:, gol_idx);

if iscell(data_to_plot(1,:))
    rf_values = vertcat(rf_values{:});
    rf_values(rf_values == 0) = [];

    gol_values = vertcat(gol_values{:});
    gol_values(gol_values == 0) = [];
end


n_values = [numel(rf_values(:)), numel(gol_values(:)) ];
[~, maxidx] = max(n_values);

% Add NaNs
if maxidx == 1
    new_val = [gol_values(:);  NaN(n_values(1)-n_values(2), 1  )];
    old_val = rf_values(:);
    labels = {'', 'RF task', 'GOL task', ''};

    [ranksum_p, rank_sum_z]  = ranksum(old_val, new_val);
    [ttest_h, ttest_p] = ttest(old_val, new_val);
elseif maxidx == 2
    new_val = [rf_values(:);  NaN(n_values(2)-n_values(1), 1  )];
    old_val = gol_values(:);
    labels = {'', 'GOL task','RF task', ''};

    [ranksum_p, rank_sum_z] = ranksum(old_val, new_val);
    [ttest_h, ttest_p] = ttest(old_val, new_val);
end

figure(10),clf,
hold on 
h(2) = subplot(1,1,1);
x_vec1 = ones(length(old_val),1);
x_vec2 = ones(length(new_val),1)*2;

b1 = boxchart(x_vec1, old_val, 'BoxFaceColor','r', 'Notch','off');
b2 = boxchart(x_vec2, new_val, 'BoxFaceColor','b', 'Notch','off');

b1.MarkerColor = 'white';
b2.MarkerColor = 'white';

jitter_factor1 = (randn( size(x_vec1))*0.05) + 1;
jitter_factor2 = (randn( size(x_vec2))*0.05) + 2;

plot(1, mean(old_val, 'omitnan'), 'kd', 'LineWidth',2)
plot(2, mean(new_val, 'omitnan'), 'kd', 'LineWidth',2)

s1 = scatter(jitter_factor1, old_val, 'filled',  'MarkerFaceColor','r');
s2 = scatter(jitter_factor2,new_val, 'filled',  'MarkerFaceColor','b');

s1.MarkerFaceAlpha = .3;
s2.MarkerFaceAlpha = .2;

ylabel(ylabel_text)
set(gca, 'FontSize', 12)    

h(2).XTick = [0 1 2 3];
title(['ranksum p = ', num2str(ranksum_p), ', t-test p = ', num2str(ttest_p)]);
xticklabels(labels);

% [p, tbl, stats] = signrank(new_data)

%% Across day comparison
% pre = merged_data_sorted(1:2:end,:);
% post = merged_data_sorted(2:2:end,:);

[p, tbl, stats] = kruskalwallis(merged_data_sorted)
c1 = multcompare(stats)

tbl1 = array2table(c1,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])