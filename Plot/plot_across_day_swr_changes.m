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
end


% if iscell(data_to_plot(1,:))
    % data_pre  =  data_to_plot(1,:) ;
    % data_post =  data_to_plot(2,:) ;
    % data_pre = reshape(data_pre, 7,[]);
    % data_post = reshape(data_post, 7,[]);
% else
    data_pre  = data_to_plot(1,:);
    data_post = data_to_plot(2,:); 

    data_pre = reshape(data_pre, 7,[]);
    data_post = reshape(data_post, 7,[]);

    % dim = size(data_pre, 2);
% end

% n_abs_swr_combined = data1;
% n_rel_swr_combined = data2;
% 
% dim = 2;
% % Format the data to plot SWRs per day
% pre_data = n_abs_swr_combined(1,:);
% post_data = n_abs_swr_combined(2, :);
% 
% pr = reshape(pre_data, 7,[])';
% po = reshape(post_data, 7,[])';

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

% Compare means between days
merged_data =  [data_pre, data_post];
merged_data_sorted = zeros( size(merged_data));

idx1 = [1 3 5 7 9 11];
idx2 = [2 4 6 8 10 12];

merged_data_sorted(:, idx1) = data_pre;
merged_data_sorted(:, idx2) = data_post;

% Mean across days
if iscell(data_to_plot(1,:))
    for i = 1:7
    merged_data_sorted_2{i} =  vertcat(merged_data_sorted{i,:});
    end

    test = cellfun(@vertcat, data_pre)
else
    data_mean = mean(merged_data_sorted,2);
    merged_data_sorted = merged_data_sorted';
end


[p, tbl, stats] = friedman(merged_data_sorted)
c1 = multcompare(stats)

tbl1 = array2table(c1,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])