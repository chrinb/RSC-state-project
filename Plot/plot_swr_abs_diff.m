function plot_swr_abs_diff(data, RF_idx, GOL_idx, session_nr, swr_type)

%{ 
Plot the difference in absolute or relative nr of SWRs in pre- vs post-RUN sessions for
(1) all sessions, (2) random foraging task sessions only, and (3)
goal-oriented task sessions only
%}

%% Select datatype
if strcmp(swr_type, 'absolute')
    figure(1); clf
    sgtitle('Absolute difference in SWR nr pre- vs- post RUN')
    ylabel_text = 'Absolute nr of SWRs';

elseif strcmp(swr_type, 'relative')
    figure(2); clf
    sgtitle('SWR occurrence (Hz) pre- vs- post RUN')
    ylabel_text = 'SWR occurrence (Hz)';

elseif strcmp(swr_type, 'amplitude')
    figure(3); clf
    sgtitle('SWR z-scored amplitude pre- vs- post RUN')
    ylabel_text = 'SWR amplitude (z-scored)';

elseif strcmp(swr_type, 'duration')
    figure(4); clf
    sgtitle('SWR durations pre- vs- post RUN')
    ylabel_text = 'SWR duration (s)';
end

%% Plot
h(1) = subplot(1,3,1);
hold on   

if iscell(data(1,:))
    data_pre  = vertcat( data{1,:});
    data_post = vertcat( data{2,:});
    dim = 1;
else
    data_pre  = data(1,:);
    data_post = data(2,:); 
    dim = 2;
end

data_pre(data_pre == 0) = [];
data_post(data_post == 0) = [];

x_vec1 = ones(1, size(data_pre, dim));
x_vec2 = ones(1, size(data_post, dim))+1;

b1 = boxchart(x_vec1 , data_pre, 'BoxFaceColor','r', 'Notch','off');
b2 = boxchart(x_vec2 , data_post, 'BoxFaceColor','b', 'Notch','off');

b1.MarkerColor = 'white';
b2.MarkerColor = 'white';

jitter_factor1 = (randn( size(x_vec1))*0.05) + 1;
jitter_factor2 = (randn( size(x_vec2))*0.05) + 2;

s1 = scatter(jitter_factor1, data_pre, 'filled', 'MarkerFaceColor','r');
s2 = scatter(jitter_factor2, data_post, 'filled', 'MarkerFaceColor','b');
s1.MarkerFaceAlpha = .2;
s2.MarkerFaceAlpha = .2;

h(1).XTick = [0 1 2 3];    
xticklabels({'', 'Pre-RUN','Post-RUN', ''});
set(gca, 'FontSize', 12)    
ylabel(ylabel_text,'FontSize',12)

if strcmp(swr_type, 'amplitude') || strcmp(swr_type, 'duration')
    [~,p]  = ttest2(data_pre, data_post );

else
    p  = signrank(data_pre, data_post );
end
title(['All sessions, p = ', num2str(p)])

%% Plot the PRE- to POST-RUN difference in SWR nr for RF task only
% plot_RF_idx = reshape(RF_idx,[2, session_nr/2]);
plot_RF_idx = reshape(RF_idx,[2, session_nr/2]);

if iscell(data(1,:))
    data_pre_RF = vertcat( data{1, plot_RF_idx(1,:) } );
    data_post_RF = vertcat( data{2, plot_RF_idx(2,:) } );
    data_pre_RF(data_pre_RF==0) = [];
    data_post_RF(data_post_RF==0) = [];

else
    data_pre_RF = data_pre(plot_RF_idx(1,:));
    data_post_RF = data_post(plot_RF_idx(2,:));
end

h(2) = subplot(1,3,2) ;
hold on   
x_vec1 = ones(1, numel(data_pre_RF) );
x_vec2 = ones(1, numel(data_post_RF) ) + 1;

b1 = boxchart(x_vec1 , data_pre_RF, 'BoxFaceColor','r', 'Notch','off');
b2 = boxchart(x_vec2 , data_post_RF, 'BoxFaceColor','b', 'Notch','off');

b1.MarkerColor = 'white';
b2.MarkerColor = 'white';

jitter_factor1 = (randn( size(x_vec1))*0.05) + 1;
jitter_factor2 = (randn( size(x_vec2))*0.05) + 2;

s1 = scatter(jitter_factor1, data_pre_RF, 'filled', 'MarkerFaceColor','r');
s2 = scatter(jitter_factor2, data_post_RF, 'filled', 'MarkerFaceColor','b');
s1.MarkerFaceAlpha = .2;
s2.MarkerFaceAlpha = .2;

h(2).XTick = [0 1 2 3];
xticklabels({'', 'Pre-RUN','Post-RUN', ''});
set(gca, 'FontSize', 12)    
% ylabel('Absolute nr of SWRs')
if strcmp(swr_type, 'amplitude') || strcmp(swr_type, 'duration')
    [~,p]  = ttest2(data_pre_RF, data_post_RF );
else
    p  = signrank(data_pre_RF , data_post_RF);
end

title(['RF, p = ', num2str(p)])

%% Test 
% plot_RF_idx = reshape(RF_idx,[2, session_nr/2]);
plot_GOL_idx = reshape(GOL_idx,[2, session_nr/2]);

if iscell(data(1,:))
    data_pre_GOL = vertcat( data{1, plot_GOL_idx(1,:) } );
    data_post_GOL = vertcat( data{2, plot_GOL_idx(2,:) } );
    data_pre_GOL(data_pre_GOL==0) = [];
    data_post_GOL(data_post_GOL==0) = [];

else
    data_pre_GOL = data_pre(plot_GOL_idx(1,:));
    data_post_GOL = data_post(plot_GOL_idx(2,:));
end

h(3) = subplot(1,3,3) ;
hold on   

x_vec1 = ones(1, numel(data_pre_GOL) );
x_vec2 = ones(1, numel(data_post_GOL) ) + 1;

b1 = boxchart(x_vec1 , data_pre_GOL, 'BoxFaceColor','r', 'Notch','off');
b2 = boxchart(x_vec2 , data_post_GOL, 'BoxFaceColor','b', 'Notch','off');

b1.MarkerColor = 'white';
b2.MarkerColor = 'white';

jitter_factor1 = (randn( size(x_vec1))*0.05) + 1;
jitter_factor2 = (randn( size(x_vec2))*0.05) + 2;

s1 = scatter(jitter_factor1, data_pre_GOL, 'filled', 'MarkerFaceColor','r');
s2 = scatter(jitter_factor2, data_post_GOL, 'filled', 'MarkerFaceColor','b');
s1.MarkerFaceAlpha = .2;
s2.MarkerFaceAlpha = .2;

h(3).XTick = [0 1 2 3];
xticklabels({'', 'Pre-RUN','Post-RUN', ''});
set(gca, 'FontSize', 12)    
% ylabel('Absolute nr of SWRs')
if strcmp(swr_type, 'amplitude') || strcmp(swr_type, 'duration')
    [~,p]  = ttest2(data_pre_GOL, data_post_GOL );
else
    p  = signrank(data_pre_GOL , data_post_GOL);
end

title(['GOL, p = ', num2str(p)])

%% Plot the PRE- to POST-RUN difference in SWR nr for GOL task only
% plot_GOL_idx = reshape(GOL_idx,[2, session_nr/2]);
% nr_swr_GOL = data_pre;
% nr_swr_GOL(~plot_GOL_idx) = NaN;
% 
% h(3) = subplot(1,3,3) ;
% hold on   
% x_vec1 = ones(1, sum(~isnan(nr_swr_GOL(1,:))) );
% x_vec2 = ones(1, sum(~isnan(nr_swr_GOL(2,:))) )+1;
% 
% b1 = boxchart(x_vec1 , data_pre(1, plot_GOL_idx(1,:)), 'BoxFaceColor','r', 'Notch','off');
% b2 = boxchart(x_vec2 , data_post(2, plot_GOL_idx(2,:)), 'BoxFaceColor','b', 'Notch','off');
% 
% b1.MarkerColor = 'white';
% b2.MarkerColor = 'white';
% 
% jitter_factor1 = (randn( size(x_vec1))*0.02) + 1;
% jitter_factor2 = (randn( size(x_vec2))*0.02) + 2;
% 
% scatter(jitter_factor1, data_pre(1, plot_GOL_idx(1,:)), 'filled', 'MarkerFaceColor','r')
% scatter(jitter_factor2, data_post(2, plot_GOL_idx(2,:)), 'filled', 'MarkerFaceColor','b')
% 
% h(2).XTick = [0 1 2 3];
% xticklabels({'', 'Pre-RUN','Post-RUN', ''});
% set(gca, 'FontSize', 12)    
% p       = ranksum(nr_swr_GOL(1,:) ,nr_swr_GOL(2,:));
% title(['GOL task, p = ', num2str(p)])