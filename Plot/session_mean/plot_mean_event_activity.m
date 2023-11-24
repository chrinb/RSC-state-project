function sort_idx = plot_mean_event_activity(mean_session_data, mouseID)

% Written by Christoffer Berge | Vervaeke Lab

%{
Plot the mean peri-event activity. If more data from more than one data is
inputted the function puts the data together.
%}

if size(mean_session_data,1) > 1
    all_session_data_cat      = horzcat(mean_session_data{:});
    mean_peri_event_activity  = vertcat(all_session_data_cat{1,:});
    params                    = mean_session_data{1, 1}{2, 1};  
    time                      = mean_session_data{1, 1}{3, 1};  
    xlabel_text               = mean_session_data{1, 1}{7, 1};
else
    sessionID                        = mean_session_data{1, 1}{4,1};
    params                           = mean_session_data{1, 1}{2,1};
    time                             = mean_session_data{1, 1}{3,1};
    mean_peri_event_activity         = mean_session_data{1, 1}{1,1};
    event_idx                        = mean_session_data{1, 1}{5,1};
    all_data                         = mean_session_data{1, 1}{6,1};
    xlabel_text                      = mean_session_data{1, 1}{7, 1};
end

% Calculate standard error
signal_SE             = std(mean_peri_event_activity, 'omitnan') ./ sqrt(size(mean_peri_event_activity,1));

if strcmp(params.signal_type, 'deconv')
    c_lim        = [-.05 .5];
elseif strcmp(params.signal_type,  'Dff')
    c_lim        = [-.5 .5];
elseif strcmp(params.signal_type,  'transients')
    c_lim        = [-.05 .5];
end

% Sort ROIs according to mean z-scored activity in the -0.5 to + 0.5 interval
% surrounding SWR peak
frameshift    = round(.5/(1/31));
interval_mean = mean(mean_peri_event_activity(:, ( median(1:187)-frameshift:median(1:187)+frameshift)),2);
[max_val,~]   = max(interval_mean,[],2);
[~, sort_idx] = sortrows(max_val);

% Find unique mouse IDs
n_mice = numel(unique(mouseID));
%% Plot results
x1 = [time(1), time(end)];
y1 = [1 size(mean_peri_event_activity,1)];

if strcmp(params.zscore, 'yes')
    title_text = 'z-score ';
else
    title_text = '';
end

figure,
sgtitle(['Peri-event activity, ' num2str(size(mean_session_data,1)) ' sessions, ', num2str(n_mice), ' mice, ', params.cell_type], 'Interpreter', 'none') 

h(1) = subplot(3,3,[2,5]);
imagesc(x1, y1, mean_peri_event_activity(sort_idx,:)) 
ylabel('ROI #', FontSize=16)
set(gca, 'xtick',[])
cbar1          = colorbar;
caxis(c_lim)
cbar1.FontSize = 10;
cbar1.Position =[0.6500 0.3900 0.0100 0.4960];
title(['Mean ',title_text params.signal_type], FontSize=12)

h(2) = subplot(3,3,8);
shadedErrorBar(time, mean(mean_peri_event_activity, 'omitnan'),  signal_SE);
xline(0, '--', 'linew', 1) % Mark SWR peak time
xlabel(['Time from', xlabel_text, '(s)'],  FontSize=16)
set(gca, 'xlim',[time(1) time(end)])

linkaxes(h, 'x');