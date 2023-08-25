function results = plot_results_nansen_sleep(varargin)

% Written by Christoffer Berge || Vervaeke lab

% Plot event-averaged results per mouse and average across all mice. 

% Select particular dataset
prompt = sprintf('Specify substring for particular dataset? ');
str   = input(prompt, 's');
if ~isempty(str)
    varargin = [varargin, str]; 

end
% Load data and mouse IDs (user types in e.g. "m6120")
[data, miceIDs] = load_multi_session_data(varargin);

% Find nr of mice in data
n_mice = numel(fieldnames(data));

mouse_data    = cell(n_mice,3);
mouse_data_SE = cell(n_mice,3);
% Loop over nr of mice
for i = 1:n_mice
    
    % Get all data for current mouse
    temp_data = data.(miceIDs{i}).mean_swr_sleep;
    
    % Loop over nr of sessions for current mouse
    temp_data_awake   = [];
    temp_data_single  = [];
    temp_data_coupled = [];
    for ii = 1:size(temp_data,1)
        % Get awake SWRs
        temp_data2 = temp_data{ii, 2}; % get z-scored data 
        temp_data_awake = [temp_data_awake; temp_data2];
        
        % Get spindle-uncoupled SWRs
        temp_data2 = temp_data{ii, 4}; % get z-scored data 
        temp_data_single = [temp_data_single; temp_data2];
    
        % Get spindle-coupled SWRs
        temp_data2 = temp_data{ii, 6}; % get z-scored data 
        temp_data_coupled = [temp_data_coupled; temp_data2];
    end
    mouse_data{i,1}    = temp_data_awake;
    mouse_data{i,2}    = temp_data_single;
    mouse_data{i,3}    = temp_data_coupled;

    mouse_data_SE{i,1} = std(temp_data_awake, 'omitnan')./sqrt( size(temp_data_awake, 1));
    mouse_data_SE{i,2} = std(temp_data_single, 'omitnan')./sqrt( size(temp_data_single, 1));
    mouse_data_SE{i,3} = std(temp_data_coupled, 'omitnan')./sqrt( size(temp_data_coupled, 1));

end

% Get plot labels
if exist('varargin')
    [y_label, y_label2, x_label] = plot_options(varargin);
end

% Specify x axis length
event_center   = round(size(mouse_data{1, 1}  ,2)/2);
nr_of_seconds  = 3; % nr of seconds before/after spindle onset
nr_of_frames   = (nr_of_seconds*31*2)+1; % nr of recording frames in spindle window
frames_to_keep = event_center-floor(nr_of_frames/2):event_center+floor(nr_of_frames/2);
time           = linspace(-nr_of_seconds,nr_of_seconds, length(frames_to_keep));
center         = round(size(time,2)/2);
frameshift     = round(31/2);

%% Plot results
all_data_unclassified = [];
all_data_single       = [];
all_data_coupled      = [];

if sum(contains(varargin, 'bulk'))
    c_lim = [-1.5 1.5];
elseif sum(contains(varargin, 'deconv'))
    c_lim = [-0.05 .5];
else
    c_lim = [-.5 .5];
end


for i = 1:n_mice
    figure(i),
    sgtitle(miceIDs{i})
    signal_swr_unclassified = mouse_data{i,1};
    signal_swr_single       = mouse_data{i,2};
    signal_swr_coupled      = mouse_data{i,3};
    
    signal_SE_unclassified = mouse_data_SE{i,1};
    signal_SE_single       = mouse_data_SE{i,2};
    signal_SE_coupled      = mouse_data_SE{i,3};
    % If analyzing bulk data, sort each SWR-type according to max activity in a
    % certain time window around SWR
    if sum(contains(varargin, 'bulk'))
        Interval_mean = mean( signal_swr_unclassified  (:,center:(center+frameshift)),2);
        [max_val,~]   = max(Interval_mean, [], 2);
        [~, sort_idx_unclassified] = sortrows(max_val);
        
        Interval_mean = mean( signal_swr_single  (:,center:(center+frameshift)),2);
        [max_val,~]   = max(Interval_mean, [], 2);
        [~, sort_idx_single] = sortrows(max_val);
        
        Interval_mean = mean( signal_swr_coupled  (:,center:(center+frameshift)),2);
        [max_val,~]   = max(Interval_mean, [], 2);
        [~, sort_idx_coupled] = sortrows(max_val);
    % If analyzing population data, sort cells after their activity in the
    % SWR-spindle events only
    else
         Interval_mean = mean( signal_swr_coupled  (:,center:(center+frameshift)),2);
        [max_val,~]   = max(Interval_mean, [], 2);
        [~, sort_idx_unclassified] = sortrows(max_val);
        [~, sort_idx_single] = sortrows(max_val);
        [~, sort_idx_coupled] = sortrows(max_val);
    end


    % Set imagesc X and Y 
    x1             = [time(1) time(end)];
    y_unclassified = [1 size(signal_swr_unclassified,1) ];
    y_single       = [1 size(signal_swr_single,1) ];
    y_coupled      = [1 size(signal_swr_coupled,1) ];

    subplot(3, 3, [1 4])
    imagesc(x1,y_unclassified, signal_swr_unclassified(sort_idx_unclassified,:)) %  colorplot of average individual roi activity during ripples
    ylabel(y_label2)
    xlabel(x_label)
    c(1) = colorbar;
    caxis(c_lim)
    c(1).Position(1) = 0.35;
    c(1).Position(2) = 0.41;
    c(1).Position(3) = 0.006;

    subplot(3,3, 7)
    shadedErrorBar(time, mean(signal_swr_unclassified, 'omitnan'), signal_SE_unclassified,'lineprops', 'r');
    xline(0, '--', 'linew', 1) % Mark SWR peak time
    xlabel(x_label)
    ylabel(y_label)
    legend('Unclassified')
    
    % Spindle-uncoupled SWRs
    subplot(3,3, [2 5])
    imagesc(x1, y_single, signal_swr_single(sort_idx_single,:)) %  colorplot of average individual roi activity during ripples
    ylabel(y_label2)
    xlabel(x_label)
    c(2) = colorbar;
    caxis(c_lim)
    c(2).Position(1) = 0.63;
    c(2).Position(2) = 0.41;
    c(2).Position(3) = 0.006;

    subplot(3,3,8)
    shadedErrorBar(time, mean(signal_swr_single,'omitnan'), signal_SE_single,'lineprops', 'b');
    xline(0, '--', 'linew', 1)

    ylabel(y_label)
    xlabel(x_label)
    legend('Spindle-uncoupled')

    subplot(3,3,[3 6])
    imagesc(x1,y_coupled, signal_swr_coupled(sort_idx_coupled,:)) %  colorplot of average individual roi activity during ripples
    xlabel(x_label)
    ylabel(y_label2)
    caxis(c_lim)
    c(3) = colorbar;
    c(3).Position(1) = 0.91;
    c(3).Position(2) = 0.41;
    c(3).Position(3) = 0.006;

    subplot(3,3,9)
    shadedErrorBar(time, mean(signal_swr_coupled, 'omitnan'), signal_SE_coupled,'lineprops', 'k');
    xline(0, '--', 'linew', 1)
    xlabel(x_label)
    ylabel(y_label)
    legend('Spindle-coupled')
    
    
    % Get data from all mice
    all_data_unclassified = [all_data_unclassified; signal_swr_unclassified];
    all_data_single       = [all_data_single; signal_swr_single];
    all_data_coupled      = [all_data_coupled; signal_swr_coupled];
end


%% Plot average of all mice
figure(n_mice+1),

% If analyzing bulk data, sort each SWR-type according to max activity in a
% certain time window around SWR
if sum(contains(varargin, 'bulk'))
    Interval_mean = mean( all_data_unclassified  (:,center:(center+frameshift)),2);
    [max_val,~]   = max(Interval_mean, [], 2);
    [~, sort_idx_unclassified] = sortrows(max_val);
    
    Interval_mean = mean( all_data_single  (:,center:(center+frameshift)),2);
    [max_val,~]   = max(Interval_mean, [], 2);
    [~, sort_idx_single] = sortrows(max_val);
    
    Interval_mean = mean( all_data_coupled  (:,center:(center+frameshift)),2);
    [max_val,~]   = max(Interval_mean, [], 2);
    [~, sort_idx_coupled] = sortrows(max_val);

    % If analyzing population data, sort cells after their activity in the
    % SWR-spindle events only
else
    Interval_mean = mean( all_data_coupled  (:,center:(center+frameshift)),2);
    [max_val,~]   = max(Interval_mean, [], 2);
    [~, sort_idx_unclassified] = sortrows(max_val);
    [~, sort_idx_single] = sortrows(max_val);
    [~, sort_idx_coupled] = sortrows(max_val);

end

all_data_unclassified_SE = std(all_data_unclassified, 'omitnan')./ sqrt(size(all_data_unclassified,1));
all_data_single_SE       = std(all_data_single, 'omitnan')./ sqrt(size(all_data_single,1));
all_data_coupled_SE      = std(all_data_coupled, 'omitnan')./ sqrt(size(all_data_coupled,1));

y_unclassified = [1 size(all_data_unclassified,1) ];
y_single       = [1 size(all_data_single,1) ];
y_coupled      = [1 size(all_data_coupled,1) ];

subplot(3, 3, [1 4])
imagesc(x1,y_unclassified, all_data_unclassified(sort_idx_unclassified,:)) %  colorplot of average individual roi activity during ripples
ylabel(y_label2)
xlabel(x_label)
caxis(c_lim)
c(1) = colorbar;
c(1).Position(1) = 0.35;
c(1).Position(2) = 0.41;
c(1).Position(3) = 0.006;
title('Awake SWR')

subplot(3,3, 7)
shadedErrorBar(time, mean(all_data_unclassified, 'omitnan'), all_data_unclassified_SE,'lineprops', 'r');
xline(0, '--', 'linew', 1) % Mark SWR peak time
xlabel(x_label)
ylabel(y_label)
legend('Unclassified')

% Spindle-uncoupled SWRs
subplot(3,3, [2 5])
imagesc(x1, y_single, all_data_single(sort_idx_single,:)) %  colorplot of average individual roi activity during ripples
ylabel(y_label2)
xlabel(x_label)
caxis(c_lim)
c(2) = colorbar;
c(2).Position(1) = 0.63;
c(2).Position(2) = 0.41;
c(2).Position(3) = 0.006;
title('Spindle-uncoupled SWR')

subplot(3,3,8)
shadedErrorBar(time, mean(all_data_single,'omitnan'), all_data_single_SE,'lineprops', 'b');
xline(0, '--', 'linew', 1)
ylabel(y_label)
xlabel(x_label)
legend('Spindle-uncoupled')

subplot(3,3,[3 6])
imagesc(x1,y_coupled, all_data_coupled(sort_idx_coupled,:)) %  colorplot of average individual roi activity during ripples
xlabel(x_label)
ylabel(y_label2)
caxis(c_lim)
c(3) = colorbar;
c(3).Position(1) = 0.91;
c(3).Position(2) = 0.41;
c(3).Position(3) = 0.006;
title('Spindle-coupled SWR')

subplot(3,3,9)
shadedErrorBar(time, mean(all_data_coupled, 'omitnan'), all_data_coupled_SE,'lineprops', 'k');
xline(0, '--', 'linew', 1)
xlabel(x_label)
ylabel(y_label)
legend('Spindle-coupled')

%% Plot average of all mice
figure(n_mice+2),

shadedErrorBar(time, mean(all_data_unclassified, 'omitnan'), all_data_unclassified_SE,'lineprops', 'r'); hold on
shadedErrorBar(time, mean(all_data_single, 'omitnan'), all_data_single_SE,'lineprops', 'b');
shadedErrorBar(time,  mean(all_data_coupled, 'omitnan'), all_data_coupled_SE,'lineprops', 'k');
xline(0, '--', 'linew', 1)
xlabel(x_label)
ylabel(y_label)
legend('Awake SWR', 'Spindle-uncoupled SWR', 'Spindle-coupled SWR')
%% Store data
results = struct();

results.miceIDs                 = miceIDs';
results.data                    = mouse_data;
results.data_SE                 = mouse_data_SE;
results.all_data_unclassified   = all_data_unclassified;
results.all_data_unclassifiedSE = all_data_unclassified_SE;
results.all_data_single         = all_data_single;
results.all_data_singleSE       = all_data_single_SE;
results.all_data_coupled        = all_data_coupled;
results.all_data_coupledSE      = all_data_coupled_SE;

