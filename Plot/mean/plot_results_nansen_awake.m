function results = plot_results_nansen_awake(varargin)

% Written by Christoffer Berge || Vervaeke lab

% Plot results obtained from using Nansen code. User inputs strings
% specifying 
% (1) signal type ('dff' or 'deconv')
% (2) event type ('swr', 'swa' or 'spindle')
% (3) experiment type ('bulk' or 'population') 
% (4)
% (5) extract the event average or significantly modulated cells ('avg' or
% 'mod')

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

mouse_data    = cell(n_mice,1);
mouse_data_SE = cell(n_mice,1);
% Loop over nr of mice
for i = 1:n_mice
    
    % Get all data for current mouse
    switch varargin{1, 2}  
        case 'swr'
            temp_data = data.(miceIDs{i}).mean_swr_awake;
            nr_of_seconds  = 3; 
        case 'spindle'
            temp_data = data.(miceIDs{i}).mean_spindle;
            nr_of_seconds  = 3; 
        case 'swa'
            temp_data = data.(miceIDs{i}).mean_swa;
            nr_of_seconds  = 1; 
    end
    
    % Loop over nr of sessions for current mouse
    temp_data1 = [];
    for ii = 1:size(temp_data,1)
%         if strcmp(varargin{1,2}, 'swa')
        temp_data2 = temp_data{ii, 2}; % get z-scored data 
        temp_data1 = [temp_data1; temp_data2];
    end
    mouse_data{i}    = temp_data1;
    mouse_data_SE{i} = std(temp_data1, 'omitnan')./sqrt( size(temp_data1, 1));
end

% Get plot labels
% if exist('varargin')
%     [y_label, y_label2, x_label] = plot_options(varargin);
% end
[x_label, y_label2, y_label] = get_xy_labels(varargin);

% Specify x axis length
event_center   = round(size(mouse_data{1, 1}  ,2)/2);
% nr_of_seconds  = 3; % nr of seconds before/after spindle onset
nr_of_frames   = (nr_of_seconds*31*2)+1; % nr of recording frames in spindle window
frames_to_keep = event_center-floor(nr_of_frames/2):event_center+floor(nr_of_frames/2);
time           = linspace(-nr_of_seconds,nr_of_seconds, length(frames_to_keep));
center         = round(size(time,2)/2);
frameshift     = round(31/2);

%% Plot results
all_data = [];
switch varargin{1, 1}  
    case 'dff'
    c_lim    = [-.3 .3];
    case 'deconv'
    c_lim    = [-.05 .5];
end

figure,
for i = 1:n_mice
%     figure(i),
    data2plot = mouse_data{i,1};
    dataSE    = mouse_data_SE{i,1};
    Interval_mean = mean( data2plot  (:,center:(center+frameshift)),2);
    [max_val,~]   = max(Interval_mean, [], 2);
    [~, sort_idx] = sortrows(max_val);


    % Set up X and Y limits of imagesc function
    x1 = [time(1) time(end)];
    y1 = [1 size(data2plot, 1)]; 

    h(1) = subplot(3,5,[0+i,5+i]);
    title(miceIDs{i}),
    imagesc(x1, y1, data2plot(sort_idx,:))
    font = gca;
    font.FontSize = 16;
    set(gca, 'xtick',[])
    if i > 1
        set(gca, 'ytick',[])
    end
%     if i == 1
%         ylabel(y_label2, fontSize=16)
%         c = colorbar('westoutside');
%         c.Position(3) = 0.02;
%     end
    caxis(c_lim)

    title(miceIDs{i});
    
    h(2) = subplot(3,5,10+i);
    shadedErrorBar(time, mean(data2plot,'omitnan'), dataSE,'lineprops', 'b');
    set(gca, 'xlim', [time(1) time(end)])
    if i == 1
        ylabel(y_label, fontSize=10)
    end
    if i == round(n_mice/2)
        xlabel(x_label, fontSize=14)
    end
    xline(0, '--', 'linew', 1)

    %     correct_dim = h(2).Position(3);
    % 
    %     h(1).Position(3) = correct_dim;
    
    % Get data from all mice
    all_data = [all_data; data2plot];
end

%% Plot average of all mice
figure(n_mice+1),

Interval_mean = mean( all_data  (:,center:(center+frameshift)),2);
[max_val,~]   = max(Interval_mean, [], 2);
[~, sort_idx] = sortrows(max_val);

all_data_SE = std(all_data, 'omitnan')./ sqrt(size(all_data,1));
y2          =  [1 size(all_data, 1)];

h1 = subplot(3,1,[1,2]);
imagesc(x1, y2, all_data(sort_idx,:))
font = gca;
font.FontSize = 16;
set(gca, 'xtick',[])
ylabel('Neuron #', fontSize=16)
caxis(c_lim)
c = colorbar;
c.Position(3) = .02;
c.Position(1) = 0.65;
% title(['n mice = ', num2str(n_mice)]);
title('Awake SWRs')
axis square

h2 = subplot(3,1,3);
shadedErrorBar(time, mean(all_data,'omitnan'), all_data_SE,'lineprops', 'b');
set(gca, 'xlim', [time(1) time(end)])
font = gca;
font.FontSize = 16;
xline(0, '--', 'linew', 1)
ylabel(y_label, fontSize=16)
xlabel(x_label, fontSize=16)
h2.Position([1 3]) = [0.395 0.245];
%% Store data
results = struct();

results.miceIDs = miceIDs;
results.signals = mouse_data;
results.SE      = mouse_data_SE;




