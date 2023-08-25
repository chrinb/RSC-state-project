function results = plot_mean_ephys(varargin)

% Written by Christoffe Berge | Vervaeke lab

% Load SWR-aligned ECoG, EMG, or running data data for all mice, and plot per-mouse data
% and grand averages for all mice
% 
% User specifies three input arguments: (1) behavioral state ('awake' or 'sleep'),
% (2) signal type ('emg' or 'ecog'), (3) mouse ID (e.g., 'm6134')


% Load data
[data_all_mice, miceIDs, ephys_type] = get_data(varargin);

n_mice = size( fieldnames(data_all_mice),1);

all_mice_data_mat        = [];
all_mice_data_mat_zscore = [];
% Loop over nr of mice
for i = 1:n_mice

    % Load data from current mouse
    mouse_data = data_all_mice.(miceIDs{1,i}).multi_session_swr_ephys;
    
    single_mouse_data        = [];
    single_mouse_data_zscore = [];
    % Loop over sessions for current mouse
    for k = 1:size(mouse_data,1)-1
        
        temp_data        = mouse_data{k,1}{1,1};
        temp_data_zscore = mouse_data{k,1}{2,1};

        single_mouse_data        = [single_mouse_data; temp_data ];
        single_mouse_data_zscore = [single_mouse_data_zscore; temp_data_zscore];
    end

    % Add data for all mice together in one matrix
    single_mouse_data_store{i}        = single_mouse_data;
    single_mouse_data_store_zscore{i} = single_mouse_data_zscore;
    all_mice_data_mat        = [all_mice_data_mat; single_mouse_data];
    all_mice_data_mat_zscore = [all_mice_data_mat_zscore; single_mouse_data_zscore];
end

%% Plot per mouse data
srate = 2500;
time  = (-(srate*3):(srate*3))./srate;

if strcmp(ephys_type, 'ecog')
    c_lim = [-3 3];
elseif strcmp(ephys_type, 'emg')
    c_lim = [-1 1];
elseif strcmp(ephys_type, 'run')
    c_lim = [0 .5];
end


figure, 

for i = 1:n_mice

     hAx(i) = subplot(3, n_mice, [i,n_mice+i]);
     temp = single_mouse_data_store_zscore{1,i};
     n_swr = [1, size(temp,1)];
     imagesc(time, n_swr, temp)
     set(gca,'xtick',[])
     title(miceIDs(i))
     
     hAx(n_mice+i) = subplot(3, n_mice, (n_mice*2+i) );
     SE = std(temp,'omitnan')./sqrt( size(temp,1));
     shadedErrorBar(time, mean(temp,1, 'omitnan'), SE)
     xline(0, 'r--', 'linew', 1) % Mark SWR peak time
     set(gca, 'xlim', [time(1) time(end)])
end
hAx(round(median(n_mice+1:n_mice*2))).XLabel.String = 'Time from SWR peak (s)';
hAx(round(median(n_mice+1:n_mice*2))).XLabel.FontSize = 16;
hAx(1).XLabel.FontSize = 16;
hAx(1).YLabel.String = 'SWR #';
hAx(n_mice+1).YLabel.String = 'Mean \muv';
hAx(1).YLabel.FontSize = 16;
hAx(n_mice+1).YLabel.FontSize = 16;

set(hAx, 'clim', [c_lim])
subplot(3, n_mice, [i,n_mice+i]);
c = colorbar;
c.Position(1) = 0.93;
c.Position(2) = 0.39;

sgtitle([varargin{1,1},' SWR-aligned ' ephys_type])

linkaxes(hAx, 'x')

set(gca, 'xlim', [-3 3])
%% Plot average for all mice
figure,

all_data_SE = std(all_mice_data_mat_zscore, 'omitnan')./ sqrt(size(all_mice_data_mat_zscore,1));
n_swr_all   =  [1 size(all_mice_data_mat_zscore, 1)];

h1 = subplot(3,1,[1,2]);
imagesc(time, n_swr_all, all_mice_data_mat_zscore)
font = gca;
font.FontSize = 16;
set(gca, 'xtick',[])
ylabel('SWR #', fontSize=16)
caxis([-3 3])
c = colorbar;
c.Position(3) = .02;
c.Position(1) = 0.65;
title([varargin{1,1}, ' mean SWR-aligned ' ephys_type])
axis square

h2 = subplot(3,1,3);
shadedErrorBar(time, mean(all_mice_data_mat_zscore,'omitnan'), all_data_SE,'lineprops', 'b');
set(gca, 'xlim', [time(1) time(end)])
xline(0, 'r--', 'linew', 1) % Mark SWR peak time
xlabel('Time from SWR peak (s)')
h2.Position([1 3]) = [0.395 0.245];

h1.XLabel.FontSize = 16;
h2.XLabel.FontSize = 16;
h1.YLabel.String = 'SWR #';
h2.YLabel.String = 'Mean \muv';
h1.YLabel.FontSize = 16;
h2.YLabel.FontSize = 16;
set(h1, 'clim', [c_lim])

linkaxes([h1, h2], 'x')

% set(gca, 'xlim', [-1 1])
%% Store data
results = struct();

results.miceIDs = miceIDs;
results.signals = all_mice_data_mat_zscore;
results.SE      = all_data_SE;

end

%% Load mouse data

function [data_all_mice, miceIDs, ephys_type] = get_data(varargin)

    % Select particular dataset
    prompt     = sprintf('Specify ephys signal: ');
    ephys_type = input(prompt, 's');

    prompt = sprintf('Type mouse ID: '); % user types in mouse ID, e.g. "m6120"
    
    miceIDs = input(prompt, 's');
    miceIDs = split(miceIDs, ' ')';
    
    % Find absolute path to data
    path = cell(size(miceIDs,2), 1 );
    
    % Loop over nr of mice
    for i = 1:size(miceIDs,2)
        full_miceIDs = replaceBetween(miceIDs{1,i}, 'm', '6', 'ouse');
        
        % Find hard drive name
        if isfolder('D:\Data\')
            hd_path = 'D:\Data\';
        elseif isfolder('E:\Data\')
            hd_path = 'E:\Data\';
        end

        path{i}      = strcat(hd_path, convertCharsToStrings(full_miceIDs), '\Results\');
    end
    
    % Load data
    if strcmp(ephys_type, 'ecog')
        data_str = ['_avg_', varargin{1,1},'_swr_ECoG.mat'];
    elseif strcmp(ephys_type, 'emg')
        data_str = ['_avg_', varargin{1,1},'_swr_EMG.mat'];
    elseif strcmp(ephys_type, 'run')
        data_str = ['_avg_', varargin{1,1},'_swr_Running speed.mat'];
    end
    
    % Loop over nr of mice
    for i = 1:numel(miceIDs)
        % Select absolute path to data
        data_path = [ path{i} miceIDs(i) data_str ];
        data_path = append(data_path(1), data_path(2), data_path(3), data_path(4), data_path(5) );
        % Load data into struct
        data_all_mice.(miceIDs{i}) = load(data_path );
    end
end
