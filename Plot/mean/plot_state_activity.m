function results = plot_state_activity(varargin)

[data_struct, miceIDs] = get_data;

n_mice = size( fieldnames(data_struct),1);

% Loop over nr of mice
for i = 1:n_mice
    data_all_mice{i,1} = data_struct.(miceIDs{1,i}).state_activity;
end

% data_all_mice = { data_all_mice.m6134.state_activity; data_all_mice.m6137.state_activity  ;...
%     data_all_mice.m6139.state_activity  ; data_all_mice.m6140.state_activity  };

%% Loop NREM 
nrem_nr_in_list = 1;
nrem_dec_list   = 3;
all_nrem_dff    = [];
all_nrem_dec    = [];
mouse_data_nrem = [];
for i = 1:n_mice
   
    nrem_dff = [];
    nrem_dec = [];
    for j = 1:size(data_all_mice{i, 1}  ,2)
        if isempty(data_all_mice{i,1}{nrem_nr_in_list,j})
            data_all_mice{i,1}{nrem_nr_in_list,j} = NaN;
            data_all_mice{i,1}{nrem_dec_list,j} = NaN;
        end

        nrem_dff = [nrem_dff, data_all_mice{i,1}{nrem_nr_in_list,j} ];
        nrem_dec = [nrem_dec, mean(data_all_mice{i,1}{nrem_dec_list,j}) ]; 
    end
    temp_nrem_dff      = nrem_dff;
    mouse_data_nrem{i} = mean(nrem_dff,'omitnan');
    per_mouse_nrem{i}  = nrem_dff;
    all_nrem_dff       = [all_nrem_dff, temp_nrem_dff];

    temp_nrem_dec          = nrem_dec;
    mouse_data_nrem_dec{i} = mean(nrem_dec, 'omitnan');
    per_mouse_nrem_dec{i}  = nrem_dff;
    all_nrem_dec           = [all_nrem_dec, temp_nrem_dec];
end

%% Loop REM

all_rem_dff    = [];
all_rem_dec    = [];
rem_nr_in_list = 7;
rem_dec_list   = 9;

mouse_data_rem = [];
for i = 1:n_mice
    
    rem_dff = [];
    rem_dec = [];
    for j = 1:size(data_all_mice{i, 1}  ,2)

        if isempty(data_all_mice{i,1}{rem_nr_in_list,j})
            data_all_mice{i,1}{rem_nr_in_list,j} = NaN;
            data_all_mice{i,1}{rem_dec_list,j} = NaN;
        end

        rem_dff = [rem_dff, data_all_mice{i,1}{rem_nr_in_list,j} ];
        rem_dec = [rem_dec, mean( data_all_mice{i,1}{rem_dec_list,j}) ];
    end
    temp_rem          = rem_dff;
    mouse_data_rem{i} = mean(rem_dff,'omitnan');
    per_mouse_rem{i}  = rem_dff;
    all_rem_dff       = [all_rem_dff, temp_rem];

    temp_rem_dec          = rem_dec;
    mouse_data_rem_dec{i} = mean(rem_dec, 'omitnan');
    per_mouse_rem_dec{i}  = rem_dff;
    all_rem_dec           = [all_nrem_dec, temp_rem_dec];
end

%% Loop AW

all_active_dff        = [];
all_active_dec        = [];
active_nr_in_list     = 13;
active_dec_list       = 15;
mouse_data_active     = [];
mouse_data_active_dec = [];
for i = 1:n_mice
    
    active_dff = [];
    active_dec = [];
    for j = 1:size(data_all_mice{i, 1}  ,2)

        if isempty(data_all_mice{i,1}{active_nr_in_list,j})
            data_all_mice{i,1}{active_nr_in_list,j} = NaN;
            data_all_mice{i,1}{active_dec_list,j} = NaN;
        end
        active_dff = [active_dff, data_all_mice{i,1}{active_nr_in_list,j} ];
        active_dec = [active_dec, mean(data_all_mice{i,1}{active_dec_list,j},'omitnan') ];
    end
    temp_active          = active_dff;
    mouse_data_active{i} = mean(active_dff,'omitnan');
    per_mouse_active{i}  = active_dff;
    all_active_dff       = [all_active_dff, temp_active];

    temp_active_dec          = active_dec;
    mouse_data_active_dec{i} = mean(active_dec, 'omitnan');
    per_mouse_active_dec{i}  = active_dff;
    all_active_dec           = [all_active_dec, temp_active_dec];
end


%% Loop QW
all_quiet_dff    = [];
all_quiet_dec       =[];
quiet_nr_in_list = 19;
quiet_dec_list   = 21;
mouse_data_quiet = [];
mouse_data_quiet_dec = [];
for i = 1:n_mice
    
    quiet_dff = [];
    quiet_dec = [];
    for j = 1:size(data_all_mice{i, 1}  ,2)

         if isempty(data_all_mice{i,1}{quiet_nr_in_list,j})
            data_all_mice{i,1}{quiet_nr_in_list,j} = NaN;
            data_all_mice{i,1}{quiet_dec_list,j} = NaN;
        end
       quiet_dff = [quiet_dff, data_all_mice{i,1}{quiet_nr_in_list,j} ];
       quiet_dec = [quiet_dec, mean( data_all_mice{i,1}{quiet_dec_list,j}) ];

    end
    temp_quiet           = quiet_dff;
    mouse_data_quiet{i} = mean(quiet_dff,'omitnan');
    per_mouse_quiet{i}  = quiet_dff;
    all_quiet_dff        = [all_quiet_dff, temp_quiet];

    temp_quiet_dec          = quiet_dec;
    mouse_data_quiet_dec{i} = mean(quiet_dec,'omitnan');
    per_mouse_quiet_dec{i}  = quiet_dec;
    all_quiet_dec           = [all_quiet_dec, temp_quiet_dec];
end


%% Create overall averages
state_averages = [mean(all_active_dff, 'omitnan'), mean(all_quiet_dff,'omitnan'), ...
    mean(all_nrem_dff, 'omitnan'), mean(all_rem_dff, 'omitnan')];

%% Plot per mouse average
per_mouse_avg = cell2mat( [mouse_data_active; mouse_data_quiet; mouse_data_nrem; mouse_data_rem] );

figure, 
plot(per_mouse_avg, 'Color',[.5, .5, .5]);
hold on

SE_all_avg = std(per_mouse_avg,'omitnan')./sqrt( size(per_mouse_avg,1));
shadedErrorBar([], mean(per_mouse_avg,2), SE_all_avg, 'lineProps', 'r')

xticks(1:4)
xticklabels({'AW', 'QW', 'NREM','REM'})
ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;
ylabel('Mean DF/F', 'FontSize',16)
%% Boxplot
% mean_dff1 = [all_active_dff, all_quiet_dff, all_nrem_dff, all_rem_dff];
% g1 = repmat({'AW'},37,1);
% g2 = repmat({'QW'},37,1);
% g3 = repmat({'NREM'},37,1);
% g4 = repmat({'REM'},37,1);
% g = [g1; g2; g3; g4];
% figure, 
% boxplot(mean_dff1, g);
% ylabel('Mean DF/F', 'FontSize',16)
%% Violin
mean_dff = {all_active_dff', all_quiet_dff', all_nrem_dff', all_rem_dff'};
categories = {'AW', 'QW', 'NREM', 'REM'};
figure, violin(mean_dff, 'xlabel', categories, 'facecolor', [.9 .9 .9; 0 0 0; .6 .6 .6; .3 .3 .3]);
ylabel('Mean DF/F', 'FontSize',16)

mean_dec =  {all_active_dec', all_quiet_dec', all_nrem_dec', all_rem_dec'};
figure, violin(mean_dec, 'xlabel', categories);
ylabel('Mean deconvolved ', 'FontSize',16)

%% Compute stats
mean_dff2 =  [all_active_dff; all_quiet_dff; all_nrem_dff; all_rem_dff];

% Kruskal-Wallis non-parametric test
[p, tbl, stats] = kruskalwallis(mean_dff2');

c = multcompare(stats);

tbl2 = array2table(c,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);

end

%% 

function [data_all_mice, miceIDs] = get_data
    % Select particular dataset
    prompt = sprintf('Specify celltype ');
    str   = input(prompt, 's');
%     if ~isempty(str)
%         varargin = [varargin, str]; 
%     
% %     end
%     % Load data and mouse IDs (user types in e.g. "m6120")    
%     [data_all_mice, miceIDs] = load_multi_session_swr(str);
    
    prompt = sprintf('Type mouse ID: '); % user types in mouse ID, e.g. "m6120"
    
    miceIDs = input(prompt, 's');
    miceIDs = split(miceIDs, ' ')';
    
    % Find absolute path to data
    path = cell(size(miceIDs,2), 1 );
    
    % Loop over nr of mice
    for i = 1:size(miceIDs,2)
        full_miceIDs = replaceBetween(miceIDs{1,i}, 'm', '6', 'ouse');
        path{i}      = strcat('D:\Data\', convertCharsToStrings(full_miceIDs), '\Results\');
    end
    
    % Load data
    if strcmp(str, 'in')
        data_str = '_state_activity_in.mat';
    elseif strcmp(str, 'pc')
        data_str = '_state_activity_pc.mat';
    elseif strcmp(str, 'axon')
        data_str = '_state_activity_axon.mat';
    end
    
    % Loop over nr of mice
    for i = 1:numel(miceIDs)
        % Select absolute path to data
        data_path = [ path{i} miceIDs(i) data_str ];
        data_path = append(data_path(1), data_path(2), data_path(3) );
        % Load data into struct
        data_all_mice.(miceIDs{i}) = load(data_path );
    end
end
