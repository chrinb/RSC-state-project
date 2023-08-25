function results = plot_mean_state_dynamics(varargin)

% Written by Christoffe Berge | Vervaeke lab

% Load neural activity matrices for 
% 
% User first specifies three input arguments: (1) cell type ('pc', 'in', 'axon'),
% (2) signal type ('dff'), and (3) absolute or normalized time(3) 'norm', 'abs')

% Then user specifies mouse ID (e.g., m6139 m6140)


% Load data
[data_all_mice, miceIDs] = get_data(varargin);

n_mice = size( fieldnames(data_all_mice),1);

if strcmp(varargin{1,3}, 'norm')
    data_text = 'state_mean_all_sessions_cell_norm';
elseif strcmp(varargin{1,3}, 'abs')
    data_text = 'state_mean_all_sessions_cell_abs';
end
%% Plot per mouse data

figure,
sgtitle(['Mean DF/F (z-score) ', varargin{1,1}], FontSize=16);
cmap2 = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.4660, 0.6740, 0.1880], [0.6350, 0.0780, 0.1840]};
labels = {'AW', 'QW', 'NREM', 'REM'};
counter = 1;
% Loop over mice
for mouse_nr = 1:n_mice
    % Get data and SE from current mouse
     mouse_data = data_all_mice.(miceIDs{1,mouse_nr}).(data_text);
     data_to_plot = mouse_data{1,1};

     % Loop over states and plot the mean and SE of each state in separate
     % subplot
     for state_nr = 1:4
         hAx(counter) = subplot(n_mice,4, state_nr+4*(mouse_nr-1));

         if strcmp(varargin{1, 3}, 'norm')
            data_SE  = mouse_data{1,2};
            time_vec = linspace(0, 1, numel(data_to_plot{state_nr}));
         elseif strcmp(varargin{1, 3}, 'abs')
            time_vec = ( 0:size(data_to_plot{state_nr},2)-1 )./31;
            data_SE{state_nr}  = NaN(1, size(data_to_plot{state_nr},2));
         end

         shadedErrorBar(time_vec, data_to_plot{state_nr}, data_SE{state_nr}, 'lineProps',{'color',cmap2{state_nr}});
         set(gca, 'xlim', [time_vec(1) time_vec(end)])
         if mouse_nr == 1
            title(labels{state_nr}, FontSize=14)
         end

         if state_nr == 1 && strcmp(varargin{1, 3}, 'norm')
             text(-.5, 0, 0, miceIDs{1,mouse_nr});
%          elseif state_nr == 1 && strcmp(varargin{1, 3}, 'abs')
%              text(-10, median( [ min(data_to_plot{state_nr}) max(data_to_plot{state_nr}) ]), 0, miceIDs{1,mouse_nr});
         end
         
         counter = counter + 1;
%          axis square
     end
end

% Plot label
if strcmp(varargin{1, 3}, 'norm')
    hAx(2).XLabel.String = 'Normalized time';
elseif strcmp(varargin{1, 3}, 'abs')
    hAx(2).XLabel.String = 'Absolute time';
end

hAx(2).XLabel.FontSize = 16;
% hAx(1).YLabel.String = 'DF/F (z-score)';
% hAx(1).YLabel.FontSize = 16;
%% Plot average for all mice



% Get all data (NOTE: first column is mean state data across sessions,
% second column is associated S.E.M, third is mean for individual sessions
mouse_data_cat = [];
for mouse_nr = 1:n_mice
     mouse_data_cat = [mouse_data_cat; data_all_mice.(miceIDs{1,mouse_nr}).(data_text)];
end
        

% Get state averages per mouse
mouse_state_avg = mouse_data_cat(:,1);

% Loop over nr of mice
state_size = [];
for mouse_nr = 1:n_mice
    % Get data from one mouse
    temp_mouse_data = mouse_state_avg{mouse_nr};
    % Concatenate state vector lengths from each mouse into one matrix
    state_size = [state_size; cell2mat( cellfun(@(x) size(x,2), temp_mouse_data, 'UniformOutput', false))]; 
end

% Find max value and indicies in matrix
[max_val, max_id] = max( state_size, [], 1, 'linear');

% Get session average per mouse
mouse_session_avg = mouse_data_cat(:,3);

mouse_data_to_plot = [];
for mouse_nr = 1:n_mice

    temp_mouse_data = mouse_session_avg{mouse_nr};
    
    % Loop over states
    for state_nr = 1:4

        session_data = temp_mouse_data{state_nr};
        
        % If analyzing time normalized data, interpolate so all vectors
        % have equal lengths.
        if strcmp(varargin{1,3}, 'norm')
    
            if ~( size(session_data,2) == max_val(state_nr) )

                % Interpolate
                session_data_interp    = interp1( 1:size(session_data,2), session_data', linspace(1, size(session_data,2), max_val(state_nr)) )';
                
                % If data that is interpolated is a single vector it might be
                % flipped. If so, transpose it to get it back to normal
                % dimensions.
                if size(session_data_interp,1) > size(session_data_interp,2)
                    session_data_interp = session_data_interp';
                end
            
                temp_mouse_data{state_nr} = session_data_interp;
            end

        % If analyzing absolute time data, pad vectors with NaNs so they
        % get the same lengths
        elseif strcmp(varargin{1,3}, 'abs')

            session_data_padded = session_data;
            session_state_size  = size(session_data,2);
            session_data_padded(:,  session_state_size:max_val(state_nr)) = NaN;
            temp_mouse_data{state_nr} = session_data_padded;

        end
    end
    
    % Store interpolated/padded data in new cell array
    mouse_data_to_plot = [mouse_data_to_plot; temp_mouse_data];
    clear temp_mouse_data
end


%% Plot    
figure,
 cmap2 = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.4660, 0.6740, 0.1880], [0.6350, 0.0780, 0.1840]};
labels = {'AW', 'QW', 'NREM', 'REM'};
for state_nr = 1:4
    data_to_plot    = cell2mat( mouse_data_to_plot(:, state_nr) );
%     data_SE         = std(data_to_plot, 'omitnan')./ sqrt( size(data_to_plot,1));

    if strcmp(varargin{1, 3}, 'norm')
        data_SE  = std(data_to_plot, 'omitnan')./ sqrt( size(data_to_plot,1));
        time_vec = linspace(0, 1, size(data_SE,2));
    elseif strcmp(varargin{1, 3}, 'abs')
        time_vec = ( 0:size(data_to_plot,2)-1 )./31;
        data_SE  = NaN(1, size(data_to_plot,2));
    end
     
    hAx(state_nr) = subplot(1,4, state_nr);
    shadedErrorBar(time_vec, mean(data_to_plot,'omitnan'), data_SE,  'lineProps',{'color',cmap2{state_nr}});
    set(gca, 'xlim', [time_vec(1) time_vec(end)])
    title([labels{state_nr}, ' (n bouts = ', num2str(size(data_to_plot,1)),')'], FontSize=14)
    axis square

end

% Plot label
if strcmp(varargin{1, 3}, 'norm')
    hAx(2).XLabel.String = 'Normalized time';
elseif strcmp(varargin{1, 3}, 'abs')
    hAx(2).XLabel.String = 'Absolute time';
end
hAx(2).XLabel.FontSize = 16;
hAx(1).YLabel.String = 'DF/F (z-score) ';
hAx(1).YLabel.FontSize = 16;
sgtitle(varargin{1,1});
%% Store data
% results = struct();
% 
% results.miceIDs = miceIDs;
% results.signals = all_mice_data_mat_zscore;
% results.SE      = all_data_SE;

end

%% Load mouse data

function [data_all_mice, miceIDs] = get_data(varargin)

    % Select particular dataset
%     prompt     = sprintf('Specify ephys signal: ');
%     ephys_type = input(prompt, 's');

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
    
      %Load data
      data_str = ['_state_dynamics_',[varargin{1,1}{1,3}], '_', [varargin{1,1}{1,1},'_', varargin{1,1}{1,2}], '.mat' ];

    % Loop over nr of mice
    for i = 1:numel(miceIDs)
        % Select absolute path to data
        data_path = [ path{i} miceIDs(i) data_str ];
        data_path = append(data_path(1), data_path(2), data_path(3) );
        % Load data into struct
        data_all_mice.(miceIDs{i}) = load(data_path );
    end
end
