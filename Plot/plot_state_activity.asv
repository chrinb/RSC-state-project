function output = plot_state_activity(varargin)

% Written by Christoffer Berge | Vervaeke lab

% Function that plots the mean and standard error of the mean of 
% cross-correlation of different e-phys events. 

% Load data
[data, miceIDs] = load_multi_session_state_activity(varargin);


% Find nr of mice in data
n_mice = numel(fieldnames(data));

mouse_data    = cell(n_mice,12);
mouse_data_SE = cell(n_mice,1);
% Loop over nr of mice
for i = 1:n_mice

    % Get all data for current mouse
    temp_data = data.(miceIDs{i}).state_activity;
    
%     % Select cross-correlation
%     switch varargin{1, 1}  
%         case 'swr_spindle'
%             idx = 1;
%             x_label = 'Time from SWR peak (s)';
%         case 'swr_delta'
%             idx = 2;
%             x_label = 'Time from SWR peak (s)';
%         case 'delta_spindle'
%             idx = 3;
%             x_label = 'Time from SWA trough (s)';
%         case 'delta_swr'
%             idx = 4;
%             x_label = 'Time from SWA trough (s)';
%         case 'spindle_swr'
%             idx = 5;
%             x_label = 'Time from spindle centre (s)';
%         case 'spindle_delta'
%             idx = 6;
%             x_label = 'Time from spindle centre (s)';
%         case 'swr_emg'
%             idx = 7;
%             x_label = 'Time from SWR peak (s)';
%     end
    
    % Loop over nr of sessions for current mouse
    temp_data1 = [];
    mean_NREM_dff = [];
    mean_REM_dff = [];
    mean_quiet_dff = [];
    mean_active_dff = [];
    NREM_transient_rate= [];
    REM_transient_rate= [];
    quiet_transient_rate= [];
    active_transient_rate= [];
    NREM_deconv_rate= [];
    REM_deconv_rate= [];
    quiet_deconv_rate= [];
    active_deconv_rate= [];
    NREM_all_transients= [];
    REM_all_transients= [];
    quiet_all_transients= [];
    active_all_transients= [];

    for ii = 1:size(temp_data,2)
    
        % Compute mean state DF/F
        temp = temp_data{1, ii};
        mean_NREM_dff = [mean_NREM_dff; temp];

        temp1 = temp_data{7, ii};
        mean_REM_dff = [mean_REM_dff; temp1];

        temp2 = temp_data{13, ii};
        mean_active_dff = [mean_active_dff; temp2];

        temp3 = temp_data{19, ii};
        mean_quiet_dff = [mean_quiet_dff; temp3];
       
        % Compute transient rate as a function of state
        temp = temp_data{2, ii};
        NREM_transient_rate = [NREM_transient_rate; temp];

        temp1 = temp_data{8, ii};
        REM_transient_rate = [REM_transient_rate; temp1];

        temp2 = temp_data{14, ii};
        active_transient_rate = [active_transient_rate; temp2];

        temp3 = temp_data{20, ii};
        quiet_transient_rate = [quiet_transient_rate; temp3];

        % Compute /deconvolved rate as a function of state
        temp = temp_data{3, ii};
        NREM_deconv_rate = [NREM_deconv_rate; temp];

        temp1 = temp_data{9, ii};
        REM_deconv_rate = [REM_deconv_rate; temp1];

        temp2 = temp_data{15, ii};
        active_deconv_rate = [active_deconv_rate; temp2];

        temp3 = temp_data{21, ii};
        quiet_deconv_rate = [quiet_deconv_rate; temp3];
        
        % Find percentage of transients per state
        temp = temp_data{4, ii};
        NREM_all_transients = [NREM_all_transients; temp];

        temp1 = temp_data{10, ii};
        REM_all_transients = [REM_all_transients; temp1];

        temp2 = temp_data{16, ii};
        active_all_transients = [active_all_transients; temp2];

        temp3 = temp_data{22, ii};
        quiet_all_transients = [quiet_all_transients; temp3];

        % Find max DF/F
    end

    

    mouse_data{i,1}    = mean_NREM_dff;
    mouse_data{i,2}    = mean_REM_dff;
    mouse_data{i,3}    = mean_active_dff;
    mouse_data{i,4}    = mean_quiet_dff;

    mouse_data{i,5}    = NREM_transient_rate;
    mouse_data{i,6}    = REM_transient_rate;
    mouse_data{i,7}    = quiet_transient_rate;
    mouse_data{i,8}    = active_transient_rate;

    mouse_data{i,9}    = NREM_deconv_rate;
    mouse_data{i,10}    = REM_deconv_rate;
    mouse_data{i,11}    = quiet_deconv_rate;
    mouse_data{i,12}    = active_deconv_rate;

    mouse_data{i,13}    = NREM_all_transients;
    mouse_data{i,14}    = REM_all_transients;
    mouse_data{i,15}    = quiet_all_transients;
    mouse_data{i,16}    = active_all_transients;

%     mouse_data_SE{i} = std(temp_data1, 'omitnan')./sqrt( size(temp_data1, 1));
end


%% Plot

% Plot mean DF/Fyy
all_mean_NREM_dff = [];
all_mean_REM_dff = [];
all_mean_active_dff = [];
all_mean_quiet_dff = [];

figure,
% grayColor = [.6 .6 .6];
x_labels = {'AW', 'QW', 'NREM', 'REM' };
for i = 1:n_mice

    temp =  mouse_data{i, 1};
    all_mean_NREM_dff = [temp;all_mean_NREM_dff];
    mouse_NREM_dff = mean(mouse_data{i, 1}, 'omitnan' );

    temp =  mouse_data{i, 2};
    all_mean_REM_dff = [temp;all_mean_REM_dff];
    mouse_REM_dff = mean(mouse_data{i, 2} , 'omitnan');

    temp =  mouse_data{i, 3};
    all_mean_active_dff = [temp;all_mean_active_dff];
    mouse_active_dff = mean(mouse_data{i, 3} ,'omitnan' );

    temp =  mouse_data{i, 4};
    all_mean_quiet_dff = [temp;all_mean_quiet_dff];
    mouse_quiet_dff = mean(mouse_data{i, 4},'omitnan'  );
% shadedErrorBar([1 2 3 4],...
%     grand_mean,[NaN NaN NaN NaN] ,'lineprops', 'r')
    plot([mouse_active_dff mouse_quiet_dff mouse_NREM_dff mouse_REM_dff], 'LineWidth',1, 'Color',[.6 .6 .6])
    xticks([1 2 3 4])
    xticklabels(x_labels)
        ylabel('Mean DF/F')
    font = gca;
    font.FontSize = 14;
    hold on

    if i == n_mice
        grand_mean    = [mean(all_mean_active_dff,'omitnan'), mean(all_mean_quiet_dff,'omitnan'),...
                         mean(all_mean_NREM_dff, 'omitnan'), mean(all_mean_REM_dff,'omitnan')];
        grand_mean_SE = [ ( std(all_mean_active_dff,'omitnan')./sqrt(size(all_mean_active_dff,1))), ...
                          ( std(all_mean_quiet_dff,'omitnan')./sqrt(size(all_mean_quiet_dff,1))),...
                          ( std(all_mean_NREM_dff,'omitnan')./sqrt(size(all_mean_NREM_dff,1))),...
                          ( std(all_mean_REM_dff,'omitnan')./sqrt(size(all_mean_REM_dff,1))) ];
        shadedErrorBar([1 2 3 4], grand_mean,grand_mean_SE ,'lineprops', 'r')
    end
end

test = 1;
% Get plot labels
% if exist('varargin')
%     [y_label, y_label2, x_label] = plot_options(varargin);
% end
% [x_label, y_label2, y_label] = get_xy_labels(varargin);

% Specify x axis length
% event_center   = round(size(mouse_data{1, 1}  ,2)/2);
% nr_of_seconds  = 3; % nr of seconds before/after spindle onset
% nr_of_frames   = (nr_of_seconds*31*2)+1; % nr of recording frames in spindle window
% frames_to_keep = event_center-floor(nr_of_frames/2):event_center+floor(nr_of_frames/2);
% time           = linspace(-nr_of_seconds,nr_of_seconds, length(frames_to_keep));
% center         = round(size(time,2)/2);
% frameshift     = round(31/2);
% time = -8:1/2500:8;
% plot_lim = dsearchn(time', [-4 4]')
%% Plot results
all_data = [];

y_label = 'Correlation coefficient';


figure,
for i = 1:n_mice
%     figure(i),
    data2plot = mouse_data{i,1};
    dataSE    = mouse_data_SE{i,1};
%     Interval_mean = mean( data2plot  (:,center:(center+frameshift)),2);
%     [max_val,~]   = max(Interval_mean, [], 2);
%     [~, sort_idx] = sortrows(max_val);


    % Set up X and Y limits of imagesc function
%     x1 = [time(1) time(end)];
%     y1 = [1 size(data2plot, 1)]; 

%     h(1) = subplot(3,5,[0+i,5+i]);
%     title(miceIDs{i}),
%     imagesc(x1, y1, data2plot(sort_idx,:))
%     font = gca;
%     font.FontSize = 16;
%     set(gca, 'xtick',[])
%     if i > 1
%         set(gca, 'ytick',[])
%     end
% %     if i == 1
%         ylabel(y_label2, fontSize=16)
%         c = colorbar('westoutside');
%         c.Position(3) = 0.02;
%     end
%     caxis(c_lim)

    
    h(i) = subplot(1,n_mice,i);
    shadedErrorBar(time, mean(data2plot,'omitnan'), dataSE,'lineprops', 'b');
    set(gca, 'xlim', [time(plot_lim(1)) time(plot_lim(2))])
    title(miceIDs{i});
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
linkaxes(h,'x')
%% Plot average of all mice
figure,

% Interval_mean = mean( all_data  (:,center:(center+frameshift)),2);
% [max_val,~]   = max(Interval_mean, [], 2);
% [~, sort_idx] = sortrows(max_val);

all_data_SE = std(all_data, 'omitnan')./ sqrt(size(all_data,1));
% y2          =  [1 size(all_data, 1)];
% 
% h1 = subplot(3,1,[1,2]);
% imagesc(x1, y2, all_data(sort_idx,:))
% font = gca;
% font.FontSize = 16;
% set(gca, 'xtick',[])
% ylabel(y_label2, fontSize=16)
% caxis(c_lim)
% c = colorbar;
% c.Position(3) = .02;
% c.Position(1) = 0.65;
% title(['n mice = ', num2str(n_mice)]);
% axis square

% h2 = subplot(3,1,3);
shadedErrorBar(time, mean(all_data,'omitnan'), all_data_SE,'lineprops', 'b');
set(gca, 'xlim',[time(plot_lim(1)) time(plot_lim(2))])
font = gca;
font.FontSize = 16;
xline(0, '--', 'linew', 1)
ylabel(y_label, fontSize=16)
xlabel(x_label, fontSize=16)
% h2.Position([1 3]) = [0.395 0.245];
%% Store data
results = struct();

results.miceIDs = miceIDs;
results.signals = mouse_data;
results.SE      = mouse_data_SE;




