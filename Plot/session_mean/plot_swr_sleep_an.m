function sort_idx = plot_swr_sleep_an(varargin)

%Written by Christoffer Berge | Vervaeke Lab

% Function that plots results from the function "swr_sleep_an". The
% variables specificed below must be inputted for the function to return
% plots. 

checkParameter = @(param, n, str) (isnumeric(param) && param==n) || strcmp(param, str);

% Misc variables
sessionID        = varargin{1,1};
opts             = varargin{1,2};
time             = varargin{1,3};
plot_z           = varargin{1,17};
celltype         = varargin{1,18};
cells_to_exclude = varargin{1,19};
label3           = varargin{1,20};
text             = varargin{1,21};

%% Data variables
% Bulk variables
signal_swr_bulk_unclassified = varargin{1,4};
signal_swr_bulk_single       = varargin{1,5};
signal_swr_bulk_coupled      = varargin{1,6};

signal_swr_bulkZ_unclassified = varargin{1,7};
signal_swr_bulkZ_single       = varargin{1,8};
signal_swr_bulkZ_coupled      = varargin{1,9};

% Population variables
signal_swr_pop_mean_unclassified = varargin{1,10};
signal_swr_pop_mean_single       = varargin{1,11};
signal_swr_pop_mean_coupled      = varargin{1,12};

signal_swr_pop_meanZ_unclassified  = varargin{1,13};
signal_swr_pop_meanZ_single        = varargin{1,14};
signal_swr_pop_meanZ_coupled       = varargin{1,15};

%% Select signals for plot
if checkParameter(opts.exp_type, 1, 'bulk') && plot_z == 0
    signal_swr_unclassified = signal_swr_bulk_unclassified;
    signal_swr_single       = signal_swr_bulk_single;
    signal_swr_coupled      = signal_swr_bulk_coupled;
elseif checkParameter(opts.exp_type, 1, 'bulk') && plot_z == 1
    signal_swr_unclassified = signal_swr_bulkZ_unclassified;
    signal_swr_single       = signal_swr_bulkZ_single;
    signal_swr_coupled      = signal_swr_bulkZ_coupled;
elseif plot_z == 0
    signal_swr_unclassified = signal_swr_pop_mean_unclassified;
    signal_swr_single       = signal_swr_pop_mean_single;
    signal_swr_coupled      = signal_swr_pop_mean_coupled;
elseif plot_z == 1
    signal_swr_unclassified = signal_swr_pop_meanZ_unclassified;
    signal_swr_single       = signal_swr_pop_meanZ_single;
    signal_swr_coupled      = signal_swr_pop_meanZ_coupled;
end

% Remove cell population not analyzed (e.g, if analyzing PCs remove INs
% from plot ( they have zero values anyway) ).
if ~isempty(cells_to_exclude)
    signal_swr_unclassified(cells_to_exclude,:) = [];
    signal_swr_single(cells_to_exclude,:)       = [];
    signal_swr_coupled(cells_to_exclude,:)      = [];
end
%% Calculate standard error
signal_SE_unclassified       = std(signal_swr_unclassified, 'omitnan') ./ ...
    sqrt(size(signal_swr_unclassified,1));
signal_SE_single             = std(signal_swr_single, 'omitnan') ./ ...
    sqrt(size(signal_swr_single,1));
signal_SE_coupled             = std(signal_swr_coupled, 'omitnan') ./ ...
    sqrt(size(signal_swr_coupled,1));
%% Select data to plot

if checkParameter(opts.exp_type, 1, 'bulk')
    scale_factor = .3;
    y_label2     = 'SWR #';
    c_lim        = [-2 2];
elseif checkParameter(opts.signal_type, 1, 'deconv')
    y_label2      = 'ROI #';
    scale_factor = 0.02;
    c_lim        = [-.05 .5];
elseif checkParameter(opts.signal_type, 2, 'dff')
    y_label2      = 'ROI #';
    scale_factor = 0.02;
    c_lim        = [-.5 .5];
end

%% Determine max-min values of Y limits in average plots

min_unclassified = min( mean(signal_swr_unclassified, 'omitnan'));
min_single       = min( mean(signal_swr_single, 'omitnan'));
min_coupled      = min( mean(signal_swr_coupled, 'omitnan'));

max_unclassified = max( mean(signal_swr_unclassified, 'omitnan'));
max_single       = max( mean(signal_swr_single, 'omitnan'));
max_coupled      = max( mean(signal_swr_coupled, 'omitnan'));

Y_lim_max = max([max_unclassified, max_single,max_coupled]) + scale_factor;
Y_lim_min = min([min_unclassified, min_single,min_coupled]) - scale_factor;

% Sort ROIs according to mean z-scored activity in the -0.5 to + 0.5 interval
% surrounding SWR peak in spindle-coupled SWRs
% sort_idx = ':';
% if ~(opts.exp_type == 1)
    frameshift    = round(.5/(1/31));
    interval_mean = mean(signal_swr_coupled(:, ( median(1:187)-frameshift:median(1:187)+frameshift)),2);
    [max_val,~]   = max(interval_mean,[],2);
    [~, sort_idx] = sortrows(max_val);
% end
%% Plot results

% Set labels

if plot_z == 1
    y_label = ['Mean z-scored ' text, ' DF/F'];
else
    y_label  = ['Mean ' text, ' DF/F'];
end
x_label = 'Time from SWR peak (sec)';

figure,
sgtitle([ sessionID, ', ', celltype], 'Interpreter', 'none') 

% Set imagesc X and Y 
x1             = [time(1) time(end)];
y_unclassified = [1 size(signal_swr_unclassified,1) ];
y_single       = [1 size(signal_swr_single,1) ];
y_coupled      = [1 size(signal_swr_coupled,1) ];

%% unclassified SWRs
if length(mean(signal_swr_unclassified, 'omitnan')) > 1

% only sort PETHs in population analyses (because of uneven nr of SWRs per
% category...)
if checkParameter(opts.exp_type, 1, 'bulk')
    sig1 = signal_swr_unclassified;
else
    sig1 = signal_swr_unclassified(sort_idx,:);
end
subplot(3, 3, [1 4])
imagesc(x1,y_unclassified, sig1) %  colorplot of average individual roi activity during ripples
ylabel(y_label2)
xlabel(x_label)
c(1) = colorbar;
c(1).Position(1) = 0.35;
c(1).Position(2) = 0.39;
c(1).Position(3) = 0.006;
if plot_z == 1
    caxis(c_lim);
end
title(['SWR n = ', num2str(size(signal_swr_bulk_unclassified,1)) ])


subplot(3,3, 7)
shadedErrorBar(time, mean(signal_swr_unclassified, 'omitnan'), signal_SE_unclassified,'lineprops', 'r');
xline(0, '--', 'linew', 1) % Mark SWR peak time
if plot_z == 1
    set(gca, 'ylim', [Y_lim_min, Y_lim_max])
end
xlabel(x_label)
ylabel(y_label)
legend('Unclassified')
end

%% Spindle-uncoupled SWRs
subplot(3,3, [2 5])
if checkParameter(opts.exp_type, 1, 'bulk')
    sig2 = signal_swr_single;
else
    sig2 = signal_swr_single(sort_idx,:);
end

imagesc(x1, y_single, sig2) %  colorplot of average individual roi activity during ripples
ylabel(y_label2)
xlabel(x_label)
c(2) = colorbar;
c(2).Position(1) = 0.63;
c(2).Position(2) = 0.39;
c(2).Position(3) = 0.006;
if plot_z == 1
    caxis(c_lim);
end
title(['SWR n = ', num2str(size(signal_swr_bulk_single,1)) ])


subplot(3,3,8)
shadedErrorBar(time, mean(signal_swr_single,'omitnan'), signal_SE_single,'lineprops', 'b');
xline(0, '--', 'linew', 1)
if plot_z == 1
    set(gca, 'ylim', [Y_lim_min, Y_lim_max])
end
ylabel(y_label)
xlabel(x_label)
legend('Spindle-uncoupled')


%% Spindle-coupled SWRs
subplot(3,3,[3 6])
% if opts.exp_type == 1
%     sig3 = signal_swr_coupled;
% else
    sig3 = signal_swr_coupled(sort_idx,:);
% end

imagesc(x1,y_coupled, sig3) %  colorplot of average individual roi activity during ripples
xlabel(x_label)
ylabel(y_label2)
c(3) = colorbar;
c(3).Position(1) = 0.91;
c(3).Position(2) = 0.39;
c(3).Position(3) = 0.006;
if plot_z == 1
    caxis(c_lim);
end
title(['SWR n = ', num2str(size(signal_swr_bulk_coupled,1)) ])


if size(signal_swr_coupled,1) > 1
    subplot(3,3,9)
    shadedErrorBar(time, mean(signal_swr_coupled, 'omitnan'), signal_SE_coupled,'lineprops', 'k');
    xline(0, '--', 'linew', 1)
    if plot_z == 1
        set(gca, 'ylim', [Y_lim_min, Y_lim_max])
    end
    xlabel(x_label)
    ylabel(y_label)
    legend('Spindle-coupled')
end
