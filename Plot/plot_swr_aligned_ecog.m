function plot_swr_aligned_ecog(varargin)

% Written by Christoffer Berge | Vervaeke Lab

% Function that plots individual & mean SWR-aligned cortical ECoG/LFP 
% traces based on input from "swr_aligned_ecog" function. 

x_label  = 'Time from SWR peak (sec)';
y_label1 = '# SWR';
y_label2 = 'Mean mV';

sessionID   = varargin{1,1};
time_vec    = varargin{1,2};
sig_to_plot = varargin{1,3};
mean_signal = varargin{1,4};
baseSub     = varargin{1,5};

% Plot "raw" or z-scored normalized data depending on input
if sig_to_plot == 1
    plot_unclassified = mean_signal{1,1};
    plot_single       = mean_signal{1,3};
    plot_coupled      = mean_signal{1,5};
elseif sig_to_plot == 2
    plot_unclassified = mean_signal{1,2};
    plot_single       = mean_signal{1,4};
    plot_coupled      = mean_signal{1,6};
end

% Compute standard error of the mean
SE_unclassified = std(plot_unclassified, 'omitnan') ./ sqrt(size(plot_unclassified,1));
SE_single       = std(plot_single, 'omitnan') ./ sqrt(size(plot_single,1));
SE_coupled      = std(plot_coupled, 'omitnan') ./ sqrt(size(plot_coupled,1));


%% Determine max-min values of Y limits in average plots

max_unclassified = max( nanmean(plot_unclassified));
max_single       = max( nanmean(plot_single));
max_coupled      = max( nanmean(plot_coupled));

min_unclassified = min( nanmean(plot_unclassified));
min_single       = min( nanmean(plot_single));
min_coupled      = min( nanmean(plot_coupled));

Y_lim_max = max([max_unclassified, max_single,max_coupled]) + 0.05;
Y_lim_min = min([min_unclassified, min_single,min_coupled]) - 0.05;
%% Plot results

% Imagesc x and y limits
x1             = [time_vec(1) time_vec(end)];
y_unclassified = [1 size(plot_unclassified,1) ];
y_single       = [1 size(plot_single,1) ];
y_coupled      = [1 size(plot_coupled,1) ];


figure;
sgtitle(sessionID),

% Unclassified SWRs
ax1 = subplot(3, 3, [1 4]);
imagesc(x1, y_unclassified, plot_unclassified)
ylabel(y_label1)
xlabel(x_label)
colorbar

subplot(337)
shadedErrorBar(time_vec, mean(plot_unclassified,'omitnan'),SE_unclassified,'lineprops', 'r');
set(gca, 'ylim', [Y_lim_min, Y_lim_max])
xlabel(x_label)
ylabel(y_label2)
legend('Unclassified')

% Spindle-uncoupled (single) SWRs
ax2 = subplot(3, 3, [2 5]);
imagesc(x1, y_single, plot_single)
ylabel(y_label1)
xlabel(x_label)
colorbar

subplot(3,3,8)
shadedErrorBar(time_vec, mean(plot_single,'omitnan'),SE_single,'lineprops', 'b');
set(gca, 'ylim', [Y_lim_min, Y_lim_max])
xlabel(x_label)
ylabel(y_label2)
legend('Spindle-uncoupled')

% Spindle-coupled SWRs
ax3 = subplot(3, 3, [3 6]);
imagesc(x1, y_coupled, plot_coupled);
ylabel(y_label1)
xlabel(x_label)
colorbar

subplot(3,3,9)
shadedErrorBar(time_vec, mean(plot_coupled,'omitnan'),SE_coupled,'lineprops', 'k');
set(gca, 'ylim', [Y_lim_min, Y_lim_max])
xlabel(x_label)
ylabel(y_label2)
legend('Spindle-coupled')

% Adjust color axis limits 
if sig_to_plot == 1 && ~isempty(baseSub)
    ax1.CLim = [-.4 .4];
    ax2.CLim = [-.4 .4];
    ax3.CLim = [-.4 .4];
elseif sig_to_plot == 2
    ax1.CLim = [-3 3];
    ax2.CLim = [-3 3];
    ax3.CLim = [-3 3];
end