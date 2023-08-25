function plot_results(varargin)


% Written by Christoffer Berge || Vervaeke lab

% Function that plots color plots and mean activity of roi x time matrices

signal    = varargin{1,1}; % roi x time matrix
signal_SE = varargin{1,2}; % 1 x time vector containg standard error during each time point

if size(varargin,2) > 2
    [y_label, x_label] = plot_options(varargin);
end

% Specify x axis length
event_center  = round(size(signal,2)/2);
nr_of_seconds  = 3; % nr of seconds before/after spindle onset
nr_of_frames  = (nr_of_seconds*31*2)+1; % nr of recording frames in spindle window
frames_to_keep = event_center-floor(nr_of_frames/2):event_center+floor(nr_of_frames/2);
time          = linspace(-nr_of_seconds,nr_of_seconds, length(frames_to_keep));

% trim data
signal = signal(:, frames_to_keep);
signal_SE = signal_SE(frames_to_keep);
%% Sort according to peak value
center        = round(size(time,2)/2);
frameshift    = round(31/2);
interval_mean = mean( signal(:,center:(center+frameshift)),2);
[max_val,~]   = max(interval_mean,[],2);
[~, sort_idx] = sortrows(max_val);
%% Plot results

% Set up X and Y limits of imagesc function
x1 = [time(1) time(end)];
y1 = [1 size(signal,1)]; %

figure,

h(1) = subplot(3,1,[1,2]);
imagesc(x1, y1, signal(sort_idx,:))
ylabel('# ROI')
caxis([-.05 .5])
c = colorbar;
c.Position(3) = .02;
c.Position(1) = 0.912;

h(2) = subplot(3,1,3);
shadedErrorBar(time, mean(signal,'omitnan'), signal_SE,'lineprops', 'b');
ylabel(y_label)
xlabel(x_label)

correct_dim = h(2).Position(3);

h(1).Position(3) = correct_dim;
% colormap viridis