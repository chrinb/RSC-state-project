function sort_idx = plot_so_sig_an(varargin)

% Written by Christoffer Berge | Vervaeke Lab

% Function that plots results from the function "so_sig_an". The
% variables specificed below must be inputted for the function to return
% plots. 

checkParameter = @(param, n, str) (isnumeric(param) && param==n) || strcmp(param, str);

% Misc variables
sessionID        = varargin{1,1};
opts             = varargin{1,2};
time             = varargin{1,3};

celltype         = varargin{1,12};
label3           = varargin{1,13};
text             = varargin{1,14};
cells_to_exclude = varargin{1,15};
delta_idx        = varargin{1,16};
SO_idx           = varargin{1,17};
nr_of_frames     = varargin{1,18};
% Data variables
signal_SO_activity        = varargin{1,4};
signal_SO_activity_zscore = varargin{1,5};
signal_SO_activity_mean   = varargin{1,6};
signal_SO_activity_meanZ  = varargin{1,7};

signal_delta_activity        = varargin{1,8};
signal_delta_activity_zscore = varargin{1,9};
signal_delta_activity_mean   = varargin{1,10};
signal_delta_activity_meanZ  = varargin{1,11};

% Remove cell population not analyzed (e.g, if analyzing PCs remove INs
% from plot ( they have zero values anyway) ).
if ~isempty(cells_to_exclude) && ~(checkParameter(opts.exp_type, 3, 'axon'))
    signal_delta_activity_mean(cells_to_exclude,:)  = [];
    signal_delta_activity_meanZ(cells_to_exclude,:) = [];

    signal_SO_activity_mean(cells_to_exclude,:)  = [];
    signal_SO_activity_meanZ(cells_to_exclude,:) = [];
end

% Select signal for plot
if checkParameter(opts.exp_type, 2, 'default') || checkParameter(opts.exp_type, 3, 'axon')
    signal_SO_activity        = signal_SO_activity_mean;
    signal_SO_activity_zscore = signal_SO_activity_meanZ;

    signal_delta_activity        = signal_delta_activity_mean;
    signal_delta_activity_zscore = signal_delta_activity_meanZ;
end

% Calculate standard error
signal_SO_SE             = std(signal_SO_activity, 'omitnan') ./ sqrt(size(signal_SO_activity,1));
signal_SO_SE_zscore      = std(signal_SO_activity_zscore, 'omitnan') ./ sqrt(size(signal_SO_activity_zscore,1));

signal_delta_SE         = std(signal_delta_activity, 'omitnan') ./ sqrt(size(signal_delta_activity,1));
signal_delta_SE_zscore  = std(signal_delta_activity_zscore, 'omitnan') ./ sqrt(size(signal_delta_activity_zscore,1));

% 
% % for plotting bulk analysis
% signal_SO_bulk_SE    = std(signal_SO_activity, 'omitnan') ./ sqrt(size(signal_SO_activity,1));
% signal_delta_bulk_SE = std(signal_delta_activity, 'omitnan') ./ sqrt(size(signal_delta_activity,1));
% 
% signal_SO_bulk_SE_zscore    = std(signal_SO_activity_zscore, 'omitnan') ./ sqrt(size(signal_SO_activity_zscore,1));
% signal_delta_bulk_SE_zscore = std(signal_delta_activity_zscore, 'omitnan') ./ sqrt(size(signal_delta_activity_zscore,1));
% 
% % for plotting population analysis
% signal_SO_pop_SE    = std(signal_SO_activity_mean, 'omitnan') ./ sqrt(size(signal_SO_activity_mean,1));
% signal_delta_pop_SE = std(signal_delta_activity_mean, 'omitnan') ./ sqrt(size(signal_delta_activity_mean,1));
% 
% signal_SO_pop_SE_zscore    = std(signal_SO_activity_meanZ, 'omitnan') ./ sqrt(size(signal_SO_activity_meanZ,1));
% signal_delta_pop_SE_zscore = std(signal_delta_activity_meanZ, 'omitnan') ./ sqrt(size(signal_delta_activity_meanZ,1));

if checkParameter(opts.exp_type, 1, 'bulk')
    c_lim        = [-1.5 1.5];
elseif checkParameter(opts.signal_type, 1, 'deconv')
    c_lim        = [-.05 .5];
elseif checkParameter(opts.signal_type, 2, 'dff')
    c_lim        = [-.5 .5];
end


% Sort ROIs according to mean z-scored activity in the -0.5 to + 0.5 interval
% surrounding SWR peak
frameshift    = round(.5/(1/31));
interval_mean = mean(signal_delta_activity_zscore(:, ( median(1:nr_of_frames)-frameshift:median(1:nr_of_frames)+frameshift)),2);
[max_val,~]   = max(interval_mean,[],2);
[~, sort_idx] = sortrows(max_val);

%% Plot results

% Set labels
label1        = ['Mean ' text, ' DF/F'];
label1z       = ['Mean z-scored ' text, ' DF/F'];
x_label_SO    = 'Time from SO trough';
x_label_delta = 'Time from delta trough';

figure,
sgtitle([ sessionID, ', Delta n = ' num2str(length(delta_idx)), ', SO n = ' ...
    num2str(length(SO_idx)), ', ' celltype], 'Interpreter', 'none', 'Fontsize',14) 

x1 = [time(1), time(end)];

subplot(3,4, [1, 5])
y1 = [1 size(signal_SO_activity,1)];
imagesc(x1, y1, signal_SO_activity) 
font = gca;
font.FontSize = 14;
ylabel(label3, 'FontSize',14)
set(gca, 'xtick',[])
c1 = colorbar('southoutside');
c1.Position(4) = 0.015;
c1.Position(2) = 0.38;
    
subplot(3,4,[2,6])
imagesc(x1, y1, signal_SO_activity_zscore) 
font = gca;
font.FontSize = 14;
% ylabel(label3, 'FontSize',14)
set(gca, 'xtick',[], 'YTick',[])
c2 = colorbar('southoutside');
c2.Position(4) = 0.015;
c2.Position(2) = 0.38;
caxis(c_lim)

subplot(3,4,[3,7])
imagesc(x1, y1, signal_delta_activity) 
font = gca;
font.FontSize = 14;
% ylabel(label3, 'FontSize',14)
set(gca, 'xtick',[],'YTick',[])
c3 = colorbar('southoutside');
c3.Position(4) = 0.015;
c3.Position(2) = 0.38;

subplot(3,4,[4,8])
imagesc(x1, y1, signal_delta_activity_zscore) 
font = gca;
font.FontSize = 14;
% ylabel(label3, 'FontSize',14)
set(gca, 'xtick',[],'ytick',[])
c4 = colorbar('southoutside');
c4.Position(4) = 0.015;
c4.Position(2) = 0.38;
caxis(c_lim)


subplot(3,4,9)
shadedErrorBar(time, mean(signal_SO_activity, 'omitnan'),  signal_SO_SE);
font = gca;
font.FontSize = 14;
xline(0, '--', 'linew', 1) % Mark SO trough
xlabel(x_label_SO,  'FontSize',14)
ylabel(label1)
set(gca, 'xlim',[time(1) time(end)])

subplot(3,4,10)
shadedErrorBar(time, mean(signal_SO_activity_zscore, 'omitnan'), signal_SO_SE_zscore);
font = gca;
font.FontSize = 14;
xline(0, '--', 'linew', 1) % Mark SO trough
xlabel(x_label_SO,  'FontSize',14)
% ylabel(label1z)
set(gca, 'xlim',[time(1) time(end)])
    
subplot(3,4,11)
shadedErrorBar(time, mean(signal_delta_activity, 'omitnan'), signal_delta_SE);
font = gca;
font.FontSize = 14;
xline(0, '--', 'linew', 1) % Mark SO trough
xlabel(x_label_delta,  'FontSize',14)
% ylabel(label1)
set(gca, 'xlim',[time(1) time(end)])

subplot(3,4,12)
shadedErrorBar(time, mean(signal_delta_activity_zscore, 'omitnan'), signal_delta_SE_zscore);
font = gca;
font.FontSize = 14;
xline(0, '--', 'linew', 1) % Mark SO trough
xlabel(x_label_delta,  'FontSize',14)
% ylabel(label1z)
set(gca, 'xlim',[time(1) time(end)])
