function sort_idx = plot_swr_awake_an(varargin)

% Written by Christoffer Berge | Vervaeke Lab

% Function that plots results from the function "swr_awake_an". The
% variables specificed below must be inputted for the function to return
% plots. 

checkParameter = @(param, n, str) (isnumeric(param) && param==n) || strcmp(param, str);

% Misc variables
sessionID        = varargin{1,1};
opts             = varargin{1,2};
time             = varargin{1,3};
celltype         = varargin{1,8};
RippleIdx        = varargin{1,9};
text             = varargin{1,10};
label3           = varargin{1,11};
cells_to_exclude = varargin{1,12};


%% Data variables
signal_swr_activity              = varargin{1,4};
signal_swr_activity_zscore       = varargin{1,5};
signal_swr_activity_mean         = varargin{1,6};
signal_swr_activity_mean_zscore  = varargin{1,7};

% Remove cell population not analyzed (e.g, if analyzing PCs remove INs
% from plot ( they have zero values anyway) ).
if ~isempty(cells_to_exclude) && ~(checkParameter(opts.exp_type, 3, 'axon'))
    signal_swr_activity_mean(cells_to_exclude,:)  = [];
    signal_swr_activity_mean_zscore(cells_to_exclude,:) = [];
end

% Select signal for plot
if checkParameter(opts.exp_type, 2, 'default') || checkParameter(opts.exp_type, 3, 'axon')
    signal_swr_activity        = signal_swr_activity_mean;
    signal_swr_activity_zscore = signal_swr_activity_mean_zscore;
end

% Calculate standard error
signal_SE             = std(signal_swr_activity, 'omitnan') ./ sqrt(size(signal_swr_activity,1));
signal_SE_zscore      = std(signal_swr_activity_zscore, 'omitnan') ./ sqrt(size(signal_swr_activity_zscore,1));

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
interval_mean = mean(signal_swr_activity_zscore(:, ( median(1:187)-frameshift:median(1:187)+frameshift)),2);
[max_val,~]   = max(interval_mean,[],2);
[~, sort_idx] = sortrows(max_val);

%% Plot results

% Set labels
label1  = ['Mean ' text, ' DF/F'];
label1z = ['Mean z-scored ' text, ' DF/F'];
x_label = 'Time from SWR peak (sec)';

figure,
sgtitle([ sessionID, ', SWR n = ' num2str(length(RippleIdx)), ', ', celltype], 'Interpreter', 'none') 

x1 = [time(1), time(end)];

if checkParameter(opts.exp_type, 1, 'bulk')

    h1 = subplot(3,2,[1,3]);
    y1 = [1 size(signal_swr_activity,1)];
    imagesc(x1, y1, signal_swr_activity(sort_idx,:)) 
    ylabel(label3)
    xlabel(x_label)
    c(1) = colorbar;
    c(1).Position(1) = 0.47;
%     c(1).Position(2) = 0.39;
    c(1).Position(3) = 0.006;
%         caxis([-0.02 .1])
%         set(h1, 'AlphaData', 1-isnan(signal_swr_activity))

    subplot(3,2,[2,4]);
    imagesc(x1, y1, signal_swr_activity_zscore(sort_idx,:)) 
    ylabel(label3)
    xlabel(x_label)
    c(2) = colorbar;
    c(2).Position(1) = 0.91;
%     c(2).Position(2) = 0.39;
    c(2).Position(3) = 0.006;
    caxis(c_lim)

    subplot(3,2,5)
    shadedErrorBar(time,mean(signal_swr_activity,'omitnan'),  signal_SE);
    font = gca;
    font.FontSize = 14;
    xline(0, '--', 'linew', 1) % Mark SWR peak time
    xlabel(x_label)
    ylabel(label1)
    set(gca, 'xlim',[time(1) time(end)])

    subplot(3,2,6)
    shadedErrorBar(time,mean(signal_swr_activity_zscore, 'omitnan'), signal_SE_zscore);
    font = gca;
    font.FontSize = 14;
    xline(0, '--', 'linew', 1) % Mark SWR peak time
    xlabel(x_label)
    ylabel(label1z)
    set(gca, 'xlim',[time(1) time(end)])

else

    subplot(3,2,[1,3]);
 
    y1 = [1 size(signal_swr_activity,1)];
    imagesc(x1, y1, signal_swr_activity(sort_idx,:))  
    font = gca;
    font.FontSize = 14;
    ylabel(label3, 'FontSize',14)
    set(gca, 'xtick',[])
    c1 = colorbar;
%     caxis([-0.02 .1])
    c1.Position(3) = 0.01;
    c1.Position(1) = 0.468;
    
    subplot(3,2,[2,4]);
    imagesc(x1, y1, signal_swr_activity_zscore(sort_idx,:)) 
    font = gca;
    font.FontSize = 14;
    ylabel(label3, 'FontSize',14)
    set(gca, 'xtick',[])
    c2 = colorbar; 
    caxis(c_lim)
    c2.Position(3) = 0.01;
    c2.Position(1) = 0.91;
    
    subplot(3,2,5)
    shadedErrorBar(time,mean(signal_swr_activity, 'omitnan'),  signal_SE);
    font = gca;
    font.FontSize = 14;
    xline(0, '--', 'linew', 1) % Mark SWR peak time
    xlabel(x_label,  'FontSize',14)
    ylabel(label1)
    set(gca, 'xlim',[time(1) time(end)])
    
    subplot(3,2,6)
    shadedErrorBar(time,mean(signal_swr_activity_zscore,'omitnan'), signal_SE_zscore);
    font = gca;
    font.FontSize = 14;
    xline(0, '--', 'linew', 1) % Mark SWR peak time
    xlabel(x_label,  'FontSize',14)
    ylabel(label1z)
    set(gca, 'xlim',[time(1) time(end)])
end
c1.Position(3) = 0.01;
c1.Position(1) = 0.468;