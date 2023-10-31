function sort_idx = plot_mean_swr(varargin)

% Written by Christoffer Berge | Vervaeke Lab

% Function that plots results from the function "swr_awake_an". The
% variables specificed below must be inputted for the function to return
% plots. 

checkParameter = @(param, n, str) (isnumeric(param) && param==n) || strcmp(param, str);

% Misc variables
sessionID                        = varargin{1,1};
params                           = varargin{1,2};
time                             = varargin{1,3};
signal_swr_activity_mean         = varargin{1,4};
RippleIdx                        = varargin{1,5};

% Calculate standard error
signal_SE             = std(signal_swr_activity_mean, 'omitnan') ./ sqrt(size(signal_swr_activity_mean,1));

if strcmp(params.signal_type, 'deconv')
    c_lim        = [-.05 .5];
elseif strcmp(params.signal_type,  'Dff')
    c_lim        = [-.5 .5];
end

% Sort ROIs according to mean z-scored activity in the -0.5 to + 0.5 interval
% surrounding SWR peak
frameshift    = round(.5/(1/31));
interval_mean = mean(signal_swr_activity_mean(:, ( median(1:187)-frameshift:median(1:187)+frameshift)),2);
[max_val,~]   = max(interval_mean,[],2);
[~, sort_idx] = sortrows(max_val);

%% Plot results
x1 = [time(1), time(end)];
y1 = [1 size(signal_swr_activity_mean,1)];

figure,
sgtitle([ sessionID, ', SWR n = ' num2str(length(RippleIdx)), ', ', params.cell_type], 'Interpreter', 'none') 

h(1) = subplot(3,3,[2,5]);
imagesc(x1, y1, signal_swr_activity_mean(sort_idx,:)) 
ylabel('ROI #', FontSize=16)
set(gca, 'xtick',[])
cbar1          = colorbar;
cbar1.FontSize = 10;
cbar1.Position =[0.6500 0.3900 0.0100 0.4960];
title(['Mean ', params.signal_type], FontSize=12)

h(2) = subplot(3,3,8);
shadedErrorBar(time, mean(signal_swr_activity_mean, 'omitnan'),  signal_SE);
xline(0, '--', 'linew', 1) % Mark SWR peak time
xlabel('Time from SWR peak (s)',  FontSize=16)
set(gca, 'xlim',[time(1) time(end)])

linkaxes(h, 'x');