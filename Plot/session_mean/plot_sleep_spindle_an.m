function sort_idx = plot_sleep_spindle_an(varargin)

% Written by Christoffer Berge | Vervaeke Lab

% Function that plots results from the function "sleep spindle_an". The
% variables specificed below must be inputted for the function to return
% plots. 

checkParameter = @(param, n, str) (isnumeric(param) && param==n) || strcmp(param, str);

% Misc variables
sessionID        = varargin{1,1};
opts             = varargin{1,2};
time             = varargin{1,3};
celltype         = varargin{1,8};
spindle_idx      = varargin{1,9};
text             = varargin{1,10};
label3           = varargin{1,11};
cells_to_exclude = varargin{1,12};

signal_spindle_activity        = varargin{1,4};
signal_spindle_activity_zscore = varargin{1,5};
signal_spindle_activity_mean   = varargin{1,6};
signal_spindle_activity_meanZ  = varargin{1,7};

% Remove cell population not analyzed (e.g, if analyzing PCs remove INs
% from plot ( they have zero values anyway) ).
if ~isempty(cells_to_exclude)
    signal_spindle_activity_mean(cells_to_exclude,:)  = [];
    signal_spindle_activity_meanZ(cells_to_exclude,:) = [];
end

% Select signal for plot
if checkParameter(opts.exp_type, 2, 'default') || checkParameter(opts.exp_type, 3, 'axon')
    signal_spindle_activity        = signal_spindle_activity_mean;
    signal_spindle_activity_zscore = signal_spindle_activity_meanZ;
end

% Calculate standard error
signal_SE             = std(signal_spindle_activity, 'omitnan') ./ sqrt(size(signal_spindle_activity,1));
signal_SE_zscore      = std(signal_spindle_activity_zscore, 'omitnan') ./ sqrt(size(signal_spindle_activity_zscore,1));


% % Bulk data
% signal_SE_spindle_bulk_start        = std(signal_spindle_activity_meanZ,'omitnan') ./ ...
%     sqrt(size(signal_spindle_activity,1));
% signal_SE_spindle_zscore_bulk_start = std(signal_spindle_activity_zscore,'omitnan') ./ ...
%     sqrt(size(signal_spindle_activity_zscore,1));
% % Population data
% signal_SE_spindle_pop_start         = std(signal_spindle_activity_mean,'omitnan') ./ ...
%     sqrt(size(signal_spindle_activity_mean,1));
% signal_SE_spindle_zscore_pop_start  = std(signal_spindle_activity_meanZ,'omitnan') ./ ...
%     sqrt(size(signal_spindle_activity_meanZ,1));


if checkParameter(opts.exp_type, 1, 'bulk')
    c_lim        = [-1.5 1.5];
elseif checkParameter(opts.signal_type, 1, 'deconv')
    c_lim        = [-.05 .5];
elseif checkParameter(opts.signal_type, 2, 'dff')
    c_lim        = [-.5 .5];
end

% Sort ROIs according to mean z-scored activity from 0-0.5 sec post spindle onset
frameshift    = round(31/2);
center        = round(size(time,2)/2);
interval_mean = mean( signal_spindle_activity(:,center:(center+frameshift)),2);
[max_val,~]   = max(interval_mean,[],2);
[~, sort_idx] = sortrows(max_val);

%% Plot results

% Set labels
label1  = ['Mean ' text, ' DF/F'];
label1z = ['Mean z-scored ' text, ' DF/F'];
x_label = 'Time from spindle onset (sec)';

x1 = [time(1) time(end)];
y1 = [1 size(signal_spindle_activity_mean,1) ];

figure,
sgtitle([ sessionID, ', Spindle n = ' num2str(length(spindle_idx)), ', ', celltype], 'Interpreter', 'none') 

% Plot bulk results
% if opts.exp_type == 1

subplot(3,2, [1 3])
imagesc(x1, y1, signal_spindle_activity(sort_idx,: )) 
font = gca;
font.FontSize = 14;
ylabel(label3, 'FontSize',14)
set(gca, 'xtick',[])
% xlabel(label_x)
% title(['Spindle n = ', num2str(spindle_nr)])
c1 = colorbar;
c1.Position(3) = 0.01;
c1.Position(1) = 0.468;

subplot(3,2, [2 4])
imagesc(x1, y1, signal_spindle_activity_zscore(sort_idx,:))  
font = gca;
font.FontSize = 14;
set(gca, 'xtick',[])
ylabel(label1, FontSize=14)
% xlabel(x_label)
c2 = colorbar;
c2.Position(3) = 0.01;
c2.Position(1) = 0.91;
caxis(c_lim)

subplot(3,2, 5)
shadedErrorBar(time, mean(signal_spindle_activity, 'omitnan'),...
    signal_SE,'lineprops', 'r');
xlabel(x_label, FontSize=14)
ylabel(label1, FontSize=14)
set(gca, 'xlim',[min(time) max(time)])
    
subplot(3,2,6)
shadedErrorBar(time, mean(signal_spindle_activity_zscore,'omitnan'),...
    signal_SE_zscore,'lineprops', 'b');
xlabel(x_label, FontSize=14)
ylabel(label1, FontSize=14)
set(gca, 'xlim',[min(time) max(time)])

% Plot population results
% else
%     
%     subplot(3,2, [1 3])
%     imagesc(x1, y1, signal_spindle_activity_mean(sort_idx,:)) 
%     font = gca;
%     font.FontSize = 14;
%     ylabel(label3, 'FontSize',14)
%     set(gca, 'xtick',[])
%     c1 = colorbar;
%     c1.Position(3) = 0.01;
%     c1.Position(1) = 0.468;
%     
%     subplot(3,2, 5)
%     shadedErrorBar(time, mean(signal_spindle_activity_mean,'omitnan'),...
%         signal_SE_spindle_pop_start,'lineprops', 'r');
%     font = gca;
%     font.FontSize = 14;
%     xlabel(x_label)
%     ylabel(label1)
%     set(gca, 'xlim',[min(time) max(time)])
% 
%     subplot(3,2, [2 4])
%     imagesc(x1, y1, signal_spindle_activity_meanZ(sort_idx,:)) %  colorplot of average individual roi activity during ripples
%     font = gca;
%     font.FontSize = 14;
%     ylabel(label3, 'FontSize',14)
%     set(gca, 'xtick',[])
% %     ylabel(label3)
% %     xlabel(x_label)
%     c2 = colorbar;
%     c2.Position(3) = 0.01;
%     c2.Position(1) = 0.91;
%     caxis(c_lim)
%     
%     subplot(3,2,6)
%     shadedErrorBar(time, mean(signal_spindle_activity_meanZ,'omitnan'),...
%         signal_SE_spindle_zscore_pop_start,'lineprops', 'b');
%     font = gca;
%     font.FontSize = 14;
%     xlabel(x_label)
%     ylabel(label1z)
%     set(gca, 'xlim',[min(time) max(time)])
% end
% 
%     
% 
% 
% 
% 
