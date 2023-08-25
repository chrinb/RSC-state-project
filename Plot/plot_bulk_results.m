function plot_bulk_results(varargin)

% Written by Christoffer Berge || Vervaeke lab

% Plot colormaps and mean/median e-phys event-aligned 2P bulk imaging data. 


label  = 'DF/F';
label2 = 'z-score ';

prompt = sprintf('SWR (1), sleep spindle (2), slow waves (3)? ');
signalSelect = input(prompt);

if signalSelect == 1
    text = 'Time from SWR peak (sec)';
    ylabel_name = ('# SWR');
elseif signalSelect == 2
    text = 'Time from random time (sec)';
    ylabel_name = ('Random #');
elseif signalSelect == 3
    text = 'Time from slow wave trough (sec)';
    ylabel_name = ('# SO/delta wave');
end

% prompt = sprintf('Bulk imaging (1) or single ROIs (2) ? ');
% exp_type = input(prompt);
% if exp_type == 2
%     ylabel_name = ('# ROI');
% end

% create new function for this

% prompt = sprintf('Plot over each other? (1)');
% plot_ontop = input(prompt);
% if plot_ontop == 1
%     prompt = sprintf('Label 1 = ');
%     label1 = input(prompt,  's');
%     prompt = sprintf('Label 2 = ');
%     label2 = input(prompt,  's');
%     
%     prompt = sprintf('z-score (1) or not (2)? ');
%     zscore_label = input(prompt);
%     if zscore_label == 1
%         text4 = 'z-score ';
%     else
%         text4 = [];
%     end
% end

% 
% if plot_ontop == 1
%     signal1     = varargin{1,1};
%     signal2     = varargin{1,2};
%     SE_signal1  = std(rmmissing(signal1))./sqrt(numel(rmmissing(signal1(:, 1))));
%     SE_signal2  = std(rmmissing(signal2))./sqrt(numel(rmmissing(signal2(:, 1))));
%     nrOfSeconds = round( (size(signal1,2) /31) /2 );
%     time        = (-(31*nrOfSeconds):(31*nrOfSeconds))./31;
% 
%     figure,
%     shadedErrorBar(time, nanmean(signal1), SE_signal1,'lineprops', 'b');hold on
%     shadedErrorBar(time, nanmean(signal2), SE_signal2,'lineprops' , 'r');
%     xline(0, '--', 'linew',1)
%     set(gca, 'xlim',[time(1) time(end)])
%     ylabel(['Mean ' text4 text3])
%     legend(label1, label2)
%     xlabel(text)
% end


data        = varargin{1,1};
data_zscore = varargin{1,2};
SE_data     = std(data, 'omitnan')./ sqrt(size(data,1));

SE_zscore_signal = std(data_zscore, 'omitnan') ./ sqrt(size(data_zscore,1));
frame_rate    = 31;
nr_of_seconds = round( (size(data,2) /frame_rate) /2 );
time          = (-(frame_rate*nr_of_seconds):(frame_rate*nr_of_seconds))./frame_rate;

%% Sort according to peak value
% for i = 1:size(data_zscore,1)
% [~, valueId] = max(data_zscore(i,:));
% distance(i) = valueId- size(data_zscore,2)/2;
% end
% test = [distance; 1:length(distance)]';
% test2 = sortrows(test,1);
%% Plot results
% Set up X and Y limits of imagesc function
x1 = [time(1) time(end)];
y1 = [1 size(data,1)]; % Y limits for data 1
y2 = [1, size(data_zscore,1)]; % Y limits for data 2

% Plot data
figure,
subplot(3,2,[1,3]),
imagesc(x1, y1, data)
ylabel(ylabel_name)
colorbar
% title(['Mean ' text3]) 

subplot(3,2,[2,4]),
imagesc(x1, y2, data_zscore)
ylabel(ylabel_name)
colormap(jet)
colorbar
caxis([-3 3])

% title(['Mean ' text4 text3]) 

subplot(3,2,5),
hold on
shadedErrorBar(time, nanmean(data), SE_data,'lineprops', 'b');
plot(time, nanmedian(data), 'linew',1)
set(gca, 'xlim',[time(1) time(end)])
ylabel('Mean DF/F')
xlabel(text)
legend('Mean', 'Median')

subplot(3,2,6),
hold on
shadedErrorBar(time, nanmean(data_zscore), SE_zscore_signal,'lineprops', 'b');
plot(time, nanmedian(data_zscore), 'linew',1)
set(gca, 'xlim',[time(1) time(end)])
ylabel('Mean z-score DF/F')
legend('Mean', 'Median')
xlabel(text)
