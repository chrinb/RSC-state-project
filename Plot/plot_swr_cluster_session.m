function plot_swr_cluster_session(varargin)

% Written by Christoffer Berge | Vervaeke Lab

% Function that plots a histogram for each mouse of the nr of 
% sleep sessions that contained different nr of SWR clusters. Use the data
% in "all_data_sleep" cell array from the "all_swr_features" dataset. Second
% argument specifies sleep (1) or awake (2).

data    = varargin{1,1};
type_rec = varargin{1,2};
if type_rec == 1 
    cell_to_choose = 9;
elseif type_rec == 2
    cell_to_choose = 6;
end


m6120_dat = cell2mat(data(2,cell_to_choose)); 
m6121_dat = cell2mat(data(3,cell_to_choose));
m6122_dat = cell2mat(data(4,cell_to_choose));
m6123_dat = cell2mat(data(5,cell_to_choose));

cat = {'Doublets','Triplets','Quartets', 'Quintets', 'Sextets' };

binE = [1.5000 2.5000 3.5000 4.5000 5.5000 6.5000];

figure,
sgtitle('Nr of sessions containing different SWR clusters')
ylabel_text = '# of sessions';

subplot(221)
histogram(m6120_dat, binE)
xticklabels(cat)
ylabel(ylabel_text)
set(gca,'ylim',[0 15])
title(['m6120 (n sessions = ' num2str(size(m6120_dat,1)) ')' ] )

subplot(222)
histogram(m6121_dat, binE)
xticklabels(cat)
ylabel(ylabel_text)
set(gca,'ylim',[0 15])
title(['m6121 (n sessions = ' num2str(size(m6121_dat,1)) ')' ] )

subplot(223)
histogram(m6122_dat, binE)
xticklabels(cat)
ylabel(ylabel_text)
set(gca,'ylim',[0 15])
title(['m6122 (n sessions = ' num2str(size(m6122_dat,1)) ')' ] )

subplot(224)
histogram(m6123_dat, binE)
xticklabels(cat)
ylabel(ylabel_text) 
set(gca,'ylim',[0 15])
title(['m6123 (n sessions = ' num2str(size(m6123_dat,1)) ')' ] )