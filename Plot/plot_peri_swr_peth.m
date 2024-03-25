function plot_peri_swr_peth(pre_data, post_data, pre_data_SE, post_data_SE)

%{
Plot ECoG/LFP peri-SWR activity as SWR x time and the mean + SEM, 
split into pre and post RUN or awake vs sleep
%}


 %% Plot
time_vec = -1:1/2500:1;

figure(111), clf
% subplot(3,2,1:4),
subplot(2,2,1)
x1       = [time_vec(1) time_vec(end)];
y1 = [1, size(pre_data,1)];
imagesc(x1, y1, pre_data)
colormap inferno
ylabel('SWR #', FontSize=16)
% set(gca,'xticklabel',[]) 
colorbar
caxis([-.5 .5])
title('Pre-RUN')

subplot(2,2,2)
y1 = [1, size(post_data,1)];
imagesc(x1, y1, post_data)
colormap inferno
% set(gca,'xticklabel',[])   
colorbar
caxis([-.5 .5])
title('Post-RUN')
% subplot(3,2, 5:6)
subplot(2,2,3:4)
shadedErrorBar(time_vec, mean(pre_data,'omitnan'),pre_data_SE, 'lineprops', 'r'), hold on
shadedErrorBar(time_vec, mean(post_data,'omitnan'),post_data_SE,'lineprops','b')
xlabel('Time from SWR peak (s)')
ylabel('\muV (z-score)')
set(gca, "FontSize", 12)
xline(0, 'r--', 'LineWidth',1)
legend({'Pre-RUN', 'Post-RUN'})
set(gca, 'xlim', [-.3 .3])