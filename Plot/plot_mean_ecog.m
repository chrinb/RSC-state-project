% prompt = sprintf('Mouse id? ');
% mouse_select = input(prompt, 's');
% 
% awake_sig   = strcat(mouse_select,' _awakeZ');
% single_sig  = strcat(mouse_select, '_singleZ');
% coupled_sig = strcat(mouse_select, '_coupledZ');
% 
% plot_awake = strcat('y1_', mouse_select, '_awakeZ');
% plot_awake = str
% plot_single = strcat('y1_', mouse_select, '_singleZ');
% plot_coupled = strcat('y1_', mouse_select, '_coupledZ');
% 
% SE_awake   = strcat('SE_', mouse_select, '_awakeZ');
% SE_single  = strcat('SE_', mouse_select, '_singleZ');
% SE_coupled = strcat('SE_', mouse_select, '_coupledZ');

figure,
subplot(231)
imagesc(x1, y1_m6123_awakeZ, m6123_awakeZ)
ylabel('# SWR')

subplot(232)
imagesc(x1, y1_m6123_singleZ, m6123_singleZ)
ylabel('# SWR')

subplot(233)
imagesc(x1,y1_m6123_coupledZ, m6123_coupledZ)
ylabel('# SWR')

subplot(234)
shadedErrorBar(time, mean(m6123_awakeZ,'omitnan'),SE_m6123_awakeZ,'lineprops', 'b');
ylabel('z-score mV')
xlabel('Time from SWR peak (sec)')

subplot(235)
shadedErrorBar(time, mean(m6123_singleZ,'omitnan'),SE_m6123_singleZ,'lineprops', 'k');
ylabel('z-score mV')
xlabel('Time from SWR peak (sec)')
subplot(236)
shadedErrorBar(time, mean(m6123_coupledZ,'omitnan'),SE_m6123_coupledZ,'lineprops', 'r');
ylabel('z-score mV')
xlabel('Time from SWR peak (sec)')