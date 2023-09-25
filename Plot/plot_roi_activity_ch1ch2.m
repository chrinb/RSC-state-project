function plot_roi_activity_ch1ch2(roi_nr, mean_ch1, mean_ch2)


figure, 
h(1) = subplot(311);
plot(mean_ch1(roi_nr,:))

h(2) = subplot(312);
plot(mean_ch2(roi_nr,:))

h(3) = subplot(313);
plot(mean_ch2(roi_nr,:)-mean_ch1(roi_nr,:))

sgtitle(['ROI ', num2str(roi_nr)])

linkaxes(h, 'x');