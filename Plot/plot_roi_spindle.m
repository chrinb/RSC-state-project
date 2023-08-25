function plot_roi_spindle(sData, x)

% Written by Christoffer Berge | Vervaeke Lab

nr_of_seconds = 3;

deconv_sig = sData.imdata.roiSignals(2).ciaDeconvolved(x,:);
dff_sig    = sData.imdata.roiSignals(2).newdff(x,:);
frames     = sData.daqdata.frame_onset_reference_frame;
spindleIdxStart = frames(round(sData.ephysdata2.NREMspindleStartEnd(:,1)));
win_length = (nr_of_seconds*31)*2+1;
time = linspace(-nr_of_seconds,nr_of_seconds,win_length);

for spindle_nr = 1:length(spindleIdxStart)
    
    spindle_window_start = spindleIdxStart(spindle_nr) - (nr_of_seconds*31); 
    spindle_window_end = spindleIdxStart(spindle_nr) + (nr_of_seconds*31); 

    % skip spindles at the beginning or end with a window shorter than
    % length specified in seconds by user above
    if length(spindle_window_start:spindle_window_end) == win_length
        dff_spindle_activity(spindle_nr, :) = ...
            dff_sig(spindle_window_start:spindle_window_end); 
        deconv_spindle_activity(spindle_nr, :) = ...
            deconv_sig(spindle_window_start:spindle_window_end); 
    end
end

dff_SE    = std(dff_spindle_activity) ./ sqrt(size(dff_spindle_activity,1));
deconv_SE = std(deconv_spindle_activity) ./ sqrt(size(deconv_spindle_activity,1));
x1 = [time(1) time(end)];
y1 = [1 size(spindleIdxStart,1)];

figure, 
subplot(221)
imagesc(x1, y1, dff_spindle_activity)
xlabel('Time from ripple peak (sec)')
ylabel('Spindle #')
% colorbar
title('dF/F')


subplot(223)
shadedErrorBar(time,nanmean(dff_spindle_activity),...
    dff_SE,'lineprops', 'b');
set(gca, 'xlim', [time(1), time(end)])
ylabel('Mean dF/F')

subplot(222)
imagesc(x1, y1, deconv_spindle_activity)
xlabel('Time from spindle center (sec)')
ylabel('Spindle #')
% colorbar
title('deconvolved dF/F')

subplot(224)
shadedErrorBar(time,nanmean(deconv_spindle_activity),...
    deconv_SE,'lineprops', 'b');
set(gca, 'xlim', [time(1), time(end)])
ylabel('Mean deconvolved dF/F')
