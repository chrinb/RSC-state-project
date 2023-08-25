function vis_spindle_resp(sData)


% Function that allows user to visualize the activity of ROI's across
% individual sleep spindles

nr_of_seconds = 3;
win_length    = (nr_of_seconds*31)*2+1;
time          = linspace(-nr_of_seconds,nr_of_seconds,win_length);
frames        = sData.daqdata.frame_onset_reference_frame;
%% Select 2P signal
signal_dff         = sData.imdata.roiSignals(2).newdff;
signal_deconvolved = sData.imdata.roiSignals(2).ciaDeconvolved;

%% Select spindle freq
prompt = sprintf('8-18Hz (1) or 10-16hz (2)? ');
spindle_band_select = input(prompt);

if spindle_band_select == 1
    spindle_select = [];
elseif spindle_band_select == 2
    spindle_select = '1016';
end

spin_start_end_str = strcat('spindleStartEnd', spindle_select);

% Convert e-phys SWR times to 2P frame rate times
spindleIdxStart = frames(round(sData.ephysdata2.(spin_start_end_str)(:,1)));


t = 1;
x1 = [time(1) time(end)];
y1 = [1 size(signal_dff,1) ];
spindle_nr = 1;
while spindle_nr <= length(spindleIdxStart)
    spindle_window_start = spindleIdxStart(spindle_nr) - (nr_of_seconds*31); 
    spindle_window_end   = spindleIdxStart(spindle_nr) + (nr_of_seconds*31); 
    if spindle_window_start > 1 && spindle_window_end < length(signal_dff)
        figure,
        subplot(211)
        imagesc(x1, y1, signal_dff(:,spindle_window_start:spindle_window_end));
        colormap jet
        
        subplot(212)
        imagesc(x1, y1, signal_deconvolved(:,spindle_window_start:spindle_window_end));
        
        txt = ['Spindle #', num2str(spindle_nr)];
        sgtitle(txt)
        
        prompt = sprintf('Store spindle index?');
        x = input(prompt,'s');
    
        if strcmp(x,'y')
            interesting_spindle_idx(t) = spindle_nr;
            t = t+1;
            close
            spindle_nr = spindle_nr + 1;
        elseif strcmp(x,'b')
            close
            spindle_nr = spindle_nr - 1;
        else
            close
            spindle_nr = spindle_nr + 1;
        end
    else
        spindle_nr = spindle_nr + 1;
    end
end

interesting_spindle_idx;
    