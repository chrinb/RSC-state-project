function [NoRunRipples,RunSpeedDuringRipples] = RipRun(sData)

%Written by Christoffer Berge | Vervaeke Lab

% Function to extract only those ripples where there is no movement of the
% running wheel at the ripple peak and 2s before and after ripple peak

% 125 corresponds to number of frames, 1 for ripple peak and 62 frames
% before and after

RunSpeedDuringRipples = zeros(length(sData.ephysdata.absRipIdx),10001); 

RippleIdx = sData.ephysdata.absRipIdx; 

v = sData.daqdata.runSpeed;
for j = 1:length(RippleIdx) % creates vector of nr of ripples during recording
        rippleNr = j;
        x = (RippleIdx(rippleNr) - (2*2500)); % 2s multiplied by LFP sampling frequency in Hz
        if x < 0
            x = 1; %to avoid negative x values
        end
        y = (RippleIdx(rippleNr) + (2*2500)); % 2s multiplied by LFP sampling frequency in Hz
        if y > length(sData.daqdata.runSpeed)
            y = length(sData.daqdata.runSpeed); %to avoid y values above length of recording
        end
        
        if length(x:y) == (10001) % total number of LFP sampling frequency 'frames' 
            RunSpeedDuringRipples(rippleNr, :) = v(x:y); %index every individual ripple frame segment into a new row for that particular roi(i) )
        elseif length(x:y) < 10001 
            % do nothing
        end
        
        if mean(RunSpeedDuringRipples(rippleNr, :)) > 0.01 % threshold for running is set to 0, try changing this if to conservative (e.g. too few ripples)
            % do nothing
        else
            AllRipples(rippleNr) = rippleNr;
        end
end

NoRunRipples = sData.ephysdata.frameRipIdx(AllRipples > 0);


