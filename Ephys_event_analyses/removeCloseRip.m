function [adjusted_awakeRippleIdx, adjusted_NREMspindleUncoupledSWRidx, adjusted_NREMspindleCoupledSWRidx] = removeCloseRip(sData)

%Written by Christoffer Berge | Vervaeke Lab

%Function that removes successive ripples with an inter-ripple interval
%smaller than a user-defined threshold. Modified version of "RemoveRip"
%function that does the same for different type of SWRs.

%% Adjust interval between ripples

%prompt = sprintf('How many seconds between successive ripples? ');
%secBetweenRipples = input(prompt);
secBetweenRipples = 2;
samplesBetweenRipples = secBetweenRipples*2500;

adjusted_awakeRippleIdx = sData.ephysdata.awake_swr; 
adjusted_NREMspindleUncoupledSWRidx = sData.ephysdata.NREM_spindle_uncoupled_swr;
adjusted_NREMspindleCoupledSWRidx = sData.ephysdata.spindle_coupled_swr;
    


Y = diff(adjusted_awakeRippleIdx); %find duration between ripples (measured in nr of frames)
riplocs = Y < samplesBetweenRipples; 
adjustedRiplocs = [false(1,1), riplocs]; %adds a logical 0 to make the vector equal in length to adjustedrippleIdx
adjusted_awakeRippleIdx(adjustedRiplocs) = [];

X = diff(adjusted_NREMspindleUncoupledSWRidx); %find duration between ripples (measured in nr of frames)
riplocs2 = X < samplesBetweenRipples; 
adjustedRiplocs2 = [false(1,1), riplocs2]; %adds a logical 0 to make the vector equal in length to adjustedrippleIdx
adjusted_NREMspindleUncoupledSWRidx(adjustedRiplocs2) = [];

Z = diff(adjusted_NREMspindleCoupledSWRidx); %find duration between ripples (measured in nr of frames)
riplocs3 = Z < samplesBetweenRipples; 
adjustedRiplocs3 = [false(1,1), riplocs3]; %adds a logical 0 to make the vector equal in length to adjustedrippleIdx
adjusted_NREMspindleCoupledSWRidx(adjustedRiplocs3) = [];





