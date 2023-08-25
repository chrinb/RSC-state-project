function [varargout]  = swrAlignedEcog(sData)

% Written by Christoffer Berge | Vervaeke Lab

prompt = sprintf('How many seconds before/after ripple peak? ');
nrOfSeconds = input(prompt);

prompt = sprintf('Awake rec.? ');
awake = input(prompt,'s');

prompt = sprintf('All ripples? (y = yes | everything else = no) ');
allrip = input(prompt,'s');

if ~strcmp(allrip,'y')
    prompt = sprintf('Remove locomotion SWR? (y = yes | everything else = no) ');
    riprun = input(prompt, 's');
    
    prompt = sprintf('Remove temporally close SWR? (y = yes | everything else = no) ');
    removerip = input(prompt, 's');
end

% select SWRs for awake recordings 
if strcmp(awake,'y') && strcmp(allrip,'y') 
    awakeSWRidx = sData.ephysdata.absRipIdx;
    
    % if remove locomotion SWR but not temporally close SWR
    elseif strcmp(awake,'y') && strcmp(riprun, 'y') && ~strcmp(removerip, 'y')
    awakeSWRidx = riprun2(sData);
    
    % if remove temporally close SWR but not locomotion SWR
    elseif strcmp(awake,'y') && strcmp(removerip, 'y') && ~strcmp(riprun, 'y')
    awakeSWRidx = RemoveRip(sData);
    
    % if remove both temporally close and locotion SWR
    elseif strcmp(awake,'y') && strcmp(removerip, 'y') && strcmp(riprun, 'y')
    sData.ephysdata.absRipIdx = riprun2(sData);
    awakeSWRidx = RemoveRip(sData);

%for sleep recordings

    % if keep all SWRs
    elseif ~strcmp(awake,'y') && strcmp(allrip,'y') 
    awakeSWRidx = sData.ephysdata.unclassified_swr;
    NREMspindleUncoupledSWRidx = sData.ephysdata.NREM_spindle_uncoupled_swr;
    NREMspindleCoupledSWRidx = sData.ephysdata.spindle_coupled_swr;

    % if remove locomotion SWR but not temporally close SWR
    elseif ~strcmp(awake,'y') && strcmp(riprun, 'y') && ~strcmp(removerip, 'y')
            [awakeSWRidx, NREMspindleUncoupledSWRidx, NREMspindleCoupledSWRidx ] = ripRunAn(sData); 

    % if remove temporally close SWR but not locomotion SWR
    elseif ~strcmp(awake,'y') && strcmp(removerip, 'y') && ~strcmp(riprun, 'y')
            [awakeSWRidx, NREMspindleUncoupledSWRidx, NREMspindleCoupledSWRidx] = removeCloseRip(sData);
   
    % if remove both temporally close and locotion SWR
    elseif ~strcmp(awake,'y') && strcmp(removerip, 'y') && strcmp(riprun, 'y')
       [awakeSWRidx, NREMspindleUncoupledSWRidx, NREMspindleCoupledSWRidx ] = ripRunAn(sData,1);
end

ECoG = sData.ephysdata2.lfp;

if strcmp(awake,'y')
    % Preallocate
    awakeSWR_ECoG_signals = zeros(length(awakeSWRidx), (nrOfSeconds*2500*2)+1);
else
    % Preallocate 
    awakeSWR_ECoG_signals = zeros(length(awakeSWRidx), (nrOfSeconds*2500*2)+1);
    spindle_uncoupled_ECoG_signals = zeros(length(NREMspindleUncoupledSWRidx), (nrOfSeconds*2500*2)+1);
    spindle_coupled_ECoG_signals = zeros(length(NREMspindleCoupledSWRidx), (nrOfSeconds*2500*2)+1);
end

time = (-(nrOfSeconds*2500):(nrOfSeconds*2500))./2500;

%% awake SWR only 
if strcmp(awake,'y')
    for swrNr = 1:length(awakeSWRidx)
    
    % find beginning and end of ECoG snippet centered on SWR
    ECoG_signal_start = awakeSWRidx(swrNr) - (nrOfSeconds * 2500);
    ECoG_signal_end = awakeSWRidx(swrNr) + (nrOfSeconds * 2500);
    
    % if ECoG snippet is starting before or outlasts the actual ECoG signal
    % from which it is taken (because SWR occurs at the very beginning or
    % end), exclude that ECoG snippet
    if ECoG_signal_start <= 0 || ECoG_signal_end > length(ECoG)
        % do nothing
    else
        awakeSWR_ECoG_signals(swrNr,:) = ECoG(ECoG_signal_start:ECoG_signal_end);
    end
    end



    zscore_awakeSWR_ECoG_signals = zscore(awakeSWR_ECoG_signals,0,2);

    avg_awakeSWR_ECoG_signal = nanmean(awakeSWR_ECoG_signals);
    avg_zscore_awakeSWR_ECoG_signals = nanmean(zscore_awakeSWR_ECoG_signals);
    varargout{1} = avg_awakeSWR_ECoG_signal;
    varargout{2} = avg_zscore_awakeSWR_ECoG_signals;
    varargout{3} = awakeSWR_ECoG_signals;
    
    N1 = length(awakeSWRidx);
    SE_avg_awakeSWR        = std(awakeSWR_ECoG_signals)./sqrt(N1);
    SE_zscore_avg_awakeSWR = std(zscore_awakeSWR_ECoG_signals)./sqrt(N1);
    
figure, 
subplot(211),
shadedErrorBar(time,avg_awakeSWR_ECoG_signal, SE_avg_awakeSWR ,'lineProps', 'r');
title(['Mean RSC ECoG during awake SWRs n = ' num2str(N1)])
xlabel('Time from ripple peak (sec)')

subplot(212)
shadedErrorBar(time, avg_zscore_awakeSWR_ECoG_signals, SE_zscore_avg_awakeSWR,'lineProps', 'r')
title(['zscore mean RSC ECoG awake SWR n = ' num2str(N1)])
xlabel('Time from ripple peak (sec)')
%% awake SWR 
else
    for swrNr = 1:length(awakeSWRidx)
    
    % find beginning and end of ECoG snippet centered on SWR
    ECoG_signal_start = awakeSWRidx(swrNr) - (nrOfSeconds * 2500);
    ECoG_signal_end = awakeSWRidx(swrNr) + (nrOfSeconds * 2500);
    
    % if ECoG snippet is starting before or outlasts the actual ECoG signal
    % from which it is taken (because SWR occurs at the very beginning or
    % end), exclude that ECoG snippet
    if ECoG_signal_start < 0 || ECoG_signal_end > length(ECoG)
        % do nothing
    else
        awakeSWR_ECoG_signals(swrNr,:) = ECoG(ECoG_signal_start:ECoG_signal_end);
    end
    
    end

    zscore_awakeSWR_ECoG_signals = zscore(awakeSWR_ECoG_signals,0,2);

    avg_awakeSWR_ECoG_signal = nanmean(awakeSWR_ECoG_signals);
    avg_zscore_awakeSWR_ECoG_signals = nanmean(zscore_awakeSWR_ECoG_signals);
    varargout{1} = avg_awakeSWR_ECoG_signal;
    varargout{2} = avg_zscore_awakeSWR_ECoG_signals;
    varargout{3} = awakeSWR_ECoG_signals;
%% NREM spindle-uncoupled SWR
for swrNr = 1:length(NREMspindleUncoupledSWRidx)
    
    % find beginning and end of ECoG snippet centered on SWR
    ECoG_signal_start = NREMspindleUncoupledSWRidx(swrNr) - (nrOfSeconds * 2500);
    ECoG_signal_end = NREMspindleUncoupledSWRidx(swrNr) + (nrOfSeconds * 2500);
    
    % if ECoG snippet is starting before or outlasts the actual ECoG signal
    % from which it is taken (because SWR occurs at the very beginning or
    % end), exclude that ECoG snippet
    if ECoG_signal_start < 0 || ECoG_signal_end > length(ECoG)
        % do nothing
    else
        spindle_uncoupled_ECoG_signals(swrNr,:) = ECoG(ECoG_signal_start:ECoG_signal_end);
    end
    
end

zscore_spindle_uncoupled_ECoG_signals = zscore(spindle_uncoupled_ECoG_signals,0,2);

avg_spindle_uncoupled_ECoG_signal = nanmean(spindle_uncoupled_ECoG_signals);
avg_zscore_spindle_uncoupled_ECoG_signals = nanmean(zscore_spindle_uncoupled_ECoG_signals);
    varargout{4} = avg_spindle_uncoupled_ECoG_signal;
    varargout{5} = avg_zscore_spindle_uncoupled_ECoG_signals;
    varargout{6} = spindle_uncoupled_ECoG_signals;
%% NREM spindle-coupled SWR
for swrNr = 1:length(NREMspindleCoupledSWRidx)
    
    % find beginning and end of ECoG snippet centered on SWR
    ECoG_signal_start = NREMspindleCoupledSWRidx(swrNr) - (nrOfSeconds * 2500);
    ECoG_signal_end = NREMspindleCoupledSWRidx(swrNr) + (nrOfSeconds * 2500);
    
    % if ECoG snippet is starting before or outlasts the actual ECoG signal
    % from which it is taken (because SWR occurs at the very beginning or
    % end), exclude that ECoG snippet
    if ECoG_signal_start < 0 || ECoG_signal_end > length(ECoG)
        % do nothing
    else
        spindle_coupled_ECoG_signals(swrNr,:) = ECoG(ECoG_signal_start:ECoG_signal_end);
    end
    
end

zscore_spindle_coupled_ECoG_signals = zscore(spindle_coupled_ECoG_signals,0, 2);

avg_spindle_coupled_ECoG_signal = nanmean(spindle_coupled_ECoG_signals);
avg_zscore_spindle_coupled_ECoG_signals = nanmean(zscore_spindle_coupled_ECoG_signals);

 varargout{7} = avg_spindle_coupled_ECoG_signal;
 varargout{8} = avg_zscore_spindle_coupled_ECoG_signals;
 varargout{9} = spindle_coupled_ECoG_signals;

%% Create plots
%determine nr of SWR in each category
N1 = length(awakeSWRidx);
N2 = length(NREMspindleUncoupledSWRidx);
N3 = length(NREMspindleCoupledSWRidx);

%calculate SE of each SWR average
SE_avg_awakeSWR        = std(awakeSWR_ECoG_signals)./sqrt(N1);
SE_zscore_avg_awakeSWR = std(zscore_awakeSWR_ECoG_signals)./sqrt(N1);

SE_avg_spindle_uncoupled_ECoG_signal        = std(spindle_uncoupled_ECoG_signals)./sqrt(N2);
SE_zscore_avg_spindle_uncoupled_ECoG_signal = std(zscore_spindle_uncoupled_ECoG_signals)./sqrt(N2);

SE_spindle_coupled_ECoG_signals        = std(spindle_coupled_ECoG_signals)./sqrt(N3);
SE_zscore_spindle_coupled_ECoG_signals = std(zscore_spindle_coupled_ECoG_signals)./sqrt(N3);

figure, 
subplot(321)
shadedErrorBar(time,avg_awakeSWR_ECoG_signal, SE_avg_awakeSWR ,'lineProps', 'r');
title(['Mean RSC ECoG during awake SWRs n = ' num2str(N1)])
xlabel('Time from ripple peak (sec)')

subplot(323)
shadedErrorBar(time, avg_spindle_uncoupled_ECoG_signal,...
    SE_avg_spindle_uncoupled_ECoG_signal,'lineProps', 'b');
title(['Mean RSC ECoG NREM spindle-uncoupled SWR n = ' num2str(N2)])
xlabel('Time from ripple peak (sec)')

subplot(325)
shadedErrorBar(time, avg_spindle_coupled_ECoG_signal, ...
    SE_spindle_coupled_ECoG_signals,'lineProps', 'k')
title(['Mean RSC ECoG during NREM spindle-coupled SWR n = ' num2str(N3)])
xlabel('Time from ripple peak (sec)')

subplot(322)
shadedErrorBar(time, avg_zscore_awakeSWR_ECoG_signals, SE_zscore_avg_awakeSWR,'lineProps', 'r')
title(['zscore mean RSC ECoG awake SWR n = ' num2str(N1)])
xlabel('Time from ripple peak (sec)')

subplot(324)
shadedErrorBar(time, avg_zscore_spindle_uncoupled_ECoG_signals,...
    SE_zscore_avg_spindle_uncoupled_ECoG_signal,'lineProps', 'b');
title(['zscore mean RSC ECoG spindle-uncoupled SWRs n = ' num2str(N2)])
xlabel('Time from ripple peak (sec)')

subplot(326)
shadedErrorBar(time, avg_zscore_spindle_coupled_ECoG_signals, ...
    SE_zscore_spindle_coupled_ECoG_signals,'lineProps', 'k')
title(['zscore mean RSC ECoG NREM spindle-coupled SWR n = ' num2str(N3)])
xlabel('Time from ripple peak (sec)')

end

