function [final_spindleLFP,final_spindleLocs, final_spindle_snip_loc] = inspectSpindles(spindleSnips,spindleCentralTimePoint,lfp,spindleHalfLength, tempListSpindles, makeSpectrogramPlot)

%use this code to manually go through each spindle in the LFP and approve or
%delete the spindle

nspindles_start = length(spindleSnips);
outOfAllSpindles = zeros(1,nspindles_start);
t = 1;
msecPerSample = 1000/2500;
extraSpindles = [];
i = 1;
while i <= nspindles_start
    if length(spindleSnips(i).lfp) < 20000
       lfpSnip = spindleSnips(i).lfp;
       snipTime = length(spindleSnips(i).lfp)*msecPerSample;
       time = linspace(0, snipTime/2500 ,length(spindleSnips(i).lfp));
       
    else 
        lfpSnip = spindleSnips(i).lfp;
        time = linspace(0,8,20001);
  
    end
    
    % calculate approximate spindle start/end for visualization
    spindleStart(i) = (median(time)-spindleHalfLength(i))/2500;
    spindleEnd(i) = (median(time)+spindleHalfLength(i))/2500;
    
    % spectrogram needs some work...
    if makeSpectrogramPlot
        figure, 
        subplot(211),
        plot(time, lfpSnip), 
        xlabel('Time (s)');
        xline( median(time) + spindleStart(i),'r--', 'LineWidth',1  ),
        xline( median(time) + spindleEnd(i),'r--', 'LineWidth',1 );

        box off
        set(gca,'TickDir','out')
        xlim([0 max(time)])
        ylim([-1 1])

        window  = 2500;              % Window size for computing the spectrogram (FFT) [# samples]
        overlap = 2400;              % Overlap of the windows for computing the spectrogram [# samples]
        nFFT    = 5:.1:20;          % Vector defining the frequencies for computing the FFT
        fs = 2500;
        x = lfpSnip;
        [s,F,T,P] = spectrogram(x,window,overlap,nFFT,fs);

        subplot(212),
        surf(T,F, (abs(P)),'edgecolor','none');
        colormap(jet)
        view([0 90])
        
        %frex = linspace(8, 20,30);
        %waves = 2*(linspace(3, 10, length(frex))./(2*pi*frex)).^2;
        %wavet = -2:1/2500:2;
        %halfW = floor(length(wavet)/2)+1;
        %nConv = length(lfpSnip) + length(wavet) - 1;
        %tf = zeros(length(frex), length(lfpSnip));
        %dataX = fft(lfpSnip, nConv);
        
        %for fi = 1:length(frex)
        %    waveX = fft( exp(2*1i*pi*frex(fi)*wavet).*exp(-wavet.^2/waves(fi)), nConv);
        %    waveX = waveX./max(waveX);
        %    as = ifft( waveX .* dataX');
        %    %as = reshape(as(halfW:end-halfW+1), [frex, lfpSnip]);
        %    tf(fi,:) = mean(abs(as),2);
        %end
        %subplot(212),
        %contourf(time, frex, tf, 30, 'linecolor', 'none')
    else
        figure,
        plot(time, lfpSnip),
        xline( median(time) + spindleStart(i) ); xline( median(time)+spindleEnd(i) );

        box off
        set(gca,'TickDir','out')
        xlim([1 max(time)])
    end
%display prompt and wait for response
    prompt = sprintf('Keep spindle? %d of %d',i,nspindles_start);
    x = input(prompt,'s');
    
    if strcmp(x,'y')
        keepSpindles(t) = i;
        outOfAllSpindles(i) = 1;
        t = t + 1;
        
        clear lfpSnip
        close
        i = i + 1;
    elseif strcmp(x,'b')
        close
        i = i - 1;
        
    elseif strcmp(x,'m')
        [manualLocs,~] = ginput;
        %calculate distance to center of plot, which is at x = 751 samples
        for k = 1:length(manualLocs)
            temp = find(round(time)==round(manualLocs(k)));
            manualLocs(k) = temp(1);
            clear temp
        end
        manualLocs = manualLocs - 10000;
        manualLocs = spindleIdx(i) + manualLocs;
        extraSpindles = [extraSpindles; manualLocs];
        clear lfpSnip
        close
    else
        clear lfpSnip
        close
        i = i + 1;
    end
end

final_spindle_snip_loc = tempListSpindles(keepSpindles,:);
final_spindleLocs = spindleCentralTimePoint(keepSpindles);
final_spindleLocs = [final_spindleLocs extraSpindles'];
final_spindleLocs = sort(final_spindleLocs);
for i = 1:length(final_spindleLocs)
    spinWindow = final_spindleLocs(i)-10000:final_spindleLocs(i)+10000;
    % delete spindle window indicies occurring prior or after length of
    % recording
    if spinWindow(1,1) < 1
        spinWindow(spinWindow < 1) = [];
    end
    
    if spinWindow(1,end) > length(lfp)
        spinWindow(spinWindow > length(lfp)) = [];
    end
    
    % zero pad spindle snippets > 20001 samples
    final_spindleLFP(i).lfp = lfp(round(spinWindow));
    if length(final_spindleLFP(i).lfp) < 20001
        final_spindleLFP(i).lfp = [final_spindleLFP(i).lfp', zeros(1, 20001-length(final_spindleLFP(i).lfp))];
    end
    
end