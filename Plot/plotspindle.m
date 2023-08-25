function plotspindle(sData, x)

prompt = sprintf('Choose LFP signal ( 2 = LFP2, 3 = LFP3, 4 = LFP4): ');
signal = input(prompt);
if signal == 2
    aa = sData.ephysdata2.absSpindleIdx;
    bb = sData.ephysdata2.spindleStartEnd;
    cc = sData.ephysdata2.spindleSnips(x).lfp;
elseif signal == 3
    aa = sData.ephysdata3.absSpindleIdx;
    bb = sData.ephysdata3.spindleStartEnd;
    cc = sData.ephysdata3.spindleSnips(x).lfp;

elseif signal == 4
    aa = sData.ephysdata4.absSpindleIdx;
    bb = sData.ephysdata4.spindleStartEnd;
    cc = sData.ephysdata4.spindleSnips(x).lfp;

end

%aa = tempListSpindles(x,:);
%Lia = ismember(sData.ephysdata2.spindleStartEnd, aa);

%testoo = Lia(:,1);
%testu = 1:length(testoo);
%testu(testoo);
time = (0:20001-1)/ 2500;

% find spindle start/end
spinStart = 4 - ( (aa(x)-bb(x,1)) /2500 );
spinEnd = 4 + ( (bb(x,2)-aa(x)) /2500 );


figure, subplot(211),
    plot(time ,cc)
    xline(spinStart); xline(spinEnd); 
    set(gca, 'ylim', [-.5 .5])
    xlabel('Time (s)')
    ylabel('mV')
    title(['Spindle nr ' num2str(x)])

subplot(212), 
    window  = 2500;              
    overlap = 2400;              
    nFFT    = 5:.1:20;         
    fs = 2500;
    x = cc;
        [s,F,T,P] = spectrogram(x,window,overlap,nFFT,fs);

        subplot(212);contourf(T,F, (abs(P)),200,'edgecolor','none');
          colormap(jet)
          view([0 90])
          xlabel('Time (s)')
          ylabel('Frequency (Hz)')
                
end