function plotNREMspindles(sData, varargin)

% plots NREM spindles in the different channels. If sData input is followed
% by 1, spindles that overlap temporally among channels are marked.

kk = cell2mat(varargin);
if kk == 1
   [~,~,~,chan2overlapSpin,chan2nonoverlapSpin,chan3overlapSpin,chan3nonoverlapSpin, ...
   chan4overlapSpin,chan4nonoverlapSpin] = spindleoverlap(sData);

RSCoverlapspindles = chan2overlapSpin/2500;
RSCsolitaryspindles = chan2nonoverlapSpin/2500;
PFCoverlapspindles = chan3overlapSpin/2500;
PFCsolitaryspindles = chan3nonoverlapSpin/2500;
S1overlapspindles = chan4overlapSpin/2500;
S1solitaryspindles = chan4nonoverlapSpin/2500;

end


NREMstartend = sData.NREMepisodes./2500;
time = (0:length(sData.ephysdata.lfp)-1) / 2500;

RSCspin = sData.ephysdata2.NREMspindleStartEnd/2500;
PFCspin = sData.ephysdata3.NREMspindleStartEnd/2500;
S1spin = sData.ephysdata4.NREMspindleStartEnd/2500;

figure, subplot(311)
plot(time, sData.ephysdata2.lfp), set(gca, 'ylim',[-.6 .5], 'xlim', [0 max(time)]), 
xlabel('Time (s)')
ylabel('mV')
title('RSC')
hold on

for i = 1:length(NREMstartend)
    x = [NREMstartend(i,1) NREMstartend(i,1) NREMstartend(i,2) NREMstartend(i,2)];
    y = [-.6 .6 .6 -.6];
    patch(x, y, 'green', 'edgecolor', 'none', 'FaceAlpha', .3);
end

if kk == 1
    for k = 1:length(RSCsolitaryspindles)
        z = [RSCsolitaryspindles(k,1) RSCsolitaryspindles(k,1) RSCsolitaryspindles(k,2) RSCsolitaryspindles(k,2) ];
        v = [-.6 .6 .6 -.6];
        patch(z,v, 'magenta', 'edgecolor', 'none', 'FaceAlpha', .5);
    end
    for p = 1:length(RSCoverlapspindles)
        g = [RSCoverlapspindles(p,1) RSCoverlapspindles(p,1) RSCoverlapspindles(p,2) RSCoverlapspindles(p,2) ];
        d = [-.6 .6 .6 -.6];
        patch(g,d, 'red', 'edgecolor', 'none', 'FaceAlpha', .5);
    end
else
    for k = 1:length(sData.ephysdata2.NREMspindleStartEnd)
    z = [RSCspin(k,1) RSCspin(k,1) RSCspin(k,2) RSCspin(k,2) ];
    v = [-.6 .6 .6 -.6];
    patch(z,v, 'magenta', 'edgecolor', 'none', 'FaceAlpha', .5);
    end
end

hold off

subplot(312)
plot(time,sData.ephysdata3.lfp), set(gca, 'ylim',[-.6 .5], 'xlim', [0 max(time)]), 
xlabel('Time (s)')
ylabel('mV')
title('PFC')
hold on

for i = 1:length(NREMstartend)
    x = [NREMstartend(i,1) NREMstartend(i,1) NREMstartend(i,2) NREMstartend(i,2)];
    y = [-.6 .6 .6 -.6];
    patch(x, y, 'green', 'edgecolor', 'none', 'FaceAlpha', .3)
end

if kk == 1
    for k = 1:length(PFCsolitaryspindles)
        z = [PFCsolitaryspindles(k,1) PFCsolitaryspindles(k,1) PFCsolitaryspindles(k,2) PFCsolitaryspindles(k,2) ];
        v = [-.6 .6 .6 -.6];
        patch(z,v, 'magenta', 'edgecolor', 'none', 'FaceAlpha', .5);
    end
    for p = 1:length(PFCoverlapspindles)
        g = [PFCoverlapspindles(p,1) PFCoverlapspindles(p,1) PFCoverlapspindles(p,2) PFCoverlapspindles(p,2) ];
        d = [-.6 .6 .6 -.6];
        patch(g,d, 'red', 'edgecolor', 'none', 'FaceAlpha', .5);
    end
else
    for k = 1:length(sData.ephysdata3.NREMspindleStartEnd)
    z = [PFCspin(k,1) PFCspin(k,1) PFCspin(k,2) PFCspin(k,2) ];
    v = [-.6 .6 .6 -.6];
    patch(z,v, 'magenta', 'edgecolor', 'none', 'FaceAlpha', .5);
    end
end
hold off, 

subplot(313)
plot(time,sData.ephysdata4.lfp), set(gca, 'ylim',[-.6 .5], 'xlim', [0 max(time)]), 
xlabel('Time (s)')
ylabel('mV')
title('S1')
hold on

for i = 1:length(NREMstartend)
    x = [NREMstartend(i,1) NREMstartend(i,1) NREMstartend(i,2) NREMstartend(i,2)];
    y = [-.6 .6 .6 -.6];
    patch(x, y, 'green', 'edgecolor', 'none', 'FaceAlpha', .3)
end

if kk == 1
    for k = 1:length(S1solitaryspindles)
        z = [S1solitaryspindles(k,1) S1solitaryspindles(k,1) S1solitaryspindles(k,2) S1solitaryspindles(k,2) ];
        v = [-.6 .6 .6 -.6];
        patch(z,v, 'magenta', 'edgecolor', 'none', 'FaceAlpha', .5);
    end
    for p = 1:length(S1overlapspindles)
        g = [S1overlapspindles(p,1) S1overlapspindles(p,1) S1overlapspindles(p,2) S1overlapspindles(p,2) ];
        d = [-.6 .6 .6 -.6];
        patch(g,d, 'red', 'edgecolor', 'none', 'FaceAlpha', .5);
    end
else
    for k = 1:length(sData.ephysdata4.NREMspindleStartEnd)
    z = [S1spin(k,1) S1spin(k,1) S1spin(k,2) S1spin(k,2) ];
    v = [-.6 .6 .6 -.6];
    patch(z,v, 'magenta', 'edgecolor', 'none', 'FaceAlpha', .5);
    end
end