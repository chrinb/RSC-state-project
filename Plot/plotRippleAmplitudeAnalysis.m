function [sortedAmps] = plotRippleAmplitudeAnalysis(rippleAmpResp,window)

%define field to collect dff response
for i = 1:4
    sortedAmps(i).respDf = [];
end


for i = 1:length(rippleAmpResp)
    amps = rippleAmpResp(i).zScoreAmp;
    %sort all ripple amplitudes, separate into quartiles
    amps(:,2) = 1:size(amps,1);
    amps = sortrows(amps,1); %1 is the smallest amplitude
    quartiles = round(linspace(1,size(amps,1),5));
    quartiles(1) = 0;
    
    %get the population dff response from the set of ripples in each
    %quartile
    for j = 1:length(quartiles)-1
        idx = quartiles(j)+1:quartiles(j+1);
        getResp = amps(idx,2);
        sortedAmps(j).respDf = [sortedAmps(j).respDf; ...
            rippleAmpResp(i).dff(getResp,:)];
        clear idx getResp
    end
end

figure;
time = linspace(-window,window,size(sortedAmps(1).respDf,2));
for i = 1:4
    
     for j = 1:size(sortedAmps(i).respDf,1)
        sortedAmps(i).respDf(j,:) = smooth(zscore(sortedAmps(i).respDf(j,:)));
        
        %subtract baseline (first second of response)
        baselineDf = nanmean(sortedAmps(i).respDf(j,1:31));
        sortedAmps(i).respDf(j,:) = sortedAmps(i).respDf(j,:)-baselineDf;
        clear baselineDf
    end
    hold on
    errorbar(time,smooth(nanmean(sortedAmps(i).respDf)),...
        nansem(sortedAmps(i).respDf,1));
    
end

xlim([-window window])
box off
set(gca,'TickDir','out')
ylabel('zScore DF/F')
xlabel('Time (s)')
