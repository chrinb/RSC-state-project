function SNR  = calculateSNR(sData)

     % peakSNR=Δpeak/ 2*σn
     % Δpeak is the difference between the biggest spike value and the baseline value, 
     % and σn is the noise SD calculated from nonactive intervals of traces.    ​
     dff = sData.imdata.roiSignals(2).newdff;
     [noROIs,noSample] = size(dff);
    
    % baseline subtraction for each ROI:
    prctl_baseline = 10; 
    DFF = NaN(noROIs,noSample);
    for i = 1:noROIs
        baseline = prctile(dff(i,:),prctl_baseline);
        if baseline < 0
           DFF(i,:) = (dff(i,:) - baseline)/(-1*baseline); 
        else
           DFF(i,:) = (dff(i,:) - baseline)/baseline;
        end
    end
    
    % Low-pass filter fluorescence signals to find noise and actual signal
    
    DFF_lowpass = NaN(noROIs,noSample);
    
    % paramteres of filter: 
    % Fp = Pass frequency (below this freq signal remains), 
    % Fst = stop frequency (above this freq signal discarded), between Fp and FST  
    % there is a transition zone, effect depends on filter type. Ap = Amplitude pass, 
    % Number of ripples to pass (don't understand what is it for), Ast = Amplitude stop (attenuation in decibel)
    d = fdesign.lowpass('Fp,Fst,Ap,Ast',0.08,0.9,1,60); % Potential changes for stronger filtering: decrease Fp (0.01) and Fst (0.08). original stting: 0.1,1,1,60
    Hd = design(d,'cheby1');   % or use 'equiripple' filter to have less amplitude filtering and more noise
    for i = 1:noROIs
        DFF_lowpass(i,:) = filter(Hd,DFF(i,:)); 
    end
    noise = DFF - DFF_lowpass;  % subtract the lowpassed data from baseline-corrected signal to get the noise (PreNoise)
    
    % calculate SNR 
    noiseStd = NaN(noROIs,1); % whole session std for each ROI
    ROIPeak  = NaN(noROIs,1); % the highest peaks in that ROI 
    SNR      = NaN(noROIs,1); % signal to noise ratio for each ROI
    for i = 1:noROIs
        noiseStd(i) = 2*std(noise(i,:),[],2); % can be adjusted to 1*std(PreNoise(i,:),[],2);
        ROIPeak(i) = prctile(DFF(i,:),99); % can be adjusted to 95%
        SNR(i) = ROIPeak(i) / noiseStd(i);
    end
    
end