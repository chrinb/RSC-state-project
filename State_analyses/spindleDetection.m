function spindle = spindleDetection(inputSignal, sampRate)
%Spindle detection and characterization using an optimal wavelet function
%   This is an improved spindle detection algorithm using an optimal wavelet
%   function and additional criteria (minimum/maximum length/cycles,
%   thresholding, central frequency, symmetry, peak to peak, ...).
%
%   Mojtaba Bandarabadi, Inselspital, Bern, June 2016
%   m.bandarabadi@gmail.com
%
%   Inputs:
%     inputSignal: matrix of EEG or LFP signal/s (single channel or multi)
%     sampRate: sampling rate of signals
%
%   Output:
%     spindle: a structure containing characterization of individual
%     detected spindles
%

%% default parameters
param.freq = [8 10 14 18]; % frequency band of spindles [9 10 14 16]
param.thrF = [1 3 20];     % threshold factor [Low High Outlier] ([1 3 20])
param.minDur = 0.5;        % min duration in second (0.4 s)
param.maxDur = 3;          % max duration in second (2 s)
param.minCycle = 5;        % min number of cycles (5)
param.maxCycle = 30;       % max number of cycles (30)


%% preprocessing of input signal/s
if size(inputSignal,2) > size(inputSignal,1)
    inputSignal = inputSignal';
end

inputSignal = detrend(inputSignal);                 % remove DC and linear trends from input signal/s
inputSignal = resample(inputSignal, 200, sampRate); % downsample to 200 Hz


%% loop for each channel
for iCh = 1:size(inputSignal,2)
    disp(['Detecting spindles of channel ' ,num2str(iCh)])
    chanSignal = inputSignal(:,iCh);
    
    spindle = struct('start',[],'end',[],'length',[],'centFreq',[],...
        'negPeak',[],'posPeak',[],'peak2peak',[],'noCycle',[],'symmetry',[]); % structure to save the results
    
    
    %% extract CWT coefficients in the spindle range
    waveName = 'fbsp2-1-2';                              % wavelet function
    wavCentFreq = centfrq(waveName);                     % central pseudo-frequency of the B-spline wavelet function
    freqCent = param.freq(2):0.5:param.freq(3);          % pseudo-frequencies of scales to extract CWT coefficients
    scales = wavCentFreq./(freqCent./sampRate);          % scales to extract cwt coefficients
    cwtCoef = cwt(chanSignal, scales, waveName);         % cwt coefficients in the spindle frequency range
    cwtCoef = abs(cwtCoef.^2);                           % cwt power
    cwtPower = freqCent*cwtCoef;                         % 1/f correction
    cwtPower = conv(cwtPower, hann(sampRate/5), 'same'); % smoothing using 200ms Hanning window
    
    % filtering in the spindle range to count number of cycles
    order = round(3*(sampRate/param.freq(1)));
    fir1Coef = fir1(order, [param.freq(1), param.freq(4)]./(sampRate/2));
    filtData = filtfilt(fir1Coef, 1, chanSignal);
    
    
    %% find upper and lower thresholds
    cwtPowerT = cwtPower(cwtPower<10*nanstd(cwtPower));                % removing outliers for a proper estimation of std-based thresholds
    param.thrL = nanmean(cwtPowerT) + param.thrF(1)*nanstd(cwtPowerT); % lower threshold
    param.thrH = nanmean(cwtPowerT) + param.thrF(2)*nanstd(cwtPowerT); % upper threshold
    param.thrM = nanmean(cwtPowerT) + param.thrF(3)*nanstd(cwtPowerT); % outlier threshold
    clear cwtPowerT
    
    
    %% find spindle events
    % detection using the upper thr, finding start/end using the lower thr
    strSpin = find(diff(sign(cwtPower-param.thrH))==2);
    endSpin = find(diff(sign(cwtPower-param.thrH))==-2);
    if endSpin(1)<strSpin(1)
        endSpin(1) = [];
    end
    if length(strSpin)>length(endSpin)
        strSpin(end) = [];
    end
    for i = 1:length(strSpin)
        strTemp = find((cwtPower(1:strSpin(i))-param.thrL)<0,1,'last');
        endTemp = find((cwtPower(endSpin(i):end)-param.thrL)<0,1,'first')+endSpin(i);
        if isempty(strTemp)||isempty(endTemp)
            strSpin(i) = NaN;
            endSpin(i) = NaN;
        else
            strSpin(i) = strTemp;
            endSpin(i) = endTemp;
        end
    end
    strSpin(isnan(strSpin)) = [];
    endSpin(isnan(endSpin)) = [];
    strSpin = unique(strSpin);
    endSpin = unique(endSpin);
    spinLen = (endSpin-strSpin)./sampRate; % spindle duration in sec
    
    
    %% go through each detected event and check spindle conditions
    nSp = 0; % spindle counter
    for i = 1:length(strSpin)
        
        % check min/max duration
        if spinLen(i)<param.minDur  || spinLen(i)>param.maxDur
            continue
        end
        
        % check min/max number of cycles
        [peakPos,locPos] = findpeaks(filtData(strSpin(i):endSpin(i)));
        [peakNeg,locNeg] = findpeaks(-filtData(strSpin(i):endSpin(i)));
        if length(peakPos)<param.minCycle || length(peakPos)>param.maxCycle
            continue
        end
        [~,idxPeakPos] = nanmax(peakPos);
        
        % check maximum power for outlier removal
        if nanmax(cwtPower(strSpin(i):endSpin(i)))>param.thrM
            continue
        end
        
        % estimate central frequency of spindle
        if i==1 || i==length(strSpin)
            segSpin = chanSignal(strSpin(i):endSpin(i));
        else
            segSpin = chanSignal(strSpin(i)-floor(sampRate/4): ...
                endSpin(i)+floor(sampRate/4));
        end
        [pxx,f] = pwelch(segSpin, hann(length(segSpin)), 0, 2*sampRate, sampRate);
        pxx = f.*pxx;       % 1/f correction
        [~,idxFreq] = max(pxx(f>=param.freq(1)&f<=param.freq(4)));
        centFreq = param.freq(1) + (idxFreq-1)*0.5;
        
        % check if power increase is spindle specific
        fTot = [6:8.5, 16.5:20];     % neighbouring frequency bands
        [~,loc] = ismember(fTot, f);
        pxxMrgn = pxx(loc);
        pxxSpin = pxx(f>param.freq(1)&f<param.freq(4));
        if max(pxxSpin)<max(pxxMrgn)
            continue
        end
        
        % save the results in the structure
        nSp = nSp+1;
        spindle(nSp).start     = strSpin(i);
        spindle(nSp).end       = endSpin(i);
        spindle(nSp).length    = spinLen(i);
        spindle(nSp).centFreq  = centFreq;
        spindle(nSp).negPeak   = -max(peakNeg);
        spindle(nSp).posPeak   = max(peakPos);
        spindle(nSp).peak2peak = max(peakPos) + max(peakNeg);
        spindle(nSp).noCycle  = length(peakPos);
        spindle(nSp).symmetry  = locPos(idxPeakPos)./(sampRate.*spinLen(i));
    end
    
    
    %% save the results
    save(['spindle_',num2str(iCh),'.mat'], 'spindle', 'param')
    fprintf('\b');
    disp(['  ' num2str(nSp), ' spindles found'])
    clear spindle cwtPower chanSignal
end

end