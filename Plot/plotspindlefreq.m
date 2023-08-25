function plotspindlefreq(sData)

% Function that plots the raw LFP signal, the signal filtered in the
% spindle frequency (sigma) band, and all detected, putative spindles
% before manual sorting. This is mainly for inspection and evaluating
% whether the threshold criteria are appropriate. 

% Written by Christoffer Berge || Vervaeke lab

%% choose signal
prompt = sprintf('Select channel: ');
signal = input(prompt);

t1 = num2str(signal);
t2 = 'ephysdata';
channelID = strcat(t2,t1);

LFP = sData.(channelID).lfp;
LFPsigma = sData.(channelID).sigmaband;

%% create signal, calculate amplitude envelope from filtered
% signal via Hilbert transform, smooth amplitude envelope, find upper/lower
% threshold, peaks in data

time = (0:length(LFP)-1) / 2500;

% find amplitude envelope
ampl_env = abs(hilbert(LFPsigma));

% smooth amplitude envelope with a Gaussian filter of 200ms
smooth_ampl_env = smoothdata(ampl_env, 'Gaussian', 500);

% determine upper and lower threshold of the amplitude envelope. Based on 
% procedure outlied in Nir et al., (2011), Sela et al., (2016), and 
% Kim et al., (2019). Parameters may need adjustment. 
upThres = mean(smooth_ampl_env) + 2*std(smooth_ampl_env);
lowThres = mean(smooth_ampl_env) + 0.2*std(smooth_ampl_env);

% Find peaks above upper threshold
[~,locs, ~,~] = findpeaks(smooth_ampl_env, 'MinPeakHeight', upThres, 'WidthReference', 'halfprom','Annotate', 'Extents');

% set all values below lower threshold to zero
smooth_ampl_env(smooth_ampl_env < lowThres) = 0;
ampl_vec_length = 1:length(smooth_ampl_env);

PeaksAboveThres = locs;

%% Go through data and find spindle start/end times, discard too short/long
% spindles, merge spindles if time interval is sufficiently short (0.5s),
% discard peaks that occur inside the previous spindle

% pre-allocate matrix for start/end times for putative spindles
putativeSpindles = zeros(length(PeaksAboveThres),2);

% loop over nr of peaks above threshold
for nrPeaks = 1:length(PeaksAboveThres)
    h = 1;
    l = 1;
    m = 1;
    p = 1;
    o = 1;
    % for all peaks > peak nr 1, check if the peak value is contained 
    %within the lower boundaries of the prior peak, and if so, skip the current peak. 
    while nrPeaks >= 2 && h == 1
        if PeaksAboveThres(nrPeaks) < max(putativeSpindles(:,2))
           l = 0; 
           m = 0;
           KK = 0;
           HH = 0;
           p = 0;
           %putativeSpindles(nrPeaks,:) = [];
        else
            
        end
        h = 0;
    end
    
            
        TimeOfPeak = PeaksAboveThres(nrPeaks);
        % go backwards and forwards from the current peak value to find the
        % nearest zero crossings with the lower threshold to determine
        % spindle length. 
        k = -1;
        while l > 0
                aa = ampl_vec_length(TimeOfPeak) + k;
            if aa == 0 || smooth_ampl_env(aa) == 0 
                KK = ampl_vec_length(aa+1);
                l = 0;
            else
                k = k -1;
            end
        end
      
        t = 1;
        while m > 0
            bb = ampl_vec_length(TimeOfPeak) + t;
            if bb > length(smooth_ampl_env) || smooth_ampl_env(bb) == 0 
                HH = ampl_vec_length(bb-1);
                m = 0;
            else
                t = t + 1;
            end
        end        
    
    % assign each peak with its crossing of the lower threshold to
    % determine spindle length. 
    putativeSpindles(nrPeaks,:) = [KK, HH];
    
    %calculate time window for first spindle
    %while nrPeaks == 1 && o == 1
    %    TimeBetweenSpindleStartStop = putativeSpindles(nrPeaks,2)-putativeSpindles(nrPeaks,1);
    %    if TimeBetweenSpindleStartStop < 1250
    %        putativeSpindles(nrPeaks,:) = [0;0];
    %    end
    %    o = 0;
    %end
    
    % calculate time between the onset of current spindle and end of
    % previous spindle (by their zero-crossing of the lower threshold). If
    % time window is < 1s (2500 sample points), merge these spindles
    while nrPeaks > 1 && p == 1 
        %find difference between onset of current spindle and end of
        %previous spindle
        TimeBetweenSpindleLowThres = putativeSpindles(nrPeaks,1) - max(putativeSpindles(1:(nrPeaks-1),2));
        %check if difference is >= 1250 frames (0.5s). If not, merge the two
        %spindles. 
        if TimeBetweenSpindleLowThres <= 500 && max(putativeSpindles(1:(nrPeaks-1),2)) > 0
            cc = putativeSpindles(:,:) == max(putativeSpindles(1:(nrPeaks-1),2));
            putativeSpindles(cc) = putativeSpindles(nrPeaks,2);
            putativeSpindles(nrPeaks,:) = [0; 0];
            p = 0;
        else
            p = 0;
        end
        
        % calculate miminum duration of spindles. Set to 
        %TimeBetweenSpindleStartStop = putativeSpindles(nrPeaks,2)-putativeSpindles(nrPeaks,1);
        %if TimeBetweenSpindleStartStop < 1250
        %    putativeSpindles(nrPeaks,:) = [0;0];
        %end
    end
end

%loop over spindles and delete putative spindles with duration < 0.5s or >
% 5s. 
for i = 1:length(putativeSpindles)
    if  (putativeSpindles(i,2)-putativeSpindles(i,1) < (2500/2) || putativeSpindles(i,2)-putativeSpindles(i,1) > (2500*5))
    putativeSpindles(i,:) = [0,0];
    end   
end

%delete empty rows and columns
putativeSpindles(putativeSpindles == 0) = [];
tempListSpindles = reshape(putativeSpindles, length(putativeSpindles)/2,2);

%% determine spindle centre
spindleHalfLength = zeros(length(tempListSpindles),1);
spindleCentralTimePoint = zeros(length(tempListSpindles),1);

for i = 1:length(tempListSpindles)
    spindleCentralTimePoint(i) = median( linspace(tempListSpindles(i,1), tempListSpindles(i,2)));
    spindleHalfLength(i) = length(tempListSpindles(i,1):spindleCentralTimePoint(i));
end

spindleSnips = struct();
spindleLocs = spindleCentralTimePoint;
spindleIdx = zeros(1,length(spindleLocs));

%% create window around spindle centre for visual inspection of spindles 

% (spindle centre +- 4s)
for i = 1:length(spindleLocs)
    lfpStartIndex = spindleLocs(i) - 4*2500;
    lfpEndIndex = spindleLocs(i) + 4*2500;
    if lfpStartIndex < 0
        lfpStartIndex = 1; 
    end
    if lfpEndIndex > length(LFP)
       lfpEndIndex = length(LFP);
    end
    
    spindleSnips(i).lfp = LFP(round(lfpStartIndex):round(lfpEndIndex));
    spindleIdx(i) = spindleLocs(i);
end

time_ephys = (0:length(sData.ephysdata2.lfp)-1)./2500;
test = zeros(1, length(LFP));
test(round(spindleLocs)) = 0.2;

figure,
hAx(1) = subplot(311);
plot(time_ephys, LFP), hold on,
set(gca, 'xlim',[time_ephys(1), time_ephys(end)])

hAx(2) = subplot(312);
plot(time_ephys, LFPsigma), hold on,
plot(time_ephys, smooth_ampl_env, 'linewidth',1),
plot(time_ephys, test, 'linewidth',1)
yline(upThres), yline(lowThres)
set(gca, 'xlim',[time_ephys(1), time_ephys(end)])

hAx(3) = subplot(313);
plot(time_ephys, sData.ephysdata2.deltaband, 'linewidth',1), hold on
plot(time_ephys, sData.ephysdata2.lfp)

linkaxes(hAx, 'x');

