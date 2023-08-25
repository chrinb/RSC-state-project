function [sData, upThres, lowThres, tempListSpindles, spindleHalfLength] = mark_spindle(sData, makeSpectrogramPlot)

% Written by Christoffer Berge || Vervaeke lab

% Function that (1) selects sleep spindle/sigma frequency band signal, (2)
% obtains the amplitude envelope of that signal using the Hilbert
% transform, (3) creates threshold and finds peaks in signal above
% threshold, (4) identifies putative sleep spindles according to user
% specified criteria, (5) opens manual curation of putative spindle events,
% (6) stores output in sData struct. 

%% choose signal
prompt = sprintf('Select channel: ');
signal = input(prompt);

prompt = sprintf('8-18Hz (1) or 10-16hz (2)? ');
spindle_band_select = input(prompt);

if spindle_band_select == 1
    ephys_select = [];
    spindle_select = [];
elseif spindle_band_select == 2
    ephys_select = num2str(2);
    spindle_select = '1016';
end


t1 = num2str(signal);
t2 = 'ephysdata';
channelID = strcat(t2,t1);

LFP = sData.(channelID).lfp;
% Use 10-16 Hz filtered signal
temp_str = strcat('sigmaband', num2str(ephys_select) );
LFPsigma = sData.(channelID).(temp_str);
LFPfilt = sData.ephysdata2.lfpFilt;

%% create signal, calculate amplitude envelope from filtered
% signal via Hilbert transform, smooth amplitude envelope, find upper/lower
% threshold, peaks in data

time = (0:length(LFP)-1) / 2500;

% find amplitude envelope
ampl_env = abs(hilbert(LFPsigma));

% smooth amplitude envelope with a Gaussian filter of 200ms
smooth_ampl_env = smoothdata(ampl_env, 'Gaussian', 500);

%% Test
% ampl_envZ = abs(hilbert(zscore(LFPsigma)));
% smooth_ampl_env = smoothdata(ampl_envZ, 'Gaussian', 500);
% upThres = 2.5;
% LowThres = 2.5;
% LFPsigma = zscore(LFPsigma);

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
        %check if difference is >= 1250 samples (0.5s). If not, merge the two
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
    if  putativeSpindles(i,2)-putativeSpindles(i,1) < (2500/2) || putativeSpindles(i,2)-putativeSpindles(i,1) > (2500*5)
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

%% go to manual inspection of spindles and assign them to sData struct
[final_spindleLFP, final_spindleLocs, final_spindle_snip_loc] = inspectSpindles(spindleSnips,spindleCentralTimePoint,LFP, spindleHalfLength,tempListSpindles ,makeSpectrogramPlot);

% assign data to correct sData struct
sData.(channelID).absSpindleIdx   = final_spindleLocs;
sData.(channelID).spindleSnips    = final_spindleLFP;
sData.(channelID).spindleStartEnd = final_spindle_snip_loc;

%% calculate spindle parameters 

% initialize vectors
%cyclesValue = zeros(length(final_spindleLFP),30);
nrOfCycles = zeros(length(final_spindleLFP),1);
intraSpindleFreq = zeros(length(final_spindleLFP),1);
spindleDur = zeros(length(final_spindleLFP),1);
spindleAmp = zeros(length(final_spindleLFP),1);
for i = 1:length(final_spindleLFP)
    
    % find nr of cycles per spindle
    cyclesValue = findpeaks(LFPsigma(final_spindle_snip_loc(i,1):final_spindle_snip_loc(i,2)));
    nrOfCycles(i) = numel(cyclesValue);
    
    % find spindle duration in seconds
    spindleDur(i) = length( (final_spindle_snip_loc(i,1):final_spindle_snip_loc(i,2)) ) / 2500 ;
    
    % calculate intra-spindle frequency
    intraSpindleFreq(i) = nrOfCycles(i)/ spindleDur(i);
    
    % find spindle amplitude
    [S,L] = bounds(LFPsigma(final_spindle_snip_loc(i,1):final_spindle_snip_loc(i,2)));
    spindleAmp(i) = abs(S)+L;
    % find spindle density (nr of spindles/minute)
end

spindleCycl_str = strcat('spindleCycl', spindle_select);
spindleDur_str = strcat('spindleDur', spindle_select);
spindleFreq_str = strcat('spindleFreq', spindle_select);
spindleAmp_str = strcat('spindleAmp', spindle_select);

sData.(channelID).(spindleCycl_str) = nrOfCycles;
sData.(channelID).(spindleDur_str)  = spindleDur;
sData.(channelID).(spindleFreq_str) = intraSpindleFreq;
sData.(channelID).(spindleAmp_str)  = spindleAmp;


