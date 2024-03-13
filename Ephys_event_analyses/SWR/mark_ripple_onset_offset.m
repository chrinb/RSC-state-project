function [swr_start_stop, swr_length, upper_thresh] = mark_ripple_onset_offset(varargin)

% Written by Christoffer Berge | Vervaeke Lab

% Function that calculates SWR length by finding the nearest threshold
% crossing points before and after each SWR peak (SWRs must be manually
% scored before applying this code). The threshold is the same as used for
% SWR detection in the "markRipples" code, i.e. a moving mean + 3*moving
% std (window length = sample rate * 5). 
sData = varargin{1,1};

fs         = 2500;
freqFilter = [100 300];
lfpSignal  = sData.ephysdata.lfp;
lfp        = lfpSignal;
time       = linspace(0,length(lfp),length(lfp))/(fs);

window_size = 5;
window_size_index = window_size * fs;

freqL = freqFilter(1);
freqU = freqFilter(2);

nyquistFs = fs/2;

% Threshold detecting SWR onset/offset
U_threshold = 1;  % in standard deviations

% Create filter and apply to LFP data
filter_kernel = fir1(600,[freqL freqU]./nyquistFs); % Different filters can also be tested her e.g. butter and firls
filtered_lfp = filtfilt(filter_kernel,1,lfp); % Filter LFP using the above created filter kernel

% Hilbert transform LFP to calculate envelope
lfp_hil_tf = hilbert(filtered_lfp);
lfp_envelop = abs(lfp_hil_tf);

% Smooth envelop using code from 
% https://se.mathworks.com/matlabcentral/fileexchange/43182-gaussian-smoothing-filter?focused=3839183&tab=function 
smoothed_envelop = gaussfilt_2017(time,lfp_envelop,.004);
moving_mean = movmean(smoothed_envelop, window_size_index);
moving_std = movstd(smoothed_envelop, window_size_index);

% Find upper/lower threshold values of the LFP
upper_thresh = moving_mean + U_threshold*moving_std;
    

% Find peaks of envelop. NB: The parameters of this function have to be properly
% chosen for best result.
% [~,locs,~,~] = findpeaks(smoothed_envelop-upper_thresh, fs,'MinPeakHeight',0,'MinPeakDistance',0.025,'MinPeakWidth',0.010,'WidthReference','halfhprom','Annotate','extents','WidthReference','halfprom');
% rippleLocs = round(locs,3);   

if nargin > 1
    RipIdx = sData.ephysdata.absRipIdx(varargin{1,2});
else
    RipIdx = sData.ephysdata.absRipIdx;
end

if isempty(RipIdx)
    [swr_start_stop, swr_length, upper_thresh] = deal([]);
else

    % set all values below lower threshold to zero
    smoothed_envelop(smoothed_envelop < upper_thresh) = 0;
    ampl_vec_length = 1:length(smoothed_envelop);
    
    for i = 1:length(RipIdx)
        
        TimeOfPeak = RipIdx(i);
        l = 1;
        m = 1;
        k = -1;
        while l > 0
            aa = ampl_vec_length(TimeOfPeak) + k;
            if aa == 0 || smoothed_envelop(aa) == 0 
            KK = ampl_vec_length(aa+1);
            l = 0;
            else
            k = k -1;
            end
         end
          
        t = 1;
        while m > 0
            bb = ampl_vec_length(TimeOfPeak) + t;
            if bb > length(smoothed_envelop) || smoothed_envelop(bb) == 0 
                HH = ampl_vec_length(bb-1);
                m = 0; 
            else
                t = t + 1;
            end
        end
        swr_start_stop(i,:) = [KK, HH];
    end
    
    swr_length = swr_start_stop(:,2)-swr_start_stop(:,1);
    swr_start_stop;

end
