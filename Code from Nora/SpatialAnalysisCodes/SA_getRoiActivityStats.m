function roiStat = SA_getRoiActivityStats(dff)
%getRoiActivityStats Return struct of different roi activity statistics
%
%   roiStat = getRoiActivityStats(sData) Returns struct with statistics
%   about the roi signal for all rois in sData. 
%
%   roiStat = getRoiActivityStats(sData, CH) calculates the roiStat on the
%   specified channel number CH. CH is an integer. The default value is 2.
%
%   Output, roiStat contains the following fields:
%       peakDff         : Peak dff of the roi time series. 
%       signalToNoise   : Signal to noise 
%       activityLevel   : Fraction of time where activity is above noise
%   
    
% %     pstr = getFilePath(sData.sessionID, 'roiStat');
% %     if exist(pstr, 'file')
% %         roiStat = loaddata(sData.sessionID, 'roiStat');
% %         return
% %     end
    
    if nargin < 2
        CH = 2; % Most of us are imaging on ch2?
    end
    
    % dff = double(squeeze(sData.imdata.roiSignals(CH).dff));
    [nRois, nSamples] = size(dff);
    

    % Get max DFF of all rois.
    peakDff = max(dff, [], 2);
    
    
    % Get SNR of all Rois.
    signalToNoise = zeros(nRois, 1);
    noiseLevel = zeros(nRois, 1);
    for i = 1:nRois
        if isnan(sum(dff(i, :)))
           continue 
        end
        noiseLevel(i) = real(GetSn(dff(i, :)));
        signalToNoise(i) = snr(dff(i, :), ones(1,nSamples) * noiseLevel(i));
    end
    
    
    % Get fraction of time above noise level
    dffSmooth = okada(dff, 2); % Smooth with okada filter
    dffSmooth = smoothdata(dffSmooth, 2, 'movmean', 7); % Smooth again with movmean
    
    isHigh = dffSmooth > noiseLevel;
    
    activityLevel = sum(isHigh, 2) ./ nSamples;
    
    roiStat = struct('peakDff', peakDff, ...
                     'signalToNoise', signalToNoise, ....
                     'activityLevel', activityLevel);
               
    
%     savedata(sData.sessionID, struct('roiStat', roiStat))
    
end


function sn = GetSn(Y, range_ff, method)
%% Estimate noise standard deviation

%% inputs:
%   Y: N X T matrix, fluorescence trace
%   range_ff : 1 x 2 vector, nonnegative, max value <= 0.5, range of frequency (x Nyquist rate) over which the spectrum is averaged
%   method: string, method of averaging: Mean, median, exponentiated mean of logvalues (default)

%% outputs:
%   sn: scalar, std of the noise

%% Authors: Pengcheng Zhou, Carnegie Mellon University, 2016
% adapted from the MATLAB implemention by Eftychios Pnevmatikakis and the
% Python implementation from Johannes Friedrich

%% References
% Pnevmatikakis E. et.al., Neuron 2016, Simultaneous Denoising, Deconvolution, and Demixing of Calcium Imaging Data

%% input arguments
if ~exist('range_ff', 'var') || isempty(range_ff)
    range_ff = [.25, .5];
end
if ~exist('method', 'var') || isempty(method)
    method = 'logmexp';
end
if any(size(Y)==1)
    Y = reshape(Y, [], 1);
else
    Y = Y';
end

%% estimate the noise
[psdx, ff] = pwelch(Y, [],[],[], 1);
indf = and(ff>=range_ff(1), ff<=range_ff(2));
switch method
    case 'mean'
        sn=sqrt(mean(psdx(indf, :)/2));
    case 'median'
        sn=sqrt(median(psdx(indf,:)/2));
    case 'logmexp'
        sn = sqrt(exp(mean(log(psdx(indf,:)/2))));    
    otherwise
        fprintf('wrong method! use logmexp instead.\n'); 
        sn = sqrt(exp(mean(log(psdx(indf,:)/2))));
end
sn = sn';
end