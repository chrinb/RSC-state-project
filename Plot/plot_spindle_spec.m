function plot_spindle_spec(sData)

LFP = sData.ephysdata2.absSpindleIdx;
nSpindles = length(LFP);
CX = 'RSC';
for i = 1:nSpindles
    try
        allSpindles(i,:) = sData.ephysdata2.spindleSnips(i).lfp;
    end
end

spindle_segments = zeros(nSpindles, 20001);

for i = 1:nSpindles
    spindle_segments(i,:) = sData.ephysdata2.spindleSnips(i).lfp;  
end
mean_spindle_waveform = mean(spindle_segments);


% frequency parameters
min_freq =  5;
max_freq = 20;
num_frex = 40;
frex = linspace(min_freq,max_freq,num_frex);

range_cycles = [ 8 14 ];
srate = 2500;
time = -4:1/srate:4;
times_to_save = -4:.10:4;
tidx = dsearchn(time', times_to_save');


s = logspace(log10(range_cycles(1)),log10(range_cycles(end)),num_frex) ./ (2*pi*frex);
wavtime = -2:1/srate:2;
half_wave = (length(wavtime)-1)/2;

% FFT parameters
nWave = length(wavtime);
nData = size(spindle_segments,1)*size(spindle_segments,2);
nConv = nWave + nData - 1;

% now compute the FFT of all trials concatenated
alldata = reshape( spindle_segments',1,[]);
dataX   = fft( alldata ,nConv );


% initialize output time-frequency data
% initialize output time-frequency data
tffull = zeros(num_frex,length(time));
tfdwns = zeros(num_frex,length(tidx));
% loop over frequencies
for fi=1:length(frex)
    
    % create wavelet and get its FFT
    wavelet  = exp(2*1i*pi*frex(fi).*wavtime) .* exp(-wavtime.^2./(2*s(fi)^2));
    waveletX = fft(wavelet,nConv);
    waveletX = waveletX ./ max(waveletX);
    
    % now run convolution in one step
    as = ifft(waveletX .* dataX);
    as = as(half_wave+1:end-half_wave);
    
    % and reshape back to time X trials
    as = reshape( as, size(spindle_segments,2), size(spindle_segments,1) );
    
    % compute power and average over trials
    powTS = mean( abs(as).^2 ,2);
    tfdwns(fi,:) =  powTS(tidx);
    tffull(fi,:) = powTS;
end

figure
contourf(times_to_save,frex,tfdwns,100,'linecolor','none'), hold on
plot(time, (mean_spindle_waveform*50)+15, 'color', 'w', 'linew',1)
colormap jet
xlabel('Time from spindle center (sec)'), 
ylabel('Frequency (Hz)')