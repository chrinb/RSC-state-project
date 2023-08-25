function [] = plotSpindleSpectrogram(sData)

% Written by Christoffer Berge | Vervaeke Lab

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
sigma_str = strcat('sigmaband', num2str(ephys_select) );
LFPsigma = sData.(channelID).(sigma_str);
spindle_idx_str = strcat('NREMAbsSpindleIdx', spindle_select);
spindle_snip_str = strcat('NREMspindleSnips', spindle_select);
spindle_idx = sData.(channelID).(spindle_idx_str);

allSpindles = zeros(length(spindle_idx), 20001);
for i = 1:length(spindle_idx)
        try
        allSpindles(i,:) = sData.(channelID).(spindle_snip_str)(i).lfp;
        end
end

% spindle_segments = zeros(length(spindle_idx), 20001);
mean_spindle_waveform = nanmean(allSpindles);

window  = 2500;              % Window size for computing the spectrogram (FFT) [# samples]
overlap = 2400;              % Overlap of the windows for computing the spectrogram [# samples]
nFFT    = 5:0.3:17;        % Vector defining the frequencies for computing the FFT
% nFFT    = linspace(5,20,50);
Fs      = 2500;             % Signal sampling frequency.


for z = 1:size(allSpindles,1)
    fprintf('Calculating spectrogram...spindle %d of %d\n',z,size(allSpindles,1));
        
    % skip any snippets containing NaNs
    if isnan(allSpindles(z,1))
        % do nothing
    else
    x = allSpindles(z,:);
    [S,~,~,cP] = spectrogram(x,window,overlap,nFFT,Fs);
    P(:,:,z) = cP;
    S(:,:,z) = S;
    end
end
P = nanmean(P,3);
S = nanmean(S,3);
[~,F,T,~] = spectrogram(x,window,overlap,nFFT,Fs);

srate = 2500;
time = -4:1/srate:4;    
spec_time = linspace(-4,4,length(T));
% spec_time = (-(length(T)/2):(length(T)/2)-1)./2500;

% Scale mean spindle waveform for plotting and center it appropriately
% (this waveform plot is only for visualization)
scaled_spindle_wav = (mean_spindle_waveform*50);
y_to_plot = median(F);
scale_factor = y_to_plot-scaled_spindle_wav(1);
spindle_wav = scaled_spindle_wav + scale_factor;
% Plot spectrogram
figure
contourf(spec_time,F,(abs(P)),200, 'edgecolor','none');
colormap(jet)
hold on
plot(time,spindle_wav, 'color', 'w', 'linew',1)
%view([0 90])
colorbar

title(['Average spindle spectrogram', ' (nr of spindles = ' num2str(size(allSpindles,1)),')'])
ylabel('Frequency (Hz)')
xlabel('Time (s)')

