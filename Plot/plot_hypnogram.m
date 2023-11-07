function plot_hypnogram(sData)

% Written by Christoffer Berge | Vervaeke la

% Function that plots session hypnogram

% TO DO: need to make scale bars!

% Initialize empty vector
hypnogram_vector = zeros(1, length(sData.behavior.NREM_vector));

% frames         = sData.daqdata.frame_onset_reference_frame;
srate = 2500;
ecog  = sData.ephysdata2.lfp;
emg   = sData.ephysdata3.lfp;
run_speed = sData.daqdata.runSpeed;

% Create hypnogram vector
hypnogram_vector(sData.behavior.quiet_wakefulness == 1)  = 0;
hypnogram_vector(sData.behavior.active_wakefulness == 1) = 1;
hypnogram_vector(sData.behavior.NREM_vector == 1)        = 2;
hypnogram_vector(sData.behavior.REM_vector  == 1)        = 3;

% Find minimum and maximum values in signals for centering plot label
ECoG_min_max_values          = [min(ecog), max(ecog)];
EMG_min_max_values           = [min(emg), max(emg)];
% Running_speed_min_max_values = [min(sData.daqdata.runSpeed), max(sData.daqdata.runSpeed)];

% Get various behavioral state snippets

% Offset the extracted signal snippet somewhat (REM theta is typically better 
% a few seconds into REM states rather than from the beginning of the
% state)
snippet_length = srate*10;
% snippet_offset = srate*5;

[eventStartIdx_QW, eventStopIdx_QW ]  = findTransitions( sData.behavior.quiet_wakefulness );
% Check if any state epochs are > 10 seconds
if sum(eventStopIdx_QW-eventStartIdx_QW > snippet_length)
    % Find index of first state epoch > 10 seconds
    idx_QW = find( (eventStopIdx_QW-eventStartIdx_QW > snippet_length),1);

    QW_bout_middle = median( eventStartIdx_QW(idx_QW): eventStopIdx_QW(idx_QW));
    QW_bout_start = QW_bout_middle - snippet_length/2;
    QW_bout_end   = QW_bout_middle + snippet_length/2 - 1;
    % Extract signal snippet: Start from state start idx by a certain
    % offset, then go to state start idx + offset + snippet length, and
    % subtract by 1 to get equal length as time vector for plotting
%     ECoG_QW_snippet = sData.ephysdata2.lfp( eventStartIdx_QW(idx_QW)+snippet_offset:eventStartIdx_QW(idx_QW) + (snippet_length-1) + snippet_offset);
    ECoG_QW_snippet = sData.ephysdata2.lfp( QW_bout_start:QW_bout_end);

end

% Repeat for REM and NREM   
[eventStartIdx_REM, eventStopIdx_REM ]  = findTransitions( sData.behavior.REM_vector);
if sum(eventStopIdx_REM-eventStartIdx_REM > snippet_length)

    idx_REM = find( (eventStopIdx_REM-eventStartIdx_REM > snippet_length),1);

    REM_bout_middle = median( eventStartIdx_REM(idx_REM): eventStopIdx_REM(idx_REM));
    REM_bout_start = REM_bout_middle - snippet_length/2;
    REM_bout_end   = REM_bout_middle + snippet_length/2 - 1;
    ECoG_REM_snippet = sData.ephysdata2.lfp(REM_bout_start:REM_bout_end);
else
    ECoG_REM_snippet = NaN;
end

[eventStartIdx_NREM, eventStopIdx_NREM ]  = findTransitions( sData.behavior.NREM_vector );
if sum(eventStopIdx_NREM-eventStartIdx_NREM > snippet_length)
    
    idx_NREM               = find( (eventStopIdx_NREM-eventStartIdx_NREM > snippet_length),1);

    NREM_bout_middle = median( eventStartIdx_NREM(idx_NREM): eventStopIdx_NREM(idx_NREM));
    NREM_bout_start = NREM_bout_middle - snippet_length/2;
    NREM_bout_end   = NREM_bout_middle + snippet_length/2 - 1;
    ECoG_NREM_snippet = sData.ephysdata2.lfp(NREM_bout_start:NREM_bout_end);
%     ECoG_NREM_snippet = sData.ephysdata2.lfp( eventStartIdx_NREM(idx_NREM):eventStartIdx_NREM(idx_NREM) + ( snippet_length)-1);
else 
    ECoG_NREM_snippet = NaN;
end

% Compute time-frequency values
[freq_ecog, time_ecog, pow_ecog, freq_lfp, time_lfp, pow_lfp] = plot_state_time_frequency(sData);
 
%% Plot hypnogram, ECoG, and EMG

% Create time vector for plotting
time_vector  = linspace(0, length(ecog), length(ecog))/2500;

%%  Downsample x 10
time_vector = downsample( time_vector, 10);
ecog        = downsample( ecog, 10);
emg         = downsample( emg, 10);
run_speed   = downsample( run_speed, 10);
hypnogram_vector = downsample( hypnogram_vector, 10);
% freq_ecog   = downsample(freq_ecog, 10);
%  time_ecog = downsample(time_ecog, 10);
%  pow_ecog = downsample(pow_ecog, 10);
% freq_lfp = downsample(freq_lfp, 10);
%  time_lfp = downsample(time_lfp, 10);
%  pow_lfp  = downsample(pow_lfp, 10);

 %% Plot
time_snippet = linspace(0, snippet_length, snippet_length)/2500;
font_size    = 10;
y_lim = [-.2 .6];

figure, 
if ~isnan(ECoG_NREM_snippet)
    subplot(7,3,1)
    plot(time_snippet, ECoG_NREM_snippet, 'Color', [0.8500 0.3250 0.0980])
    set(gca, 'xlim', [time_snippet(1) time_snippet(end)], 'ylim',y_lim)
    title('NREM')
    axis off
end

if ~isnan(ECoG_REM_snippet)
    subplot(7,3,2)
    plot(time_snippet, ECoG_REM_snippet, 'color',[0.6350 0.0780 0.1840])
    set(gca, 'xlim', [time_snippet(1) time_snippet(end)], 'ylim',y_lim)
    title('REM')
    axis off
end

subplot(7,3,3)
plot(time_snippet, ECoG_QW_snippet, 'color',[0.3010 0.7450 0.9330])
set(gca, 'xlim', [time_snippet(1) time_snippet(end)], 'ylim', y_lim)
title('QW')
axis off,

% ECoG
h(1) = subplot(7,3,(4:6));
plot(time_vector, ecog, 'k'),
set(gca, 'xlim', [time_vector(1) time_vector(end)])
axis on
ytickvalue = median(ECoG_min_max_values);
x = zeros(size(ytickvalue));
str = {'ECoG'};
text(x, ytickvalue, str, 'HorizontalAlignment', 'right','FontSize',font_size, 'Position',[-5 ytickvalue 0]);
set(gca,'xtick',[])

% Plot QW snippet limits
xline( QW_bout_start/srate, '--', 'LineWidth',2, 'color', [0.3010 0.7450 0.9330])
xline( QW_bout_end/srate, '--', 'LineWidth',2, 'color', [0.3010 0.7450 0.9330])

% Plot REM snippet limits
xline( REM_bout_start/srate, '--', 'LineWidth',2, 'color',[0.6350 0.0780 0.1840])
xline( REM_bout_end/srate, '--', 'LineWidth',2, 'color', [0.6350 0.0780 0.1840])

% Plot NREM snippet limites
xline( NREM_bout_start/srate, '--', 'LineWidth',2, 'color', [0.8500 0.3250 0.0980])
xline( NREM_bout_end/srate, '--', 'LineWidth',2, 'color', [0.8500 0.3250 0.0980])


% EMG
h(2) = subplot(7,3,(7:9));
plot(time_vector, emg, 'k'),
set(gca, 'xlim', [time_vector(1) time_vector(end)])
axis on
ytickvalue = median(EMG_min_max_values);
x = zeros(size(ytickvalue));
str = {'EMG'};
text(x, ytickvalue, str, 'HorizontalAlignment', 'right','FontSize',font_size, 'Position',[-5 ytickvalue 0]);
set(gca,'xtick',[])

% Run speed
h(3) = subplot(7,3,(10:12));
plot(time_vector, run_speed, 'k')
set(gca, 'xlim', [time_vector(1) time_vector(end)], 'ylim',[0 10])
axis on
ytickvalue = 5;
x = zeros(size(ytickvalue));
str = {'Running speed'};
text(x, ytickvalue, str, 'HorizontalAlignment', 'right','FontSize',font_size, 'Position',[-5 ytickvalue 0]);
set(gca,'xtick',[])

% Hypnogram
h(4) = subplot(7,3, (13:15));
plot(time_vector, hypnogram_vector, 'k', 'LineWidth',1);
set(gca, 'xlim', [time_vector(1) time_vector(end)], 'ylim', [0 3])
axis off
ytickvalues = 0:3;
x = zeros(size(ytickvalues));
str = {'QW','AW','NREM', 'REM'};
set(gca,'xtick',[])

for i = ytickvalues+1
    text(x(i), ytickvalues(i), str(i), 'HorizontalAlignment', 'right','FontSize',font_size, 'Position',[-5 i-1 0]);
end

% Time-frequency plot
h(6) = subplot(7,3,(16:18));
% contourf(time_ecog,freq_ecog,10*log10( pow_ecog),300,'linecolor','none'), 
imagesc(time_ecog,flipud(freq_ecog), 10*log10(pow_ecog) )

clim([-45 -20]), colormap jet
text(-60, 15,0, 'S1 ECoG')
ylabel('(Hz)', FontSize = 10);
h(6).FontSize = 10;
set(gca,'xtick',[])
cb = colorbar;
cb.Position(1) = 0.91;
cb.Position(3) = 0.01;
ylabel(cb,'Power (dB)')
y_ticks = [ 20 10 1 ];
yticklabels(y_ticks)

h(7) = subplot(7,3,(19:21));
% contourf(time_lfp,freq_lfp,10*log10( pow_lfp),300,'linecolor','none')
imagesc(time_lfp,flipud(freq_lfp), 10*log10(pow_lfp) ), colorbar
yticklabels(y_ticks)
clim([-45 -20]), colormap jet
text(-57, 15,0, 'CA1 LFP')
h(7).FontSize = 10;
xlabel('Time (s)', FontSize=14)
ylabel('(Hz)', FontSize=9)
cb = colorbar;
cb.Position(1) = 0.91;
cb.Position(3) = 0.01;
ylabel(cb,'Power (dB)')

linkaxes(h, 'x');




