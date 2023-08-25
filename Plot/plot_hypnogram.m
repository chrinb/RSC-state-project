function plot_hypnogram(sData)

% Written by Christoffer Berge | Vervaeke la

% Function that plots session hypnogram

% TO DO: need to make scale bars!

% Initialize empty vector
hypnogram_vector = zeros(1, length(sData.behavior.NREM_vector));

% Create hypnogram vector
hypnogram_vector(sData.behavior.quiet_wakefulness == 1)  = 0;
hypnogram_vector(sData.behavior.active_wakefulness == 1) = 1;
hypnogram_vector(sData.behavior.NREM_vector == 1)        = 2;
hypnogram_vector(sData.behavior.REM_vector  == 1)        = 3;

% Find minimum and maximum values in signals for centering plot label
ECoG_min_max_values          = [min(sData.ephysdata2.lfp), max(sData.ephysdata2.lfp)];
EMG_min_max_values           = [min(sData.ephysdata3.lfp), max(sData.ephysdata3.lfp)];
Running_speed_min_max_values = [min(sData.daqdata.runSpeed), max(sData.daqdata.runSpeed)];

% Get various behavioral state snippets

% Snippet length = 10s
snippet_length = 2500*10;

% Offset the extracted signal snippet somewhat (REM theta is typically better 
% a few seconds into REM states rather than from the beginning of the
% state)
snippet_offset = 25000;

[eventStartIdx, eventStopIdx ]  = findTransitions( sData.behavior.quiet_wakefulness );
% Check if any state epochs are > 10 seconds
if sum(eventStopIdx-eventStartIdx > snippet_length)
    % Find index of first state epoch > 10 seconds
    idx = find( (eventStopIdx-eventStartIdx > snippet_length),1);
    % Extract signal snippet: Start from state start idx by a certain
    % offset, then go to state start idx + offset + snippet length, and
    % subtract by 1 to get equal length as time vector for plotting
    ECoG_QW_snippet = sData.ephysdata2.lfp( eventStartIdx(idx):eventStartIdx(idx) + ( snippet_length)-1);
end

% Repeat for REM and NREM
[eventStartIdx, eventStopIdx ]  = findTransitions( sData.behavior.REM_vector);
if sum(eventStopIdx-eventStartIdx > snippet_length)

    idx = find( (eventStopIdx-eventStartIdx > snippet_length),1);
    ECoG_REM_snippet = sData.ephysdata2.lfp( eventStartIdx(idx)+snippet_offset:eventStartIdx(idx) + ( snippet_offset+snippet_length)-1);
else
    ECoG_REM_snippet = NaN;
end

[eventStartIdx, eventStopIdx ]  = findTransitions( sData.behavior.NREM_vector );
if sum(eventStopIdx-eventStartIdx > snippet_length)
    
    idx               = find( (eventStopIdx-eventStartIdx > snippet_length),1);
    ECoG_NREM_snippet = sData.ephysdata2.lfp( eventStartIdx(idx):eventStartIdx(idx) + ( snippet_length)-1);
else 
    ECoG_NREM_snippet = NaN;
end

% Compute time-frequency values
[freq_ecog, time_ecog, pow_ecog, freq_lfp, time_lfp, pow_lfp] = plot_state_time_frequency(sData);
 
%% Plot hypnogram, ECoG, and EMG

% Create time vector for plotting
time_vector  = linspace(0, length(sData.ephysdata.lfp), length(sData.ephysdata.lfp))/2500;
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

h1 = subplot(7,3,[4:6]);
plot(time_vector, sData.ephysdata2.lfp, 'k'),
set(gca, 'xlim', [time_vector(1) time_vector(end)])

% Create y label
axis off
ytickvalue = median(ECoG_min_max_values);
x = zeros(size(ytickvalue));
str = {'ECoG'};
text(x, ytickvalue, str, 'HorizontalAlignment', 'right','FontSize',font_size, 'Position',[-5 ytickvalue 0]);

h2 = subplot(7,3,[7:9]);
plot(time_vector, sData.ephysdata3.lfp, 'k'),
set(gca, 'xlim', [time_vector(1) time_vector(end)])
% Create y label
axis off
ytickvalue = median(EMG_min_max_values);
x = zeros(size(ytickvalue));
str = {'EMG'};
text(x, ytickvalue, str, 'HorizontalAlignment', 'right','FontSize',font_size, 'Position',[-5 ytickvalue 0]);

h3 = subplot(7,3,[10:12]);
plot(time_vector, sData.daqdata.runSpeed, 'k')
set(gca, 'xlim', [time_vector(1) time_vector(end)], 'ylim',[0 10])

% Create y label
axis off
% ytickvalue = median(Running_speed_min_max_values);
ytickvalue = 5;
x = zeros(size(ytickvalue));
str = {'Running speed'};
text(x, ytickvalue, str, 'HorizontalAlignment', 'right','FontSize',font_size, 'Position',[-5 ytickvalue 0]);

h4 = subplot(7,3, [13:15]);
plot(time_vector, hypnogram_vector, 'k', 'LineWidth',1);
set(gca, 'xlim', [time_vector(1) time_vector(end)], 'ylim', [0 3])
% Create y label
axis off
ytickvalues = 0:3;
x = zeros(size(ytickvalues));
str = {'QW','AW','NREM', 'REM'};
for i = ytickvalues+1
    text(x(i), ytickvalues(i), str(i), 'HorizontalAlignment', 'right','FontSize',font_size, 'Position',[-5 i-1 0]);
end

% Time-frequency plot
h5 = subplot(7,3,[16:18]);
contourf(time_ecog,freq_ecog,10*log10( pow_ecog),300,'linecolor','none'), caxis([-50 -20]), colormap jet
text(-55, 15,0, 'S1 ECoG')
% xlabel('Time (sec)')
ylabel('(Hz)')
% cb = colorbar;
% cb.Position(1) = 0.91;
% cb.Position(3) = 0.01;
% ylabel(cb,'Power (dB)')

h6 = subplot(7,3,[19:21]);
contourf(time_lfp,freq_lfp,10*log10( pow_lfp),300,'linecolor','none'), caxis([-50 -20]), colormap jet
text(-55, 15,0, 'CA1 LFP')
xlabel('Time (sec)')
ylabel('(Hz)')
cb = colorbar;
cb.Position(1) = 0.91;
cb.Position(3) = 0.01;
ylabel(cb,'Power (dB)')

linkaxes([h1, h2, h3, h4], 'x');




