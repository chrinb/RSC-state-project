function mean_sleep_features = sleepfeatures(sData)

% find length of recording in hours
rec_length = length(sData.ephysdata.lfp)/(2500*3600);

% determine NREM, REM, and IS bouts per hour
nr_of_NREMepisodes = sData.episodes.state == 'NREM';
boutsNREM_perH = sum( nr_of_NREMepisodes)/rec_length;

nr_of_REMepisodes = sData.episodes.state == 'REM';
boutsREM_perH = sum( nr_of_REMepisodes)/rec_length;

nr_of_ISepisodes = sData.episodes.state == 'IS';
boutsIS_perH = sum( nr_of_ISepisodes)/rec_length;

% determine percentage sleep time in each stage



% determine bout duration
sleep_durations = sData.episodes.state_duration;

% NREM
boutNREM_dur = sleep_durations(nr_of_NREMepisodes);
mean_NREMbout_duration = mean(boutNREM_dur);

% REM
boutREM_dur = sleep_durations(nr_of_REMepisodes);
mean_REMbout_duration = mean(boutREM_dur);

% IS
boutIS_dur = sleep_durations(nr_of_ISepisodes);
mean_ISbout_duration = mean(boutIS_dur);

% calculate total NREM sleep minutes 
NREMep = sData.episodes.state == 'NREM';
NREMarr = table2array(sData.episodes(:,4));
NREMepisodes = NREMarr(NREMep);
if length(NREMepisodes) > 1
    totalNREMmin = sum( vertcat(NREMepisodes))/60;
else
    totalNREMmin = NREMepisodes/60;
end


% calculate total REM sleep minutes 
REMep = sData.episodes.state == 'REM';
REMarr = table2array(sData.episodes(:,4));
REMepisodes = REMarr(REMep);
if length(REMepisodes) > 1
    totalREMmin = sum( vertcat(REMepisodes))/60;
else
    totalREMmin = REMepisodes/60;
end

% calculate total IS sleep minutes 
ISep = sData.episodes.state == 'IS';
ISarr = table2array(sData.episodes(:,4));
ISepisodes = ISarr(ISep);
if length(ISepisodes) > 1
    totalISmin = sum( vertcat(ISepisodes))/60;
else
    totalISmin = ISepisodes/60;
end

mean_sleep_features = {boutsNREM_perH, boutsREM_perH, boutsIS_perH, ...
            boutNREM_dur, boutREM_dur, boutIS_dur, ...
            totalNREMmin, totalREMmin, totalISmin};
