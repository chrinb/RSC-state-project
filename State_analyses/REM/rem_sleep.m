function REM_episodes = rem_sleep(sData)

% Written by Christoffer Berge | Vervaeke Lab

% Function that returns an array consisting of the nr of REM episodes along 
% along the row(s) and their start and end time stamps in ephys time 
% in column 1 and 2, respectively. (Depends on user having scored data beforehand 
% using Begonia). 

REMep = sData.episodes.state == 'REM';
% check that there is REM episodes in recording
if ~isempty(REMep)
    % copy array in 2nd column
    logidx = repmat(REMep, 1, 2);

    % convert table to array
    REMarr = table2array(sData.episodes(:,2:3));

    % select NREM bouts 
    REM = REMarr(logidx);

    % reshape NREM start/end vector into original matrix and multiply by sample
    % rate (2500)
    REM_episodes = reshape(REM, length(REM)/2, 2).*2500;
end
