function value = REM_check(sData)

% Written by Christoffer Berge | Vervaeke lab

% Function that looks for REM episode(s) in sData.episode struct and
% returns true if session contains any REM episodes, false if not.

if isfield(sData, 'episodes')

    value = sum( ismember(sData.episodes{:,1}, 'REM')) > 0;

end
