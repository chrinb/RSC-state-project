function test_window = select_mod_win(varargin)

% Written by Christoffer Berge | Vervaeke lab

% Function that selects a baseline time window for which to compare
% significant modulation responses to various events.

win_length = varargin{1,1};
keyword    = varargin{1,2};

frames_in_1sec = 31;

switch keyword
    case 'SWR'
        % window (0 to +500ms)
        test_window    = round( (win_length/2) : (win_length/2) + frames_in_1sec/2 );
    case 'Spindle'
        % window 0ms - 1s
        test_window    = round( (win_length/2) :  (win_length/2) + frames_in_1sec);
end

               
