function [x_label, y_label, y_label2] = get_xy_labels(varargin)

% Get labels for plot

keyword = varargin{1,1};

switch keyword{1,1}
    case 'deconv'
        y_label2_part2 = 'dec.';
    case 'dff'
        y_label2_part2 = 'DF/F';
end

switch keyword{1,2}
    case 'spindle'
        x_label =  'Time from spindle onset (s)';
    case 'swr'
        x_label =  'Time from SWR peak (s)';
    case 'swa'
        x_label =  'Time from SWA trough (s)';     
end

switch keyword{1,3}
    case 'population'
        y_label = '# ROI';
    case 'bulk'
        y_label = ['# ', keyword{1,2}];
end

switch keyword{1,5}
    case 'avg'
        y_label2 = ['Mean ', y_label2_part2, '(z-score)'];
    case 'mod'
        y_label2 = ['# ', keyword{1,2}];
end