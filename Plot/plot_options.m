function [y_label, y_label2, x_label] = plot_options(keyword)

% Written by Christoffer Berge | Vervaeke Lab

% Check if user specified signal type as third input argument
try
%     switch keyword{1,3}
%         case 'dff'
%             y_label = 'Mean DF/F (z-score)';
%         case 'deconv'
%             y_label = 'Mean deconv. DF/F (z-score)';
%     end
     if sum( cellfun(@(c) contains('dff',c), keyword))
            y_label = 'Mean DF/F (z-score)';
     elseif sum( cellfun(@(c) contains('deconv',c), keyword))
            y_label = 'Mean deconv. DF/F (z-score)';
    end
catch
end

% Check if user specified signal type as third input argument
try
     if sum( cellfun(@(c) contains('bulk',c), keyword))
            y_label2 = '# SWR';
     else
            y_label2 = '# ROI';
    end
catch
end
% Check if user specified event type as fourth input argument
try
%     switch keyword{1,4}
%         case 'swr'
%             x_label = 'Time from SWR peak (s)';
%         case 'spindle'
%             x_label = 'Time from spindle onset (s)';
%         case 'so'
%             x_label = 'Time from so trough (s)';
%     end
   if sum( cellfun(@(c) contains('swr',c), keyword))
            x_label = 'Time from SWR peak (s)';
   elseif sum( cellfun(@(c) contains('spindle',c), keyword))
            x_label = 'Time from spindle onset (s)';
   elseif sum( cellfun(@(c) contains('so',c), keyword))
            x_label = 'Time from so trough (s)';
    end
catch
end


