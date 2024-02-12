function all_data_modified = pad_short_roi_stats(all_data)

%{
Pad shorter ROI stats vectors (because of uneven nr of cells in each
sextile) with NaNs.
%}

% Find max length
 max_length = max( length(vertcat( all_data{:, 9})), length(vertcat( all_data{:, 10})));

 % Pad the shorter vector with NaN values
for i = 9:20

    tmp_data{i} = vertcat( all_data{:, i});

    if length(tmp_data{i}) < max_length
        tmp_data{i}(end+1:max_length) = NaN;
    end

end

all_data_modified = tmp_data;

