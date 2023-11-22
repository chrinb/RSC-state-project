function plot_grid_rois_on_fov(sData, sessionObjects, rem)

image_folder = 'C:\Users\chrinb\UIO Physiology Dropbox Dropbox\Lab Data\Christoffer Berge\Conferences, summer schools, grants\Conferences\SfN 2023\Figures\Blood vessel analysis\axon_fov\';

% uiopen(image_folder);
% 
% cd(image_folder);
fov_image_name = uigetfile(image_folder, '*.*');

fov_image_path = strcat(image_folder, fov_image_name);
fov_image = imread(fov_image_path);
%% 
grids_to_keep = sData.imdata.grids_to_keep_idx;  

roi_arr = sessionObjects.Data.RoiArray(1,1);

n_grids = 64;

image_dim = roi_arr.FovImageSize;

grid_size = 8;
obj.imHeight = image_dim(1);
obj.imWidth  = image_dim(2);
x0 = linspace(1, obj.imWidth, grid_size+1);
y0 = linspace(1, obj.imHeight, grid_size+1);
x0 = round(x0(1:    end-1));
y0 = round(y0(1:end-1));

tileSize = ([obj.imHeight, obj.imWidth] - 1) ./ grid_size;
X = arrayfun(@(x) (x-1) + (1:tileSize(2)), x0, 'uni', 0);
Y = arrayfun(@(y) (y-1) + (1:tileSize(1)), y0, 'uni', 0);

B = imresize(fov_image,image_dim);
%% Plot 
figure,
imagesc(B)
colormap gray
hold on
for i = 1:n_grids
    
    % test = roi_arr.roiArray(1,i).mask;
    temp_boundary = roi_arr.roiArray(1, i).boundary{1, 1};
    size_X = size(X{1,1},2);
    size_Y = size(Y{1,1},2);

    tmp  = size_X+size_Y-1;
    tmp2 = tmp+size_X-1;

    tmp3 = size_X+1;

    % length(117:172)

    x_bottom = temp_boundary( 1:size_X, 2);    
    y_bottom = temp_boundary( tmp:tmp2, 1);
    y_top = temp_boundary( 1:size_X, 1);
    y_left = temp_boundary( tmp3:tmp, 1);
    x_right = temp_boundary( tmp3:tmp, 2);
    % y3 = temp_boundary( 120:175, 2);
    % x4 = temp_boundary( 176:end, 1);
    x_left = temp_boundary( tmp2+1:end, 2);

    % x_bottom = temp_boundary( 1:57, 2);    
    % y_bottom = temp_boundary( 120:176, 1);
    % y_top = temp_boundary( 1:57, 1);
    % y_left = temp_boundary( 58:119, 1);
    % x_right = temp_boundary( 58:119, 2);
    % % y3 = temp_boundary( 120:175, 2);
    % % x4 = temp_boundary( 176:end, 1);
    % x_left = temp_boundary( 176:end, 2);
    

    
    if ismember(i, grids_to_keep) && rem == 1
        col = 'red';
    else
        col = 'white';
    end
    linew = 2;
    % Bottom
    plot(x_bottom, y_bottom, 'Color', col, 'LineWidth',linew)
    % Top
    plot(x_bottom, y_top, 'Color', col, 'LineWidth',linew)

    % Left side
    plot(x_left, y_left, 'Color', col, 'LineWidth',linew)

    plot(x_right, y_left,'Color', col, 'LineWidth',linew)

    if ismember(i, grids_to_keep) && rem == 1
        center = roi_arr.roiArray(1, i).center;
        % txt    = ['ROI ', num2str(i)];
        text(center(1)-10, center(2), 0, num2str(i), 'FontSize',15, 'Color','Red' );
    end

end
axis square
axis off
