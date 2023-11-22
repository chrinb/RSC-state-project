function locate_ch1_grid_dff_signals

% Define the folder where you want to start the search.
rootFolder = 'D:\Data\';  % Replace with the path to your root folder.

% rootFolder = 'D:\Data\m6134';  % Replace with the path to your root folder.
% fileList = dir(fullfile(rootFolder, '**', ['*' desiredFileName '*' desiredFileType]));

relevant_subfolders = {'mouse6134', 'mouse6136', 'mouse6137', 'mouse6139', 'mouse6140'};

% Define the file name and file type you're looking for.
desiredFileName = 'mean_f_ch_1';      % Replace with the desired file name.
desiredFileType = '.mat';          % Replace with the desired file type (e.g., '.txt').

% Loop through the relevant subfolders.
for i = 1:numel(relevant_subfolders)

    currentSubfolder = relevant_subfolders{i};
    subfolderPath = fullfile(rootFolder, currentSubfolder);
    
    % Recursively search for files in the current subfolder and its subfolders.
    fileList = dir(fullfile(subfolderPath, '**', ['*' desiredFileName '*' desiredFileType]));
    
    % Loop through the list of found files and load them into the workspace or perform other actions.
    for j = 1:numel(fileList)
        fullFilePath = fullfile(fileList(j).folder, fileList(j).name);
        
        % You can load the file into the workspace or perform other actions here.
        % For example, to load a text file into a variable, you can use 'load':
        data{i,j} = load(fullFilePath);
        
        % Or, if you want to display the file paths:
        disp(['Found file: ', fullFilePath]);
    end
end