% Get the directory where this startup.m file is located (project root)
projectRoot = fileparts(mfilename('fullpath'));

% Add specific folders and their subfolders to path
folders_to_add = {'core', 'parameters', 'tests', 'scripts', 'scenarios'};
for i = 1:length(folders_to_add)
    folder_path = fullfile(projectRoot, folders_to_add{i});
    if exist(folder_path, 'dir')
        addpath(genpath(folder_path));
    end
end

format short g
