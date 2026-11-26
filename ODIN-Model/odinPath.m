function odinPath()
    % odinPath: adds the ODIN-Model folders to the MATLAB path for this session.
    
    root = fileparts(mfilename('fullpath')); % root: repository root folder
    folderlist = dir(fullfile(root,"**"));
    folderlist = folderlist([folderlist.isdir] & ~ismember({folderlist.name},[".",".."]));
    
    addpath(root);
    for i=1:length(folderlist)
        if ~contains(path,fullfile(folderlist(i).folder,folderlist(i).name))
            addpath(fullfile(folderlist(i).folder,folderlist(i).name));
        end
    end
end
