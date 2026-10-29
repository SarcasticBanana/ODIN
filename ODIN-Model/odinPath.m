function odinPath()
% odinPath: adds the ODIN-Model folders to the MATLAB path for this session.

root = fileparts(mfilename('fullpath')); % root: repository root folder
addpath(root, fullfile(root, 'model'), fullfile(root, 'external', 'j77'), fullfile(root, 'demos'), ...
        fullfile(root, 'power'), fullfile(root, 'studies', 'orbitTrade'), fullfile(root, 'validation'), ...
        fullfile(root, 'studies', 'debrisDetection'));

end
