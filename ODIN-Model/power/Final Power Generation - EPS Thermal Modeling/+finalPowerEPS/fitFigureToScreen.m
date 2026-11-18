function outerPosition = fitFigureToScreen(fig, preferredSize_px)
% FITFIGURETOSCREEN Keep the complete undocked window inside one monitor.
% preferredSize_px is [width height] for the OUTER window, including MATLAB
% menus, toolbars and title bar. Smaller screens reduce it automatically.
% The monitor under the mouse is used; no figure-size global defaults change.
if nargin < 2, preferredSize_px = [1100 740]; end
validateattributes(preferredSize_px, {'numeric'}, ...
    {'real','finite','positive','vector','numel',2});
root = groot;
oldUnits = get(root,'Units');
restoreUnits = onCleanup(@() set(root,'Units',oldUnits));
set(root,'Units','pixels');
monitors = get(root,'MonitorPositions');
pointer = get(root,'PointerLocation');
clear restoreUnits;
index = find(pointer(1)>=monitors(:,1) & pointer(1)<sum(monitors(:,[1 3]),2) & ...
    pointer(2)>=monitors(:,2) & pointer(2)<sum(monitors(:,[2 4]),2),1);
if isempty(index), index = 1; end
monitor = monitors(index,:);
% Leave room for the OS taskbar/dock, window borders and display scaling.
availableSize = [0.92 0.84].*monitor(3:4);
windowSize = min(preferredSize_px(:).',availableSize);
usableOrigin = monitor(1:2)+[0.04 0.10].*monitor(3:4);
origin = usableOrigin+(availableSize-windowSize)/2;
outerPosition = round([origin windowSize]);
set(fig,'WindowStyle','normal','Units','pixels');
if isprop(fig,'WindowState'), set(fig,'WindowState','normal'); end
set(fig,'OuterPosition',outerPosition,'PaperPositionMode','auto');
end
