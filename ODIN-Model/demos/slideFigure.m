function [f, ax] = slideFigure(name, w, h, fs)
% slideFigure: white figure of w x h inches with one axes for a slide, Arial, grid on, light theme where MATLAB has themes.
% name: figure name
% w: width [in]
% h: height [in]
% fs: font size, 10 when omitted [pt]
% f: figure handle
% ax: axes handle

if nargin < 4
    fs = 10;
end

f = figure('Color', 'w', 'Name', name, 'Units', 'inches', 'Position', [1 1 w h]);
try
    theme(f, 'light'); % R2025a and newer follow the desktop theme otherwise
catch
end

ax = axes(f, 'FontName', 'Arial', 'FontSize', fs);
hold(ax, 'on'); grid(ax, 'on'); box(ax, 'off');

end
