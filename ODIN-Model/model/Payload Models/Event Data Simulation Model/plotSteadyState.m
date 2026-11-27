clear; clc; close all;
odinPath();

% needs to match dimensions from eventDataStreamFun.m 
nx = 50;
ny = 50;

% faster timestep for earlier, then slower for later once presumed steady state has occured, this can be changed to
% any desired timestep
tEnd = [(0.1:0.1:5)'; (10:10:2000)'];
time_step = diff([0; tEnd]);
nSamples = numel(tEnd);

stats = zeros(nSamples,6);
state = [];

for k = 1:nSamples
    % Each window starts from the previous ones ending state
    [events, state] = eventDataStreamFun( ...
        nx/2, ny/2, 0, 0, 0, time_step(k), 0, state, false);

    V = state.V(:);
    Vref = state.V_ref(:);

    stats(k,1) = mean(V);
    stats(k,2) = mean(Vref);
    stats(k,3) = std(V);
    stats(k,4) = std(Vref);
    stats(k,5) = std(V - Vref);

    fprintf('%.1f seconds\n',tEnd(k));
end


% The time are plotted at each windows midpoint because they describe the mean of the window.
rateTime = tEnd - time_step/2;
panelTitles = {'Mean V', 'Mean V_{ref}', 'Std V', 'Std V_{ref}', ...
               'Std (V - V_{ref})'};

figure;
tiledlayout(3,2);
for j = 1:5
    nexttile;
    plot(tEnd,stats(:,j),'.-');
    ylabel('Voltage');
    title(panelTitles{j});
    xlabel('Time (s)');
    grid on;
end
sgtitle('Event Camera Startup');

