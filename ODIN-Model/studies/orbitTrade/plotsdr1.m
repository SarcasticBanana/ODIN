% plotsdr1: Model 1 slide plots, outputs fig_dr1_ltan.png (longest shadow and passes in shadow against plane LTAN) and fig_dr1_shadow_bars.png (battery energy used in the worst shadow of the year).

clear; close all; clc;

L = 0:2:24; % L: plane LTAN, local time of the northbound equator crossing [h]
sh550 = [35.80 34.74 29.57 22.14 29.79 34.82 35.80 34.67 29.60 22.19 29.37 34.59 35.80]; % sh550: longest shadow in one year at 550 km [min]
sh650 = [35.56 34.41 28.66 20.22 28.91 34.49 35.56 34.34 28.69 20.25 28.43 34.25 35.56]; % sh650: longest shadow in one year at 650 km [min]
ps550 = [9.76 9.42 6.09 0.58 6.05 9.36 9.79 8.43 3.96 0.00 4.03 8.36 9.76]; % ps550: ground station pass time inside shadow at 550 km, yearly average [min/day]
ps650 = [11.72 11.24 7.10 0.30 7.08 11.18 11.73 9.91 4.20 0.00 4.34 9.84 11.72]; % ps650: ground station pass time inside shadow at 650 km, yearly average [min/day]

grp = {'Plane at 06:00 or 18:00', '4 h from noon or midnight', '2 h from noon or midnight', 'Plane at 00:00 or 12:00'}; % grp: plane LTAN groups, bottom to top
tSh = [22.08 28.90 34.48 35.55]; % tSh: worst shadow of the year in each group [min]
tPs = [2.38 8.96 8.99 9.03]; % tPs: ground station pass time inside that shadow [min]
Pcam = 4.5; % Pcam: event camera power [W]
Pbus = 5.0; % Pbus: bus power [W]
Ptx = 15; % Ptx: radio transmit power [W]
etaD = 0.95; % etaD: battery discharge efficiency
Ebat = 45; % Ebat: battery capacity [Wh]
soh = 0.80; % soh: end-of-life capacity fraction
dod = 0.25; % dod: depth-of-discharge limit

% Battery discharge energy, observing and pass in shadow:
Eobs = (Pcam + Pbus)*(tSh - tPs)/60/etaD; % Eobs: battery energy while observing [Wh]
Eps = (Pbus + Ptx)*tPs/60/etaD; % Eps: battery energy during the pass in shadow [Wh]

% Battery energy over the whole shadow:
Etot = Eobs + Eps; % Etot: battery energy used in the worst shadow [Wh]

% Depth-of-discharge limit on the end-of-life capacity:
Eref = dod*soh*Ebat; % Eref: battery energy at the depth-of-discharge limit [Wh]

band = [0 2; 10 14; 22 24]; % band: LTAN ranges within 2 h of noon or midnight [h]
cBand = [0.92 0.92 0.92]; % cBand: band color
tickLbl = arrayfun(@(x) sprintf('%02d:00', mod(x, 24)), L, 'UniformOutput', false); % tickLbl: LTAN tick labels
ytl = cellfun(@(g, t) sprintf('%s\\newline%.0f min shadow', g, t), grp, num2cell(tSh), 'UniformOutput', false); % ytl: two-line group labels

f1 = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 7 4.5]); % f1: LTAN sweep figure

subplot(2, 1, 1); hold on; grid on; box off;
for k = 1:size(band, 1) % k: band index
    patch(band(k, [1 2 2 1]), [15 15 40 40], cBand, 'EdgeColor', 'none');
end
p1 = plot(L, sh550, '-o', 'LineWidth', 1.5); % p1: 550 km line
p2 = plot(L, sh650, '-o', 'LineWidth', 1.5); % p2: 650 km line
text(0.2, 15.5, 'midnight', 'VerticalAlignment', 'bottom');
text(12, 15.5, 'within 2 h of noon', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(23.8, 15.5, 'midnight', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
xlim([0 24]); ylim([15 40]); xticks(L); xticklabels({});
ylabel({'longest shadow', '[min]'});
legend([p1 p2], {'550 km', '650 km'}, 'Location', 'northoutside', 'Orientation', 'horizontal');

subplot(2, 1, 2); hold on; grid on; box off;
for k = 1:size(band, 1)
    patch(band(k, [1 2 2 1]), [0 0 13 13], cBand, 'EdgeColor', 'none');
end
plot(L, ps550, '-o', 'LineWidth', 1.5);
plot(L, ps650, '-o', 'LineWidth', 1.5);
xlim([0 24]); ylim([0 13]); xticks(L); xticklabels(tickLbl);
ylabel({'passes in shadow', '[min/day]'});
xlabel('plane LTAN (local time at the northbound equator crossing)');

exportgraphics(f1, 'fig_dr1_ltan.png', 'Resolution', 300);

f2 = figure('Color', 'w', 'Units', 'inches', 'Position', [1 1 7 4.5]); % f2: battery energy figure
hold on; box off; set(gca, 'XGrid', 'on', 'TickLabelInterpreter', 'tex');
b = barh(1:4, [Eobs(:) Eps(:)], 0.5, 'stacked'); % b: observing and pass bar series
pRef = plot([Eref Eref], [0.4 4.6], 'k--', 'LineWidth', 1.2); % pRef: depth-of-discharge limit line
text(Etot + 0.15, 1:4, arrayfun(@(e) sprintf('%.1f Wh', e), Etot, 'UniformOutput', false));
xlim([0 10]); ylim([0.4 4.6]); yticks(1:4); yticklabels(ytl);
xlabel('battery energy used in the worst shadow of the year [Wh]');
legend([b(:); pRef], {sprintf('observing: camera + bus, %.1f W', Pcam + Pbus), sprintf('pass in shadow: radio + bus, %.0f W', Pbus + Ptx), ...
    sprintf('%.0f%% of the %.0f Wh battery at end of life (%.0f Wh)', 100*dod, Ebat, Eref)}, 'Location', 'northoutside');

exportgraphics(f2, 'fig_dr1_shadow_bars.png', 'Resolution', 300);