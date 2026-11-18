function figures=plotModelDiagnostics(out,cfg,array,components,iv)
% PLOTMODELDIAGNOSTICS Additional figures only; existing figures are untouched.
% Panel terminal power and delivered bus load remain on separately labeled axes.
s=finalPowerEPS.modelSettings('diagnostics'); figures=gobjects(0);
if ~s.makeAdditionalPlots, return; end
wiring=sprintf('%gS%gP assumed',iv.nominal.params.seriesCells,iv.nominal.params.parallelStrings);
caption=sprintf('Area %.4g m^2 | altitude %.4g km | %s', ...
    array.params.selectedArrayArea_m2,cfg.h,wiring);
fig=figure('WindowStyle','normal','Color','w','Name','Panel front-minus-rear temperature');
finalPowerEPS.fitFigureToScreen(fig,s.figureSize_px); figures(end+1)=fig;
layout=tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
title(layout,'Temperature difference across panel'); subtitle(layout,caption);
a=nexttile(layout);
plot(a,array.time_s/60,array.cellTemperature_C-array.backTemperature_C,'LineWidth',1.4);
grid(a,'on'); ylabel(a,'Front - rear [K]'); title(a,'Available-power area study');
b=nexttile(layout); hold(b,'on');
plot(b,components.time_s/60,components.nominal.solar.cellTemperature_C-components.nominal.solar.backTemperature_C,'LineWidth',1.4);
plot(b,components.time_s/60,components.peak.solar.cellTemperature_C-components.peak.solar.backTemperature_C,'--','LineWidth',1.4);
grid(b,'on'); ylabel(b,'Front - rear [K]'); xlabel(b,'Time [min]');
title(b,'Coupled extraction cases'); legend(b,{'Nominal','Peak'},'Location','southoutside','Orientation','horizontal');
linkaxes([a b],'x');

fig=figure('WindowStyle','normal','Color','w','Name','Coupled solar power and battery operation');
finalPowerEPS.fitFigureToScreen(fig,s.figureSize_px); figures(end+1)=fig;
layout=tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');
title(layout,'Available power, extraction, load and storage'); subtitle(layout,caption);
colors=lines(2); names={'nominal','peak'};
a=nexttile(layout); hold(a,'on'); b=nexttile(layout); hold(b,'on'); c=nexttile(layout); hold(c,'on');
for k=1:2
    v=components.(names{k}); x=components.timeEdges_s/60;
    stairs(a,x,[v.solar.potentialElectricalPower_W v.solar.potentialElectricalPower_W(end)], ...
        'Color',colors(k,:),'LineStyle','--','LineWidth',1.2,'DisplayName',[names{k} ' available']);
    stairs(a,x,[v.solar.electricalPower_W v.solar.electricalPower_W(end)], ...
        'Color',colors(k,:),'LineWidth',1.4,'DisplayName',[names{k} ' extracted']);
    stairs(b,x,[v.load_W v.load_W(end)],'Color',colors(k,:),'LineStyle','--', ...
        'LineWidth',1.2,'DisplayName',[names{k} ' demanded']);
    stairs(b,x,[v.servedLoad_W v.servedLoad_W(end)],'Color',colors(k,:), ...
        'LineWidth',1.4,'DisplayName',[names{k} ' served']);
    plot(c,x,100*v.batterySOC,'Color',colors(k,:),'LineWidth',1.4,'DisplayName',names{k});
end
ylabel(a,'Panel terminals [W]'); title(a,'Available versus actual solar extraction');
ylabel(b,'Delivered bus load [W]'); title(b,'Demanded versus served load');
ylabel(c,'Battery SOC [%]'); title(c,'State of charge'); xlabel(c,'Time [min]');
for ax=[a b c]
    grid(ax,'on'); legend(ax,'Location','southoutside','Orientation','horizontal');
end
linkaxes([a b c],'x'); xlim(c,[out.t(1) out.t(end)]/60);
end
