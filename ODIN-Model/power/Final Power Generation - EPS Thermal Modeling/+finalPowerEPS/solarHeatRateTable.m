function values=solarHeatRateTable(r)
Quantity={'Qout';'QradiationToSpace';'QtoStructure';'QheatInput';'Qstorage'};
summaries={r.outflowSummary;r.spaceRadiationSummary;r.structureHeatSummary; ...
    r.totalHeatSummary;r.storageSummary};
Mean_W=zeros(5,1); Min_W=Mean_W; Max_W=Mean_W;
for k=1:5
    Mean_W(k)=summaries{k}.mean_W; Min_W(k)=summaries{k}.min_W; Max_W(k)=summaries{k}.max_W;
end
values=table(Quantity,Mean_W,Min_W,Max_W,'VariableNames',{'Quantity','Mean_W','Min_W','Max_W'});
end
