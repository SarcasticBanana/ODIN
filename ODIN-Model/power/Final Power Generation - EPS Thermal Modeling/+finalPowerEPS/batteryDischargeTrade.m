function trade=batteryDischargeTrade(t,solar_W,load_W,s)
% BATTERYDISCHARGETRADE Sweep an allowed discharge ceiling with demand-following dispatch.
% Each row starts with the same SOC. Solar supplies load first, surplus recharges,
% and battery discharge supplies only the deficit. Capacity is an energy/Wh proxy.
% Headroom is the extra supply CAPABILITY for this interval at its initial SOC;
% it is not continuously exported power and not a sustained extra-load guarantee.
t=t(:).'; solar_W=solar_W(:).'; load_W=load_W(:).'; dt=diff(t);
assert(numel(dt)==numel(solar_W)&&numel(dt)==numel(load_W)&&all(dt>0),'Invalid time basis.');
assert(all(isfinite(solar_W)&solar_W>=0)&all(isfinite(load_W)&load_W>0),'Invalid solar/load power.');
trade.command_percent=s.batteryTradeDischarge_percent(:); % retained API: now a CEILING, not constant output
trade.commandPack_W=trade.command_percent*s.batteryMaxDischarge_W/100;
trade.time_s=t; trade.targetMargin_percent=s.supplementalMargin_percent;
n=numel(trade.command_percent); m=numel(dt);
capacity=s.batteryCapacityBOL_Wh*s.batteryCapacityRetained;
floorE=capacity*s.batteryMinimumSOC; ceilE=capacity*s.batteryMaximumSOC;
trade.storedEnergy_Wh=zeros(n,m+1); trade.storedEnergy_Wh(:,1)=capacity*s.batteryInitialSOC;
trade.packDischarge_W=zeros(n,m); trade.packCharge_W=zeros(n,m);
trade.availableBatteryBus_W=zeros(n,m);
deficit=max(load_W-solar_W,0); surplus=max(solar_W-load_W,0);
% The trade allows battery assistance whenever a deficit exists, including sunlight.
for j=1:m
    E=trade.storedEnergy_Wh(:,j);
    availablePack=min(trade.commandPack_W,max(E-floorE,0)*3600/dt(j));
    trade.availableBatteryBus_W(:,j)=availablePack*s.batteryToBusEfficiency;
    discharge=min(availablePack,deficit(j)/s.batteryToBusEfficiency);
    charge=min(min(surplus(j)*s.batteryChargeEfficiency,s.batteryMaxCharge_W),max(ceilE-E,0)*3600/dt(j));
    trade.packDischarge_W(:,j)=discharge; trade.packCharge_W(:,j)=charge;
    trade.storedEnergy_Wh(:,j+1)=E+(charge-discharge)*dt(j)/3600;
end
trade.SOC=trade.storedEnergy_Wh/capacity;
trade.busDischarge_W=trade.packDischarge_W*s.batteryToBusEfficiency;
trade.servedLoad_W=min(load_W,solar_W+trade.busDischarge_W);
trade.unserved_W=max(load_W-trade.servedLoad_W,0);
trade.supply_W=solar_W+trade.availableBatteryBus_W;
trade.balance_percent=100*(trade.supply_W-load_W)./load_W; % capability, not delivered excess
trade.minimumMargin_percent=min(trade.balance_percent,[],2);
trade.unservedEnergy_Wh=sum(trade.unserved_W.*dt,2)/3600;
trade.dischargeEnergy_Wh=sum(trade.packDischarge_W.*dt,2)/3600;
trade.chargeEnergy_Wh=sum(trade.packCharge_W.*dt,2)/3600;
trade.minimumSOC_percent=100*min(trade.SOC,[],2);
trade.finalSOC_percent=100*trade.SOC(:,end);
trade.meetsLoadAllOrbit=trade.unservedEnergy_Wh<=s.energyTolerance_Wh;
atTarget=trade.balance_percent>=trade.targetMargin_percent-s.batteryTradeTolerance_percent;
trade.meetsTargetAllOrbit=all(atTarget,2)&trade.meetsLoadAllOrbit;
trade.restoresInitialSOC=trade.storedEnergy_Wh(:,end)>=trade.storedEnergy_Wh(:,1)-s.energyTolerance_Wh;
trade.repeatableAtInitialSOC=trade.meetsLoadAllOrbit & trade.restoresInitialSOC;
trade.timeAtTarget_min=sum(atTarget.*dt,2)/60;
trade.bestIndex=find(trade.meetsTargetAllOrbit & trade.restoresInitialSOC,1,'first');
trade.bestLoadIndex=find(trade.repeatableAtInitialSOC,1,'first');
trade.bestCommand_percent=[];
if ~isempty(trade.bestIndex), trade.bestCommand_percent=trade.command_percent(trade.bestIndex); end
trade.energyResidual_Wh=diff(trade.storedEnergy_Wh,1,2)-(trade.packCharge_W-trade.packDischarge_W).*dt/3600;
trade.note=['Independent power-limit cases; demand-following discharge and solar-surplus recharge. ' ...
 'Margin is interval-start available supply headroom, not exported excess. Repeated-orbit feasibility ' ...
 'requires load coverage and restoration of initial SOC under the same periodic solar/load history. ' ...
 'The +margin test does not propagate a continuously increased load.'];
trade.summary=table(trade.command_percent,trade.commandPack_W,trade.minimumMargin_percent, ...
 trade.unservedEnergy_Wh,trade.dischargeEnergy_Wh,trade.chargeEnergy_Wh, ...
 trade.minimumSOC_percent,trade.finalSOC_percent,trade.meetsTargetAllOrbit,trade.repeatableAtInitialSOC, ...
 'VariableNames',{'Limit_percent','Limit_W','MinHeadroom_percent','Unserved_Wh','Discharge_Wh', ...
 'Charge_Wh','MinSOC_percent','FinalSOC_percent','MeetsHeadroom','RepeatableLoad'});
trade.summary.Properties.VariableDescriptions={ ...
 'Discharge ceiling as percent of configured pack power limit; not SOC or energy consumed.'; ...
 'Maximum allowed pack discharge power [W].'; ...
 'Minimum interval-start available supply margin relative to load [%].'; ...
 'Unserved load energy over the orbit [Wh].'; ...
 'Actual battery discharge energy [Wh].'; ...
 'Actual stored charging energy [Wh], after charge efficiency.'; ...
 'Minimum state of charge [% of retained energy capacity].'; ...
 'Final state of charge [% of retained energy capacity].'; ...
 'Available headroom meets target at every sampled interval and base load is served.'; ...
 'Base load served and initial SOC restored under this repeated periodic envelope.'};
end
