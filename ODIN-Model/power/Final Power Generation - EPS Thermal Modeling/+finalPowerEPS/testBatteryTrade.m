function testBatteryTrade()
% TESTBATTERYTRADE Analytic dispatch/energy checks independent of the solar solver.
s=struct('batteryTradeDischarge_percent',[0 50 55 100], ...
 'batteryMaxDischarge_W',50,'batteryCapacityBOL_Wh',45,'batteryCapacityRetained',.8, ...
 'batteryMinimumSOC',.4,'batteryMaximumSOC',1,'batteryInitialSOC',1, ...
 'batteryToBusEfficiency',.95,'batteryChargeEfficiency',.9,'batteryMaxCharge_W',20, ...
 'supplementalMargin_percent',10,'batteryTradeTolerance_percent',1e-8,'energyTolerance_Wh',1e-6);
t=0:60:6000; S=[zeros(1,30) 50*ones(1,70)]; L=23.65*ones(size(S));
r=finalPowerEPS.batteryDischargeTrade(t,S,L,s);
assert(max(abs(r.energyResidual_Wh(:)))<1e-12);
assert(all(all(r.packDischarge_W(:,31:end)==0))); % no pointless sunlit drain
assert(all(all(r.packCharge_W.*r.packDischarge_W==0)));
assert(all(all(r.SOC>=.4-1e-12 & r.SOC<=1+1e-12)));
assert(abs(r.dischargeEnergy_Wh(2)-23.65/.95*.5)<1e-10);
assert(all(abs(r.packDischarge_W(2,:)-r.packDischarge_W(4,:))<1e-10)); % more capacity does not force more use
assert(abs(r.finalSOC_percent(2)-100)<1e-9);
assert(r.bestLoadIndex==2 && r.bestIndex==3);
assert(~r.meetsLoadAllOrbit(1)&&r.meetsLoadAllOrbit(2));
assert(all(all(abs(r.servedLoad_W+r.unserved_W-L)<1e-10)));
assert(abs(r.minimumMargin_percent(3)-100*(27.5*.95/23.65-1))<1e-10);
% Low capacity cannot be made feasible by a higher discharge ceiling.
s.batteryCapacityBOL_Wh=5; low=finalPowerEPS.batteryDischargeTrade(t,S,L,s);
assert(isempty(low.bestIndex)&&all(~low.meetsLoadAllOrbit));
% A covered first orbit without recharge must not be labeled repeatable.
s.batteryCapacityBOL_Wh=45; s.batteryChargeEfficiency=0;
norecharge=finalPowerEPS.batteryDischargeTrade(t,S,L,s);
assert(norecharge.meetsLoadAllOrbit(2)&&~norecharge.repeatableAtInitialSOC(2)&&isempty(norecharge.bestIndex));
% Time refinement leaves integrated dispatch invariant for aligned discontinuities.
s.batteryChargeEfficiency=.9; tf=0:30:6000; Sf=[zeros(1,60) 50*ones(1,140)];
fine=finalPowerEPS.batteryDischargeTrade(tf,Sf,23.65*ones(size(Sf)),s);
assert(max(abs(fine.dischargeEnergy_Wh-r.dischargeEnergy_Wh))<1e-10);
assert(max(abs(fine.finalSOC_percent-r.finalSOC_percent))<1e-10);
fprintf('PASS: 15 battery trade conservation, dispatch, feasibility and refinement checks.\n');
end
