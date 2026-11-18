function out = propagate(cfg, nOrbits) % integrates the equations of motion with ode45.
% Verification-only standalone dynamics. ODIN integration consumes the team's
% already propagated orbit via finalPowerEPS.run and never calls this helper.
% cfg: parameter struct from config.m
% nOrbits: number of orbital periods
% out: struct holding the state history and derived quantities

if nargin < 2
    nOrbits = 3;
end
if (isfield(cfg.use,'drag') && cfg.use.drag) || (isfield(cfg.use,'srp') && cfg.use.srp)
    error('finalPowerEPS:standaloneDynamics', ...
        'Use the unchanged ODIN propagate for drag/SRP, then pass cfg/out to finalPowerEPS.run.');
end

a = cfg.Re + cfg.h; % a: semi-major axis [km]
T = 2*pi*sqrt(a^3/cfg.mu); % T: orbital period [s]

[r0, v0] = finalPowerEPS.oe2rv(a, cfg.ecc, cfg.inc, cfg.raan0, cfg.argp0, cfg.nu0, cfg.mu);
% r0: initial position vector in ECI [km]
% v0: initial velocity vector in ECI [km/s]

validateattributes(nOrbits,{'numeric'},{'real','finite','scalar','positive'});
validateattributes(cfg.ode.dtOut,{'numeric'},{'real','finite','scalar','positive'});
finish=nOrbits*T;
tspan=unique([0:cfg.ode.dtOut:finish,finish]);

opts = odeset('RelTol', cfg.ode.relTol, 'AbsTol', cfg.ode.absTol);
% opts: ode45 options structure

solution=ode45(@(t,x) finalPowerEPS.eom(t,x,cfg),[0 finish],[r0;v0],opts);
eventTimes=struct();
if cfg.use.eclipse
    searchStep=10; transitionStep=1;
    if isfield(cfg.eclipse,'searchStep_s'), searchStep=cfg.eclipse.searchStep_s; end
    if isfield(cfg.eclipse,'transitionStep_s'), transitionStep=cfg.eclipse.transitionStep_s; end
    validateattributes(searchStep,{'numeric'},{'real','finite','scalar','positive'});
    validateattributes(transitionStep,{'numeric'},{'real','finite','scalar','positive'});
    survey.t=unique([0:min(searchStep,cfg.ode.dtOut):finish,finish]);
    surveyState=deval(solution,survey.t); survey.r=surveyState(1:3,:);
    survey.rsun=finalPowerEPS.sunVector(survey.t,cfg); survey.orbitSolution=solution;
    eventTimes=finalPowerEPS.eclipseTimes(survey,cfg);
    contacts=sort([eventTimes.umbraIn,eventTimes.umbraOut,eventTimes.penumbraIn,eventTimes.penumbraOut]);
    transitionTimes=[];
    % Refine penumbra and any grazing passage, independently of figure cadence.
    passageEdges=unique([0,contacts,finish]);
    for j=1:numel(passageEdges)-1
        passageStart=passageEdges(j); passageEnd=passageEdges(j+1);
        state=deval(solution,(passageStart+passageEnd)/2);
        fraction=finalPowerEPS.shadow(state(1:3),finalPowerEPS.sunVector((passageStart+passageEnd)/2,cfg),cfg);
        if fraction>0 && fraction<1
            transitionTimes=[transitionTimes,passageStart:transitionStep:passageEnd,passageEnd]; %#ok<AGROW>
        end
    end
    tspan=unique([tspan,contacts,transitionTimes]);
end
t=tspan(:); X=deval(solution,t).';
% t: column vector of output times [s]
% X: N x 6 state history, one row per output time [km, km/s]

out.t = t(:)'; % out.t: 1 x N output times since sim start [s]
out.r = X(:,1:3)'; % out.r: 3 x N position history in ECI [km]
out.v = X(:,4:6)'; % out.v: 3 x N velocity history in ECI [km/s]
out.eclipseTimes=eventTimes;
out.requestedDuration_s=finish;
out.T = T; % out.T: orbital period [s]
out.a = a; % out.a: semi-major axis [km]
out.cfg = cfg; % out.cfg: copy of the parameter struct used for this run

% Specific angular momentum h = r x v
out.hvec = cross(out.r, out.v);   % out.hvec: 3 x N specific angular momentum [km^2/s]

if cfg.use.eclipse
    [out.rsun, out.sunhat] = finalPowerEPS.sunVector(out.t, cfg);
    % out.rsun: 3 x N Sun position in ECI [km]
    % out.sunhat: 3 x N unit vector toward the Sun

    out.illum = finalPowerEPS.shadow(out.r, out.rsun, cfg);
    % out.illum: 1 x N visible fraction of the solar disk, 0 in umbra and 1 in full sunlight
end

end
