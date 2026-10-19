function [sat, R] = odinToScenario(sc, out, name)
% odinToScenario: ephemeris satellite in a satelliteScenario from a propagate.m state history.
% sc: satelliteScenario object that receives the satellite
% out: state history struct from propagate.m
% name: satellite name, 'ODIN' when omitted
% sat: Satellite object carrying the GCRF position and velocity history
% R: 3x3 rotation matrix from the mean equator and equinox of cfg.epoch to GCRF

if nargin < 3
    name = 'ODIN';
end

t0 = datetime(out.cfg.epoch, 'TimeZone', sc.StartTime.TimeZone); % t0: UTC instant of cfg.epoch
tt = t0 + seconds(out.t(:)); % tt: N x 1 UTC sample times

if sc.StartTime < tt(1) || sc.StopTime > tt(end)
    error('odinToScenario:span', 'scenario window must lie inside the propagated span');
end

R = gcrfRotation(out.cfg);

rG = R * out.r * 1e3; % rG: 3 x N position in GCRF [m]
vG = R * out.v * 1e3; % vG: 3 x N velocity in GCRF [m/s]

pos = timetable(tt, rG.'); % pos: position ephemeris [m]
vel = timetable(tt, vG.'); % vel: velocity ephemeris [m/s]

sat = satellite(sc, pos, vel, 'CoordinateFrame', 'inertial', 'Name', name);

end
