function S = solarFlux(out, cfg)
% solarFlux: solar irradiance at the spacecraft, independent of any panel property.
% out: state history struct from propagate.m
% cfg: parameter struct from config.m
% S: 1 x N irradiance on a surface held normal to the Sun [W/m^2]

if ~isfield(out, 'rsun')
    error('solarFlux:noSun','cfg.use.eclipse must be true for propagate to store the Sun history');
end

% Inverse square falloff of the solar constant, attenuated by the Earth shadow:
S = cfg.S0 * (cfg.AU ./ vecnorm(out.rsun)).^2 .* out.illum;

end
