# ODIN-Model

Orbit, sunlight and drag model for ODIN, ASEN 4018 Power. Propagates the orbit (two-body + J2, drag optional), works out the Earth shadow and gives the solar irradiance at the spacecraft.

## How to run

1. Open MATLAB in this folder
2. Run `odinPath`, adds all the folders to the path (once per session)
3. Run `main_demo` for the standard plots
4. To change anything (orbit, epoch, spacecraft, solver settings), edit the Inputs sheet in `assumptions.xlsx`, save it, run again

Own runs:

```matlab
odinPath
cfg = config(); % cfg: inputs from assumptions.xlsx
out = propagate(cfg, 15); % out: state history over 15 orbits
S = solarFlux(out, cfg); % S: irradiance at the spacecraft [W/m^2]
ecl = eclipseTimes(out, cfg); % ecl: umbra and penumbra entry and exit times [s]
```

- `config('myCase.xlsx')` runs off your own copy of the sheet instead
- one-off change in a script: `cfg = config(); cfg.h = 600;`
- drag is off by default, turn it on with `cfg.use.drag` in the sheet
- use `eclipseTimes` for eclipse lengths, the 60 s output grid overstates each one by about a minute
- after changing the model run `dragCheck(30)`, takes about a minute, should end with `checks passed: 6 of 6` (5 of 6 without `numericalCrossCheck_365d.mat`, not on GitHub as it's too big)

Needs MATLAB R2019a or newer. Demos and checks need the Aerospace Toolbox.

## Where everything is

- `assumptions.xlsx`: all inputs (Inputs sheet), drag assumptions, check results, sources
- `odinPath.m`: path setup
- `model/`: the model itself (config, propagate, eom (gravity, J2, drag), oe2rv, sunVector, shadow, shadowGeometry, eclipseTimes, solarFlux, dragTable, gcrfRotation, odinToScenario)
- `demos/`: main_demo, scenarioDemo (3D viewer), runSatViewer (one day with a ground station), plotScale
- `power/`: arrayPower, array sizing and incidence sweeps
- `studies/orbitTrade/`: pre-DR1 orbit trade (cases A to F), still used as test orbits by the checks
- `validation/`: checks against independent references, logs and plots in `validation/results` (.mat files stay local, too big)
- `external/j77/`: Jacchia 1977 atmosphere code, third party, licence in `license.txt`
