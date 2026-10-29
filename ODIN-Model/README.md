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
- `constellation(cfg)` gives one cfg per satellite plus names, plane by plane, from the Constellation section of the sheet
- main_demo runs arrayPower at the panel area `cfg.power.panelArea`, Power section of the sheet
- `scenarioDemo([], [], 24*7, [], true)`: a week instead of the default day, plus each satellite's previous and next orbit as ground tracks (off by default)
- `RUN`: all satellites over a week, orbit, shadow, array power, ground track and downlink figures, every contact window to the millisecond in `results/contacts.csv`, then the 3D viewer, run length and viewer length at the top
- ground stations: Stations sheet of `assumptions.xlsx`, one row per station, Use set to FALSE leaves one out, `contactTimes(out, cfg, st(1))` gives one station's windows from `st = stations()`
- `dl = downlinkTimes(out, cfg)`: downlink start and end times in seconds from sim start, one row per window, any station, RUN leaves one per satellite in `dl{j}`
- use `eclipseTimes` for eclipse lengths, the 60 s output grid overstates each one by about a minute
- after changing the model run `dragCheck(30)`, takes about a minute, should end with `checks passed: 6 of 6` (5 of 6 without `numericalCrossCheck_365d.mat`, not on GitHub as it's too big)

Needs MATLAB R2019a or newer. Demos and checks need the Aerospace Toolbox.

## Where everything is

- `assumptions.xlsx`: all inputs (Inputs sheet), ground stations (Stations sheet), drag assumptions, check results, sources
- `odinPath.m`: path setup
- `RUN.m`: one run of everything, contact table in `results/` (kept local)
- `model/`: the model itself (config, propagate, eom (gravity, J2, drag), oe2rv, sunVector, shadow, shadowGeometry, eclipseTimes, solarFlux, dragTable, gcrfRotation, odinToScenario, constellation, stations, earthRotation, stationVector, contactTimes, downlinkTimes)
- `demos/`: main_demo, scenarioDemo (3D viewer), runSatViewer (one day with the stations in the Stations sheet), plotScale
- `power/`: arrayPower, array sizing and incidence sweeps
- `studies/orbitTrade/`: pre-DR1 orbit trade (cases A to F), still used as test orbits by the checks
- `validation/`: checks against independent references, logs and plots in `validation/results` (.mat files stay local, too big)
- `external/j77/`: Jacchia 1977 atmosphere code, third party, licence in `license.txt`
