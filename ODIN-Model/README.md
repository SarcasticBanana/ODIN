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
b = batteryAging(cfg); % b: battery capacity, resistance and lowest charge, day by day over cfg.power.years
```

- `config('myCase.xlsx')` runs off your own copy of the sheet instead
- one-off change in a script: `cfg = config(); cfg.h = 600;`
- drag is off by default, turn it on with `cfg.use.drag` in the sheet, density from `cfg.drag.model`: NRLMSISE-00 (day and night, latitude, season, daily Ap) or J77 (one global profile)
- `constellation(cfg)` gives one cfg per satellite plus names, plane by plane, from the Constellation section of the sheet
- main_demo runs arrayPower at the panel area `cfg.power.panelArea`, Power section of the sheet
- arrayPower's aging factor is the product of the five end-of-life factors in the Power section, pessimistic set by default (0.8895), nominal set in the notes (0.9505)
- batteryAging runs off the Battery section, one 35.8 min shadow every orbit at 550 km by default, `plot(b.t/365.25, b.socMin)` for the lowest charge against years
- `scenarioDemo([], [], 24*7, [], true)`: a week instead of the default day, plus each satellite's previous and next orbit as ground tracks (off by default)
- `RUN`: all satellites over a week, orbit, shadow, array power, battery and operating mode, observation slew, ground track and downlink figures, every contact window to the millisecond in `results/contacts.csv`, drag decay of each plane over the mission at ECSS low, moderate and high solar activity with the day and night density, then the 3D viewer, run length and viewer length at the top
- ground stations: Stations sheet of `assumptions.xlsx`, one row per station, Use set to FALSE leaves one out, `contactTimes(out, cfg, st(1))` gives one station's windows from `st = stations()`
- `dl = downlinkTimes(out, cfg)`: downlink start and end times in seconds from sim start, one row per window, any station, RUN leaves one per satellite in `dl{j}`
- `pm = powerModes(out, cfg)`: battery charge and operating mode through a run, starting full in Nominal, loads and thresholds in the Operating modes section, `pm.mode` 1 Nominal, 2 Raven, 3 Safe, a panel tilt history as third input
- `fr = slewFrequency(out, cfg, [30 45 60], [63.3 90])`: observation windows tilted from the Sun, each started once the battery is full again, windows per day, recovery times and lowest charge per tilt and length, RUN prints them for ODIN 1.1 and the 550 km battery sizing orbit
- `dec = dragDecay(cfg, 365)`: drag decay of a circular orbit, density averaged around one orbit every quarter day so the day and night sides both count, `dec.h` the semi-major axis less `cfg.Re` against `dec.t`, stops at `cfg.decay.hFloor`, step and the ECSS activity levels in the Mission drag decay section of the sheet
- use `eclipseTimes` for eclipse lengths, the 60 s output grid overstates each one by about a minute
- after changing the model run `dragCheck(30)`, takes about a minute, should end with `checks passed: 6 of 6` (5 of 6 without `numericalCrossCheck_365d.mat`, not on GitHub as it's too big)
- `decayCheck(30)`: atmosnrlmsise00 against pymsis densities, then dragDecay against ode45 with NRLMSISE-00 at every step for the three activity levels, slow as ode45 calls NRLMSISE-00 every step, should end with `checks passed: 4 of 4`

Needs MATLAB R2019a or newer. Demos, checks and NRLMSISE-00 drag need the Aerospace Toolbox.

## Where everything is

- `assumptions.xlsx`: all inputs (Inputs sheet), ground stations (Stations sheet), drag and power assumptions, check results, sources
- `odinPath.m`: path setup
- `RUN.m`: one run of everything, contact table in `results/` (kept local)
- `model/`: the model itself (config, propagate, eom (gravity, J2, drag), oe2rv, sunVector, shadow, shadowGeometry, eclipseTimes, solarFlux, dragTable, dragDensity, dragDecay, gcrfRotation, odinToScenario, constellation, stations, earthRotation, stationVector, contactTimes, downlinkTimes)
- `demos/`: main_demo, scenarioDemo (3D viewer), runSatViewer (one day with the stations in the Stations sheet), plotScale
- `power/`: arrayPower, array sizing and incidence sweeps, batteryAging, battery capacity fade and resistance growth over the mission (NREL BLAST-Lite cell models), powerModes and modeStep, battery charge and operating mode, slewFrequency, observation windows the battery can afford, panelStep, arrayPower's plate and cell equations one step at a time (change both together)
- `studies/orbitTrade/`: pre-DR1 orbit trade (cases A to F), still used as test orbits by the checks
- `validation/`: checks against independent references, logs and plots in `validation/results` (.mat files stay local, too big)
- `external/j77/`: Jacchia 1977 atmosphere code, third party, licence in `license.txt`
