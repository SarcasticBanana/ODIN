# ODIN-Model

Orbit, solar irradiance and drag model for the ODIN space debris science constellation, CU Boulder ASEN 4018 Senior Design.

The model propagates an orbit under two-body gravity plus J2 and, when switched on, atmospheric drag in the Jacchia 1977 atmosphere. It computes the Earth shadow and returns the solar irradiance at the spacecraft, finds exact eclipse entry and exit times, and can hand the orbit to MATLAB's `satelliteScenario` for ground contacts and a 3D viewer. Irradiance is panel independent: it is the flux on a surface held normal to the Sun, attenuated by the shadow. Turning it into electrical power needs an attitude model, which is not yet included.

Every input comes from one spreadsheet, `assumptions.xlsx`, so changing the orbit, the spacecraft or a solver setting does not mean editing code.

## Quick start

Open MATLAB in this folder, then:

```matlab
odinPath                         % adds the model folders to the path, once per session
cfg = config();                  % every input, read from assumptions.xlsx
out = propagate(cfg, 15);        % 15 orbital periods
S = solarFlux(out, cfg);         % irradiance at the spacecraft [W/m^2]
ecl = eclipseTimes(out, cfg);    % umbra and penumbra entry and exit times and durations [s]
main_demo                        % example plots
```

Requirements: MATLAB R2019a or later for the model (`readcell`), no toolboxes. The 3D viewer demos and the validation scripts use the Aerospace Toolbox; `numericalCrossCheck` also uses the Parallel Computing Toolbox when it is installed. Everything was last run on R2026b.

## Inputs: `assumptions.xlsx`

The **Inputs** sheet holds every model input, one row per field of the parameter struct `cfg`. `config.m` reads every row whose `config.m field` column starts with `cfg.` and builds the struct from them.

- To change an input, change its Value and save the workbook before running. MATLAB reads the saved file.
- Values must be in the units of the Unit column. The code does not convert units.
- Values can be numbers, `TRUE`/`FALSE`, or text. The epoch is text in the form `yyyy-mm-dd HH:MM:SS` (UTC).
- To add an input, add a row with a new `cfg.` field name; it appears in `cfg` with no code change. Rows under **Derived** are formulas for reference and are not read.
- `cfg = config('myCase.xlsx')` reads another copy of the workbook, so a team can keep its own case without changing the shared one.
- For a one-off change in a script, overwrite the field after reading it, for example `cfg = config(); cfg.h = 600;`.

The other sheets document the drag model: **Assumptions** lists 25 assumptions with their value, basis, effect on results and how to relax each one; **Checks** holds the drag validation numbers; **Sources** lists the references.

## Layout

| Path | Contents |
|---|---|
| `assumptions.xlsx` | Model inputs and the drag assumption list |
| `odinPath.m` | Adds the folders below to the MATLAB path |
| `model/` | The model itself, listed below |
| `external/j77/` | Jacchia 1977 atmosphere code `j77sri.m` with its reference output and licence, unchanged |
| `demos/` | `main_demo` (plots and the array power example), `scenarioDemo` (3D viewer of the trade orbits for meetings), `runSatViewer` (one day with a ground station), `plotScale` (plot helper for `main_demo`) |
| `power/` | `arrayPower`: array area sizing and incidence sweeps with a one-sided transient plate temperature, Power subteam |
| `studies/orbitTrade/` | Pre-DR1 orbit trade, cases A to F: `tradeLighting`, `plotTradeLighting`, `run_trade`, `plotsdr1` and its DR1 figures. The validation scripts use these cases as test orbits. |
| `validation/` | Checks of the model against independent references, with logs and figures in `validation/results` |

### `model/`

| File | Purpose |
|---|---|
| `config.m` | Reads the Inputs sheet of `assumptions.xlsx` and returns `cfg` |
| `propagate.m` | `ode45` integration on a uniform output grid; returns the state history struct `out` |
| `eom.m` | State derivative: two-body gravity, J2 (`accelJ2`) and drag (`accelDrag`) |
| `oe2rv.m` | Classical orbital elements to ECI position and velocity |
| `dragTable.m` | J77 mass density against height for the exospheric temperature set by `cfg.drag.F107` |
| `sunVector.m` | Low precision Sun ephemeris, mean equator and equinox of date |
| `shadowGeometry.m` | Apparent Sun and Earth disk radii and their separation, shared by the shadow routines |
| `shadow.m` | Visible fraction of the solar disk, conical or cylindrical shadow |
| `eclipseTimes.m` | Exact umbra and penumbra entry and exit times, and durations, by root finding |
| `solarFlux.m` | Irradiance in W/m^2, panel independent |
| `gcrfRotation.m` | Rotation from the mean equator and equinox of the epoch to GCRF |
| `odinToScenario.m` | Adds a `propagate.m` trajectory to a `satelliteScenario` as an ephemeris satellite in GCRF |

### Output struct `out`

| Field | Contents |
|---|---|
| `t` | 1 x N output times since epoch [s] |
| `r`, `v` | 3 x N ECI position [km] and velocity [km/s] |
| `T`, `a` | Orbital period [s] and semi-major axis [km] from the initial elements |
| `hvec` | 3 x N specific angular momentum [km^2/s] |
| `rsun`, `sunhat` | 3 x N Sun position [km] and unit vector (only when `cfg.use.eclipse` is true) |
| `illum` | 1 x N visible fraction of the solar disk, 0 to 1 (only when `cfg.use.eclipse` is true) |
| `cfg` | Copy of the parameter struct used for the run |

## Drag

With `cfg.use.drag` set to `TRUE`, `propagate.m` builds a table of J77 density from 90 to 2500 km once per run and `eom.m` adds the drag acceleration. The exospheric temperature comes from `cfg.drag.F107`, held constant for the run, through J77 equation 20 at Kp = 0. The atmosphere rotates rigidly with the Earth and the height is measured above the WGS84 ellipsoid. The spacecraft is described by `cfg.drag.mass`, `cfg.drag.area` and `cfg.drag.Cd`, which are placeholders until Structures has numbers. Drag is off by default, so results without it are unchanged. Below 90 km the run stops with an error rather than extrapolating.

The Assumptions sheet of `assumptions.xlsx` lists what this leaves out and how much each omission matters.

## Validation

Each script writes a log to `validation/results`. Large outputs (`*.mat`, `SW-All.csv`) stay local and are not pushed.

| Script | What it checks | Result |
|---|---|---|
| `dragCheck` | J77 exospheric temperature against J77 Table 5; `j77sri.m` against its reference output; J77 density against the NRLMSISE-00 global mean; drag acceleration units; drag-free run against the earlier J2 run; decay against an orbit-averaged rate and the toolbox high-fidelity propagator | 6 of 6 pass. One year at 650 km, F10.7 150, m/(Cd A) = 30 kg/m^2: drag lowers the orbit-averaged altitude by 4.671 km, against 4.649 km from the orbit-averaged rate and 4.864 km from the high-fidelity run with NRLMSISE-00 and EGM2008 |
| `eclipseCrossCheck` | `eclipseTimes` against the Aerospace Toolbox eclipse intervals, trade case B | 29 of 29 shadows matched; entry within 0.92 s and exit within 0.74 s at 1 s scenario sampling |
| `numericalCrossCheck` | The J2 model against the toolbox numerical propagator with EGM2008 120 x 120, with and without drag, one year, six trade cases | Time in full Sun moves by at most 0.062 percentage points and the longest shadow by at most 0.27 min |
| `numericalSettings` | One day runtime and accuracy of toolbox propagator settings | Toolbox J2 against the model J2: 209 to 305 m after one day |
| `frameCorrection` | Size of the mean of epoch to GCRF rotation applied by `odinToScenario` | 0.373 deg rotation; up to 48 km in position and 41.5 km on the ground |
| `groundContacts` | Ground station contacts, a first UHF link budget with a placeholder radio, and contact time inside shadow | 3.7 to 4.5 contacts and 26 to 42 min per day from Boulder at a 10 deg mask |

`dragCheck(30)` runs in about a minute and is a quick regression check after a change to the model; the last line should read `checks passed: 6 of 6`. `dragCheck` without an argument runs the full year in about ten minutes. `numericalCrossCheck` takes several hours, saves checkpoints in the temporary folder and resumes if interrupted; its CelesTrak runs need `SW-All.csv` from CelesTrak in `validation/results` and are skipped without it.

## Notes on use

- Use `eclipseTimes` for any eclipse duration, including worst case battery sizing. Durations counted from the 60 s output grid are biased by about one sample (+59 s, +3.2 percent). The mean illuminated fraction from the grid is accurate to 0.01 percent.
- `eclipseTimes` always uses the conical umbra and penumbra boundaries, regardless of `cfg.eclipse.model`.
- For accuracy checks under J2 without drag, the conserved quantities are total energy including the J2 potential term, and the z component of angular momentum `hz`. Two-body energy and `|h|` are not conserved under J2, and nothing is conserved with drag on.
- `ode45` global error scales as tol^0.975. At 1e-11 the position error is about 1.3 km after six months, about 0.2 s of along-track timing.
- The integrated nodal rate at 98.6 deg and 800 km is 0.98949 deg/day against the Sun's 0.98565 deg/day, so the beta angle drifts about 1.4 deg per year. An inclination of about 98.57 deg nulls this drift in this model.
- Drag slows the integration: a year of case A took 489 s with drag and 67 s without on R2026b.

## Status, October 2026

Working, with the checks in Validation:

- Orbit propagation with two-body gravity and J2, Sun position, conical and cylindrical shadow, exact eclipse times and irradiance.
- Atmospheric drag in a static J77 atmosphere, off by default.
- Export to `satelliteScenario` in GCRF for ground contacts and the 3D viewer.

Also in place:

- All inputs read from `assumptions.xlsx`.
- Array power with a transient plate temperature from the Power subteam in `power/arrayPower.m`, not yet covered by the validation scripts.

Not modelled yet:

- Attitude. Irradiance is on a Sun-facing surface; electrical power needs the panel orientation.
- Earth albedo and infrared, solar radiation pressure (the `cfg.use.srp` switch does nothing yet), third-body gravity and gravity harmonics beyond J2.
- Battery state of charge through the orbit.
- Time-varying solar and geomagnetic activity for drag.

Inputs that are not final:

- Spacecraft mass, area and drag coefficient are placeholders.
- The default orbit is still 800 km at 98.6 deg. DR1 settled a Sun-synchronous orbit at 550 to 650 km with RAAN still open.
- The six trade cases in `studies/orbitTrade` are the pre-DR1 candidates; they stay because the validation scripts use them as test orbits.

## Next steps

1. Attitude model, then power from the double-deployable arrays on that attitude, starting from `arrayPower`.
2. Battery state of charge over the orbit, sized against the worst eclipse at end of life.
3. Set the baseline orbit in `assumptions.xlsx` once RAAN is chosen, and add it as a test orbit next to the trade cases.
4. Replace the spacecraft placeholders with the Structures mass and area, and tie the drag area to attitude.
5. Time-varying F10.7 and Ap from the CelesTrak file for decay and lifetime, through the J77 time-dependent terms or NRLMSISE-00.
6. Speed up the drag evaluation (489 s against 67 s for a year of case A).
7. Albedo, infrared and solar radiation pressure, if the power or lifetime margins turn out to depend on them.

## Working on the model

- Put new inputs in `assumptions.xlsx`, not in code.
- MATLAB comment convention: one line per function stating what it is, then its inputs and outputs. Every variable annotated at its first appearance only, as `name: meaning [units]`, with units in brackets only where they exist and nothing for dimensionless quantities. A single space before the trailing `%`. A short line naming the relation immediately above the line implementing it, with symbols only where the code does not already show them. No explanatory or narrative comments.
- Team repository naming: branches `Subteam(s)_function_main#.workingMain#`, commits `Branchname.push#: change description`.

## Third-party code

`external/j77/j77sri.m` is the Jacchia 1977 code by D. L. Huestis, SRI International (1999, 2002), in the MATLAB version by M. Mahooti from MathWorks File Exchange 54912. The File Exchange package is under the BSD licence in `external/j77/license.txt`; the SRI header in `j77sri.m` reads "All Rights Reserved" and "provided on an as is basis". Check that redistribution is allowed before this repository is made public. The File Exchange readme is `external/j77/README.md`.
