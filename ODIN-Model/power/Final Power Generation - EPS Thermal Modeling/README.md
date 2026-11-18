# Final Power Generation / EPS Thermal Modeling

An additive MATLAB package for the updated solar-array, battery and PMAD thermal analysis. The MATLAB namespace is `finalPowerEPS` (the folder name must remain `+finalPowerEPS`). It uses ODIN orbit/configuration results without replacing the team's functions.

## Use with ODIN

Upload this entire new folder to **`ODIN-Model/power/Final Power Generation - EPS Thermal Modeling`**. Do not replace, move or edit any existing ODIN file. In MATLAB, start in `ODIN-Model`, initialize ODIN's normal path, and add this package's parent folder explicitly:

```matlab
odinPath
packageRoot = fullfile(pwd,'power','Final Power Generation - EPS Thermal Modeling');
addpath(packageRoot)

cfg = config();                 % ODIN's existing workbook-driven configuration
out = propagate(cfg,10);         % ODIN's existing dynamics, including enabled drag/SRP
result = finalPowerEPS.run(cfg,out);
```

Use the orbit output you already calculated rather than propagating it again:

```matlab
options = struct('makePlots',false,'refineOrbit',false, ...
    'includePowerBalance',false);
result = finalPowerEPS.run(cfg,out,options);
```

The supported integration entry point is `finalPowerEPS.run`. `result.array`, `result.thermal`, and `result.iv.nominal` / `.peak` contain the updated analysis. `result.orbit` contains the separate analysis ephemeris and `result.mapping` describes the configuration mapping. The source `cfg` and `out` are unchanged.

For use alongside the existing standalone MATLAB folder, add the new package parent there and supply that standalone run's `cfg` and `out`:

```matlab
addpath(fullfile(pwd,'Final Power Generation - EPS Thermal Modeling'))
result = finalPowerEPS.run(cfg,out);
```

That uses the supplied standalone scenario. Use the ODIN setup above when you want the team's workbook configuration and dynamics.

`options.arrayOverrides`, `options.componentOverrides`, `options.ivOverrides`, and `options.powerBalanceOverrides` supply explicit per-run overrides. Keep their fields consistent with the corresponding package function's parameter checks. `options.showSummary` controls printing, and `options.refinementOptions` accepts `searchStep_s` / `transitionStep_s` for the private eclipse grid. Shared defaults are in `+finalPowerEPS/modelSettings.m`.

The integration entry point disables the full incidence sweep by default to avoid silently repeating an expensive study. Enable it explicitly with `options.arrayOverrides.runIncidenceSweep=true`; with `makePlots=true`, this produces the package's original incidence-sweep figures in addition to its other figures. The direct package array function and standalone workflow retain their own default controls. No existing ODIN figure is edited or replaced.

The coupled integration run uses `power_consistent_proxy` I-V mode so its curve power remains consistent with the solved array/thermal history. `options.ivOverrides.mode='xtj_reference'` is rejected by that entry point. For a separate counterfactual Spectrolab reference-cell study, call `finalPowerEPS.solarArrayIV` directly in `xtj_reference` mode with reference irradiance/temperature inputs and **without** an `availablePower_W` field. That independent curve is not a NanoAvionics assembly prediction or the coupled EPS result.

Start from the repository folder when constructing the team cfg/out: a current MATLAB folder containing another `config.m` takes precedence over added paths. `which config -all` and `which propagate -all` should identify the ODIN copies for the calls above. Do not add the old flat standalone power-code folder or the review snapshot to ODIN's path. Add only this package's parent folder, not `+finalPowerEPS` itself; do not use `genpath`. ODIN's existing `odinPath` adds `power`, but does not automatically add this nested package parent.

## How configuration is combined

The integration maps the workbook's selected panel area (`cfg.power.panelArea`) and complete five-factor end-of-life allowance into the updated array analysis when provided. The five factors are `eolRadiation`, `eolUV`, `eolDebris`, `eolThermalCycling`, and `eolContamination`. Explicit per-run array overrides take precedence and the mapping is reported.

The retained 64-cell geometry/coverage proxy, cell efficiency and temperature coefficients, rear emissivity, laminate properties, 127 C initial temperatures, battery layout/resistance/OCV assumptions, and component schedules are not silently replaced with different workbook fields. ODIN's battery/mode models and this package's component dispatch remain different, explicitly identified analyses. Scheduled communications in the component model are assumptions, not the result of ODIN's ground-station contact calculation.

The optional `includePowerBalance` study remains the existing independent worst-shadow envelope. It is not the propagated mission orbit or a replacement for ODIN's battery aging study.

## Private eclipse and environment analysis

When refinement is enabled, `finalPowerEPS.refineOrbit` inserts eclipse contacts and penumbra samples in a private ephemeris copy. It uses cubic Hermite interpolation of the team's supplied positions and velocities; original r/v knots are retained exactly. It does not substitute dynamics, alter the team propagator, or propagate a missing endpoint.

The final analysis time is the exact **supplied** final time. If the source stops short of a requested duration and records `out.requestedDuration_s`, the refinement report exposes that shortfall rather than extrapolating. Otherwise it records that the requested duration is unknown. A finer power-analysis grid does not improve the underlying orbital dynamics accuracy. Very short grazing eclipses require adequate source/search resolution.

The separate standalone model's updated propagator can include the exact requested end while it is solving its own dynamics. This ODIN package does not substitute that propagator, re-propagate ODIN's orbit, or extend an incomplete team ephemeris. It preserves the supplied ODIN endpoint and reports its scope.

Known ephemeris fields are rebuilt consistently: time, r/v, angular momentum, Sun vector/direction, illumination and eclipse events. Supplied panel normals are linearly interpolated and normalized between original samples; antipodal/ambiguous interpolation is rejected. Unrecognized numeric sampled fields are omitted from the refined copy and listed in the refinement information, so they cannot silently retain an incompatible old time dimension. Static metadata is preserved.

Earth IR/albedo, two-face radiation, through-thickness conduction and the original EPS loss boundaries are the previously validated model. Environmental averages remain explicit assumptions, and unsupported mounting conductance remains zero. The copied physics functions use qualified package calls; no MATLAB-path switching is used inside the analysis.

## Existing ODIN behavior stays intact

No existing repository file needs replacement. In particular, do not upload our standalone `config`, `propagate`, `eom`, `shadow`, `eclipseTimes`, `sunVector`, `solarFlux`, `main_demo` or `plotScale` over their team counterparts.

ODIN's original `RUN`, `MISSION`, `UNCERTAINTY`, `powerModes`, `panelStep`, `modeStep`, `slewFrequency`, battery aging, contacts and orbit validation remain as they were. Their legacy power calculations keep their existing equations. This package's updated results are additional outputs; calling it does not redirect those original routines to the new thermal model.

The standalone flat power model and its existing figures remain available separately. Namespacing prevents this package's helpers from shadowing them or the team functions. The optional additional diagnostics are labeled with scenario context.

## Validation and upload boundary

The release's [validation report](VALIDATION_REPORT.md) records the actual native-MATLAB checks and reproduction steps. The tests are in `tests/testFinalPowerEPS.m`. The source manifest records provenance; `UPLOAD_MANIFEST.json` and `UPLOAD_ALLOWLIST.txt` list the exact new package files and their repository destinations. The team repository files are not part of the upload set.

Before uploading, compare the proposed file list with the allowlist and confirm that no protected team file appears in the diff. Retain the repository's required branch naming, sequential push-number commit messages, and reviewed integration workflow. The reviewed repository README referred to `workingMain`, but that branch was not visible during review; resolve the publishing base with the team before creating a merge request. Adding this package locally does not publish or merge anything.
