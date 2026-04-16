# Ocean Flow Lab

Julia + Oceananigans + Dash app for launching several idealized ocean-flow examples, saving each run to NetCDF, and reopening completed runs for animation.

## Included configurations

- Barotropic jet
- Baroclinic jet
- Dipole vortex (two-eddy system)
- Double gyre
- Zonal jet instability
- Vortex merger
- Isolated eddy propagation

## Run

```bash
source ~/.bashrc
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. run_app.jl
```

The dashboard opens on `http://127.0.0.1:8050` by default.

## Interface

- Dropdown to select the flow configuration
- Editable controls for `dt`, stop time, output cadence, and domain size / resolution
- GPU-first execution with CPU fallback
- Background run monitor showing which models are still executing
- Unique run folders under `runs/<flow>_<timestamp>/`
- NetCDF output containing raw model fields plus derived diagnostics
- Completed-run viewer with variable selection, timestep playback, and vertical-level control
- Diagnostics reloaded from the saved NetCDF file

## Notes

- The examples are intentionally idealized so the dashboard stays interactive.
- The app requests GPU execution by default and falls back to CPU only if CUDA is unavailable.
- The baroclinic case writes buoyancy along with the velocity fields; the other examples still save all velocity components and derived outputs such as depth-mean vorticity.
- The main viewer no longer follows live model state. It opens only completed NetCDF outputs after a run has finished.
