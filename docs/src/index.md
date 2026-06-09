# OscillatingHeatPipe.jl

`OscillatingHeatPipe.jl` is a Julia package for coupled oscillating heat pipe and plate heat-transfer simulations.

The tutorials in this site are generated from runnable scripts in `test/literate/`, so the documentation and examples stay in sync.

## Tutorials

- [Simulation](manual/OHP simulation.md): build and run a compact coupled OHP simulation.
- [OHP DIY](manual/OHP DIY.md): customize the plate, heaters, condensers, and OHP channel.
- [Post-Processing](manual/PostProcessing-oneresult.md): load saved results and create diagnostic plots.

## Local Use

```julia
using OscillatingHeatPipe
using Plots
```

Run a sample simulation from the repository root:

```bash
julia --project=@. "test/literate/OHP simulation.jl"
```

Then inspect it:

```bash
julia --project=@. "test/literate/PostProcessing-oneresult.jl"
```
