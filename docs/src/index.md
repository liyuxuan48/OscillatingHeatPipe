# OscillatingHeatPipe.jl

`OscillatingHeatPipe.jl` simulates an oscillating heat pipe (OHP) coupled to a two-dimensional
heat-conduction plate. It provides tools for defining OHP channels, heaters, and condensers,
initializing the tube and plate systems, time-marching the coupled problem, saving results, and
post-processing temperatures, pressure, slug motion, and boiling statistics.

The package builds on [`ComputationalHeatTransfer.jl`](https://github.com/JuliaIBPM/ComputationalHeatTransfer.jl)
from the [JuliaIBPM](https://github.com/JuliaIBPM) ecosystem.

The tutorials in this site are generated from runnable scripts in `test/literate/`, so the documentation and examples stay in sync.

## Tutorials

- [Simulation](manual/OHP simulation.md): build and run a compact coupled OHP simulation.
- [OHP DIY](manual/OHP DIY.md): customize the plate, heaters, condensers, and OHP channel.
- [Post-Processing](manual/PostProcessing-oneresult.md): load saved results and create diagnostic plots.

## Local Use

The fastest way to learn the package is to run the three example notebooks in Jupyter, in order.
This workflow takes you from nothing to a working simulation.

### Step 0 — Prerequisites

You need **Julia 1.11 or newer**. Download it from [julialang.org/downloads](https://julialang.org/downloads/)
and make sure `julia` runs in your terminal:

```bash
julia --version
```

### Step 1 — Install the package and notebook dependencies

Start Julia and run:

```julia
using Pkg
Pkg.add(url="https://github.com/liyuxuan48/OscillatingHeatPipe.jl.git")
Pkg.add(["Plots", "JLD2", "XLSX", "ProgressMeter", "IJulia"])
```

`IJulia` lets you run the notebooks in Jupyter. The first install compiles a lot of code and can
take several minutes — this only happens once.

### Step 2 — Set up your work folder

Make a work directory anywhere you like on your computer — for example `ohp-work` in your home folder.
You will run the notebooks from there instead of editing the package's own copies, which `Pkg.add`
tucked away deep inside Julia's package depot.

The notebooks read and write data through the relative path `../numedata/…`, so they need to live in a
subfolder next to a `numedata/` folder. The commands below build that layout and copy the example
notebooks out of the installed package for you. Run them in the same Julia session:

```julia
using OscillatingHeatPipe

work  = joinpath(homedir(), "ohp-work")          # change this to wherever you want your work folder
nbdir = joinpath(work, "notebooks")
mkpath(nbdir)
mkpath(joinpath(work, "numedata"))               # the notebooks save their results here

src = joinpath(pkgdir(OscillatingHeatPipe), "examples")
for nb in filter(f -> endswith(f, ".ipynb"), readdir(src))
    cp(joinpath(src, nb), joinpath(nbdir, nb); force = true)
end
```

This creates the following layout:

```text
ohp-work/
├── notebooks/   ← your editable notebook copies (run Jupyter here)
└── numedata/    ← simulation results are saved here
```

Re-run the loop any time you want to refresh your copies from the installed package. Up to this step, you only need to do once. In the future, you just need to start at Step 3.

### Step 3 — Launch Jupyter

```julia
using IJulia
notebook()
```

This opens Jupyter in your browser. If you run this for the first time you need to install several packages by following its instructions.

### Step 4 — Run the notebooks in order

Open each notebook and run all cells (menu: **Run ▸ Run All Cells**).

1. **`OHP simulation.ipynb`** — a compact coupled OHP simulation. It writes the result to
   `../numedata/example.jld2`.
2. **`PostProcessing-oneresult.ipynb`** — loads `../numedata/example.jld2` and produces the
   temperature, pressure, slug-motion, and boiling plots. **Run this after the simulation**, because
   it needs the saved result file.
3. **`OHP DIY.ipynb`** — optional. Shows how to build your own geometry, heaters, condensers, and
   channel; it writes its own result to `../numedata/DIY.jld2`.

That's the whole workflow. Edit the parameters in `OHP simulation.ipynb`, re-run it, then re-run
`PostProcessing-oneresult.ipynb` to see the effect.

## Citation

If you use this package, please cite the associated paper:

Li, Y., Eldredge, J. D., Lavine, A. S., Fisher, T. S., & Drolen, B. L. (2024).
A conjugate heat transfer model of oscillating heat pipe dynamics, performance, and dryout.
*International Journal of Heat and Mass Transfer, 227*, 125530.
[https://doi.org/10.1016/j.ijheatmasstransfer.2024.125530](https://doi.org/10.1016/j.ijheatmasstransfer.2024.125530)
