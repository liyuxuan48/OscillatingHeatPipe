using Documenter

makedocs(
    sitename = "OscillatingHeatPipe",
    format = Documenter.HTML(
        size_threshold = 10 * 1024 * 1024,
        size_threshold_warn = 5 * 1024 * 1024,
    ),
    # modules = [OscillatingHeatPipe],
    pages = [
        "Home" => "index.md",
        "Simulation" => "manual/OHP simulation.md",
        "Post-Processing" => "manual/PostProcessing-oneresult.md",
        "OHP-DIY" => "manual/OHP DIY.md"
    ]
)

if get(ENV, "CI", "false") == "true"
    deploydocs(
        repo = "github.com/liyuxuan48/OscillatingHeatPipe.jl.git",
        devbranch = "main",
        target = "build",
        deps = nothing,
        make = nothing,
    )
end
