##using TestSetExtensions
using OscillatingHeatPipe
using Literate
using Test

const GROUP = get(ENV, "GROUP", "All")

ENV["GKSwstype"] = "nul" # removes GKS warnings during plotting

macro mysafetestset(args...)
    name, expr = args
    quote
        ex = quote
          name_str = $$(QuoteNode(name))
          expr_str = $$(QuoteNode(expr))
          mod_name = gensym(name_str)
          ex2 = quote
              @eval module $mod_name
                      using Test
                      @testset $name_str $expr_str
                    end
              nothing
          end
          eval(ex2)
        end
        eval(ex)
    end
end

notebookdir = "../examples"
docdir = "../docs/src/manual"
litdir = "./literate"

notebookdir = "../examples"
docdir = "../docs/src/manual"
litdir = "./literate"

if GROUP == "All" || GROUP == "Auxiliary"
    include("thermomodel.jl")
    include("correlations.jl")
    include("integrator.jl")
    include("datastructures.jl")
end

if GROUP == "All" || GROUP == "Literate"
  for (root, dirs, files) in walkdir(litdir)
    for file in files
      global file_str = "$file"
      global body = :(begin include(joinpath($root,$file)) end)
      endswith(file,".jl") && startswith(file,"OHP DIY.jl") && @mysafetestset file_str body
      #endswith(file,".jl") && @mysafetestset file_str body
    end
  end
end


if GROUP == "Notebooks"
  # Literate's default sandbox is `Module(gensym())`, so each notebook executes in a
  # uniquely-named module (Main.##NNN). Functions saved into a JLD2 file (e.g. the BC/forcing
  # callbacks held by the integrator in `OHP simulation.jl`) get serialized under that module
  # name, and a later notebook (`PostProcessing-oneresult.jl`) cannot resolve them on load.
  # Override the sandbox to a fresh-but-stably-named module so each file still runs clean while
  # JLD2 always stores/resolves the callbacks under `Main.NotebookSandbox.*`.
  @eval Literate function sandbox()
      m = Core.eval(Main, :(module NotebookSandbox end))  # fresh module, stable name Main.NotebookSandbox
      Core.eval(m, :(eval(x) = Core.eval($m, x)))
      Core.eval(m, :(include(x) = Base.include($m, x)))
      return m
  end
  for (root, dirs, files) in walkdir(litdir)
    for file in sort(files)
      # endswith(file,".jl") && startswith(file,"OHP DIY.jl") && Literate.notebook(joinpath(root, file),notebookdir)
      endswith(file,".jl") && Literate.notebook(joinpath(root, file),notebookdir; execute=true)
    end
  end
end

if GROUP == "Documentation"
  for (root, dirs, files) in walkdir(litdir)
    for file in files
      endswith(file,".jl") && Literate.markdown(joinpath(root, file),docdir; execute=false)
    end
  end
end
