module OscillatingHeatPipe

export inches,gravity,expfileDict

using CoolProp
using Interpolations
using LinearAlgebra
using RecipesBase
using UnPack
using CartesianGrids
# using Plots
using Statistics
using SparseArrays 
using Reexport
using Revise


@reexport using ComputationalHeatTransfer
import RigidBodyTools: arccoord,arclength

const inches = 2.54e-2; 
const gravity = 9.8;

include("utils/Systems.jl")
include("utils/CoolProp.jl")
include("utils/Preprocessing.jl")
include("utils/Mapping.jl")
include("utils/Postprocessing.jl")
include("utils/ohp.jl")
include("utils/Heatercondenser.jl")
include("utils/Plotrecipe.jl")

include("solver/Thermomodel.jl")
include("solver/Tools.jl")
include("solver/Timemarching.jl")

include("solver/callback/boiling.jl")
include("solver/callback/vapormerging.jl")
include("solver/callback/liquidmerging.jl")
include("solver/callback/fixdx.jl")
include("solver/callback/slugbc.jl")

expfileDict = Dict([
    ("O001_H002_P010","20190607_F_PD_%23013_O001_H002_P010_expA.xlsx"),
    ("O001_H002_P020","20190608_F_PD_%23014_O001_H002_P020_expA.xlsx"),
    ("O001_H002_P030","20190614_F_PD_%23015_O001_H002_P030_expA.xlsx"),
    ("O001_H002_P040","20190617_F_PD_%23016_O001_H002_P040_expA.xlsx"),
    ("O001_H001_P010","20190604_F_PD_%23001_O001_H001_P010_expA.xlsx"),
    ("O001_H001_P020","20190606_F_PD_%23002_O001_H001_P020_expA.xlsx"),
    ("O001_H001_P030","20190612_F_PD_%23003_O001_H001_P030_expA.xlsx"),
    ("O001_H001_P040","20190613_F_PD_%23004_O001_H001_P040_expA.xlsx"),
    ("O002_H002_P010","20190607_F_PD_%23017_O002_H002_P010_expA.xlsx"),
    ("O002_H002_P020","20190608_F_PD_%23018_O002_H002_P020_expA.xlsx"),
    ("O002_H002_P030","20190614_F_PD_%23019_O002_H002_P030_expA.xlsx"),
    ("O002_H002_P040","20190617_F_PD_%23020_O002_H002_P040_expA.xlsx"),
    ("O002_H001_P010","20190604_F_PD_%23005_O002_H001_P010_expA.xlsx"),
    ("O002_H001_P020","20190606_F_PD_%23006_O002_H001_P020_expA.xlsx"),
    ("O002_H001_P030","20190612_F_PD_%23007_O002_H001_P030_expA.xlsx"),
    ("O002_H001_P040","20190613_F_PD_%23008_O002_H001_P040_expA.xlsx"),
    ("O003_H001_P040","20190613_F_PD_%23012_O003_H001_P040_expA.xlsx"),
    ("O003_H002_P010","20190607_F_PD_%23021_O003_H002_P010_expA.xlsx"),
    ("O003_H002_P020","20190608_F_PD_%23022_O003_H002_P020_expA.xlsx"),
    ("O003_H002_P030","20190614_F_PD_%23023_O003_H002_P030_expA.xlsx"),
    ("O003_H002_P040","20190617_F_PD_%23024_O003_H002_P040_expA.xlsx")
        ]);

end