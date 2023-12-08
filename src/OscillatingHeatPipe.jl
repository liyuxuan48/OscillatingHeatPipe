module OscillatingHeatPipe

export inches,gravity,expfileDict

using Reexport

@reexport using ComputationalHeatTransfer

const inches = 2.54e-2; 
const gravity = 9.8;

include("intube/CoolProp.jl")
include("intube/Systems.jl")
include("intube/Preprocessing.jl")
include("intube/Thermomodel.jl")
include("intube/Tools.jl")
include("intube/Mapping.jl")
include("intube/Timemarching.jl")
include("intube/Postprocessing.jl")
include("intube/DrawingOHP.jl")
# include("intube/FluidProperty.jl")
include("intube/Heatercondenser.jl")
include("intube/Plotrecipe.jl")
include("intube/callback/boiling.jl")
include("intube/callback/vapormerging.jl")
include("intube/callback/liquidmerging.jl")
include("intube/callback/fixdx.jl")
include("intube/callback/slugbc.jl")

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