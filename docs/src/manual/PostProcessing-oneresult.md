```@meta
EditURL = "../../../test/literate/PostProcessing-oneresult.jl"
```

  # PostProcessing
  This notebook is initially designed for PostProcessing ASETS-II numerical and experimental data.
  It can also be used for other configurations's numerical results. But if you want to compare with other experimental data with a different format than ASETS-II. You should find a way to read them to RTD and RTDt in this notebook.

````@example PostProcessing-oneresult
using ComputationalHeatTransfer
using Plots
using Interact
````

  # Read simulation data

````@example PostProcessing-oneresult
read_path = "../numedata/solution.jld2"
SimuResult = load(read_path)["SimulationResult"];
nothing #hide
````

  ### get time array

````@example PostProcessing-oneresult
t = SimuResult.tube_hist_t;
nothing #hide
````

  # Plot 2D graphs

  ### film and slug dynamics

````@example PostProcessing-oneresult
@gif for i in eachindex(t)
    plot(OHPSlug(),i,SimuResult)
end
````

  ### plate T [K]

````@example PostProcessing-oneresult
@gif for i in eachindex(t)
    plot(OHPTemp(),i,SimuResult,clim=(291.2,294.0))
end
````

  ### 2D superheat

````@example PostProcessing-oneresult
@gif for i in eachindex(t)
    plot(OHPSuper(),i,SimuResult)
end
````

  ### 2D pressure

````@example PostProcessing-oneresult
@gif for i in eachindex(t)
    plot(OHPPres(),i,SimuResult)
end
````

  # Plot 2D interpolated curves

  ### Interpolate 2D T data from the plate for fixed sensors on the plate

  place the 2D sensors

````@example PostProcessing-oneresult
x2Dsensors = [-2.75,-1.4,-0.8,0.0,0.0,0.8,1.4,2.75] .* inches
y2Dsensors = [0.0,   0.0, 0.0,0.0,0.4,0.0,0.0,0.0] .* inches
plate_sensors = (x2Dsensors,y2Dsensors);
nothing #hide
````

  get the curve

````@example PostProcessing-oneresult
t_hist,g_hist = getTcurve(plate_sensors,SimuResult);
nothing #hide
````

  ### Read experiment T data
  This part can be customized as long as you can get a matrix for RTD (sensor data) and an array for RTDt (time points)

````@example PostProcessing-oneresult
import XLSX
````

 read experiment file for ASETS-II data.

````@example PostProcessing-oneresult
expfile = expfileDict["O002_H001_P040"]
exppath = "../expdata/"
xf = XLSX.readxlsx(exppath*expfile);
nothing #hide
````

  get experiment data for ASETS-II data

````@example PostProcessing-oneresult
Onum, Hnum, power_exp = getconfig(expfile)
RTDt,RTD = getRTD(xf,Onum);
nothing #hide
````

  ### 2D interpolated temperature curve at fixed sensors

````@example PostProcessing-oneresult
RTD_for_plotting = [1,4,8];
nothing #hide
````

plot OHP

````@example PostProcessing-oneresult
plot(OHP(),SimuResult)
````

plot sensors

````@example PostProcessing-oneresult
scatter!(x2Dsensors[RTD_for_plotting],y2Dsensors[RTD_for_plotting])
annotate!(x2Dsensors[RTD_for_plotting], y2Dsensors[RTD_for_plotting].-0.005, RTD_for_plotting)
````

plot temperature curve

````@example PostProcessing-oneresult
plot(OHPTcurve(),RTD_for_plotting,(t_hist,g_hist),SimuResult)
plot!(OHPTexp() ,RTD_for_plotting,(RTDt,RTD)     ,SimuResult)
````

  ### 2D interpolated thermal conductance

````@example PostProcessing-oneresult
ihot = 4 # hot sensor  for calculating thermal conductance
icold = 8 # cold sensor  for calculating thermal conductance;
````

plot them separately

````@example PostProcessing-oneresult
plot(OHPCond(),(ihot,icold),(t_hist,g_hist),(RTDt,RTD),SimuResult)
````

  ### Liquid slug velocity statistics

fix title and ylabel and legend

````@example PostProcessing-oneresult
plot(OHPV(), SimuResult::SimulationResult,ylimit=(-2,2))
````

  # Plot 1D interpolated curves

````@example PostProcessing-oneresult
@manipulate for i in 1:1:length(t)
    plot(OHP1DT(),i,SimuResult,xlim=(1,2))
    plot!(twinx(),OHPTwall(),i,SimuResult,xlim=(1,2))
#     plot!(twinx(),OHP1DΔT(),i,SimuResult,xlim=(1,2))
#     plot!(twinx(),OHP1DP(),i,SimuResult,xlim=(1,2))
end
````

  ### 1D sensor selector

````@example PostProcessing-oneresult
L = SimuResult.integrator_tube.p.tube.L
@manipulate for ξ in 0:1e-3:L
    plot(OHP(),SimuResult) # plot the ohp layout

    xprobe,yprobe = oneDtwoDtransform(ξ,SimuResult)
    scatter!([xprobe],[yprobe])
end
````

  # Plot 1D property curve for a fixed location sensor

````@example PostProcessing-oneresult
xsensors1D = [2.097, 3.0, 4,4.1]

θhist1D,phist1D = get1DTandP(xsensors1D, SimuResult);

plot(t,θhist1D,label=string.("ξ=", xsensors1D'),xlabel="time [s]", ylabel="temperature [K]")
````

  ### get boiling data (if there are any)

````@example PostProcessing-oneresult
if length(SimuResult.boil_hist) != 0
boil_data,boil_num_x,boil_num_t,t_boil,x2D_boil,y2D_boil,boil_dt = get_boil_matrix(SimuResult::SimulationResult);
end
````

boiling frequency scatter graph (if there are any)

````@example PostProcessing-oneresult
plt = plot()
if length(SimuResult.boil_hist) != 0
plt = plot(OHP(),SimuResult)
scatter!(x2D_boil,y2D_boil,
    colorbar=true,markeralpha=delta.(boil_num_x),colorbar_title="\n boiling frequency [Hz]",right_margin=3Plots.mm,marker_z=boil_num_x./SimuResult.tube_hist_t[end],markerstrokewidth=0,markercolor=cgrad(:greys, rev = true))
end
plt
````

boiling frequency curve (if there are any)

````@example PostProcessing-oneresult
plt = plot()
if length(SimuResult.boil_hist) != 0
plt = plot(t_boil,boil_num_t./boil_dt,
color=:orange, legend=:topleft, ylabel="f [HZ]",xlabel="time [s]", label="overall boiling frequency")
end
plt
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

