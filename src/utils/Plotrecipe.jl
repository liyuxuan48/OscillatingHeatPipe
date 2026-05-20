import Plots: mm,palette

export OHP,OHPTemp,OHPSlug,OHPPres,OHPSuper,OHP2DΔT,OHPTexp,OHPTcurve,OHPCond,OHPV,OHP1DT,OHP1DP,OHP1DΔT,OHPTwall

mutable struct OHP end

mutable struct OHPTemp end
mutable struct OHPSlug end
mutable struct OHPPres end
mutable struct OHPSuper end
mutable struct OHP2DΔT end

mutable struct OHPTexp end
mutable struct OHPTcurve end
mutable struct OHPCond end
mutable struct OHPV end

mutable struct OHP1DT end
mutable struct OHP1DP end
mutable struct OHP1DΔT end
mutable struct OHPTwall end

# mutable struct OHP end

@recipe function f(::OHP, SimuResult::SimulationResult)
    grid = SimuResult.integrator_plate.p.base_cache.g
    xlim = grid.xlim[1]
    ylim = grid.xlim[2]

    # ohp = SimuResult.integrator_plate.p.qline[1].body # one body only for now:)
    ohp =  SimuResult.integrator_plate.p.forcing["heating models"][end].shape
    closedornot := SimuResult.integrator_tube.p.tube.closedornot

    aspectratio --> 1
    grid := false

    xlabel --> "x [m]"
    ylabel --> "y [m]"

    xlims := xlim
    ylims := ylim

    fillalpha := 0
    framestyle := :plain

    linecolor := palette([:blue,:red], 2)
    
    ohp
end

@recipe function f(::OHPTemp, i::Int64, SimuResult::SimulationResult; plain=false)
    grid = SimuResult.integrator_plate.p.base_cache.g
    plate_T = SimuResult.plate_T_hist[i]
    tube_hist_t = SimuResult.tube_hist_t

    time := tube_hist_t[i]
    minimal := plain

    OHPTemp(),plate_T,grid
end

@recipe function f(::OHPTemp, plate_T::Nodes, grid::PhysicalGrid; time =:none, minimal=false)

    seriestype --> :heatmap

    if minimal == false
        xlabel --> "x [m]"
        ylabel --> "y [m]"
        legend --> true
        
        xlimit --> grid.xlim[1]
        ylimit --> grid.xlim[2]
        
        colorbar_title --> "\n T[K]"
        right_margin --> 5mm

            if time != :none
            # @series begin
            annotation --> [(0.0, 0.04, string("Time = ", round(time, digits=2), "[s]"))]
            # end
            end
    else
        showaxis := false;
        legend := false;
        colorbar:=false;
    end
    
    plate_T,grid
end

@recipe function f(::OHPSlug, i::Int64, SimuResult::SimulationResult; plain=false)
    adjust = 1e-2;
    
    tube_sys = getcurrentsys(SimuResult,i)
    closedornot := tube_sys.tube.closedornot
    
    tube_hist_t = SimuResult.tube_hist_t
    
    sluginterp = slug_interp(tube_sys)
    xs =  tube_sys.wall.Xarray
    markers  = map(sluginterp, xs)
    
    grid = SimuResult.integrator_plate.p.base_cache.g

    ohp_model = SimuResult.integrator_plate.p.forcing["heating models"][end]
    ohp = ohp_model.transform(ohp_model.shape)
    
    # seriestype --> :heatmap
    xlabel --> "x [m]"
    ylabel --> "y [m]"
    fillalpha := 0
    linewidth := 2.0
    clim      :=(0,2)
    line_z    := markers
    linecolor := palette([:red,  :blue, :yellow])
    bbox_inches := "tight"
    border  := :none

    # legend --> true
    
    xlimit --> grid.xlim[1]
    ylimit --> grid.xlim[2]

    if plain == false

        annotation := [(-0.052+adjust, -0.028, "Dry vapor"),(-0.01+0.007+adjust, -0.028, "Vapor with film"),
        (0.04+0.002+adjust, -0.028, "Liquid"),(0.0, 0.03, string("Time = ", round(tube_hist_t[i], digits=4), "[s]"))]

        @series begin
            # annotation := (0, 0, "Look up!")
            seriestype := :scatter
            color := :red
            [-0.07+adjust],[-0.028]
        end

        @series begin
            seriestype := :scatter
            color := :yellow
            [-0.03+adjust],[-0.028]
        end

        @series begin
            seriestype := :scatter
            color := :blue
            [0.03+adjust],[-0.028]
        end
    else
        annotation := [(0.0, 0.03, string("Time = ", round(tube_hist_t[i], digits=4), "[s]"))]
    end
    
    ohp
end

@recipe function f(::OHPPres, i::Int64, SimuResult::SimulationResult; plain=false)
    adjust = 1e-2;

    tube_sys = getcurrentsys(SimuResult,i)
    closedornot := tube_sys.tube.closedornot
    
    tube_hist_t = SimuResult.tube_hist_t
    
    Ptmp = tube_sys.mapping.P_interp_liquidtowall[tube_sys.wall.Xarray]
    
    # grid = SimuResult.grid
    grid = SimuResult.integrator_plate.p.base_cache.g
    # ohp = SimuResult.integrator_plate.p.qline[1].body
    ohp_model = SimuResult.integrator_plate.p.forcing["heating models"][end]
    ohp = ohp_model.transform(ohp_model.shape)
    
    # seriestype --> :heatmap
    xlabel --> "x [m]"
    ylabel --> "y [m]"
    fillalpha := 0
    linewidth := 2.0
    clim      := (minimum(Ptmp),maximum(Ptmp))
    line_z    := Ptmp
    # clim         := :batlow
    linecolor := :autumn1
    bbox_inches := "tight"
    border  := :none

    # color := palette(:heat, 3)
    # colorbar_title --> "\n P[Pa]"
    right_margin --> 5mm
    # colorbar --> true
    # legend --> true
    # grid --> true
    
    xlimit --> grid.xlim[1]
    ylimit --> grid.xlim[2]

    if plain == false

        annotation := [(-0.03, -0.028, "Low pressure"),
        (0.04, -0.028, "High pressure"),(0.0, 0.03, string("Time = ", round(tube_hist_t[i], digits=4), "[s]"))]

        @series begin
            # annotation := (0, 0, "Look up!")
            seriestype := :scatter
            color := :red
            [-0.07+adjust],[-0.028]
        end

        # @series begin
        #     seriestype := :scatter
        #     color := :red
        #     [-0.03+adjust],[-0.028]
        # end

        @series begin
            seriestype := :scatter
            color := :yellow
            [0.01],[-0.028]
        end

    else
        annotation := [(0.0, 0.03, string("Time = ", round(tube_hist_t[i], digits=4), "[s]"))]

    end
    ohp
end

@recipe function f(::OHPSuper, i::Int64, SimuResult::SimulationResult)
    adjust = 1e-2;

    tube_sys = getcurrentsys(SimuResult,i)
    closedornot := tube_sys.tube.closedornot
    
    tube_hist_t = SimuResult.tube_hist_t
    
    # Twall = tube_sys.mapping.θ_interp_walltoliquid[tube_sys.wall.Xarray]
    @unpack PtoT,TtoP = tube_sys.propconvert
    @unpack P = tube_sys.vapor

    Twall = SimuResult.tube_hist_θwall[i]
    Tref = (PtoT.(maximum(P)) + PtoT.(minimum(P)))/2
    Tsat = PtoT(tube_sys.mapping.P_interp_liquidtowall[tube_sys.wall.Xarray])

    ΔT = Twall - Tsat

    ΔTmin = RntoΔT(SimuResult.integrator_tube.p.wall.Rn,Tref,SimuResult.integrator_tube.p.wall.fluid_type,SimuResult.integrator_tube.p.tube.d,TtoP)
    
    # grid = SimuResult.grid
    grid = SimuResult.integrator_plate.p.base_cache.g
    # ohp = SimuResult.integrator_plate.p.qline[1].body
    ohp_model = SimuResult.integrator_plate.p.forcing["heating models"][end]
    ohp = ohp_model.transform(ohp_model.shape)
    
    # seriestype --> :heatmap
    xlabel --> "x [m]"
    ylabel --> "y [m]"
    fillalpha := 0
    linewidth := 2.0
    clim      := (ΔTmin-1e-9,ΔTmin+1e-9)
    line_z    := ΔT
    # clim         := :batlow
    linecolor := palette([:gray,:red], 2)
    bbox_inches := "tight"
    border  := :none

    # color := palette(:heat, 3)
    # colorbar_title --> "\n P[Pa]"
    right_margin --> 5mm
    # colorbar --> true
    # legend --> true
    # grid --> true
    
    xlimit --> grid.xlim[1]
    ylimit --> grid.xlim[2]

    annotation := [(0.0, -0.032, string("ΔT > critical value ",round(ΔTmin, digits=3), "[K], \n  Tref =",round(Tref, digits=3)," [K]")),(0.0, 0.03, string("Time = ", round(tube_hist_t[i], digits=2), "[s]"))]

    # @series begin
    #     # annotation := (0, 0, "Look up!")
    #     seriestype := :scatter
    #     color := :black
    #     [-0.07+adjust],[-0.028]
    # end

    # @series begin
    #     seriestype := :scatter
    #     color := :red
    #     [-0.03+adjust],[-0.028]
    # end

    @series begin
        seriestype := :scatter
        color := :red
        [-0.06],[-0.028]
    end
    ohp
end

@recipe function f(::OHP2DΔT, i::Int64, SimuResult::SimulationResult; plain=false)
    adjust = 1e-2;

    tube_sys = getcurrentsys(SimuResult,i)
    closedornot := tube_sys.tube.closedornot
    
    tube_hist_t = SimuResult.tube_hist_t
    
    # Twall = tube_sys.mapping.θ_interp_walltoliquid[tube_sys.wall.Xarray]
    @unpack PtoT,TtoP = tube_sys.propconvert
    @unpack P = tube_sys.vapor
    @unpack θ_interp_liquidtowall = tube_sys.mapping

    Twall = SimuResult.tube_hist_θwall[i]
    Ts = θ_interp_liquidtowall(tube_sys.wall.Xarray)

    ΔT = Twall - Ts

    # grid = SimuResult.grid
    grid = SimuResult.integrator_plate.p.base_cache.g
    # ohp = SimuResult.integrator_plate.p.qline[1].body
    ohp_model = SimuResult.integrator_plate.p.forcing["heating models"][end]
    ohp = ohp_model.transform(ohp_model.shape)


    sluginterp = slug_interp(tube_sys)
    xs =  tube_sys.wall.Xarray
    markers  = heaviside.(map(sluginterp, xs) .- 1.5)
    
    # seriestype --> :heatmap
    xlabel --> "x [m]"
    ylabel --> "y [m]"

    linewidth := 2.0
    clims      := (-5.0,5.0)
    line_z    := markers .* ΔT
    # clim         := :batlow
    # linecolor := :bluesreds
    bbox_inches := "tight"
    border  := :none

    color := :bluesreds
    colorbar_title --> "\n ΔT [K]"
    right_margin --> 5mm

    xlimit --> grid.xlim[1]
    ylimit --> grid.xlim[2]

    if plain == false
            
        colorbar --> true

        annotation := [(-0.04+adjust, -0.028, "Plate hotter"), 
                        (0.03+adjust, -0.028, "Plate cooler"), 
                        (0.0, 0.03, string("Time = ", round(tube_hist_t[i], digits=2), "[s]"))]

        @series begin
            seriestype := :scatter
            color := :red
            colorbar_entry  := false
            [-0.06],[-0.028]
        end
        @series begin
            seriestype := :scatter
            color := :blue
            colorbar_entry  := false
            [0.01],[-0.028]
        end
    else
        colorbar --> false

        annotation := [(0.0, 0.03, string("Time = ", round(tube_hist_t[i], digits=2), "[s]"))]

    end


    ohp
end

@recipe function f(::OHPTexp, label_for_plotting::Vector, exp_data::Tuple, SimuResult::SimulationResult)  
    
    seriestype --> :scatter
    
    RTDt = exp_data[1]
    RTD_T = exp_data[2]

    
    T0 = SimuResult.integrator_tube.u[end]
    
    title := "temperature curve"
    linewidth := 2
#     legend := :topleft
    
    xlabel:="Time [s]"
    ylabel:="T-T₀ [K]"
    
    n = length(label_for_plotting)
    color := palette(:default)[1:n]'
    
#     palette := :tab10
#     ribbon := 1
    
#     xlimit = (0.0,SimuResult.tube_hist_t[end])
    
    
    label := string.("RTD", label_for_plotting', " exp")
    
    RTDt,RTD_T[:,label_for_plotting] .- RTD_T[1,label_for_plotting]'
end

@recipe function f(::OHPTcurve, label_for_plotting::Vector, sensor_data::Tuple, SimuResult::SimulationResult)  
    
    thist = sensor_data[1]
    ghist = sensor_data[2]
    
    T0 = SimuResult.integrator_tube.u[end]
    
    title --> "Temperature curve"
    linewidth := 2
#     legend := :topleft
    
    xlabel:="Time [s]"
    ylabel:="T-T₀ [K]"
    
    xlimit := (0.0,thist[end])
    ylimit := (0.0,maximum(ghist[:,label_for_plotting] .- T0))
    
#     ribbon := 1
    n = length(label_for_plotting)
    color := palette(:default)[1:n]'
    
    
    label --> string.("RTD", label_for_plotting')
    
    thist,ghist[:,label_for_plotting] .- T0
end

@recipe function f(::OHPCond, label_for_conductance::Tuple, sensor_data::Tuple, SimuResult::SimulationResult)  
    
    i1 = label_for_conductance[1]
    i2 = label_for_conductance[2]
     
    thist = sensor_data[1]
    ghist = sensor_data[2]
    
    power = SimuResult.integrator_tube.p.wall.power
    
    T0 = SimuResult.integrator_tube.u[end]
    
    title := "Thermal conductance"
    linewidth := 2
    # legend := :topright
    
    xlabel:="Time [s]"
    ylabel:="G [W/K]"
    
    xlimit --> (0.0,thist[end])
    ylimit --> (0.0,20.0)
    
    color := :blue
#     ribbon := 1
#     n = length(label_for_plotting)
#     color := palette(:default)[1:n]'
    
    label := string.("G model")

    thist,power./(ghist[:,i1] .-ghist[:,i2])
end

@recipe function f(::OHPCond, label_for_conductance::Tuple, sensor_data::Tuple, exp_data::Tuple, SimuResult::SimulationResult)  
    
    i1 = label_for_conductance[1]
    i2 = label_for_conductance[2]
     
    thist = sensor_data[1]
    ghist = sensor_data[2]
    
    RTDt = exp_data[1]
    RTD = exp_data[2]
    
    power = SimuResult.integrator_tube.p.wall.power
    
    T0 = SimuResult.integrator_tube.u[end]
    
    title := "Thermal conductance"
    linewidth := 2
    # legend := :topright
    
    xlabel:="Time [s]"
    ylabel:="G [W/K]"
    
    xlimit --> (0.0,thist[end])
    ylimit --> (0.0,20.0)
    
    color := :blue
#     ribbon := 1
#     n = length(label_for_plotting)
#     color := palette(:default)[1:n]'
    
    label := string.("G model")
    
    @series begin
        # annotation := (0, 0, "Look up!")
        seriestype := :scatter
        color := :darkred
        label := string.("G experiment")
        RTDt,power./(RTD[:,i1] .- RTD[:,i2])
    end
    
    thist,power./(ghist[:,i1] .-ghist[:,i2])
end

@recipe function f(::OHPV, SimuResult::SimulationResult)  
    
    thist = SimuResult.tube_hist_t
    sysfinal = getcurrentsys.([SimuResult],1:length(thist));

    
    Vavg_hist = []
    Vabsavg_hist = []
    Vmax_hist = []
    Vmin_hist = []
    for sysi in sysfinal
        V = [elem[2] for elem in sysi.liquid.dXdt]
        Vabsavg = mean(abs.(V))
        Vavg = mean(V)
        Vmax = maximum((V))
        Vmin = minimum((V))

        push!(Vavg_hist, Vavg)
        push!(Vabsavg_hist, Vabsavg)
        push!(Vmax_hist, Vmax)
        push!(Vmin_hist, Vmin)

    end

    title := "Velocity data"
    linewidth := 2
    
    xlabel:="Time [s]"
    ylabel:="V [m/s]"
    label := "|V| (avg)"
    
    xlimit --> (0.0,thist[end])
    
    @series begin
        # label := string.("RTD", i1, "-RTD", i2, " exp")
        
        fillalpha := 1
        c := :green
        fillrange := Float64.(Vmin_hist)
        fillcolor :=:green
        label := "V range"
        thist,Vmax_hist
    end
    
    @series begin
        # label := string.("RTD", i1, "-RTD", i2, " exp")
            fillalpha := 0.0
        # fillcolor :=:green
        label := false
        c := :green
        thist,Vmin_hist
    end
    
    @series begin
        # label := string.("RTD", i1, "-RTD", i2, " exp")
        c := :yellow
        label :="V (avg)"
        thist,Vavg_hist
    end
    
    c := :red
    
    thist,Vabsavg_hist
end

@recipe function f(::OHP1DT,i::Int64,SimuResult::SimulationResult)
    tube_sys = getcurrentsys(SimuResult,i)

    # title --> string(" time = ",SimuResult.tube_hist_t)
    plottype := "T"
    time := SimuResult.tube_hist_t[i]
    tube_sys
end

@recipe function f(::OHP1DP,i::Int64,SimuResult::SimulationResult)
    tube_sys = getcurrentsys(SimuResult,i)

    plottype := "P"
    time := SimuResult.tube_hist_t[i]
    tube_sys
end

@recipe function f(::OHP1DΔT,i::Int64,SimuResult::SimulationResult)
    tube_sys = getcurrentsys(SimuResult,i)
    L = tube_sys.tube.L

    @unpack TtoP = tube_sys.propconvert

    ΔTmin = RntoΔT(SimuResult.integrator_tube.p.wall.Rn,291.2,SimuResult.integrator_tube.p.wall.fluid_type,SimuResult.integrator_tube.p.tube.d,TtoP)

    @series begin
        # label := string.("RTD", i1, "-RTD", i2, " exp")
        c := :red
        linestyle := :dash
        linewidth := 2
        label :="critical ΔT (Tref = 291.2 K)"
        [0,L],[ΔTmin,ΔTmin]
    end

    plottype := "ΔT"
    time := SimuResult.tube_hist_t[i]
    tube_sys
end

@recipe function f(::OHPTwall,i::Int64,SimuResult::SimulationResult)
    tube_sys = getcurrentsys(SimuResult,i)

    c := :blue

    title := string(" time = ", round(SimuResult.tube_hist_t[i]; digits = 3))
    label := "wall T [K]"
    linewidth := 2
    ylabel := "wall T [K]"

    tube_sys.wall.Xarray,SimuResult.tube_hist_θwall[i]
end

@recipe function f(val::PHPSystem;plottype="T",time=0)

    # T0 = 295.0
    # P0 = 220337

    xlabel := "ξ [m]"
    linewidth := 2
    title := string(" time = ", round(time; digits = 3))

    # layout := (3,1)
    if plottype == "T"

            x1,y1 = stackXpTemp(val)

            popfirst!(x1)
            popfirst!(y1)

            legend := false
            # linecolor := [:red, :blue]
            # color := palette(:tab10)[2]
            color_palette  := :seaborn_dark
            # title := "OHP temperatures"
            # xlabel := "ξ [m]"
            ylabel := "OHP temperatures T [K]"
            # y1 = y1 .* T0

            return x1,y1

    elseif plottype == "ΔT"

            x2,y2 = stackXpTemp(val)

            label := "ΔT"
            color_palette  := :seaborn_dark
            color := palette(:seaborn_dark)[1]
            # title := "OHP ΔT"

            # xlabel := "ξ [m]"
            ylabel := "OHP superheat ΔT [K]"

            x2 = x2[1]
            y2 = y2[1]

            @unpack PtoT = val.propconvert
            y2 -= PtoT(val.mapping.P_interp_liquidtowall.(x2))

            # y2 -= nondi_PtoT(val.mapping.P_interp_liquidtowall.(x2))

            # y2 = y2 .* T0

            return x2,y2


    elseif plottype == "P"

            x3,y3 = stackXpTemp(val)

            legend := false
            color_palette  := :seaborn_dark
            # title := "OHP pressures"

            # xlabel := "ξ [m]"
            ylabel := "OHP pressures P [Pa]"

            popfirst!(x3)
            popfirst!(y3)

            for i = 1:length(y3)
                y3[i] = val.mapping.P_interp_liquidtowall.(x3[i])
            end

            # y3 = y3 .* P0

            return x3,y3
    end
end

function stackXpTemp(val::PHPSystem)
    Xpvapor = getXpvapor(val.liquid.Xp,val.tube.closedornot)
    # θvapor  = nondi_PtoT.(val.vapor.P)
    
    @unpack PtoT = val.propconvert
    θvapor  = PtoT.(val.vapor.P)
    Xp = val.liquid.Xp

    all_θ  = []
    all_Xp = []

    push!(all_Xp,val.wall.Xarray); push!(all_θ, val.wall.θarray)

    j=1
    while j <= length(Xp)
        if Xp[j][end] >= Xp[j][1]
            push!(all_Xp,val.liquid.Xarrays[j]); push!(all_θ, val.liquid.θarrays[j])
            else
            # find the index at the end
            index = findfirst(x->x <= val.liquid.Xarrays[j][end], val.liquid.Xarrays[j])

            push!(all_Xp,val.liquid.Xarrays[j][1:index-1]); push!(all_θ, val.liquid.θarrays[j][1:index-1])
            push!(all_Xp,val.liquid.Xarrays[j][index:end]); push!(all_θ, val.liquid.θarrays[j][index:end])
        end

        j += 1
    end


    j=1
    while j <= length(Xpvapor)
        if Xpvapor[j][end] >= Xpvapor[j][1]
            push!(all_Xp,[Xpvapor[j][1],Xpvapor[j][end]]); push!(all_θ,[θvapor[j], θvapor[j]])
            else
            push!(all_Xp,[0.0,Xpvapor[j][end]]); push!(all_θ,[θvapor[j], θvapor[j]])
            push!(all_Xp,[Xpvapor[j][1],val.tube.L]); push!(all_θ,[θvapor[j], θvapor[j]])
        end

        j += 1
    end

    return all_Xp,all_θ
end

# @recipe function f(val::BasicBody)
#     fillalpha  := 0
#     linecolor  := :black
#     framestyle := :box
#     return val
# end

# @recipe function f(val::PrescribedHeatFluxRegion)
#             fillcolor := :red
#             alpha     := 0.5
#             return val.body
# end

# @recipe function f(val::PrescribedHeatModelRegion)
#     fillcolor := :navyblue
#     alpha     := 0.5
#     return val.body
# end

@recipe function f(ohp::BasicBody;closedornot::Union{Bool, Nothing}=nothing)
    legend := false
    label := false

    # Extract user-provided line_z (if any)
    line_z_data = get(plotattributes, :line_z, nothing)

    if(closedornot==true)
        color --> :blue # change default plot settings
        x_closed = [ohp.x; ohp.x[1]]
        y_closed = [ohp.y; ohp.y[1]]

        # If user supplied line_z, close it too
        if line_z_data !== nothing
            line_z --> [line_z_data; line_z_data[1]]
        end
        return x_closed, y_closed
        
    elseif(closedornot==false)
            color --> :red
            return ohp.x, ohp.y
    else return error("need to specify closedornot")
    end   
end