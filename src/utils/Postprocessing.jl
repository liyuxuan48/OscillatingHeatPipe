export getcurrentsys!,getcurrentsys,getRTD,getconfig,getghist,getthist,getgt,getsysfinal,getwetness,getV,getδ,getHtmp_marker,
translateOHPdata,getTcurve,oneDtwoDtransform,get_boil_matrix,get1DTandP

"""
    give a new u and an old system, return a new system sysnew
"""

function getcurrentsys!(u,sys0)

        indexes = Int64[]
        θliquidrec = Array[]

        for i in eachindex(u)
            if abs(u[i]+1e10) <= 10^(-1)
                push!(indexes,i)
            end
        end


    Xp,dXdt,M,δstart,δend,Lfilm_start,Lfilm_end = vectoXMδL(u[1:indexes[1]-1])

    if sys0.tube.closedornot == true 
        modX!(Xp,sys0.tube.L)
    end
    for i = 1:length(indexes)-1
    push!(θliquidrec, u[indexes[i]+1:indexes[i+1]-1])
    end
    push!(θliquidrec, u[indexes[end]+1:end])

    Lvaporplug = XptoLvaporplug(Xp,sys0.tube.L,sys0.tube.closedornot)

    Ac = sys0.tube.Ac

    d = sys0.tube.d
    δarea_start = Ac .* (1 .- ((d .- 2*δstart) ./ d) .^ 2);
    δarea_end = Ac .* (1 .- ((d .- 2*δend) ./ d) .^ 2);

    volume_vapor = Lvaporplug .* Ac - Lfilm_start .* δarea_start - Lfilm_end .* δarea_end
    ρ = M ./ volume_vapor

    # println(Lvaporplug)
    # println(Lfilm_start)
    # println(Lfilm_end)

    @unpack DtoP = sys0.tube
    P = DtoP.(ρ)
  
    sysnew = sys0

    sysnew.liquid.Xp = Xp
    sysnew.liquid.dXdt = dXdt
    sysnew.liquid.θarrays = θliquidrec
    sysnew.liquid.Xarrays = updateXarrays(Xp,sysnew.liquid.θarrays,sysnew.tube.L)

    sysnew.vapor.P = P
    sysnew.vapor.δstart = δstart
    sysnew.vapor.δend = δend
    sysnew.vapor.Lfilm_start = Lfilm_start
    sysnew.vapor.Lfilm_end = Lfilm_end

    θ_interp_walltoliquid, θ_interp_liquidtowall, H_interp_liquidtowall, P_interp_liquidtowall = sys_interpolation(sysnew)
    heightg_interp = sysnew.mapping.heightg_interp
    sysnew.mapping = Mapping(θ_interp_walltoliquid, θ_interp_liquidtowall, H_interp_liquidtowall, P_interp_liquidtowall,heightg_interp)

    return sysnew
end

function getcurrentsys(u,sys0)
    systemp = deepcopy(sys0)
    getcurrentsys!(u,systemp)
end

function modX!(Xp,L)

    for i in eachindex(Xp)
        Xp[i] = mod.(Xp[i],L)
    end

    return Xp
end

"""
    When having a new Xp because of dynamics, Xarrays need to be updated, too.
"""

function updateXarrays(Xp,θarrays,L)

    # Xarrays = zero.(deepcopy(θarrays))
    Xarrays = 0 .* θarrays

    for i = 1:length(Xarrays)
        if Xp[i][1] < Xp[i][2]
            Xarrays[i] = range(Xp[i][1], Xp[i][2], length=length(Xarrays[i]))
        else
            Xarrays[i] = mod.(range(Xp[i][1], Xp[i][2]+L, length=length(Xarrays[i])),L)
        end
        Xarrays[i][1] = Xp[i][1]
        Xarrays[i][end] = Xp[i][end]
    end

    return Xarrays
end


function getRTD(xf,Onum)
    if Onum == 1
        RTD = xf["Raw Data!G3603:N7201"]
    elseif Onum == 2
        RTD = xf["Raw Data!U3603:AB7201"]
    elseif Onum == 3
        RTD = xf["Raw Data!AI3603:AP7201"]
    else
        return error
    end
    RTDt = (0:size(RTD,1)-1);
    
    RTDt,Float64.(RTD)
end

function getconfig(filepath)
    Pindex = findfirst("_P0",filepath)[2]
    Oindex = findfirst("_O0",filepath)[2]
    Hindex = findfirst("_O0",filepath)[2]
    power = parse(Int64,filepath[Pindex+1:Pindex+3])
    Onum = parse(Int64,filepath[Oindex+1:Oindex+3])
    Hnum = parse(Int64,filepath[Hindex+1:Hindex+3])
    
    Onum, Hnum, power
end

function getghist(g,H,plate_T_hist)
    ghist = []
    for j in eachindex(g)
        gtemp = []
        for i in eachindex(plate_T_hist)
            H(g,plate_T_hist[i])
        append!(gtemp,deepcopy(g[j]))
    end
    push!(ghist,deepcopy(gtemp))
    end
    ghist
end

function getthist(tube_hist)
    thist = []
    for i in eachindex(tube_hist)
        append!(thist,tube_hist[i].t)
    end
    thist
end

# function getgt(g,H,plate_T_hist,tube_hist)
#     ghist = getghist(g,H,plate_T_hist)
#     thist = getthist(tube_hist)

#     ghist,thist
# end

function getsysfinal(tube_hist)
    sysfinal = []
    for tube_i in tube_hist
        push!(sysfinal, deepcopy(getcurrentsys!(tube_i.u,tube_i.p)))
    end
    
    sysfinal
end

function getsysfinal(tube_hist_u,tube_hist_θwall,integrator_tube)
    sysfinal = []
    for i in eachindex(tube_hist_u)
        push!(sysfinal, deepcopy(getcurrentsys!(tube_hist_u[i],integrator_tube.p)))
        sysfinal[i].wall.θarray = tube_hist_θwall[i]
    end
    
    sysfinal
end

function getwetness(sysfinal)
    wetness = Float64[]
    for sysfinali in sysfinal
        Lvaporplug = XptoLvaporplug(sysfinali.liquid.Xp,sysfinali.tube.L,sysfinali.tube.closedornot)
        Lvaporsum = sum(Lvaporplug)
        Ldryvapor = Lvaporplug - sysfinali.vapor.Lfilm_start  - sysfinali.vapor.Lfilm_end
        Ldrysum = sum(Ldryvapor)

        push!(wetness, 1 - Ldrysum/Lvaporsum)
    end
    
    wetness
end

function getV(sysfinal)
    Vabs_avg = Float64[]
    Vabs_max = Float64[]
    Vavg = Float64[]

    for sysfinali in sysfinal
    
        V = [elem[2] for elem in sysfinali.liquid.dXdt]
        Vabs = mean(abs.(V))
        Vmax = maximum(abs.(V))
        Vrealavg = mean(V)
        
        push!(Vabs_avg, Vabs)
        push!(Vabs_max, Vmax)
        push!(Vavg, Vrealavg)
    end

    Vabs_avg,Vavg,Vabs_max
end

function getδ(sysfinal)
    δ_avg_start = Float64[]
    δ_avg_end = Float64[]

    for sysfinali in sysfinal
        δstart = sysfinali.vapor.δstart
        δend = sysfinali.vapor.δend

        push!(δ_avg_start, mean(δstart))
        push!(δ_avg_end, mean(δend))
    end

    δ_avg_start,δ_avg_end
end

function getHtmp_marker(Htmp,Hₗ,Hᵥ)
    Htmp_marker = (Htmp < Hₗ + 1e-10) && (Htmp > Hₗ - 1e-10) ? 1.0 : (Htmp > Hᵥ + 1e-10) ? 2.0 : 0.0
end

function translateOHPdata(OHPdata::Dict{String, Any})
    ohp = OHPdata["ohp"]
    sr  = OHPdata["SimulationResult"]
    @unpack boil_hist,tube_hist_t,tube_hist_θwall,tube_hist_u,plate_T_hist = sr;
    integrator_tube = OHPdata["integrator_tube"]
    integrator_plate = OHPdata["integrator_plate"];
    sysfinal = getsysfinal(tube_hist_u,tube_hist_θwall,integrator_tube);
    sys = integrator_plate.p
    plate_T_hist = sr.plate_T_hist
    
    sysfinal,ohp,integrator_plate,integrator_tube,boil_hist,plate_T_hist
end

function getTcurve(x::Vector,y::Vector,SimulationResult::SimulationResult)
    
    sys = SimulationResult.integrator_plate.p
    plate_T_hist = SimulationResult.plate_T_hist
    X =  VectorData(x,y);
    H = Regularize(X,cellsize(sys),I0=origin(sys.grid))
    g = ScalarData(X);
    ghist = getghist(g,H,plate_T_hist);
    thist = SimulationResult.tube_hist_t; 
    
    thist,hcat(ghist...)
end

function getTcurve(plate_sensors::Tuple,SimulationResult)
    getTcurve(plate_sensors[1],plate_sensors[2],SimulationResult)
end

function oneDtwoDtransform(ξ,sr::SimulationResult)
    plate_sys = sr.integrator_plate.p
    ql = plate_sys.qline[1]
    interp_linear_x = LinearInterpolation(ql.arccoord, ql.body.x,extrapolation_bc = Line());
    interp_linear_y = LinearInterpolation(ql.arccoord, ql.body.y,extrapolation_bc = Line());

    x2D = interp_linear_x[ξ]
    y2D = interp_linear_y[ξ];
    
    x2D,y2D
end

function get_boil_matrix(SimuResult::SimulationResult;boil_dt = 0.1)
    
    boil_hist = SimuResult.boil_hist;
    
    boil_station_hist = [elem[1] for elem in boil_hist]
    boil_t_hist = [elem[2] for elem in boil_hist]
    boil_data = Array(reshape([boil_station_hist;boil_t_hist], length(boil_t_hist),2));
    
    Xstations = SimuResult.integrator_tube.p.wall.Xstations;
    tend = boil_data[end,2]
#     boil_dt = 0.1;
    
    boil_matrix=spzeros(size(Xstations,1),Int64(round(tend/boil_dt)));
    for tnum = 1:size(boil_matrix,2)
         boil_index = findall(x->Int64(round(x/boil_dt))==tnum,boil_data[:,2])

        if length(boil_index) != 0
            i_station_temp = Int64.(boil_data[boil_index,1])
            for i in i_station_temp
                boil_matrix[i,tnum] = 1
            end
        end
    end
    
    boil_num_t = sum(boil_matrix, dims=1)'
    boil_num_x = sum(boil_matrix, dims=2);
    t_boil = boil_dt:boil_dt:boil_dt*length(boil_num_t)
    
    integrator_tube = SimuResult.integrator_tube
    ξ = integrator_tube.p.wall.Xstations
    x2D_boil,y2D_boil = oneDtwoDtransform(ξ::Vector,SimuResult)
    
    boil_data,boil_num_x,boil_num_t,t_boil,x2D_boil,y2D_boil,boil_dt
end

function get1DTandP(xsensors::Vector, SimuResult::SimulationResult)
    θhist = []
    phist = []
    sysfinal = get
    for i in eachindex(SimuResult.tube_hist_u)
        tube_sys = getcurrentsys!(SimuResult.tube_hist_u[i],SimuResult.integrator_tube.p)
        ptemp = tube_sys.mapping.P_interp_liquidtowall[xsensors]
        @unpack PtoT = tube_sys.tube
        θtemp = PtoT.(ptemp)
        push!(phist, ptemp)
        push!(θhist, θtemp)
    end
    phist1D = Float64.(hcat(phist...)')
    θhist1D = Float64.(hcat(θhist...)');
    
    θhist1D,phist1D
end