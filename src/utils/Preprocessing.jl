export SimulationResult,onesideXp,randomXp,initialize_ohpsys,store!,newstate

"""
(WIP)

a more general way to generate the liquid slugs
"""

function onesideXp(ohp,tube::Tube,line)
    
    A = line[1]
    B = line[2]
    C = line[3]
    sign = line[4]

    L = tube.L
    Lmin = tube.d
    
    body = ohp.body
    
    largeorsmall = A .* body.x .+ B .* body.y .+ C .< 0

    ls_label = xor.(largeorsmall, [largeorsmall[2:end];largeorsmall[1]])
    X0array_label = findall(!iszero,ls_label)
    X0array = ohp.arccoord[X0array_label]

    X0 = map(tuple,X0array[1:2:end], X0array[2:2:end])

    dXdt0 = [zero.(X) for X in X0]

    Ls = XptoLliquidslug(X0,L)
    real_ratio = sum(Ls)/L
    
    X0,dXdt0,real_ratio
end

"""
generate an array of tuple "X0" of the liquid slugs (X0 stands for location), 
"dXdt" of the same format but in zero for the liquid slugs (dXdt stands for velocity)
"real_ratio" means the volume fraction of the liquid slugs
"""

function randomXp(tube::Tube; kwargs...)
    @unpack L,d,closedornot = tube

    randomXp(L::Real,d::Real,closedornot::Bool; kwargs...)
end

function randomXp(L::Real,Lmin::Real,closedornot::Bool;numofslugs=30,chargeratio=0.46,σ_charge=0.1)


    σ_persection = σ_charge*L/sqrt(numofslugs)

    L_perslug=L/numofslugs*chargeratio
    L_persection=L/numofslugs

    Ls = abs.((rand(numofslugs) .- 0.5).*σ_persection .+ L_perslug)

    Xp1s = zeros(numofslugs);
    Xp2s = deepcopy(Xp1s);

    if minimum(Ls) > Lmin && maximum(Ls) < L_persection
        if closedornot == true

            for i in eachindex(Xp1s)
                Xp1s[i] = (i-1)*L_persection
                Xp2s[i] = Xp1s[i] + Ls[i]
            end

            displacement = L*rand()

            Xp1s = mod.(Xp1s.+displacement,L)
            Xp2s = mod.(Xp2s.+displacement,L)

        # add openloop(starting from 0 for simplicity for now)
        elseif closedornot == false && numofslugs != 1
            for i in eachindex(Xp1s[1:end-1])
                Xp1s[i] = (i-1)*L_persection
                Xp2s[i] = Xp1s[i] + Ls[i]
            end
            Xp2s[end] = L
            Xp1s[end] = Xp2s[end] - Ls[end]

            displacement = 0.0
        end

    else println("generation failed")
    end

    X0 = map(tuple,Xp1s,Xp2s)
    dXdt0 = [zero.(X) for X in X0]
    real_ratio = sum(Ls)/L

    X0,dXdt0,real_ratio
end

"""
give d and shape of tube, return perimeter and cross-section area
"""
function peri_Ac(d::Float64,tubeshape::String)
    if tubeshape=="square"
        peri = d*4
        Ac = d*d
    elseif tubeshape=="circle"
        peri = d*π
        Ac = d*d*π/4
    end

    peri,Ac
end

function initialize_ohpsys(sys,p_fluid,power;closedornot=true,boil_waiting_time=1.0,Rn_boil=3e-6,inertia_f=1.3,d=1e-3,tubeshape="square",Nu=3.6,slugnum=30,δfilm_relative=0.04,film_fraction=0.3,g = [0.0,0.0], ηplus=0.6, ηminus=0.0, nucleatenum = 250, L_newbubble = 6e-3, ch_ratio=0.46)

    # unpack CoolProp Properties
    @unpack fluid_type,Tref,kₗ,ρₗ,Cpₗ,αₗ,μₗ,σ = p_fluid  

    # PropConvert
    # an interpolation between different working fluid propeties (currently assuming saturated gas curve)
    propconvert = PropConvert(fluid_type)

    # Tube
    ohp = sys.forcing["heating models"][end] # by default, the last heating model is the ohp
    L = arccoord(ohp.shape)[end]             # total length of the pipe when streched to a 1D pipe (an approximate here)
    N=numpts(ohp.shape)                      # number of ohp points
    peri,Ac = peri_Ac(d,tubeshape)
    tube = Tube(d,peri,Ac,L,g,closedornot,N);

    # Liquid
    Hₗ = p_fluid.kₗ/d * Nu # Nusselt number given
    X0,dXdt0,realratio = randomXp(tube,numofslugs=slugnum,chargeratio=ch_ratio)
    Xarrays,θarrays = constructXarrays(X0,N,Tref,L);
    
    liquids=Liquid(Hₗ,ρₗ,Cpₗ,αₗ,μₗ,σ,X0,dXdt0,Xarrays,θarrays);

    # Vapor
    @unpack TtoP = propconvert

    Lvaporplug = XptoLvaporplug(X0,L,tube.closedornot)
    P_initial = zero(Lvaporplug) .+ TtoP(Tref);
    δfilm = δfilm_relative * d/2
    δstart_initial = zero(Lvaporplug) .+ δfilm ;
    δend_initial   = zero(Lvaporplug) .+ δfilm ;
    Lfilm_start_initial = 0.5 .* film_fraction .* Lvaporplug
    Lfilm_end_initial   = deepcopy(Lfilm_start_initial)
    vapors=Vapor(ad_fac=inertia_f,k = p_fluid.kₗ,P=P_initial,δstart=δstart_initial,δend=δend_initial,Lfilm_start=Lfilm_start_initial,Lfilm_end=Lfilm_end_initial,Eratio_plus=ηplus,Eratio_minus=ηminus);

    # Wall
    Xstations = sort(rand(nucleatenum) .* L);
    Xstation_time = zeros(nucleatenum);
    boil_type = "wall T"
    boil_interval = boil_waiting_time
    Xwallarray,θwallarray = constructwallXθarray(arccoord(ohp.shape),Tref);
    wall = Wall(boil_interval=boil_interval,fluid_type=fluid_type,boil_type=boil_type,power=power,L_newbubble=L_newbubble,Xstations=Xstations,boiltime_stations=Xstation_time,Xarray=Xwallarray,θarray=θwallarray,Rn=Rn_boil);

    # Mapping
    @unpack x,y = ohp.transform(ohp.shape)
    sys0_nomapping = PHPSystem_nomapping(tube,liquids,vapors,wall,propconvert);
    θ_interp_walltoliquid, θ_interp_liquidtowall, H_interp_liquidtowall, P_interp_liquidtowall = sys_interpolation(sys0_nomapping)
    ht = getgh(g,x,y);
    heightg_interp = LinearInterpolation(Xwallarray,ht,extrapolation_bc = Line())
    mapping = Mapping(θ_interp_walltoliquid, θ_interp_liquidtowall, H_interp_liquidtowall, P_interp_liquidtowall,heightg_interp);

    # Cache
    Mfilm_left, Mfilm_right = getMfilm(sys0_nomapping)
    totalmass = sum(getMvapor(sys0_nomapping)) + sum(getMliquid(sys0_nomapping)) + sum(Mfilm_left .+ Mfilm_right)
    boil_hist_int = []
    cache = Cache(boil_hist_int,totalmass)

    # PHPSystem
    sys0 = PHPSystem(tube,liquids,vapors,wall,propconvert,mapping,cache);

    sys0
end


function newstate(sys0::PHPSystem)
    @unpack Ac,d,L,closedornot = sys0.tube
    @unpack Xp,dXdt = sys0.liquid
    @unpack δstart,δend,Lfilm_start,Lfilm_end,P = sys0.vapor
    @unpack PtoD = sys0.propconvert

    M = systoM(sys0::PHPSystem)

    u=[XMδLtovec(Xp,dXdt,M,δstart,δend,Lfilm_start,Lfilm_end); liquidθtovec(sys0.liquid.θarrays)];
end

mutable struct SimulationResult
    tube_hist_t      ::Vector{Any}
    tube_hist_u      ::Vector{Any}
    tube_hist_θwall  ::Vector{Any}
    boil_hist        ::Vector{Any}
    plate_T_hist     ::Vector{Any}
    integrator_tube  ::Any
    integrator_plate ::Any
    integrator_tube_resume  ::Any
    integrator_plate_resume ::Any    
    # grid             ::Any
end

function SimulationResult(tube_hist_t,tube_hist_u,tube_hist_θwall,boil_hist,plate_T_hist,integrator_tube,integrator_plate)
    integrator_tube_resume = deepcopy(integrator_tube)
    integrator_plate_resume = deepcopy(integrator_plate)
    SimulationResult(tube_hist_t,tube_hist_u,tube_hist_θwall,boil_hist,plate_T_hist,integrator_tube,integrator_plate,integrator_tube_resume,integrator_plate_resume)
    # SimulationResult(tube_hist_t,tube_hist_u,tube_hist_θwall,boil_hist,plate_T_hist,integrator_tube,integrator_plate,integrator_tube_resume,integrator_plate_resume,grid)
end

function SimulationResult(int_tube,int_plate)
    
    boil_hist= []
    plate_T_hist = []
    tube_hist_u  = []
    tube_hist_t = []
    tube_hist_θwall = []
    integrator_tube = deepcopy(int_tube)
    integrator_plate = deepcopy(int_plate)

    # grid = int_plate.p.grid

    return SimulationResult(tube_hist_t,tube_hist_u,tube_hist_θwall,boil_hist,plate_T_hist,integrator_tube,integrator_plate)
    # return SimulationResult(tube_hist_t,tube_hist_u,tube_hist_θwall,boil_hist,plate_T_hist,integrator_tube,integrator_plate,grid)
end

function store!(sr,integrator_tube,integrator_plate)
        
        append!(sr.boil_hist,deepcopy(integrator_tube.p.cache.boil_hist));
        integrator_tube.p.cache.boil_hist = []

        push!(sr.plate_T_hist,deepcopy(temperature(integrator_plate)));
        push!(sr.tube_hist_θwall,deepcopy(integrator_tube.p.wall.θarray))
        push!(sr.tube_hist_u,deepcopy(integrator_tube.u));
        push!(sr.tube_hist_t,deepcopy(integrator_tube.t));
        sr.integrator_tube_resume = deepcopy(integrator_tube)
        sr.integrator_plate_resume = deepcopy(integrator_plate)

end