export SimulationResult,onesideXp,randomXp,initialize_ohpsys,store!,newstate

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

    # println(typeof(X0))

    Ls = XptoLliquidslug(X0,L)
    real_ratio = sum(Ls)/L
    
    X0,real_ratio
end

function randomXp(tube::Tube;numofslugs=30,chargeratio=0.46,σ_charge=0.1)

    L = tube.L
    Lmin = tube.d
    closedornot = tube.closedornot

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
    real_ratio = sum(Ls)/L

    X0,real_ratio
end

# function L_to_boiltime(L_newbubble,Rn,fluid_type,vapor::Vapor,tube::Tube)
#     P = vapor.P[1]

#     @unpack TtoP,PtoT,PtoD,PtoHfg = tube
#     property = SaturationFluidProperty(fluid_type,PtoT(P));

#     d = tube.d
#     ρₗ = property.ρₗ
#     Cpₗ = property.Cpₗ 
#     kₗ = property.kₗ

#     Hfg = PtoHfg(P)
#     Tref = PtoT(P)
#     ρv = PtoD(P)
#     ΔTthres = RntoΔT(Rn,Tref,fluid_type,d,TtoP)

#     A = (2 * (ΔTthres/Tref) * Hfg*ρv/ρₗ)^0.5
#     Ja = ΔTthres*Cpₗ*ρₗ/ρv/Hfg
#     B = (12*kₗ/pi/ρₗ/Cpₗ)^0.5 * Ja

#     t = 1e-2:1e-2:1e2
#     tstar = t .* A^2 ./ B^2 

#     Rplus = 2/3 .* ((tstar .+ 1).^1.5 .- (tstar).^1.5 .- 1);
#     R = Rplus .* B^2 ./ A; 
#     L_equivalent = (4/3) .* R

#     interp_Rtot = LinearInterpolation(L_equivalent, t);

#     return interp_Rtot(L_newbubble)
# end

# function initialize_ohpsys(OHPtype,fluid_type,sys,p_fluid,Tref,power)
#     initialize_ohpsys(fluid_type,sys,p_fluid,Tref,power)
# end

function initialize_ohpsys(sys,p_fluid,power;closedornot=true,boil_waiting_time=1.0,Rn_boil=3e-6,inertia_f=1.3,d=1e-3,tubeshape="square",Nu=3.6,slugnum=30,film_fraction=0.3,g = [0.0,0.0], ηplus=0.6,ηminus=0.0,nucleatenum = 250,L_newbubble = 6e-3)

    ohp = sys.forcing["heating models"][end]
    L = arccoord(ohp.shape)[end] # total length of the pipe when streched to a 1D pipe (an approximate here)
    @unpack x,y = ohp.transform(ohp.shape)
    
    N=numpts(ohp.shape)

    Tref = p_fluid.Tref

        if tubeshape=="square"
            peri = d*4
            Ac = d*d
        elseif tubeshape=="circle"
            peri = d*π
            Ac = d*d*π/4
        end

        tube = Tube(d,peri,Ac,L,g,closedornot,N);

        # liquid
        # Nu = 3.60
        Hₗ = p_fluid.kₗ/d * Nu # Nusselt number 3.60

        # line = []
        # X0,realratio = onesideXp(ohp,tube,line)
        X0,realratio = randomXp(tube,numofslugs=slugnum)
        dXdt0_l = zeros(length(X0))
        dXdt0_r = zeros(length(X0))
        dXdt0 = map(tuple,dXdt0_l,dXdt0_r);

        # N=numpts(ohp.body)

        Xarrays,θarrays = constructXarrays(X0,N,Tref,L);

        liquids=Liquid(Hₗ,p_fluid.ρₗ,p_fluid.Cpₗ,p_fluid.αₗ,p_fluid.μₗ,p_fluid.σ,X0,dXdt0,Xarrays,θarrays);

        # CoolProp
        fluid_type = p_fluid.fluid_type
        propconvert = PropConvert(fluid_type)

        # Vapor
        @unpack TtoP = propconvert
        Lvaporplug = XptoLvaporplug(X0,L,tube.closedornot)

        P_initial = 0*zeros(length(Lvaporplug)) .+ TtoP(Tref);
        δfilm = 2e-5
        δstart_initial = 0*zeros(length(Lvaporplug)) .+ δfilm ;
        δend_initial = 0*zeros(length(Lvaporplug)) .+ δfilm ;

        Lfilm_start_initial =  0.5 .* film_fraction .* Lvaporplug
        Lfilm_end_initial = 0.5 .* film_fraction .* Lvaporplug
        # δmin = 2e-6
        # vapors=Vapor(ad_fac,Hᵥ,p_fluid.kₗ,δmin,Eratio_plus,Eratio_minus,P,δfilm_deposit,δstart,δend,Lfilm_start,Lfilm_end);
        vapors=Vapor(ad_fac=inertia_f,k = p_fluid.kₗ,P=P_initial,δstart=δstart_initial,δend=δend_initial,Lfilm_start=Lfilm_start_initial,Lfilm_end=Lfilm_end_initial,Eratio_plus=ηplus,Eratio_minus=ηminus);

        # Wall
        # ΔTthres = RntoΔT(Rn,Tref,fluid_type,d)
        # nucleate_density = 0.005
        # nucleatenum = 1000
        Xstations = sort(rand(nucleatenum) .* L);
        Xstation_time = zeros(nucleatenum);

        boil_type = "wall T"
        # L_newbubble = 2tube_d
        # boil_interval = L_to_boiltime(L_newbubble,Rn,fluid_type,vapors::Vapor,tube::Tube)
        boil_interval = boil_waiting_time
        Xwallarray,θwallarray = constructXarrays(arccoord(ohp.shape),L,Tref);
        θwallarray .= Tref

        wall = Wall(boil_interval=boil_interval,fluid_type=fluid_type,boil_type=boil_type,power=power,L_newbubble=L_newbubble,Xstations=Xstations,boiltime_stations=Xstation_time,Xarray=Xwallarray,θarray=θwallarray,Rn=Rn_boil);
        # wall = Wall(fluid_type,boil_type,boil_interval,Rn,L_newbubble,Xstations,Xstation_time,Xwallarray,θwallarray);

    # end

    sys0_nomapping = PHPSystem_nomapping(tube,liquids,vapors,wall,propconvert);
    θ_interp_walltoliquid, θ_interp_liquidtowall, H_interp_liquidtowall, P_interp_liquidtowall = sys_interpolation(sys0_nomapping)

    ht = getgh(g,x,y);
    heightg_interp = LinearInterpolation(Xwallarray,ht,extrapolation_bc = Line())
    mapping = Mapping(θ_interp_walltoliquid, θ_interp_liquidtowall, H_interp_liquidtowall, P_interp_liquidtowall,heightg_interp);

    Mfilm_left, Mfilm_right = getMfilm(sys0_nomapping)
    totalmass = sum(getMvapor(sys0_nomapping)) + sum(getMliquid(sys0_nomapping)) + sum(Mfilm_left .+ Mfilm_right)
    boil_hist_int = []
    cache = Cache(boil_hist_int,totalmass)

    sys0 = PHPSystem(tube,liquids,vapors,wall,propconvert,mapping,cache);

    # Lvaporplug = XptoLvaporplug(X0,sys0.tube.L,sys0.tube.closedornot)
    # M = PtoD(P) .* Lvaporplug .* Ac

    # u=[XMδLtovec(X0,dXdt0,M,δstart,δend,Lfilm_start,Lfilm_end); liquidθtovec(sys0.liquid.θarrays)];
    sys0
end

function newstate(sys0::PHPSystem)
    Ac = sys0.tube.Ac
    d = sys0.tube.d

    X0 = sys0.liquid.Xp
    dXdt0 = sys0.liquid.dXdt

    δstart = sys0.vapor.δstart
    δend = sys0.vapor.δend
    Lfilm_start = sys0.vapor.Lfilm_start
    Lfilm_end = sys0.vapor.Lfilm_end
    P = sys0.vapor.P
    Lvaporplug = XptoLvaporplug(X0,sys0.tube.L,sys0.tube.closedornot)

    δarea_start = Ac .* (1 .- ((d .- 2*δstart) ./ d) .^ 2);
    δarea_end = Ac .* (1 .- ((d .- 2*δend) ./ d) .^ 2);

    volume_vapor = Lvaporplug .* Ac - Lfilm_start .* δarea_start - Lfilm_end .* δarea_end

    @unpack PtoD = sys0.propconvert

    M = PtoD.(P) .* volume_vapor

    u=[XMδLtovec(X0,dXdt0,M,δstart,δend,Lfilm_start,Lfilm_end); liquidθtovec(sys0.liquid.θarrays)];
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