using Statistics
using Interpolations
using LinearAlgebra

heaviside(x::AbstractFloat) = ifelse(x < 0, zero(x), ifelse(x > 0, one(x), oftype(x,0.0)))

hc = 4500 # condenser coeff
h_adiabatic = 160 # adiabatic coeff
length_ohpchannel_mm = 141 # x range of the ohp in mm
power = 30 # Watt

ρₛ = 2730; # material density [kg/m^3] Aluminum3003-O
cₛ  = 8.93e02; # material specific heat [J/kg K]
# kₛ  = 2.37e02; # material heat conductivity Aluminum6061
kₛ  = 1.93e02; # material heat conductivity Aluminum3003-O


αₛ = kₛ/ρₛ/cₛ

dₛ = 1.37e-3; # effective d (The thickness of an ideal uniform thickness plate occupying the same volume)

Tref = 291.15# reference temperature (try higher!)
fluid_type = "Butane"
p_fluid = SaturationFluidProperty(fluid_type,Tref) # This function relies on CoolProp.jl package

tube_d = 1e-3

Δx = 0.000633 # [m]

Lx = 6*INCHES*1.02; # plate size x [m]
Ly = 1.981*INCHES*1.06; # plate size y [m]
xlim = (-Lx/2,Lx/2) # plate x limits
ylim = (-Ly/2,Ly/2) # plate y limits

g = PhysicalGrid(1.03 .* xlim,1.1 .* ylim,Δx); # build a gird slightly larger than the plate

# power = 80 # Watt
areaheater_area = 50e-3 * 50e-3

phys_params = Dict( "diffusivity"              => αₛ,
                    "flux_correction"          => ρₛ*cₛ*dₛ,
                    # "angular velocity"         => 0.0,
                    "Fourier"                  => 1.0,
                    "ohp_flux"                 => [NaN], # initial value, the value here is useless
                    "areaheater_power"         => power, # total power
                    "areaheater_area"          => areaheater_area, # total area
                    "areaheater_temp"          => 0.0,   # relative temperature compared with "background temperature"
                    "areaheater_coeff"         => hc,
                    "adiabatic_coeff"          => h_adiabatic,
                    "background temperature"   => Tref
                     )
     

Δs = 1.4*cellsize(g)

trim = 0.006
cond_block = 1.1INCHES
xbound = [ -Lx/2,-Lx/2, 
            Lx/2, Lx/2];
ybound = [  Ly/2,-Ly/2, 
           -Ly/2, Ly/2];
body = Polygon(xbound,ybound,Δs)

X = MotionTransform([0,0],0)
joint = Joint(X)
m = RigidBodyMotion(joint,body)
x = zero_motion_state(body,m)
update_body!(body,x,m)

Lx_a = 6*INCHES; # plate size x [m]
Ly_a = 1.981*INCHES; # plate size y [m]
xbound_a = [ -Lx_a/2,-Lx_a/2, 
            Lx_a/2, Lx_a/2];
ybound_a = [  Ly_a/2,-Ly_a/2, 
           -Ly_a/2, Ly_a/2];
body_a = Polygon(xbound_a,ybound_a,Δs)
update_body!(body_a,x,m)

function get_qbplus(t,x,base_cache,phys_params,motions)
    nrm = normals(base_cache)
    qbplus = zeros_surface(base_cache)
    return qbplus
end

function get_qbminus(t,x,base_cache,phys_params,motions)
    nrm = normals(base_cache)
    qbminus = zeros_surface(base_cache)
    # qbminus .= nrm.u
    return qbminus
end

bcdict = Dict("exterior" => get_qbplus,"interior" => get_qbminus)

function heatermodel!(σ,T,t,fr::AreaRegionCache,phys_params)
    σ .= phys_params["areaheater_power"] / phys_params["areaheater_area"] / phys_params["flux_correction"] 
end


function condensermodel!(σ,T,t,fr::AreaRegionCache,phys_params)
    T0 = phys_params["areaheater_temp"]
    h = phys_params["areaheater_coeff"]
    corr = phys_params["flux_correction"] 

    σ .= h*(T0 - T) / corr
end

function adiabaticmodel!(σ,T,t,fr::AreaRegionCache,phys_params)
    T0 = phys_params["areaheater_temp"]
    h = phys_params["adiabatic_coeff"]
    corr = phys_params["flux_correction"] 

    σ .= h*(T0 - T)/ corr
end

fregion1_h = Rectangle(25e-3,25e-3,1.4*Δx) # H01
tr1_h = RigidTransform((0.0,-0.0),0.0)

heater1 = AreaForcingModel(fregion1_h,tr1_h,heatermodel!);

fregion1_c = Rectangle(7.5e-3,Ly_a/2,1.4*Δx)
tr1_c = RigidTransform((3INCHES-7.5e-3,-0.0),0.0)

fregion2_c = deepcopy(fregion1_c)
tr2_c = RigidTransform((-(3INCHES-7.5e-3),-0.0),0.0)

cond1 = AreaForcingModel(fregion1_c,tr1_c,condensermodel!);
cond2 = AreaForcingModel(fregion2_c,tr2_c,adiabaticmodel!);

# x, y = construct_ohp_curve("ASETS",Δx) # get x and y coordinates for the channel
    ds = 1.5Δx
    nturn = 16
    width_ohp = 1.821*INCHES - tube_d
    length_ohp = 5.34*INCHES + (length_ohpchannel_mm *1e-3 - 5.34*INCHES)
    gap = 0.11INCHES - 0.5*tube_d
    pitch = width_ohp/(2*nturn+1)

    x0, y0 = Lx_a/2 - 0.34*INCHES - 0.5*tube_d + (length_ohpchannel_mm *1e-3 - 5.34*INCHES), width_ohp/2
    x, y, xf, yf = construct_ohp_curve(nturn,pitch,length_ohp,gap,ds,x0,y0,false,false,3pi/2)

ohp = BasicBody(x,y) # build a BasicBody based on x,y

tr_ohp = RigidTransform((0.0,0.0),0.0)

# phys_params["ohp_flux"] = zero(x);

function ohpmodel!(σ,T,t,fr::LineRegionCache,phys_params)
    σ .= phys_params["ohp_flux"] ./ phys_params["flux_correction"] 
end
ohp_linesource = LineForcingModel(ohp,tr_ohp,ohpmodel!);

# forcing_dict = Dict("heating models" => [heater1,cond1,ohp_linesource])
forcing_dict = Dict("heating models" => [heater1,cond1,cond2,ohp_linesource])

tspan = (0.0, 300.0); # start time and end time
dt_record = 0.2   # saving time interval

tstep = 4e-4     # actrual time marching step

timestep_fixed(u,sys) = tstep

prob = NeumannHeatConductionProblem(g,body,scaling=GridScaling,
                                             phys_params=phys_params,
                                             bc=bcdict,
                                             motions=m,
                                             forcing=forcing_dict,
                                             # timestep_func=timestep_fourier
                                             timestep_func=timestep_fixed
    );

sys_plate = construct_system(prob);

# test 1, 30 slugs, no heat transfer, given an initial velocity, with and without gravity
@testset "dynamicsmodel test1" begin
    numofslugs = 30

    for gi in [[0.0,0.0],[-9.8,0.0]]
  
        sys_tube = initialize_ohpsys(sys_plate,p_fluid,power,g=gi,slugnum=numofslugs,ηplus=rand())
        Vint = 1.0 + rand()
        sys_tube.liquid.dXdt = [zero.(X) .+ Vint for X in sys_tube.liquid.dXdt]# initial velocity of the slugs
        sys_tube.vapor.δstart .= zeros(numofslugs) .+ 1e-5 # initial velocity of the slugs

        u_tube = newstate(sys_tube) # initialize OHP tube
        # integrator_tube = init(u_tube,tspan,sys_tube); # construct integrator_tube

        σ = sys_tube.liquid.σ
        L = sys_tube.tube.L
        ρₗ = sys_tube.liquid.ρₗ
        μₗ = sys_tube.liquid.μₗ
        ad_fac = sys_tube.vapor.ad_fac
        d = sys_tube.tube.d
        Ac = sys_tube.tube.Ac
        peri = sys_tube.tube.peri
        Xp = sys_tube.liquid.Xp
        dXdt = sys_tube.liquid.dXdt
        δend = sys_tube.vapor.δend
        # characteristic bulk velocities for each liquid slug
        V = [mean(elem) for elem in sys_tube.liquid.dXdt]
        # get a characteristic Capilarry number based on the average velocities
        Vavg = mean(abs.(V))
        Ca = getCa.(μₗ,σ,Vavg)

        δdep = Catoδ(d,Ca,adjust_factor=ad_fac)

        Lliquidslug = XptoLliquidslug(Xp,L)

        Adeposit = getAdeposit(sys_tube,δdep)
        Adeposit_left = [elem[1] for elem in Adeposit]
        Adeposit_right = [elem[2] for elem in Adeposit]
        
    # The first test is to check the case where there is no heat transfer and an initial velocity
        uu_test1 = dynamicsmodel(u_tube[1:9*numofslugs],sys_tube)

        @test all(uu_test1[1:2:2*numofslugs-1] .≈ Ac ./ (Ac .- Adeposit_left) .*  V) # velocity should remain the same
        @test all(uu_test1[2:2:2*numofslugs]   .≈ Ac ./ (Ac .- Adeposit_right) .*  V) # velocity should remain the same

        # get differential equation factors
        lhs = ρₗ*Ac .* Lliquidslug

        # analytical solution for dXdt term (friction)
        Re = ρₗ .* abs.(V) .* d ./ μₗ
        f_coefficient = f_churchill.(Re)
        dXdt_to_stress = -1/8 .* f_coefficient .* ρₗ .* V .* abs.(V)
        rhs_dXdt = peri .* Lliquidslug .* dXdt_to_stress ./ lhs

        # analytical solution for dLdt term (mass conservation)
        dVdt_nodLdt = dXdt_to_stress*peri .* Lliquidslug ./ lhs # if dLdt = 0
        dLdt = (Ac ./ (Ac .- Adeposit_right) - Ac ./ (Ac .- Adeposit_left)) .*  V
        rhs_dLdt = -ρₗ .* Ac .* V .*  dLdt ./ lhs

        # analytical solution for gravity term, gh is the gravitational energy per unit mass
        heightg_interp = sys_tube.mapping.heightg_interp
        Xp1 = [elem[1] for elem in Xp]
        Xp2 = [elem[2] for elem in Xp]
        heightg_Xp1 = heightg_interp(Xp1)
        heightg_Xp2 = heightg_interp(Xp2)
        Δgh = heightg_Xp2 .- heightg_Xp1
        rhs_g_gh = -Δgh .* Ac*ρₗ ./ lhs

        if all(isapprox.(gi,0.0,atol=1e-12))
            @test all(isapprox.(rhs_g_gh,0.0,atol=1e-12))
        else
            @test !all(isapprox.(rhs_g_gh,0.0,atol=1e-12))
        end

        @test all(isapprox.(uu_test1[2*numofslugs+1:2:4*numofslugs-1], rhs_dXdt .+ rhs_dLdt .+ rhs_g_gh, rtol=1e-7))

        @test all(uu_test1[2*numofslugs+2:2:4*numofslugs]   .== uu_test1[2*numofslugs+1:2:4*numofslugs-1])

        @test all(isapprox.(uu_test1[4*numofslugs+1:5*numofslugs],0.0,atol=1e-12)) #dMdt = 0

        @test all(isapprox.(uu_test1[5*numofslugs+1:6*numofslugs],0.0,atol=1e-12)) #dδstart/dt = 0

        F_end = ρₗ .* getδarea.(Ac,d,δend)
        F2_end = ρₗ .* getδarea.(Ac,d,δdep)

        peri_end = peri .* (d .- 2δend)/d
        C_end = ρₗ .* peri_end
        Lfilm_end = sys_tube.vapor.Lfilm_end

        @test all(uu_test1[6*numofslugs+1:7*numofslugs] .≈ (F2_end .- F_end) .* Ac ./ (Ac .- Adeposit_left) .*  V ./ C_end ./ Lfilm_end)

        @test all(uu_test1[7*numofslugs+1:8*numofslugs] .≈ -uu_test1[2:2:2*numofslugs])

        @test all(uu_test1[8*numofslugs+1:9*numofslugs] .≈ uu_test1[1:2:2*numofslugs-1])


        # lastly, just to test if the Re is in laminar range, you get a stress value based on Poiseuille flow
        V_laminar = 1e-1
        Re_laminar = ρₗ .* abs.(V_laminar) .* d ./ μₗ
        @test (Re_laminar < 2300.0)
        f_coefficient_laminar = f_churchill.(Re_laminar)
        dXdt_to_stress_laminar = -1/8 .* f_coefficient_laminar .* ρₗ .* V_laminar .* abs.(V_laminar)
        @test isapprox(dXdt_to_stress_laminar, -8 * μₗ * V_laminar / d, atol=1e-12)
    end
                                                 
end

# test 2, 2 slugs, with heat transfer, given initial temperature differences, no initial velocity
@testset "dynamicsmodel test2" begin

    ηplusvalue = rand()
    ηminusvalue = 0.0
    numofslugs = 2
    sys_tube = initialize_ohpsys(sys_plate,p_fluid,power,slugnum=numofslugs,ηplus=ηplusvalue)

    ΔT = rand()
    ΔT_array = [ΔT,-ΔT]
    sys_tube.vapor.P[1] = sys_tube.propconvert.TtoP(Tref-ΔT_array[1])
    sys_tube.vapor.P[2] = sys_tube.propconvert.TtoP(Tref-ΔT_array[2])

    closedornot = sys_tube.tube.closedornot
    P1 = sys_tube.vapor.P[1]
    P2 = sys_tube.vapor.P[2]
    L = sys_tube.tube.L
    Ac = sys_tube.tube.Ac
    d = sys_tube.tube.d
    ρₗ = sys_tube.liquid.ρₗ
    Xp = sys_tube.liquid.Xp
    Hfg = sys_tube.propconvert.PtoHfg.(sys_tube.vapor.P)
    k = sys_tube.vapor.k
    peri = sys_tube.tube.peri
    δstart = sys_tube.vapor.δstart
    δend = sys_tube.vapor.δend
    Lfilm_start = sys_tube.vapor.Lfilm_start
    Lfilm_end = sys_tube.vapor.Lfilm_end

    u_tube = newstate(sys_tube) # initialize OHP tube

    uu_test2 = dynamicsmodel(u_tube[1:9*numofslugs],sys_tube)


    @test all(isapprox.(uu_test2[1:2:2*numofslugs-1],0.0,atol=1e-12)) # velocity should remain the same
    @test all(isapprox.(uu_test2[2:2:2*numofslugs],0.0,atol=1e-12)) # velocity should remain the same

    Lliquidslug = XptoLliquidslug(Xp,L)
    lhs = ρₗ*Ac .* Lliquidslug
    @test all(uu_test2[2*numofslugs+1:2:4*numofslugs-1] .≈ [1,-1] .* (P1-P2)*Ac ./ lhs) # velocity should remain the same
    @test all(uu_test2[2*numofslugs+2:2:4*numofslugs]  .== uu_test2[2*numofslugs+1:2:4*numofslugs-1])

    # analytical solution for dMdt_latent
    dMdt_latent_start_analytical = Lfilm_start .* peri .* k ./ δstart ./ Hfg .* ΔT_array
    dMdt_latent_end_analytical = Lfilm_end .* peri .* k ./ δend ./ Hfg .* ΔT_array
    dMdt_latent_start_positive_analytical = heaviside.(dMdt_latent_start_analytical) .* dMdt_latent_start_analytical
    dMdt_latent_end_positive_analytical = heaviside.(dMdt_latent_end_analytical) .* dMdt_latent_end_analytical
    dMdt_latent_start_negative_analytical = heaviside.(-dMdt_latent_start_analytical) .* dMdt_latent_start_analytical
    dMdt_latent_end_negative_analytical = heaviside.(-dMdt_latent_end_analytical) .* dMdt_latent_end_analytical

    # run the dynamics model to get numerical dMdt_latent
    Xpvapor = getXpvapor(Xp,closedornot)
    dMdt_latent_start,dMdt_latent_end,dMdt_latent_start_positive,dMdt_latent_end_positive = dMdtdynamicsmodel(Xpvapor,sys_tube)
    dMdt_latent_start_negative = dMdt_latent_start .- dMdt_latent_start_positive
    dMdt_latent_end_negative = dMdt_latent_end .- dMdt_latent_end_positive

    @test all(isapprox.(dMdt_latent_start,dMdt_latent_start_analytical,atol=1e-12)) 
    @test all(isapprox.(dMdt_latent_end,dMdt_latent_end_analytical,atol=1e-12)) 
    @test all(isapprox.(dMdt_latent_start_positive,dMdt_latent_start_positive_analytical,atol=1e-12))
    @test all(isapprox.(dMdt_latent_end_positive,dMdt_latent_end_positive_analytical,atol=1e-12))
    @test all(isapprox.(dMdt_latent_start_negative,dMdt_latent_start_negative_analytical,atol=1e-12))
    @test all(isapprox.(dMdt_latent_end_negative,dMdt_latent_end_negative_analytical,atol=1e-12))

    dMdt_latent = dMdt_latent_start .+ dMdt_latent_end

    @test all(isapprox.(uu_test2[4*numofslugs+1:5*numofslugs],dMdt_latent,atol=1e-12)) #dMdt

    F_start = ρₗ .* getδarea.(Ac,d,δstart)
    F_end = ρₗ .* getδarea.(Ac,d,δend)
    dLdt_start = -(ηplusvalue.*dMdt_latent_start_positive .+ ηminusvalue.*dMdt_latent_start_negative) ./ F_start
    dLdt_end = -(ηplusvalue.*dMdt_latent_end_positive .+ ηminusvalue.*dMdt_latent_end_negative) ./ F_end

    @test all(isapprox.(uu_test2[7*numofslugs+1:8*numofslugs],dLdt_start,atol=1e-12))
    @test all(isapprox.(uu_test2[8*numofslugs+1:9*numofslugs],dLdt_end,atol=1e-12))

    peri_start = peri .* (d .- 2δstart)/d
    peri_end = peri .* (d .- 2δend)/d

    C_start = ρₗ .* peri_start
    C_end = ρₗ .* peri_end
    dδdt_start = -(dMdt_latent_start .+ dLdt_start.*F_start) ./ C_start ./ Lfilm_start
    dδdt_end = -(dMdt_latent_end .+ dLdt_end.*F_end) ./ C_end ./ Lfilm_end

    @test all(isapprox.(uu_test2[5*numofslugs+1:6*numofslugs], dδdt_start,atol=1e-12))
    @test all(isapprox.(uu_test2[6*numofslugs+1:7*numofslugs], dδdt_end,atol=1e-12))

end

# test 3, five slugs, with heat transfer, given different initial film lengthes and liquid velocities to represent five different film states

@testset "five possible film states" begin

    ηplusvalue = rand()
    ηminusvalue = 0.0
    numofslugs = 5


    for Xmag in [-1.0,1.0]
        sys_tube = initialize_ohpsys(sys_plate,p_fluid,power,slugnum=numofslugs,ηplus=ηplusvalue)

        sys_tube.vapor.δstart .= zeros(numofslugs) .+ 1e-5 .+ 1e-5 .* rand(numofslugs) # initial velocity of the slugs
        sys_tube.vapor.δend .= zeros(numofslugs) .+ 1e-5 .+ 1e-5 .* rand(numofslugs) # initial velocity of the slugs
        sys_tube.liquid.dXdt[1:5] = [zero.(X) .+ Xmag .* rand() for X in sys_tube.liquid.dXdt[1:5]]# initial velocity of the slugs


        closedornot = sys_tube.tube.closedornot
        σ = sys_tube.liquid.σ
        L = sys_tube.tube.L
        ρₗ = sys_tube.liquid.ρₗ
        μₗ = sys_tube.liquid.μₗ
        ad_fac = sys_tube.vapor.ad_fac
        d = sys_tube.tube.d
        Ac = sys_tube.tube.Ac
        Xp = sys_tube.liquid.Xp
        Hfg = sys_tube.propconvert.PtoHfg.(sys_tube.vapor.P)
        k = sys_tube.vapor.k
        peri = sys_tube.tube.peri
        δstart = sys_tube.vapor.δstart
        δend = sys_tube.vapor.δend
        Lfilm_start = sys_tube.vapor.Lfilm_start
        Lfilm_end = sys_tube.vapor.Lfilm_end

        Lvaporplug = XptoLvaporplug(Xp,L,closedornot)

        # set up five different film states
        V = [mean(elem) for elem in sys_tube.liquid.dXdt]
        Vavg = mean(abs.(V))
        Ca = getCa.(μₗ,σ,Vavg)

        δdep = Catoδ(d,Ca,adjust_factor=ad_fac)
        Adeposit = getAdeposit(sys_tube,δdep)
        Adeposit_left = [elem[1] for elem in Adeposit]
        Adeposit_right = [elem[2] for elem in Adeposit]
        V_normal_start = circshift(Ac ./ (Ac .- Adeposit_right) .*  V,1)
        V_normal_end = Ac ./ (Ac .- Adeposit_left) .*  V

        Astart = getδarea.(Ac,d,δstart)
        Aend = getδarea.(Ac,d,δend)
        V_start_case5 = Ac / (Ac - Aend[5]) *  V[4]
        V_end_case5   = Ac / (Ac - Astart[4]) *  V[4]


        #case 2 to 5: set different initial film thicknesses
        Lfilm_start[2] = 0.1*sys_tube.wall.L_newbubble
        Lfilm_start[3] = Lvaporplug[3]/2 .- 0.1*sys_tube.wall.L_newbubble
        Lfilm_start[4] = Lvaporplug[4] .- 0.2*sys_tube.wall.L_newbubble
        Lfilm_start[5] = 0.1*sys_tube.wall.L_newbubble

        Lfilm_end[2] = 0.1*sys_tube.wall.L_newbubble
        Lfilm_end[3] = Lvaporplug[3]/2 .- 0.1*sys_tube.wall.L_newbubble
        Lfilm_end[4] = 0.1*sys_tube.wall.L_newbubble
        Lfilm_end[5] = Lvaporplug[5] .- 0.2*sys_tube.wall.L_newbubble

        u_tube = newstate(sys_tube) # initialize OHP tube

        uu_test3 = dynamicsmodel(u_tube[1:9*numofslugs],sys_tube)
        # dX2,ddXdt2,dM2,dδstart2,dδend2,dLfilm_start2,dLfilm_end2 = ComputationalHeatTransfer.vectoXMδL(uu_test2)

        V_start_analytical = zeros(numofslugs)
        V_end_analytical = zeros(numofslugs)

        if Xmag < 0

            V_start_analytical = [V_normal_start[1],
                                V_normal_start[2],
                                V_normal_start[3],
                                V_normal_start[4],
                                V_normal_start[5]]
            V_end_analytical = [V_normal_end[1],
                                V[2],
                                V_normal_end[3],
                                V_end_case5,
                                V_normal_end[5]]
            else

            V_start_analytical = [V_normal_start[1],
                            V[1],
                            V_normal_start[3],
                            V_normal_start[4],
                            V_start_case5]
            V_end_analytical = [V_normal_end[1],
                            V_normal_end[2],
                            V_normal_end[3],
                            V_normal_end[4],
                            V_normal_end[5]]
        end

        @test all(isapprox.(uu_test3[1:2:2*numofslugs-1],V_end_analytical,atol=1e-12))
        @test all(isapprox.(uu_test3[2:2:2*numofslugs],circshift(V_start_analytical,-1),atol=1e-12))

    end

end

# test 4, 2 slugs, with heat transfer, given initial temperature differences, with initial velocity, test mass conservation δ L
@testset "mass conservation δ L" begin

    ηplusvalue = rand()
    ηminusvalue = 0.0
    numofslugs = 2
    sys_tube = initialize_ohpsys(sys_plate,p_fluid,power,slugnum=numofslugs,ηplus=ηplusvalue)
    sys_tube.liquid.dXdt = [zero.(X) .+ 1.0 .* rand() for X in sys_tube.liquid.dXdt]# initial velocity of the slugs
   
    ΔT = rand()
    ΔT_array = [ΔT,-ΔT]
    sys_tube.vapor.P[1] = sys_tube.propconvert.TtoP(Tref-ΔT_array[1])
    sys_tube.vapor.P[2] = sys_tube.propconvert.TtoP(Tref-ΔT_array[2])

    closedornot = sys_tube.tube.closedornot
    P1 = sys_tube.vapor.P[1]
    P2 = sys_tube.vapor.P[2]
    σ = sys_tube.liquid.σ
    L = sys_tube.tube.L
    ρₗ = sys_tube.liquid.ρₗ
    μₗ = sys_tube.liquid.μₗ
    ad_fac = sys_tube.vapor.ad_fac
    d = sys_tube.tube.d
    Ac = sys_tube.tube.Ac
    Xp = sys_tube.liquid.Xp
    Hfg = sys_tube.propconvert.PtoHfg.(sys_tube.vapor.P)
    k = sys_tube.vapor.k
    peri = sys_tube.tube.peri
    δstart = sys_tube.vapor.δstart
    δend = sys_tube.vapor.δend
    Lfilm_start = sys_tube.vapor.Lfilm_start
    Lfilm_end = sys_tube.vapor.Lfilm_end

    u_tube = newstate(sys_tube) # initialize OHP tube

    uu_test2 = dynamicsmodel(u_tube[1:9*numofslugs],sys_tube)


    V = [mean(elem) for elem in sys_tube.liquid.dXdt]
    # get a characteristic Capilarry number based on the average velocities
    Vavg = mean(abs.(V))
    Ca = getCa.(μₗ,σ,Vavg)

    δdep = Catoδ(d,Ca,adjust_factor=ad_fac)
    Adeposit = getAdeposit(sys_tube,δdep)
    Adeposit_left = [elem[1] for elem in Adeposit]
    Adeposit_right = [elem[2] for elem in Adeposit]
    Adeposit_start = circshift(Adeposit_right,1)
    Adeposit_end = Adeposit_left

    # run the dynamics model to get numerical dMdt_latent
    Xpvapor = getXpvapor(Xp,closedornot)
    dMdt_latent_start,dMdt_latent_end,dMdt_latent_start_positive,dMdt_latent_end_positive = dMdtdynamicsmodel(Xpvapor,sys_tube)
    dMdt_latent_start_negative = dMdt_latent_start .- dMdt_latent_start_positive
    dMdt_latent_end_negative = dMdt_latent_end .- dMdt_latent_end_positive

    dMdt_latent = dMdt_latent_start .+ dMdt_latent_end

    @test all(isapprox.(uu_test2[4*numofslugs+1:5*numofslugs],dMdt_latent,atol=1e-12)) #dMdt

    F_start = ρₗ .* getδarea.(Ac,d,δstart)
    F_end = ρₗ .* getδarea.(Ac,d,δend)
    dLdt_start_analytical = -(ηplusvalue.*dMdt_latent_start_positive .+ ηminusvalue.*dMdt_latent_start_negative) ./ F_start .- Ac ./ (Ac .- Adeposit_start) .*  circshift(V,-1)
    dLdt_end_analytical  = -(ηplusvalue.*dMdt_latent_end_positive .+ ηminusvalue.*dMdt_latent_end_negative) ./ F_end .+ Ac ./ (Ac .- Adeposit_end) .*  V

    dLdt_start = uu_test2[7*numofslugs+1:8*numofslugs]
    dLdt_end = uu_test2[8*numofslugs+1:9*numofslugs]
    @test all(isapprox.(dLdt_start_analytical,dLdt_start,atol=1e-12))
    @test all(isapprox.(dLdt_end_analytical,dLdt_end,atol=1e-12))

    peri_start = peri .* (d .- 2δstart)/d
    peri_end = peri .* (d .- 2δend)/d

    C_start = ρₗ .* peri_start
    C_end = ρₗ .* peri_end
    dδdt_start = -(dMdt_latent_start .+ dLdt_start.*F_start) ./ C_start ./ Lfilm_start
    dδdt_end = -(dMdt_latent_end .+ dLdt_end.*F_end) ./ C_end ./ Lfilm_end

    dδdt_start = uu_test2[5*numofslugs+1:6*numofslugs]
    dδdt_end = uu_test2[6*numofslugs+1:7*numofslugs]
    # @test all(isapprox.(uu_test2[5*numofslugs+1:6*numofslugs], dδdt_start,atol=1e-12))
    # @test all(isapprox.(uu_test2[6*numofslugs+1:7*numofslugs], dδdt_end,atol=1e-12))

    dMdt_start_lhs = F_start .* dLdt_start .+ C_start .* Lfilm_start .* dδdt_start
    dMdt_end_lhs = F_end .* dLdt_end .+ C_end .* Lfilm_end .* dδdt_end
    dMdt_start_rhs = -dMdt_latent_start - ρₗ .* Adeposit_start .* Ac ./ (Ac .- Adeposit_start) .*  circshift(V,-1)
    dMdt_end_rhs = -dMdt_latent_end + ρₗ .* Adeposit_end .* Ac ./ (Ac .- Adeposit_end) .*  V

    @test all(isapprox.(dMdt_start_lhs,dMdt_start_rhs,atol=1e-12))
    @test all(isapprox.(dMdt_end_lhs,dMdt_end_rhs,atol=1e-12))
end

@testset "liquidmodel" begin

    # p_fluid = SaturationFluidProperty("Butane",Tref)

    ηplusvalue = rand()
    ηminusvalue = 0.0
    numofslugs = 2
    sys_tube = initialize_ohpsys(sys_plate,p_fluid,power,slugnum=numofslugs,ηplus=ηplusvalue)

    closedornot = sys_tube.tube.closedornot

    N_varyT = 10
    dTs = [1.0,-1.0]

    sys_tube.liquid.θarrays[1][N_varyT] += dTs[1]
    sys_tube.liquid.θarrays[2][N_varyT] += dTs[2]

    # use the liquidmodel function to get the du
    dus = OscillatingHeatPipe.liquidmodel(sys_tube)


    α = sys_tube.liquid.αₗ
    # k = sys_tube.liquid.kₗ
    Cpₗ = sys_tube.liquid.Cpₗ
    Hₗ = sys_tube.liquid.Hₗ
    ρₗ = sys_tube.liquid.ρₗ
    peri = sys_tube.tube.peri
    Ac = sys_tube.tube.Ac

    @test all(isapprox.(α,p_fluid.kₗ/(Cpₗ*ρₗ),atol=1e-12))

    H_rhs = peri / (ρₗ*Cpₗ*Ac)

    # hand calculate the du
    θarrays = sys_tube.liquid.θarrays
    dx1 = mod(sys_tube.liquid.Xarrays[1][2] - sys_tube.liquid.Xarrays[1][1], sys_tube.tube.L)

    dθarray1 = zero(θarrays[1])
    dθarray1[N_varyT-1] = α .* dTs[1] ./ dx1 ./ dx1
    dθarray1[N_varyT] = α .* -2dTs[1] ./ dx1 ./ dx1 - Hₗ .* dTs[1] .* peri ./ (ρₗ*Cpₗ*Ac)
    dθarray1[N_varyT+1] = α .* dTs[1] ./ dx1 ./ dx1

    @test all(isapprox.(dθarray1,dus[1],atol=1e-12))


    dx2 = mod(sys_tube.liquid.Xarrays[2][2] - sys_tube.liquid.Xarrays[2][1], sys_tube.tube.L)

    dθarray2 = zero(θarrays[2])
    dθarray2[N_varyT-1] = α .* dTs[2] ./ dx2 ./ dx2
    dθarray2[N_varyT] = α .* -2dTs[2] ./ dx2 ./ dx2 - Hₗ .* dTs[2] .* peri ./ (ρₗ*Cpₗ*Ac)
    dθarray2[N_varyT+1] = α .* dTs[2] ./ dx2 ./ dx2

    @test all(isapprox.(dθarray2,dus[2],atol=1e-12))

    u_tube = newstate(sys_tube) # initialize OHP tube
    


    liquiddu = duliquidθtovec(dus)
    liquiddu_analytical = [0.0; dθarray1;
                           0.0; dθarray2]


    @test all(isapprox.(liquiddu,liquiddu_analytical,atol=1e-12))    
end






# test 5, still under construction still 30 slugs, with heat transfer, testing mass conservation for these five states
@testset "five possible film states" begin


    tspan = (0.0, 1.0); # start time and end time
    dt_record = 1.0   # saving time interval
    tstep = 4e-4 # one-step mass-conservation check

    ηplusvalue = 0.5
    ηminusvalue = 0.0
    numofslugs = 30
    Rn_disable = 1e-8
    sys_tube = initialize_ohpsys(sys_plate,p_fluid,power,Rn_boil=Rn_disable, slugnum=numofslugs,ηplus=ηplusvalue)

    for casei in [1,2,3,4,5]
        sys_tube.vapor.δstart .= zeros(numofslugs) .+ 1e-5 .+ 0e-5 .* rand(numofslugs) # initial velocity of the slugs
        sys_tube.vapor.δend .= zeros(numofslugs) .+ 1e-5 .+ 0e-5 .* rand(numofslugs) # initial velocity of the slugs
        # 0.0 need to be replaced by 1.0
        sys_tube.liquid.dXdt[1:numofslugs] = [zero.(X) .- 1 for X in sys_tube.liquid.dXdt[1:numofslugs]]# initial velocity of the slugs
        if casei == 1
            sys_tube.liquid.dXdt[1:numofslugs] = [zero.(X) .- 1 for X in sys_tube.liquid.dXdt[1:numofslugs]]# initial velocity of the slugs
            sys_tube.vapor.P .= fill(sys_tube.propconvert.TtoP(Tref - 5.0),numofslugs)
        elseif casei == 2
            sys_tube.liquid.dXdt[1:numofslugs] = [zero.(X) .- 1 for X in sys_tube.liquid.dXdt[1:numofslugs]]# initial velocity of the slugs
            sys_tube.vapor.P .= fill(sys_tube.propconvert.TtoP(Tref - 5.0),numofslugs)
            Lvaporplug = XptoLvaporplug(sys_tube.liquid.Xp,sys_tube.tube.L,sys_tube.tube.closedornot)
            sys_tube.vapor.Lfilm_start .= Lvaporplug ./ 2
            sys_tube.vapor.Lfilm_end .= 1e-1 * sys_tube.wall.L_newbubble
        elseif casei == 3
            sys_tube.liquid.dXdt[1:numofslugs] = [zero.(X) .- 1 for X in sys_tube.liquid.dXdt[1:numofslugs]]# initial velocity of the slugs
            sys_tube.vapor.P .= fill(sys_tube.propconvert.TtoP(Tref + 5.0),numofslugs)
            Lvaporplug = XptoLvaporplug(sys_tube.liquid.Xp,sys_tube.tube.L,sys_tube.tube.closedornot)
            sys_tube.vapor.Lfilm_start .= Lvaporplug ./ 2 .- 1e-1 * sys_tube.wall.L_newbubble
            sys_tube.vapor.Lfilm_end .= Lvaporplug ./ 2 .- 1e-1 * sys_tube.wall.L_newbubble
        elseif casei == 4
            sys_tube.liquid.dXdt[1:numofslugs] = [zero.(X) .- 1 for X in sys_tube.liquid.dXdt[1:numofslugs]]# initial velocity of the slugs
            sys_tube.vapor.P .= fill(sys_tube.propconvert.TtoP(Tref + 5.0),numofslugs)
            Lvaporplug = XptoLvaporplug(sys_tube.liquid.Xp,sys_tube.tube.L,sys_tube.tube.closedornot)
            sys_tube.vapor.Lfilm_start .= Lvaporplug .- 2e-1 * sys_tube.wall.L_newbubble
            sys_tube.vapor.Lfilm_end .= 1e-1 * sys_tube.wall.L_newbubble
        elseif casei == 5
            sys_tube.liquid.dXdt[1:numofslugs] = [zero.(X) .+ 1 for X in sys_tube.liquid.dXdt[1:numofslugs]]# initial velocity of the slugs
            sys_tube.vapor.P .= fill(sys_tube.propconvert.TtoP(Tref + 5.0),numofslugs)
            Lvaporplug = XptoLvaporplug(sys_tube.liquid.Xp,sys_tube.tube.L,sys_tube.tube.closedornot)
            sys_tube.vapor.Lfilm_start .= 1e-1 * sys_tube.wall.L_newbubble
            sys_tube.vapor.Lfilm_end .= Lvaporplug .- 2e-1 * sys_tube.wall.L_newbubble
        end
        # sys_tube.vapor.P .= fill(sys_tube.propconvert.TtoP(Tref - 5.0),numofslugs)
    
        u_tube = newstate(sys_tube) # initialize OHP tube 
        integrator_tube = init(u_tube,tspan,deepcopy(sys_tube)); # construct integrator_tube

        u_plate = init_sol(sys_plate) # initialize plate T field to uniform Tref
        integrator_plate = init(u_plate,tspan,sys_plate) # construct integrator_plate

        SimuResult = SimulationResult(integrator_tube,integrator_plate);

        p_old = deepcopy(integrator_tube.p)
        Mvapor_old = sum(getMvapor(p_old))
        Mfilm_old = sum(sum.(getMfilm(p_old)))
        Mliquid_old = sum(getMliquid(p_old))
        Mold = Mvapor_old + Mfilm_old + Mliquid_old

        uu_test1 = dynamicsmodel(u_tube[1:9*numofslugs],sys_tube)

        for t in tspan[1]:tstep:tspan[1]
    
            timemarching!(integrator_tube,integrator_plate,tstep)

        end

        p_new = deepcopy(integrator_tube.p)
        Mvapor_new = sum(getMvapor(p_new))
        Mfilm_new = sum(sum.(getMfilm(p_new)))
        Mliquid_new = sum(getMliquid(p_new))
        Mnew = Mvapor_new + Mfilm_new + Mliquid_new

        println(1-Mnew/Mold)

        @test isapprox.(Mold,Mnew,rtol=1e-9)
    end

end




@testset "liquid Nu and Hₗ" begin

    Nuₗ=3.6

    sys_tube = initialize_ohpsys(sys_plate,p_fluid,power,Nu=Nuₗ)

    Hₗ_analytical = Nuₗ * p_fluid.kₗ / sys_tube.tube.d

    @test all(isapprox.(sys_tube.liquid.Hₗ,Hₗ_analytical,atol=1e-12))    
end

@testset "wall heat flux" begin

    numofslugs = 2
    sys_tube = initialize_ohpsys(sys_plate,p_fluid,power,slugnum=numofslugs)
    sys_tube.wall.θarray .+= [rand() for i in 1:length(sys_tube.wall.θarray)]
    u_tube = newstate(sys_tube) # initialize OHP tube
    sys_tube = getcurrentsys_nowall!(u_tube,sys_tube)

    qwall = sys_to_heatflux(sys_tube)

    xs = sys_tube.wall.Xarray

    @test all(isapprox.(sys_tube.wall.θarray,sys_tube.mapping.θ_interp_walltoliquid(xs),atol=1e-12))


    qwall_analytical = sys_tube.mapping.H_interp_liquidtowall(xs) .* (sys_tube.mapping.θ_interp_walltoliquid(xs) .- sys_tube.mapping.θ_interp_liquidtowall(xs)) .* sys_tube.tube.peri

    @test all(isapprox.(qwall,qwall_analytical,atol=1e-12))    
end


@testset "heater area (for Rectangle only)" begin
    total_heater_area = (maximum(heater1.shape.x) - minimum(heater1.shape.x)) * (maximum(heater1.shape.y) - minimum(heater1.shape.y))
    @test isapprox(total_heater_area, areaheater_area, atol=1e-12)
end