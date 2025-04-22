@testset "Liquid and vapor arrays" begin

    L, Lmin, ϕ0 = 4.0 + rand(), 0.001, 0.6 + 0.1*rand()
    closedornot = true
    X, dXdt, realratio = randomXp(L,Lmin,closedornot;chargeratio=ϕ0)

    lslug = mod.(map(u -> u[2],X) - map(u -> u[1],X),L)
    @test all(lslug .> 0)

    ϕ = sum(lslug)/L
    @test abs(ϕ - ϕ0) < 0.1
    @test ϕ ≈ realratio

    lvap = XptoLvaporplug(X,L,closedornot)
    
    # Continuity test
    @test sum(lvap) + sum(lslug) ≈ L


end

@testset "Plate" begin

    ρₛ = 2730 # material density [kg/m^3]
    cₛ  = 8.93e02 # material specific heat [J/kg K]
    kₛ  = 1.93e02 # material heat conductivity
    αₛ = kₛ/ρₛ/cₛ
    
    dₛ = 1.5e-3

    Tref = 291.2 # reference temperature
    fluid_type = "Butane"
    p_fluid = SaturationFluidProperty(fluid_type,Tref) # This function relies on CoolProp.jl package

    power = 70 # [W], total power
    areaheater_area = 50e-3 * 50e-3 # [m] total area

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
    

    phys_params = Dict( "diffusivity"              => αₛ,
                    "flux_correction"          => ρₛ*cₛ*dₛ,
                    # "angular velocity"         => 0.0,
                    "Fourier"                  => 1.0,
                    "ohp_flux"                 => [NaN], # initial value, the value here is useless
                    "areaheater_power"         => power, # total power
                    "areaheater_area"          => areaheater_area, # total area
                    "areaheater_temp"          => 0.0,   # relative temperature compared with "background temperature"
                    "areaheater_coeff"         => 4000.0,
                    "background temperature"   => Tref
                     )

    Δx = 0.0007 # [m] # grid size, at the same order of 1D OHP channel node spacing ~ 0.001[m]

    Lx = 6*INCHES*1.02 # plate size x [m]
    Ly = 2*INCHES*1.05 # plate size y [m]
    xlim = (-Lx/2,Lx/2) # plate x limits
    ylim = (-Ly/2,Ly/2) # plate y limits
    
    g = PhysicalGrid(1.03 .* xlim,1.1 .* ylim,Δx)

    Δs = 1.4*cellsize(g) # 1D OHP node spacing, here it is 1.4Δx

    xbound = [ -Lx/2,-Lx/2, 
             Lx/2, Lx/2] # x coordinates of the shape

    ybound = [  Ly/2,-Ly/2, 
            -Ly/2, Ly/2] # y coordinates of the shape

    body = Polygon(xbound,ybound,Δs)
    
    X = MotionTransform([0,0],0) # move the plate or rotate the plate
    joint = Joint(X)
    m = RigidBodyMotion(joint,body)
    x = zero_motion_state(body,m)
    update_body!(body,x,m)

    function heatermodel!(σ,T,t,fr::AreaRegionCache,phys_params)
        σ .= phys_params["areaheater_power"] / phys_params["areaheater_area"] / phys_params["flux_correction"] 
    end

    function condensermodel!(σ,T,t,fr::AreaRegionCache,phys_params)
        T0 = phys_params["areaheater_temp"]
        h = phys_params["areaheater_coeff"]
        corr = phys_params["flux_correction"] 
    
        σ .= h*(T0 - T) / corr
    end

    fregion1_h = Rectangle(25e-3,25e-3,1.4*Δx)
    tr1_h = RigidTransform((0.0,-0.0),0.0)
    heater1 = AreaForcingModel(fregion1_h,tr1_h,heatermodel!)

    fregion1_c = Rectangle(15e-3,1.0INCHES,1.4*Δx)
    tr1_c = RigidTransform((2.4INCHES,-0.0),0.0)
    cond1 = AreaForcingModel(fregion1_c,tr1_c,condensermodel!)

    ds = 1.5Δx
    nturn = 13
    width_ohp = 46.25*1e-3
    length_ohp = 147.0*1e-3
    gap = 3e-3
    pitch = width_ohp/(2*nturn+1)
    x0, y0 = length_ohp/2 +2e-3, width_ohp/2

    x, y, xf, yf = construct_ohp_curve(nturn,pitch,length_ohp,gap,ds,x0,y0,false,false,3pi/2)
    ohp = BasicBody(x,y) # build a BasicBody based on x,y
    tr_ohp = RigidTransform((0.0,0.0),0.0)

    function ohpmodel!(σ,T,t,fr::LineRegionCache,phys_params)
        σ .= phys_params["ohp_flux"] ./ phys_params["flux_correction"] 
    end
    ohp_linesource = LineForcingModel(ohp,tr_ohp,ohpmodel!)

    forcing_dict = Dict("heating models" => [heater1,cond1,ohp_linesource])
    
    timestep_fixed(u,sys) = tstep

    prob = NeumannHeatConductionProblem(g,body,scaling=GridScaling,
                                             phys_params=phys_params,
                                             bc=bcdict,
                                             motions=m,
                                             forcing=forcing_dict,
                                             # timestep_func=timestep_fourier
                                             timestep_func=timestep_fixed)

    sys_plate = construct_system(prob)

    ohp = OscillatingHeatPipe._get_ohp_from_forcing_list(sys_plate)
    @test ohp == ohp_linesource

                                                 
end