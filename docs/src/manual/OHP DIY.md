```@meta
EditURL = "../../../test/literate/OHP DIY.jl"
```

  # DIY

  This notebook shows how to customize the heater/condenser and ohp configuration

  # Packages

  Firstly, let's import the necessary packages, you may need to install them
  for the first time.

````@example OHP_DIY
using OscillatingHeatPipe # our main package
using Plots # for plotting
using ProgressMeter # to have a progress bar in the calculation
````

  # Specify properties

  ### Solid Physical parameters

  The numbers below represent aluminum alloy 3003

````@example OHP_DIY
ρₛ = 2730; # material density [kg/m^3]
cₛ  = 8.93e02; # material specific heat [J/kg K]
kₛ  = 1.93e02; # material heat conductivity [W/m K]
plate_d = 1.5e-3; # effective d [m] (The thickness of an ideal uniform thickness plate occupying the same volume)
αₛ = kₛ/ρₛ/cₛ

Tref = 291.2 # reference temperature [K]

power = 70 # [W], total power
Lheater_x = 50e-3 # [m], length of heater along x axis
Lheater_y = 50e-3 # [m], length of heater along y axis
areaheater_area = Lheater_x * Lheater_y # [m^2] total area

phys_params = Dict( "diffusivity"              => αₛ,
                    "flux_correction"          => ρₛ*cₛ*plate_d,
                    "Fourier"                  => 1.0,
                    "ohp_flux"                 => [NaN], # initial value, the value here is useless
                    "areaheater_power"         => power, # total power
                    "areaheater_area"          => areaheater_area, # total area
                    "areaheater_temp"          => 0.0,   # relative temperature compared with "background temperature"
                    "areaheater_coeff"         => 4000.0,
                    "background temperature"   => Tref
                     )
````

  ### Fluid Physical parameters

  Here, we set the p_fluid contains the vapor and liquid properties at a constant reference
  temperature. Noted that the vapor pressure and the vapor density will be
  functions of temperatures during the simulation, other properties are
  extracted from pfluid as an approximate value.

````@example OHP_DIY
fluid_type = "Butane"
p_fluid = SaturationFluidProperty(fluid_type,Tref)
````

  # Set the geometries

  ### Geometry parameters

  The 2D domain is of rectangular shape (slightly different from ASETS-II). In
  the future it can be of arbitrary shape using the immersedlayers.jl package.

````@example OHP_DIY
Δx = 0.0007 # [m] # grid size, at the same order of 1D OHP channel node spacing ~ 0.001[m]
Lx = 6*INCHES*1.02; # plate size x [m]
Ly = 2*INCHES*1.05; # plate size y [m]
xlim = (-Lx/2,Lx/2) # plate x limits
ylim = (-Ly/2,Ly/2) # plate y limits

g = PhysicalGrid(1.03 .* xlim,1.1 .* ylim,Δx) # build a gird slightly larger than the plate

xbound = [ -Lx/2,-Lx/2,
            Lx/2, Lx/2] # x coordinates of the shape

ybound = [  Ly/2,-Ly/2,
           -Ly/2, Ly/2] # y coordinates of the shape

Δs = 1.4*cellsize(g)
body = Polygon(xbound,ybound,Δs)

X = MotionTransform([0,0],0) # move the plate or rotate the plate
joint = Joint(X)
m = RigidBodyMotion(joint,body)
x = zero_motion_state(body,m)
update_body!(body,x,m)

function get_qbplus(t,x,base_cache,phys_params,motions)
    nrm = normals(base_cache)
    qbplus = zeros_surface(base_cache)
    return qbplus
end

function get_qbminus(t,x,base_cache,phys_params,motions)
    nrm = normals(base_cache)
    qbminus = zeros_surface(base_cache)
    return qbminus
end

bcdict = Dict("exterior" => get_qbplus,"interior" => get_qbminus)
````

  # Set up the evaporators and condensers

  In the "OHP simulation" notebook, I use "OHPtype" to look up a preset dictionary of OHP evaporators and condensers.

  You can also customize them, following the procedure below in this notebook.

  Firstly let's give the total heater power

````@example OHP_DIY
function heatermodel!(σ,T,t,fr::AreaRegionCache,phys_params)
    σ .= phys_params["areaheater_power"] / phys_params["areaheater_area"] / phys_params["flux_correction"]
end

function condensermodel!(σ,T,t,fr::AreaRegionCache,phys_params)
    T0 = phys_params["areaheater_temp"]
    h = phys_params["areaheater_coeff"]
    corr = phys_params["flux_correction"]

    σ .= h*(T0 - T) / corr
end
````

Then let's construct a heater

````@example OHP_DIY
eb1 = Rectangle(Lheater_x/2,Lheater_x/2,1.4*Δx)
tr1_h = RigidTransform((0.0,-0.0),0.0)
heater1 = AreaForcingModel(eb1,tr1_h,heatermodel!);
nothing #hide
````

Then let's consctruct a condenser

````@example OHP_DIY
Lcondenser_x = 15e-3
Lcondenser_y = 1.0INCHES
cb1 = Rectangle(Lcondenser_x,Lcondenser_y,1.4*Δx)
tr1_c = RigidTransform((2.4INCHES,-0.0),0.0)
cond1 = AreaForcingModel(cb1,tr1_c,condensermodel!);
nothing #hide
````

# Set up OHP channel's shape

Similarly, In the "OHP simulation" notebook, I used **construct_ohp_curve("ASETS",Δx)** to look up a preset dictionary of ASETS-II OHP.

You can customize the ohp curve in either of the two ways:

 1. simply supply two arrays of x and y of the same length:

````@example OHP_DIY
a = 0.03
θ = 0:2π/1000:2π
r = a*sin.(2θ)
x = r .* cos.(θ)
y = r .* sin.(θ);

plot(x,y,aspectratio=1)
````

2. **construct_ohp_curve(nturn, pitch, height, gap, ds, x0, y0, flipx, flipy, angle)**, a built-in function to generate a closed loop multi-turn channel

````@example OHP_DIY
ds = 1.5Δx # point interval
nturn = 9 # number of turns
width_ohp = 46.25*1e-3
length_ohp = 147.0*1e-3
gap = 3e-3 # gap between the closed loop end to the channel(not the distance between each channels)
pitch = width_ohp/(2*nturn+1) # pitch between channels
rotation_angle = 3π/2
x0, y0 = length_ohp/2 + 2e-3, width_ohp/2

x,y, xf, yf = construct_ohp_curve(nturn,pitch,length_ohp,gap,ds,x0,y0,false,false,rotation_angle)

plot(x,y,aspectratio=1)
````

Now set the coordinates into a "body", which allows us to transform it in various helpful ways.
Here, we will just place it at the origin

````@example OHP_DIY
ohp = BasicBody(x,y) # build a BasicBody based on x,y
tr_ohp = RigidTransform((0.0,0.0),0.0)
````

We also need a function that will introduce the OHP heat flux into the plate

````@example OHP_DIY
function ohpmodel!(σ,T,t,fr::LineRegionCache,phys_params)
    σ .= phys_params["ohp_flux"] ./ phys_params["flux_correction"]
end
ohp_linesource = LineForcingModel(ohp,tr_ohp,ohpmodel!);
nothing #hide
````

  ### Plot what you have so far

  This is a exmaple of the compuational domain (the box) and the OHP channel
  serpentine (in blue)

````@example OHP_DIY
plot(body,fillalpha=0)
update_body!(eb1,tr1_h)
update_body!(cb1,tr1_c)
plot!(eb1)
plot!(cb1)
update_body!(ohp,tr_ohp)

plot!(ohp,fillalpha=0,closedornot=true)
````

 Assemble into a forcing list

````@example OHP_DIY
forcing_dict = Dict("heating models" => [heater1,cond1,ohp_linesource])
````

  # Construct the systems
Now we will set up the thermal conduction problem, and then set up the data structures
for the plate and OHP channels

 ### Set time step
We first set the time step size (in seconds), and a function that will supply this time step.

````@example OHP_DIY
tstep = 1e-3
timestep_fixed(u,sys) = tstep
````

### Create heat conduction system
The solid module dealing with the 2D conduction, evaporator, condenser, and
the OHP line heat source is constructed here.

````@example OHP_DIY
prob = NeumannHeatConductionProblem(g,body,phys_params=phys_params,
                                           bc=bcdict,
                                           motions=m,
                                           forcing=forcing_dict,
                                           timestep_func=timestep_fixed);
nothing #hide
````

The `sys_plate` structure contains everything about the plate

````@example OHP_DIY
sys_plate = construct_system(prob);
nothing #hide
````

 ### Create OHP channel system
 The `sys_tube` structure contains everything about the OHP channels and fluid

````@example OHP_DIY
sys_tube = initialize_ohpsys(sys_plate,p_fluid,power);
nothing #hide
````

 # Initialize the problem
We set intial conditions for the plate and channel here
For plate, the time span should be a range larger than the TOTAL time you plan to simulate (including saving and re-run),
if the range is smaller than the total time range, there will be errors in temperature interpolations

````@example OHP_DIY
tspan_init = (0.0,1e4)
u_plate = init_sol(sys_plate)# initialize plate T field to uniform Tref
integrator_plate = init(u_plate,tspan_init,sys_plate,save_on=false) # construct integrator_plate
````

Set the tubes time span for simulation and its initial condition

````@example OHP_DIY
tspan = (0.0, 1.0); # start time and end time
dt_record = 0.2   # saving time interval
u_tube = newstate(sys_tube) # initialize OHP tube
integrator_tube = init(u_tube,tspan,sys_tube); # construct integrator_tube
nothing #hide
````

  ### initialize arrays for saving

````@example OHP_DIY
SimuResult = SimulationResult(integrator_tube,integrator_plate);
nothing #hide
````

  # Solve
  Run the simulation and store data

````@example OHP_DIY
@showprogress for t in tspan[1]:tstep:tspan[2]

    timemarching!(integrator_tube,integrator_plate,tstep)

    if (mod(integrator_plate.t,dt_record) < 1e-6) || (mod(-integrator_plate.t,dt_record) < 1e-6)
        store!(SimuResult,integrator_tube,integrator_plate)
    end

end
````

  # Store data

````@example OHP_DIY
#save_path = "../numedata/solution.jld2"
#save(save_path,"SimulationResult",SimuResult)
````

### take a peek at the solution (more at the PostProcessing notebook)
First, a movie of temperature in the plate

````@example OHP_DIY
@gif for i in eachindex(SimuResult.tube_hist_t)
    plot(OHPTemp(),i,SimuResult,clim=(291.2,294.0))
    plot!(body,fillalpha=0)
end
````

Show a movie of the channels and the locations of slugs/film vapor/dry vapor

````@example OHP_DIY
@gif for i in eachindex(SimuResult.tube_hist_t)
    plot(OHPSlug(),i,SimuResult,aspectratio=1)
    plot!(body,fillalpha=0)
end
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

