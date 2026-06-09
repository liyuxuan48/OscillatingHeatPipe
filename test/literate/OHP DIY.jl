#   # DIY Geometry And Channel Design

#   This notebook shows how to build a custom OHP setup by defining the plate
#   boundary, heater/condenser regions, and OHP channel path directly.

#   ### What do we need to customize an OHP problem?

#   **specify properties** : Solid property, fluid property, heater power, and
#   condenser strength
#
#   **set the geometries** : A custom plate polygon, heater/condenser patches,
#   and an OHP channel path
#
#   **construct the systems** : Fluid system (1D) and heat-conduction system
#   (2D)
#
#   **initialize** : Build the integrators and the data structure for saving
#
#   **solve** : March the coupled tube and plate systems forward in time
#
#   **save/examine** : Save the data and preview the custom result

#   # Packages

# Load the OHP package and plotting/progress utilities used throughout the notebook.

using OscillatingHeatPipe # our main package
using Plots # for plotting
using ProgressMeter # to have a progress bar in the calculation
using JLD2 # to save and load data

#   # Specify Properties

#   ### Solid Physical Parameters

# Define aluminum plate material properties and compute the thermal diffusivity used by the immersed-layer heat-conduction model.

ρₛ = 2730; # material density [kg/m^3]
cₛ  = 8.93e02; # material specific heat [J/kg K]
kₛ  = 1.93e02; # material heat conductivity [W/m K]
plate_d = 1.5e-3; # effective plate thickness [m]
αₛ = kₛ/ρₛ/cₛ

Tref = 291.2 # reference temperature [K]
power = 4 # total heater power [W]

#   ### Fluid Physical Parameters

# Set the working fluid and construct the saturation-property helper used by the one-dimensional tube model.

fluid_type = "Butane"
p_fluid = SaturationFluidProperty(fluid_type,Tref)

#   # Set The Geometry

#   ### Grid And Plate Boundary

# Create a slightly oversized computational grid and a deliberately irregular polygonal plate that still encloses the heaters, condensers, and tilted OHP channel.

Δx = 0.0007 # [m]
Lx = 6.0INCHES
Ly = 2.2INCHES
xlim = (-0.68Lx,0.68Lx)
ylim = (-1.08Ly,0.86Ly)

g = PhysicalGrid(1.04 .* xlim,1.08 .* ylim,Δx)
Δs = 1.4*cellsize(g)

xbound = [
    -0.62Lx, -0.67Lx, -0.56Lx, -0.20Lx,
     0.10Lx,  0.45Lx,  0.66Lx,  0.61Lx,
     0.32Lx, -0.05Lx, -0.42Lx
]
ybound = [
     0.52Ly, -0.08Ly, -0.84Ly, -1.00Ly,
    -0.64Ly, -0.90Ly, -0.32Ly,  0.48Ly,
     0.76Ly,  0.62Ly,  0.70Ly
]
body = Polygon(xbound,ybound,Δs)

# Create a stationary rigid-body motion object for the irregular plate boundary.

X = MotionTransform([0,0],0)
joint = Joint(X)
m = RigidBodyMotion(joint,body)
x_motion = zero_motion_state(body,m)
update_body!(body,x_motion,m)

# Define zero-flux boundary-condition functions for the exterior and interior immersed boundaries.

function get_qbplus(t,x,base_cache,phys_params,motions)
    qbplus = zeros_surface(base_cache)
    return qbplus
end

function get_qbminus(t,x,base_cache,phys_params,motions)
    qbminus = zeros_surface(base_cache)
    return qbminus
end

bcdict = Dict("exterior" => get_qbplus,"interior" => get_qbminus)

#   # Set Up Heaters And Condensers

# Collect the physical parameters passed into the heat-conduction system, including the total heater area and condenser coefficient.

heater_size = (0.62INCHES,0.48INCHES)
areaheater_area = prod(heater_size)

phys_params = Dict( "diffusivity"              => αₛ,
                    "flux_correction"          => ρₛ*cₛ*plate_d,
                    "Fourier"                  => 1.0,
                    "ohp_flux"                 => [NaN],
                    "areaheater_power"         => power,
                    "areaheater_area"          => areaheater_area,
                    "areaheater_temp"          => 0.0,
                    "areaheater_coeff"         => 3600.0,
                    "background temperature"   => Tref
                     )

# Define the area-source models for heater input and condenser cooling.

function heatermodel!(σ,T,t,fr::AreaRegionCache,phys_params)
    σ .= phys_params["areaheater_power"] / phys_params["areaheater_area"] / phys_params["flux_correction"]
end

function condensermodel!(σ,T,t,fr::AreaRegionCache,phys_params)
    T0 = phys_params["areaheater_temp"]
    h = phys_params["areaheater_coeff"]
    corr = phys_params["flux_correction"]

    σ .= h*(T0 - T) / corr
end

# Create one tilted heater patch and place it in the lower-left part of the irregular plate.

heater_body = Rectangle(heater_size[1],heater_size[2],1.4*Δx)
heater_transform = RigidTransform((0.38INCHES,-0.30INCHES),π/10)
heater = AreaForcingModel(heater_body,heater_transform,heatermodel!)

# Create one tilted condenser patch on the opposite side of the plate.

condenser_body = Rectangle(0.36INCHES,1.02INCHES,1.4*Δx)
condenser_transform = RigidTransform((1.95INCHES,0.05INCHES),π/10)
condenser = AreaForcingModel(condenser_body,condenser_transform,condensermodel!)

#   # Set Up The OHP Channel Shape

# First, build a fully custom curve from arrays. This is a quick way to sketch an unusual channel before choosing the final simulation path.

a = 0.018
θ = 0:2π/800:2π
r = a .* (1 .+ 0.28cos.(5θ)) .* sin.(2θ)
x = r .* cos.(θ)
y = 0.72r .* sin.(θ)

plot(x,y,aspectratio=1)

# Next, use the built-in multi-turn generator, but tilt the channel so it sits diagonally inside the strange polygonal plate.

ds = 1.5Δx
nturn = 7
width_ohp = 38e-3
length_ohp = 125e-3
gap = 3.2e-3
pitch = width_ohp/(2*nturn+1)
rotation_angle = -5π/9
x0, y0 = length_ohp/2 - 4e-3, width_ohp/2 - 5e-2

x,y, xf, yf = construct_ohp_curve(nturn,pitch,length_ohp,gap,ds,x0,y0,false,true,rotation_angle)

plot(x,y,aspectratio=1)

# Convert the generated coordinates into a BasicBody and shift the tilted channel slightly off center.

ohp = BasicBody(x,y)
tr_ohp = RigidTransform((-0.10INCHES,0.05INCHES),0.0)

# Define the line-source model that injects the tube-to-plate heat flux into the solid problem.

function ohpmodel!(σ,T,t,fr::LineRegionCache,phys_params)
    σ .= phys_params["ohp_flux"] ./ phys_params["flux_correction"]
end
ohp_linesource = LineForcingModel(ohp,tr_ohp,ohpmodel!)

#   ### Plot The Full Custom Layout

# Plot the irregular plate, heater, condenser, and tilted OHP channel together for visual inspection.

plot(body,fillalpha=0)
update_body!(heater_body,heater_transform)
update_body!(condenser_body,condenser_transform)
plot!(heater_body,fillcolor=:red,fillalpha=0.5)
plot!(condenser_body,fillcolor=:blue,fillalpha=0.5)
update_body!(ohp,tr_ohp)
plot!(ohp,fillalpha=0,closedornot=true)

# Assemble all thermal forcing models into the dictionary expected by the heat-conduction problem.

forcing_dict = Dict("heating models" => [heater,condenser,ohp_linesource])

#   # Construct The Systems

# Set the coupled time-step size and define the fixed time-step callback used by the plate solver.

tstep = 4e-4
timestep_fixed(u,sys) = tstep

# Construct the Neumann heat-conduction problem using the custom geometry, boundary conditions, and forcing models.

prob = NeumannHeatConductionProblem(g,body,phys_params=phys_params,
                                           bc=bcdict,
                                           motions=m,
                                           forcing=forcing_dict,
                                           timestep_func=timestep_fixed)

# Build the concrete two-dimensional plate system.

sys_plate = construct_system(prob)

# Initialize the one-dimensional OHP tube system using the custom plate system.

sys_tube = initialize_ohpsys(sys_plate,p_fluid,power)

#   # Initialize The Problem

# Create initial states and ODE integrators for the plate and tube systems.

tspan = (0.0,0.2)
dt_record = 0.01

tspan_init = (0.0,1e4)
u_plate = init_sol(sys_plate)
integrator_plate = init(u_plate,tspan_init,sys_plate,save_on=false)

u_tube = newstate(sys_tube)
integrator_tube = init(u_tube,tspan,sys_tube)

# Create the SimulationResult container that stores snapshots for plotting and post-processing.

SimuResult = SimulationResult(integrator_tube,integrator_plate)

#   # Solve

# Advance the coupled tube and plate systems in time and store snapshots at the requested interval.

@showprogress for t in tspan[1]:tstep:tspan[2]

    timemarching!(integrator_tube,integrator_plate,tstep)

    if (mod(integrator_plate.t,dt_record) < 1e-6) || (mod(-integrator_plate.t,dt_record) < 1e-6)
        store!(SimuResult,integrator_tube,integrator_plate)
    end

end

#   # Store Data

# Save the custom-geometry simulation result for later inspection.

save_path = "../numedata/DIY.jld2"
save(save_path,"SimulationResult",SimuResult)

#   ### Preview The Solution

# Generate an animation of plate temperature over time with the strange polygon boundary overlaid.

@gif for i in eachindex(SimuResult.tube_hist_t)
    plot(OHPTemp(),i,SimuResult,clim=(291.2,294.0))
    plot!(body,fillalpha=0)
end

# Generate an animation of the slug distribution over time with the tilted OHP channel and plate boundary overlaid.

@gif for i in eachindex(SimuResult.tube_hist_t)
    plot(OHPSlug(),i,SimuResult,aspectratio=1)
    plot!(body,fillalpha=0)
end
