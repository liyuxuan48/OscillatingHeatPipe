export PHPSystem_nomapping,PHPSystem,Tube,Liquid,Vapor,Wall,PropConvert,Mapping,Cache
# ,PHPResult

using Parameters

"""
Tube is a struct containing tube geometries
    d           ::Float64   tube characteristic diameter
    peri        ::Float64   tube perimeter (may be an arbitrary shape so cannot be derived from d)
    Ac          ::Float64   tube cross-sectional area
    L           ::Float64   tube one dimensional length
    g           ::Float64   gravity vector along the plate
    closedornot ::Bool      if the tube is closed loop or not (open loop)
"""

struct Tube
    d::Float64
    peri::Float64
    Ac::Float64
    L::Float64
    g::Vector{Float64}
    closedornot::Bool
    N::Int64  
end

# function Tube(d,peri,Ac,L,gravity,closedornot,N);
#     # PtoT,TtoP,PtoD,DtoP,PtoHfg = createCoolPropinterpolation(fluid_type::String)
#     # Tube(d,peri,Ac,L,L2D,angle,gravity,closedornot,N,fluid_type,PtoT,TtoP,PtoD,DtoP,PtoHfg);
#     # Tube(d,peri,Ac,L,gravity,closedornot,N,fluid_type,PtoT,TtoP,PtoD,DtoP,PtoHfg);
#     Tube(d,peri,Ac,L,gravity,closedornot,N);
# end

"""
Liquid is a struct containing liquid properties at a ref temperature
    Hₗ::Float64                              heat transfer coefficient between the wall and the pure liquid slug
    ρ::Float64                              liquid density
    Cp::Float64                             liquid specific heat capacity
    α::Float64                              liquid heat diffusivity
    μₗ::Float64                              liquid dynamic viscosity
    σ::Float64                              liquid surface tension
    Xp::Array{Tuple{Float64,Float64},1}     interface locations for each liquid slug
    dXdt::Array{Tuple{Float64,Float64},1}   interface velocity for each liquid slug
    Xarrays::Array{Array{Float64,1},1}      finite difference location points within each liquid slug
    θarrays::Array{Array{Float64,1},1}      finite difference temperature points within each liquid slug
"""

struct PropConvert
    fluid_type::String
    PtoT    :: AbstractInterpolation
    TtoP    :: AbstractInterpolation
    PtoD    :: AbstractInterpolation
    DtoP    :: AbstractInterpolation
    PtoHfg  :: AbstractInterpolation
end

function PropConvert(fluid_type::String);
    PtoT,TtoP,PtoD,DtoP,PtoHfg = createCoolPropinterpolation(fluid_type::String)
    PropConvert(fluid_type,PtoT,TtoP,PtoD,DtoP,PtoHfg);
end

mutable struct Liquid
    Hₗ::Float64
    ρₗ::Float64
    Cpₗ::Float64
    αₗ::Float64
    μₗ::Float64
    σ::Float64
    Xp::Array{Tuple{Float64,Float64},1}
    dXdt::Array{Tuple{Float64,Float64},1}
    Xarrays::Array{Array{Float64,1},1}
    θarrays::Array{Array{Float64,1},1}

end

"""
Vapor is a struct containing vapor properties at a ref temperature
    Hᵥ::Float64             heat transfer coefficient between the wall and the pure vapor bubble
    k::Float64              heat conductivity of liquid in the film
    δmin::Float64           the delta with maximum heat transfer coefficient in the H interpolation
    P::Array{Float64,1}     pressure in each vapor
    δ::Array{Float64,1}     film thickness in each vapor
end
"""

@with_kw mutable struct Vapor
    ad_fac::Float64                 = 1.3
    Hᵥ::Float64                     = 0.0
    k::Float64
    δmin::Float64                   = 2e-6
    Eratio_plus::Float64            = 0.15
    Eratio_minus::Float64           = 0.0
    P::Array{Float64,1}
    δstart::Array{Float64,1}
    δend::Array{Float64,1}
    Lfilm_start::Array{Float64,1}
    Lfilm_end::Array{Float64,1}
end

"""
Wall is a struct containing wall properties
    ΔTthres::Float64                superheat threshold to trigger boiling
    Xstations::Array{Float64,1}     locations of boiling stations on the wall
    Xarray::Array{Float64,1}        wall discrete location points for immersed boundary method
    θarray::Array{Float64,1}        wall discrete temperature points for immersed boundary method
end
"""

@with_kw mutable struct Wall
    fluid_type::String
    boil_type::String
    power::Float64
    boil_interval::Float64     = 1.0
    Rn::Float64                = 3e-6
    L_newbubble::Float64       = 4e-3
    Xstations::Array{Float64,1}
    boiltime_stations::Array{Float64,1}
    Xarray::Array{Float64,1}
    θarray::Array{Float64,1}
end

"""
Mapping is a struct containing interpolation data
    θ_interp_walltoliquid   temperature interpolation from wall to OHP
    θ_interp_liquidtowall   temperature interpolation from OHP to wall
    H_interp_liquidtowall   heat transfer coefficient interpolation from OHP to wall
    P_interp_liquidtowall   pressure interpolation from OHP to wall
"""

# mutable struct Mapping
#     walltoliquid::Array{Tuple{Int64,Int64},1}
#     liquidtowall::Array{Array{Int64,1},1}
# end

mutable struct Mapping
    θ_interp_walltoliquid :: AbstractInterpolation
    θ_interp_liquidtowall :: AbstractInterpolation
    H_interp_liquidtowall :: AbstractInterpolation
    P_interp_liquidtowall :: AbstractInterpolation
    heightg_interp        :: AbstractInterpolation
end

"""
Cache is a struct containing
    boil_hist    ::Vector{Any}
"""

mutable struct Cache
    boil_hist    ::Vector{Any}
    mass ::Float64
end


"""
PHPSystem is a struct containing
    tube    ::Tube
    liquid  ::Liquid
    vapor   ::Vapor
    wall    ::Wall
    mapping ::Mapping
"""

mutable struct PHPSystem
    tube    ::Tube
    liquid  ::Liquid
    vapor   ::Vapor
    wall    ::Wall
    propconvert :: PropConvert
    mapping ::Mapping
    cache   ::Cache
end

"""
PHPSystem_nomapping is a struct containing
    Tube
    Liquid
    Vapor
    Wall
It is used to construct PHPSystem, no other use.
"""
mutable struct PHPSystem_nomapping
    tube    ::Tube
    liquid  ::Liquid
    vapor   ::Vapor
    wall    ::Wall
    propconvert :: PropConvert
end

function Base.show(io::IO, sys::PHPSystem)
    N = sys.tube.N
    fluidtype = sys.propconvert.fluid_type
    println(io, "$N point OHP system filled with $fluidtype")
end
