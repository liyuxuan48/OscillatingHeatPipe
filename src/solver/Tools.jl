export getgvec,getgh,delta, # get actrual heightg of the tube
XptoLvaporplug,XptoLliquidslug,getXpvapor, # transfer Xp to the length of vapors, length of liquids, and Xp for vapor.
# ifamong,constructXarrays,
duliquidθtovec,duwallθtovec,liquidθtovec,wallθtovec, # transfer temperature field to state vector for liquid and wall.
Hfilm,getδarea,getδFromδarea,getMvapor,getMfilm,getMliquid,getVolumevapor,
getCa,filmδcorr,getAdeposit,f_churchill,Catoδ,RntoΔT,systoM

function getgvec(g0::T,g_angle::T=3/2*π) where {T<:Real}
    g = g0*[cos(g_angle),sin(g_angle)]
end

function getgh(g::Vector{T},x::Vector{T},y::Vector{T}) where {T<:Real}
    xy = [x';y']
    vec(sum(-g .* xy,dims=1));
end

delta(n) = n != 0


#=
Some notes on structures of liquid and vapor positions:
- The Xp is coupled by every liquid slug. For instance, if there is one liquid slug. Xp is a one-element tuple (Xp[1][1], Xp[1][2]).
    But sometimes we need Xp to be coupled by every vapor plug. For one liquid slug, we have two vapor plugs.
    So by adding 0 and L at the beginning and the end,
    we construct a two-element tuple ((0.0,Xp[1][1]) and ((Xp[1][2],L). Generally, for every N-element Xp, we construct an N+1 element Xpvapor
        Xp    ::   the locations of all interfaces, each element means a liquid slug.
        L     ::   the length of the 1D tube
- Xp contains the liquid slug begin/end positions (in arclength). It is arranged as a
  vector of tuples, where the first tuple element is the beginning and second element is
  the end
- If the tube is periodic, then `tube.closedornot = true` and there are as many vapor plugs as liquid slugs.
  Vapor plug i is bounded by liquid slugs i-1 and i. Plug 1 is bounded by N and 1.
- If the tube is not periodic, then there is one less vapor plug than liquid slugs. Vapor plug i is bounded
  by slugs i and i+1.
=#

"""
    Xtovec(Xp,dXdt) -> Vector

Transform vectors `Xp`, `dXdt` with the coordinates and velocities of the liquid slug interfaces into
a single state vector 
"""
function Xtovec(Xp::Array{Tuple{Float64,Float64},1},dXdt::Array{Tuple{Float64,Float64},1})
        
    Np = length(Xp)
    u = zeros(4*Np)
    
    for i = 1:Np
        # input Xp
        u[2*i-1] = Xp[i][1]
        u[2*i] = Xp[i][end]
        # input dXdt
        u[2*Np + 2*i-1] = dXdt[i][1]
        u[2*Np + 2*i] = dXdt[i][end]
    end
    
    return u
end


"""
    XMδLtovec(Xp,dXdt,M,δstart,δend,Lfilm_start,Lfilm_end) -> Vector

Assemble vectors `Xp`, `dXdt` of the coordinates and velocities of the slug interfaces, and
mass `M` of the vapor regions, `δstart` and `δend` of film thicknesses, and `Lfilm_start`
and `Lfilm_end` of film lengths, to form state vector.
"""
function XMδLtovec(Xp,dXdt,M,δstart,δend,Lfilm_start,Lfilm_end)
    return ([Xtovec(Xp,dXdt);M;δstart;δend;Lfilm_start;Lfilm_end])
end

"""
    vectoXMδL(u) -> Vector{Tuple}, Vector{Tuple}, Vector, Vector, Vector, Vector, Vector, Vector, Vector

Disassemble the state vector `u` into vectors `Xp`, `dXdt` of the coordinates and velocities of the slug interfaces, and
mass `M` of the vapor regions, `δstart` and `δend` of film thicknesses, and `Lfilm_start`
and `Lfilm_end` of film lengths.
"""
function vectoXMδL(u::Array{Float64,1})

    maxindex = div(length(u),9)
    maxliquidindex = 0
    if mod(length(u),9) == 4
        maxliquidindex = maxindex + 1
    elseif mod(length(u),9) == 0
        maxliquidindex = maxindex
    else
        println("new function error")
        return error("new function error")
    end

    Xp =   map(tuple, zeros(maxliquidindex), zeros(maxliquidindex))
    dXdt = map(tuple, zeros(maxliquidindex), zeros(maxliquidindex))
    for i = 1:maxliquidindex
        # input Xp
        Xp[i] = (u[2*i-1],u[2*i])
        # input dXdt
        dXdt[i] = (u[2*maxliquidindex + 2*i-1],u[2*maxliquidindex + 2*i])
    end

    M =             u[4*maxliquidindex+1:           4*maxliquidindex+maxindex]
    δstart =        u[4*maxliquidindex+maxindex+1:  4*maxliquidindex+2maxindex]
    δend =          u[4*maxliquidindex+2maxindex+1: 4*maxliquidindex+3maxindex]
    Lfilm_start =   u[4*maxliquidindex+3maxindex+1: 4*maxliquidindex+4maxindex]
    Lfilm_end =     u[4*maxliquidindex+4maxindex+1: 4*maxliquidindex+5maxindex]
    return Xp,dXdt,M,δstart,δend,Lfilm_start,Lfilm_end
end

"""
    XptoLvaporplug(Xp,L,closedornot) -> Vector{Float64}

This function uses the set of coordinates `Xp` of every liquid/vapor interface,
and length `L` of the tube to form an array of vapor lengths.
"""
function XptoLvaporplug(Xp::Vector{Tuple{Float64, Float64}},L::Float64,closedornot::Bool)

    if closedornot == false
        maxindex = length(Xp) - 1
        Lvaporplug = zeros(maxindex)

        if maxindex > 0
            for i = 1:maxindex
                Lvaporplug[i] = Xp[i+1][1] - Xp[i][end]
            end
        else
            return error("need more than one liquid slugs")
        end
    end

    if closedornot == true
        maxindex = length(Xp)
        Lvaporplug = zeros(maxindex)

        Lvaporplug[1] = mod((Xp[1][1]-Xp[end][end]),L)
        for i = 2:maxindex
            Lvaporplug[i] = mod((Xp[i][1] - Xp[i-1][end]),L)
        end
    end
        
    return Lvaporplug
end

"""
    XptoLliquidslug(Xp,L) -> Vector{Float64}

This function uses the set of coordinates `Xp` of every liquid/vapor interface,
and length `L` of the tube to form an array of slug lengths.
"""
function XptoLliquidslug(Xp::Vector{Tuple{Float64, Float64}},L::Float64)

    Lliquidslug = zeros(length(Xp))

    for i = 1:length(Xp)
        Lliquidslug[i] = mod((Xp[i][end] - Xp[i][1]),L)
    end

    return Lliquidslug
end



"""
    getXpvapor(Xp::Vector{Tuple},closedornot) -> Vector{Tuple}

Given the vector of liquid slug coordinates `Xp`, return the vector of vapor plug coordinates. 
"""
function getXpvapor(Xp::Vector{Tuple{Float64, Float64}},closedornot::Bool)

    # this is not periodic -- it has ends.  
    if closedornot == false 
        Xpvapor = map(tuple,zeros(length(Xp)-1),zeros(length(Xp)-1))
        maxindex = length(Xp) - 1
        if maxindex > 1
            for i = 1:maxindex
                Xpvapor[i] = (Xp[i][end], Xp[i+1][1])
            end
        else
            return error("need more than one liquid slugs")
        end
    end

    # this is a periodic tube, so there are as many vapor plugs and liquid slugs
    if closedornot == true  
        Xpvapor = map(tuple,zeros(length(Xp)),zeros(length(Xp)))
        Xpvapor[1]=(Xp[end][end],Xp[1][1])
        for i = 2:(length(Xp))
            Xpvapor[i]=(Xp[i-1][end],Xp[i][1])
        end
    end

    return Xpvapor
end


"""
    ifamongone(value::Real,range::Tuple)

Determine if `value` lies between the pair of values in the tuple `range`.
This accounts for cases in closed tubes in which `range` crosses the branch cut,
so that the start point in `range` has a larger value than the end point.
"""
function ifamongone(value::Float64, range::Tuple{Float64,Float64})
    # Note: this function does not make use of the tube length, so it might have errors.
    # It should also be possible to do it faster.
    return ((value >= range[1]) && (value <= range[end])) || ((value <= range[end]) && (range[1] >= range[end])) || ((value >= range[1]) && (range[1] >= range[end])) ? true : false
end

"""
    ifamong(value::Real,X::Vector{Tuple})

Determine if `value` lies between any of the tuples in the vector `X`.
"""
function ifamong(value::Float64, X::Vector{Tuple{Float64,Float64}})
    return Bool(sum(ifamongone.(value,X)) >= 1 ? true : false)
end

"""
    constructXarrays(X::Vector{Tuple},N::Int,θi::Float64,L::Float64) -> Vector{Vector}, Vector{Vector}

Initialize the 1D arc coordinate and temperature fields for every liquid slug. The
number of field points in each slug is never smaller than 2, and set equal to the fraction
of `N` based on the slug's length relative to `L`. 
The coordinate array is interpolated between the values in each `X` element (with 
the branch cut at `L` in a closed tube respected), and the temperature array
is set to `θi` uniformly. 
"""
function constructXarrays(X0::Vector{Tuple{Float64, Float64}},N::Int,θinitial::Float64,L::Float64)
    Xarrays=Array{Array{Float64, 1}, 1}(undef, length(X0))

    Lliquid = XptoLliquidslug(X0,L)
    Nliquid = max.(ceil.(Int, N.*Lliquid./L),2)

    for i = 1:length(Xarrays)
        if X0[i][1] < X0[i][2]
            Xarrays[i] = range(X0[i][1], X0[i][2], length=Nliquid[i])
        else
            Xarrays[i] = mod.(range(X0[i][1], X0[i][2]+L, length=Nliquid[i]), L)
            # Xarrays[i] = mod.(Xarrays[i], L)
        end
    end

    θarrays = deepcopy(Xarrays)
    for i in eachindex(θarrays)
        θarrays[i][:] .= θinitial
    end

    return(Xarrays,θarrays)
end

"""
    constructoneXarray(X::Tuple,N::Int,L::Float64) -> Vector

Initialize the 1D arc coordinate and temperature fields for one liquid slug. The
number of field points in each slug is set equal to `N`. 
The coordinate array is interpolated between the values in `X` (with 
the branch cut at `L` in a closed tube respected), and the temperature array
is set to `θi` uniformly. 
"""
function constructoneXarray(X0::Tuple{Float64,Float64},Nliquid,L)
    Xarray=Array{Float64, 1}

        if X0[1] < X0[2]
            Xarray = Array(range(X0[1], X0[2], length=Nliquid))
        else
            Xarray = Array(mod.(range(X0[1], X0[2]+L, length=Nliquid),L))
        end

        Xarray[1] = X0[1]
        Xarray[end] = X0[2]

    return Xarray
end


"""
    initialize X and θ field for wall, return Array{Float64, 1} and Array{Float64, 1}

    X0       :: Array{Tuple{Float64,Float64},1}
    N        :: Int, the number of cells in the wall (ΔX for liquid equals ΔX for the wall)
    θinitial :: value
    L        :: tube length
"""

"""
    constructwallXθarray(X::Vector,θi::Float64,curv::Vector) -> Vector, Vector, Vector

Return wall vectors for arc length coordinate, temperature, and curvature, given
the vector `X` for coordinates, a temperature value `θi`, and a vector of curvatures
`curv` 
"""
function constructwallXθarray(line::Vector{Float64},θinitial::Float64,curv::Vector{Float64})
    Xwallarray = deepcopy(line)
    θwallarray = Xwallarray .* 0 .+ θinitial
    curvwallarray = deepcopy(curv)

    return(Xwallarray,θwallarray,curvwallarray)
end



"""
    A bunch of functions to transfer θ to state vector rate du
"""

function duliquidθtovec(duθarrays::Vector{Vector{Float64}})
    return vcat(map(oneduliquidθtovec, duθarrays)...)
end

function oneduliquidθtovec(duθarray::Vector{Float64})
    return [0.0; duθarray]
end

"""
    A bunch of functions to transfer θ to state vector u
"""

function liquidθtovec(θarrays::Vector{Vector{Float64}})
    return vcat(map(wallθtovec, θarrays)...)
end

function wallθtovec(θwall::Vector{Float64})
    return [SEPERATION_VAR; θwall]
end

"""
    Hfilm(δfilm,sys::PHPSystem) -> Float64

Return the heat transfer coefficient for the film with thickness `δfilm`
"""
Hfilm(δfilm,sys::AbstractPHP) = Hfilm(δfilm,sys.vapor)

function Hfilm(δfilm,vapor::Vapor)
    @unpack δmin,k,Hᵥ  = vapor
    δthreshold = FILM_THRESHOLD
    δmax = FILM_MAX_THICKNESS
    kₗ = k

    if (δfilm > δthreshold) && (δfilm < δmax)
        return kₗ/δfilm
    elseif (δfilm > δmax) && (δfilm < 2δmax)
        return  kₗ/δmax - (δfilm-δmax)*(kₗ/δmax^2) + 1e-6
    elseif δfilm > δmin
        return  Hᵥ + (δfilm-δmin)*(kₗ/δthreshold - Hᵥ)/(δthreshold-δmin) + 1e-6
    else
        # return Hᵥ  + 1e-6
        return 0.0
    end
end



"""
    getδarea(Ac,d,δ)

Given cross-sectional tube area `Ac`, diameter `d`, and film thickness `δ`,
return the cross-sectional area of the film. 
"""
getδarea(Ac,d,δ) = Ac * (1 - ((d - 2*δ ) / d) ^ 2)

"""
    getδFromδarea(Ac,d,δarea)

Given cross-sectional tube area `Ac`, diameter `d`, and film cross-sectional
area `δarea`, return the film thickness. 
"""
getδFromδarea(Ac,d,δarea) = sqrt(δarea/Ac) * d/2
 
"""
    getMvapor(sys::PHPSystem) -> Vector

Return the masses of all of the vapor regions
"""
function getMvapor(sys)

    @unpack propconvert, vapor = sys
    @unpack PtoD = propconvert
    @unpack P = vapor

    ρᵥ = PtoD.(P)

    return ρᵥ .* getVolumevapor(sys)
    
end


"""
    getVolumevapor(sys::PHPSystem)

Return the volumes of the vapor regions, given the tube system `sys`.
"""
function getVolumevapor(sys)
    @unpack tube, vapor, liquid = sys
    @unpack Ac, L, d, closedornot = tube
    @unpack δstart, δend, Lfilm_start, Lfilm_end = vapor
    @unpack Xp = liquid
    
    Lvaporplug = XptoLvaporplug(Xp,L,closedornot)
    Astart = getδarea.(Ac,d,δstart)
    Aend = getδarea.(Ac,d,δend)
    
    return _getVolumevapor.(Ac,Astart,Aend,Lvaporplug,Lfilm_start,Lfilm_end)

end


_getVolumevapor(Ac,Astart,Aend,Lvaporplug,Lfilm_start,Lfilm_end) = 
                        Ac * Lvaporplug - Astart * Lfilm_start - Aend * Lfilm_end



"""
    getMfilm(sys::PHPSystem)

Return the masses of all of the liquid films in the system.
"""
function getMfilm(sys)
    @unpack tube, vapor, liquid = sys
    @unpack Ac, d = tube
    @unpack δstart, δend, Lfilm_start, Lfilm_end = vapor
    @unpack ρₗ = liquid

    Astart = getδarea.(Ac,d,δstart)
    Aend = getδarea.(Ac,d,δend)

    return _getM.(ρₗ,Astart,Lfilm_start), _getM.(ρₗ,Aend,Lfilm_end)
end

"""
    getMliquid(sys::PHPSystem)

Return the masses of all of the liquid slugs in the system.
"""
function getMliquid(sys)
    @unpack tube, liquid = sys
    @unpack Ac, L = tube
    @unpack ρₗ, Xp = liquid

    Lliquidslug = XptoLliquidslug(Xp,L)
    return  _getM.(ρₗ,Ac,Lliquidslug)

end

_getM(ρₗ,A,L) = ρₗ*A*L


"""
    getCa(μ,σ,velocity)

Return the capillary numbers for all liquid interfaces
"""
function getCa(μ,σ,velocity)
    Ca = abs.(μ.*velocity./σ)
end

"""
    filmδcorr(Ca,d)

Return the film thickness of deposited film, based on Aussillous and Quere (2000),
using the capillary number and diameter.
"""
function filmδcorr(Ca,d)
    filmδ = d .* 0.67.*Ca.^(2/3)./(1 .+ 3.35.*Ca.^(2/3))
end

"""
    getAdeposit(sys,δdeposit) -> Vector{Tuple}

Return the cross-sectional areas of deposited films, associated with each liquid slug's
interfaces. For an open tube, the first and last slugs are stuck at the ends of the tube,
so the deposited film areas are set to zero and ignored.
"""
function getAdeposit(sys,δdeposit)
    dXdt= sys.liquid.dXdt
    Ac= sys.tube.Ac
    d = sys.tube.d
    closedornot = sys.tube.closedornot
    δstart = sys.vapor.δstart
    δend = sys.vapor.δend

    Nliquid = length(dXdt)

    # numofliquidslug = length(dXdt)

    # Areas of deposited films
    δdepositArea = getδarea.(Ac,d,δdeposit)

    # Areas of existing films
    δarea_start = getδarea.(Ac,d,δstart)
    δarea_end = getδarea.(Ac,d,δend)

# need to initialize it later on
    loop_plus_index = circshift(1:Nliquid,-1)
    Adeposit = deepcopy(dXdt)

# left and right are relative for liquid. Left = behind, right = ahead
# if advancing, then use deposited film area. otherwise, existing film area

if closedornot == true
    for i in eachindex(Adeposit)

        # film behind slug. use "end" film from vapor region i if slug not depositing.
        Adeposit_left = dXdt[i][1] > 0 ? δdepositArea : δarea_end[i]
        # film in front of slug. use "start" film from vapor region i+1 if slug not depositing.
        Adeposit_right = dXdt[i][end] < 0 ? δdepositArea : δarea_start[loop_plus_index[i]]

        Adeposit[i]  =   (Adeposit_left, Adeposit_right)
    end
else
    Adeposit[1] =   (0.0,0.0)
    Adeposit[end] = (0.0,0.0)
    for i in eachindex(Adeposit[2:end-1])
        # film behind slug. use "end" film from vapor region i if slug not depositing.
        Adeposit_left = dXdt[i][1] > 0 ? δdepositArea : δarea_end[i]
        # film in front of slug. use "start" film from vapor region i+1 if slug not depositing.
        Adeposit_right = dXdt[i][end] < 0 ? δdepositArea : δarea_start[i+1]

        Adeposit[i]  =   (Adeposit_left, Adeposit_right)
    end
end

    Adeposit
end

"""
    f_churchill(Re,ϵ)

Return the Darcy friction factor, based on the Churchill correlation.
"""
function f_churchill(Re,ϵ=0.001)
    Θ1 = (-2.457*log((7/Re)^(0.9)  +  0.27 * ϵ))^16
    Θ2 = (37530/Re)^16
    f=8*((8/Re)^12+(1/(Θ1+Θ2)^1.5))^(1/12)
    
    f
end

"""
    Catoδ(d,Ca[;adjust_factor=1,δmin=2e-6,δmax=1e-4])

Given capillary number `Ca` and `d`, return the thickness.
This is a generalization of the Aussillous and Quere formula.
It uses an adjustment factor that can be tuned empirically.
"""
function Catoδ(d,Ca;adjust_factor=1,δmin=2e-6,δmax=1e-4)

    δ = Ca .^ (2/3) ./ (1 .+ Ca .^ (2/3)) .* d ./ 2 .* adjust_factor
    if (δ < δmin)
        return δmin
    elseif (δ > δmax)
        return δmax
    else 
        return δ
    end
end


"""
    RntoΔT(Rn,Tref,fluid_type,d,TtoP) -> Vector

Given the nucleation radius `Rn`, the fluid type `fluid_type` and reference
temperature `Tref`, the channel diameter `d`, and the temperature-to-pressure
relation `TtoP`, return the superheat threshold temperature for boiling.
"""
function RntoΔT(Rn,Tref,fluid_type,d,TtoP)
    p_fluid = SaturationFluidProperty(fluid_type,Tref);

    Rkg = p_fluid.R/p_fluid.M
    Rin = d/2
    P = TtoP(Tref)

    y = Rkg .* Tref ./ (p_fluid.hᵥ-p_fluid.hₗ) .* log.(1 .+ 2 .* p_fluid.σ ./ P .* (1 ./ Rn .- 1/(2Rin)))
    ΔTref = Tref .* (1 ./ (1 .- y) .- 1)
end

"""
    systoM(sys::PHPSystem)

Given tube system `sys`, return the masses of the vapor regions.
"""
function systoM(sys0::PHPSystem)

    @unpack Ac,d,L,closedornot = sys0.tube
    @unpack Xp = sys0.liquid
    @unpack δstart,δend,Lfilm_start,Lfilm_end,P = sys0.vapor
    @unpack PtoD = sys0.propconvert

    Lvaporplug = XptoLvaporplug(Xp,L,closedornot)

    δarea_start = getδarea.(Ac,d,δstart)
    δarea_end   = getδarea.(Ac,d,δend)

    volume_vapor = Lvaporplug .* Ac .- Lfilm_start .* δarea_start .- Lfilm_end .* δarea_end

    M = PtoD.(P) .* volume_vapor
end