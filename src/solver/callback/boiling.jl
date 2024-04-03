export boiling_affect!,nucleateboiling,boiling_condition
# boiling_condition,
function boiling_condition(u,t,integrator)
    t_interval = 1e-2

    ϵ = 1e-5

    return (abs(mod(t,t_interval)-t_interval) < ϵ) || mod(t,t_interval) < ϵ
end

function boiling_affect!(integrator)

    # Δθthreshold = integrator.p.wall.ΔTthres

    p = deepcopy(getcurrentsys!(integrator.u,integrator.p))

    boil_type = p.wall.boil_type
    fluid_type = p.propconvert.fluid_type
    d = p.tube.d
    Rn = p.wall.Rn
    boil_interval = p.wall.boil_interval

    @unpack PtoT,TtoP,PtoD = p.propconvert
    Tref = (PtoT.(maximum(p.vapor.P)) + PtoT.(minimum(p.vapor.P)))/2

    Δθthreshold = RntoΔT(Rn,Tref,fluid_type,d,TtoP)

    # println(Δθthreshold)
  
    Δθ_array = getsuperheat.(p.wall.Xstations,[p])
    superheat_flag = (Δθ_array .> Δθthreshold) .* ((integrator.t .- p.wall.boiltime_stations) .> boil_interval)

    # println(p.wall.boiltime_stations)
    b_count = 0;
    boiltime_update_flags = Bool.(zero(p.wall.boiltime_stations))
    for i = 1:length(p.wall.Xstations)
        if ifamong(p.wall.Xstations[i], p.liquid.Xp) && suitable_for_boiling(p,i) && superheat_flag[i]
                push!(integrator.p.cache.boil_hist,[i,integrator.t]);
                b_count += 1;

                if boil_type == "liquid T"
                    Pinsert = p.mapping.P_interp_liquidtowall(p.wall.Xstations[i])
                elseif boil_type == "wall T"
                    Pinsert = TtoP(p.mapping.θ_interp_walltoliquid(p.wall.Xstations[i]))
                end
                # println("boiling",get_vapor_energy(p)+get_liquid_energy(p))
                p = nucleateboiling(p,(p.wall.Xstations[i]-p.wall.L_newbubble/2,p.wall.Xstations[i]+p.wall.L_newbubble/2),Pinsert) # P need to be given from energy equation
                # println("boiling",get_vapor_energy(p)+get_liquid_energy(p))
                boiltime_update_flags[i] = true
                # p.wall.boiltime_stations[i] = integrator.t
            elseif !ifamong(p.wall.Xstations[i], p.liquid.Xp)
                boiltime_update_flags[i] = true
                # p.wall.boiltime_stations[i] = integrator.t
        end
    end

    boiltime_stations = p.wall.boiltime_stations + boiltime_update_flags .* (integrator.t .- p.wall.boiltime_stations)
    integrator.p.wall.boiltime_stations = boiltime_stations
    

    Lvaporplug = XptoLvaporplug(p.liquid.Xp,p.tube.L,p.tube.closedornot)
    Ac = p.tube.Ac
    d = p.tube.d

    δstart = p.vapor.δstart
    δend = p.vapor.δend

    Lfilm_start = p.vapor.Lfilm_start
    Lfilm_end = p.vapor.Lfilm_end

    δarea_start = Ac .* (1 .- ((d .- 2*δstart) ./ d) .^ 2);
    δarea_end = Ac .* (1 .- ((d .- 2*δend) ./ d) .^ 2);

    volume_vapor = Lvaporplug .* Ac - Lfilm_start .* δarea_start - Lfilm_end .* δarea_end
    M = PtoD.(p.vapor.P) .* volume_vapor

    # println(p.wall.boiltime_stations)


    unew=[XMδLtovec(p.liquid.Xp,p.liquid.dXdt,M,p.vapor.δstart,δend,Lfilm_start,Lfilm_end);liquidθtovec(p.liquid.θarrays)];

    resize!(integrator.u,length(unew))
    integrator.u = deepcopy(unew)
end

function nucleateboiling(sys,Xvapornew,Pinsert)
    ρₗ = deepcopy(sys.liquid.ρₗ)
    Ac = sys.tube.Ac
    d = deepcopy(sys.tube.d)
    Xp = deepcopy(sys.liquid.Xp)
    dXdt = deepcopy(sys.liquid.dXdt)
    δstart = deepcopy(sys.vapor.δstart)
    δend = deepcopy(sys.vapor.δend)
    
    P = deepcopy(sys.vapor.P)
    L = sys.tube.L
    
    Xarrays = sys.liquid.Xarrays
    θarrays = sys.liquid.θarrays
    closedornot = sys.tube.closedornot
    Lfilm_start = deepcopy(sys.vapor.Lfilm_start)
    Lfilm_end = deepcopy(sys.vapor.Lfilm_end)


    index = getinsertindex(Xp,(Xvapornew[2]+Xvapornew[1])/2,sys.tube.L,sys.tube.closedornot)

    Linsert = mod(Xvapornew[end] - Xvapornew[1],L)

    """let's do constant film thickness for now!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"""
    # δdeposit = maximum([5e-5;δstart;δend]) + 1e-6
    δdeposit = 2e-5
    δstart_new = insert!(δstart,index+1,δdeposit)
    δend_new = insert!(δend,index+1,δdeposit)


    Nvapor = length(P)
    loop_plus_index = [2:Nvapor;1]
    loop_plus_index_new = [3:Nvapor+1;1:2]

    Lfilm_start_new = insert!(Lfilm_start,index+1,Linsert/8)
    Lfilm_end_new = insert!(Lfilm_end,index+1,Linsert/8)

    Lliquid_adjust = 0
    Xpnew = getnewXp(Xp,index,Xvapornew,Lliquid_adjust,L,closedornot)

    # const P in adjacent vapors
    Pnew = insert!(P,index+1,Pinsert)

    dXdtnew = deepcopy(dXdt) # momentum conservation
    insert!(dXdtnew,index+1,dXdtnew[index])

    sysnew = deepcopy(sys)

    sysnew.liquid.Xp = Xpnew
    sysnew.liquid.dXdt = dXdtnew
    sysnew.vapor.P = Pnew
    sysnew.vapor.δstart = δstart_new
    sysnew.vapor.δend = δend_new
    sysnew.vapor.Lfilm_start = Lfilm_start_new
    sysnew.vapor.Lfilm_end = Lfilm_end_new

    # Mvapor_old = sum(getMvapor(sys))
    # Mfilm_old = sum(sum.(getMfilm(sys)))
    Mvapor_new = sum(getMvapor(sysnew))
    Mfilm_new = sum(sum.(getMfilm(sysnew)))
    Mliquid_new = sum(getMliquid(sysnew))
    Mold = sysnew.cache.mass
    L_adjust = (Mvapor_new + Mfilm_new + Mliquid_new - Mold)./ (ρₗ*Ac)

    Lvaporplug = XptoLvaporplug(sysnew.liquid.Xp,sysnew.tube.L,sysnew.tube.closedornot)
    Lpurevapor = Lvaporplug .- Lfilm_start_new .- Lfilm_end_new
    Lliquidslug = XptoLliquidslug(Xpnew,sys.tube.L)



    maxindex = 0;

    if L_adjust < 0
        L_newbubble = sysnew.wall.L_newbubble
        L_adjust = (-L_adjust > L_newbubble) ? -L_newbubble : L_adjust

        maxvalueindex = findmax(Lpurevapor)
        maxvalue = maxvalueindex[1]
        maxindex = maxvalueindex[2]

        if maxvalue > L_adjust
            sysnew.liquid.Xp[maxindex] = mod.((sysnew.liquid.Xp[maxindex][1]+L_adjust,sysnew.liquid.Xp[maxindex][2]),L)
            # println(L_adjust, "-")
        else 
            maxindex = 0
            println("boiling error!")
        end
    else
        # L_newbubble = sysnew.wall.L_newbubble
        # L_adjust = (L_adjust > L_newbubble) ? L_newbubble : L_adjust

        maxvalueindex = findmax(Lliquidslug)
        maxvalue = maxvalueindex[1]
        maxindex = maxvalueindex[2]

        if maxvalue > L_adjust
            sysnew.liquid.Xp[maxindex] = mod.((sysnew.liquid.Xp[maxindex][1]+L_adjust,sysnew.liquid.Xp[maxindex][2]),L)
        else 
            maxindex = 0
            println("boiling error!")
        end
    end

    # up to now we got the correct Xpnew, next step is to get Xarraysnew, the splitted Xarrays.

    Xarraysnew,θarraysnew = getnewXθarrays(index,sysnew.liquid.Xp,Xarrays,θarrays,L,maxindex)
    # θarraysnew = getnewθarrays(index,Xp,sysnew.liquid.Xp,Xarrays,θarrays,L,closedornot)

    sysnew.liquid.Xarrays = Xarraysnew
    sysnew.liquid.θarrays = θarraysnew

    θ_interp_walltoliquid, θ_interp_liquidtowall, H_interp_liquidtowall, P_interp_liquidtowall = sys_interpolation(sysnew)
    heightg_interp = sysnew.mapping.heightg_interp
    sysnew.mapping = Mapping(θ_interp_walltoliquid, θ_interp_liquidtowall, H_interp_liquidtowall, P_interp_liquidtowall,heightg_interp)
    # sysnew.mapping = Mapping(θ_interp_walltoliquid, θ_interp_liquidtowall, H_interp_liquidtowall, P_interp_liquidtowall)

    # println(Pnew[index-2:index+5])

return sysnew
end

function getnewXθarrays(index,Xpnew,Xarrays_old,θarrays_old,L,maxindex)
    Nold= length(Xarrays_old[index])

    Xarraysnew = deepcopy(Xarrays_old)
    θarraysnew = deepcopy(θarrays_old)

    # Lliquidslug_old =       mod((Xarrays_old[index][end] - Xarrays_old[index][1]),L)
    Lliquidslug_new_left =  mod((Xpnew[index][end] - Xpnew[index][1]),L) 
    Lliquidslug_new_right = mod((Xpnew[index+1][end] - Xpnew[index+1][1]),L) 

    arrayindex = ceil.(Int, Nold * Lliquidslug_new_left / (Lliquidslug_new_left + Lliquidslug_new_right))
    arrayindex = (arrayindex != 1) ? arrayindex : 2
    arrayindex = (arrayindex != Nold) ? arrayindex : Nold-1

    Xarraysnewleft = constructoneXarray(Xpnew[index],arrayindex,L)
    Xarraysnewright = constructoneXarray(Xpnew[index+1],Nold-arrayindex+1,L)

    θarraysnewleft = θarrays_old[index][1:arrayindex]
    θarraysnewright= θarrays_old[index][arrayindex+1:end]
    insert!(θarraysnewright, 1, θarrays_old[index][arrayindex])

    splice!(Xarraysnew, index)
    insert!(Xarraysnew, index,Xarraysnewleft)
    insert!(Xarraysnew, index+1,Xarraysnewright)

    splice!(θarraysnew, index)
    insert!(θarraysnew, index,θarraysnewleft)
    insert!(θarraysnew, index+1,θarraysnewright)

    if maxindex != 0
        Xarraysnew[maxindex] = constructoneXarray(Xpnew[maxindex],length(Xarraysnew[maxindex]),L)
    end

    Xarraysnew,θarraysnew
end


function getinsertindex(Xp,Xvapornew_center,L,closedornot)
if closedornot
    for index = 1:length(Xp)
        if Xp[index][1] <= Xvapornew_center && Xp[index][2] >= Xvapornew_center && Xp[index][2] >= Xp[index][1]
            return index
        end

        if (mod(Xvapornew_center - Xp[index][1],L) < mod(Xp[index][2] - Xp[index][1],L)) && Xp[index][2] <= Xp[index][1]
            return index
        end
    end
end


if !closedornot
    for index = 1:length(Xp)
        if Xp[index][1] <= Xvapornew_center && Xp[index][2] >= Xvapornew_center
            return index
        end
    end
end
# println((Xvapornew[2]+Xvapornew[1])/2)
# println(Xp)
        return NaN
end

# function getarrayindex(X,Xarray)

# for arrayindex = 1:length(Xarray)
#     if (X >= Xarray[arrayindex] && X <= Xarray[arrayindex+1]) || (X >= Xarray[arrayindex] && Xarray[arrayindex] >= Xarray[arrayindex+1]) || (X <= Xarray[arrayindex+1] && Xarray[arrayindex] >= Xarray[arrayindex+1])
#         return (arrayindex > 1.1) ? arrayindex : 2
# end
#     end
#         return NaN
# end

function getnewXp(Xp,index,Xvapornew,Lliquid_adjust,L,closedornot)

    Xpnew = deepcopy(Xp)

    Linsert = mod(Xvapornew[end] - Xvapornew[1],L)

    insertXp1=mod.((Xp[index][1]+Lliquid_adjust/2,Xvapornew[1]),L)
    insertXp2=mod.((Xvapornew[2],Xp[index][2]-Lliquid_adjust/2),L)

    splice!(Xpnew, index)
    insert!(Xpnew, index,insertXp1)
    insert!(Xpnew, index+1,insertXp2)

    return Xpnew
end

# # simplified non mass conservation
# function getnewXp(Xp,index,Xvapornew,L)

#     Xpnew = deepcopy(Xp)

#     insertXp1=mod.((Xp[index][1],Xvapornew[1]),L)
#     insertXp2=mod.((Xvapornew[2],Xp[index][2]),L)

#     splice!(Xpnew, index)
#     insert!(Xpnew, index,insertXp1)
#     insert!(Xpnew, index+1,insertXp2)

#     return Xpnew
# end


function getnewM(M,index,Minsert,closedornot)

    Mnew = deepcopy(M)

    insert!(Mnew,index+1,Minsert)

    return Mnew
end


function getsuperheat(Xstation,sys)
    @unpack PtoT = sys.propconvert
    Δθ = sys.mapping.θ_interp_walltoliquid(Xstation) - PtoT(sys.mapping.P_interp_liquidtowall(Xstation))
    return Δθ
end


# ```
#     get the index of element of Xarray closest to X.
#     closed loop considered
# ```
# function getoneXarrayindex(X,Xarray)
#     for i = 1:length(Xarray)
#         if (Xarray[i] >= Xarray[i+1]) && ((Xarray[i] <= X) || (Xarray[i+1] >= X))
#             return i
#         end
#         if (X >= Xarray[i] && X <= Xarray[i+1])
#             return i
#         end
#     end

#     return length(Xarray)
# end

function suitable_for_boiling(p,i)
    suitable_flag =  true
    L_newbubble = p.wall.L_newbubble

        index = getinsertindex(p.liquid.Xp,p.wall.Xstations[i],p.tube.L,p.tube.closedornot)

        L_liquid_left =  mod(p.wall.Xstations[i] - p.liquid.Xp[index][1],p.tube.L)
        L_liquid_right = mod(p.liquid.Xp[index][2] - p.wall.Xstations[i],p.tube.L)

        if (2*L_newbubble> L_liquid_left) || (2*L_newbubble > L_liquid_right)
            suitable_flag = false
        end

    return suitable_flag
end

function get_vapor_energy(sys0::PHPSystem)
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

    PtoD = sys0.propconvert.PtoD
    PtoT = sys0.propconvert.PtoT

    M = PtoD.(P) .* volume_vapor
    
    U =  PtoT.(P) .* sys0.liquid.Cpₗ .* volume_vapor
    sum(U)
end

function get_liquid_energy(sys0::PHPSystem)
    Ac = sys0.tube.Ac
    ρ = sys0.liquid.ρₗ
    Cp = sys0.liquid.Cpₗ
    L = sys0.tube.L
    sys0.liquid.Xarrays;
    ds_liquid = [mod(Xarray[end]-Xarray[1],L)/length(Xarray) for Xarray in sys0.liquid.Xarrays];
    # println(sys0.liquid.Xp)
    # println(ds_liquid)
    liquid_energy = sum.(sys0.liquid.θarrays).*ds_liquid .* (Ac .* ρ .* Cp)
    # println(liquid_energy)
    sum(liquid_energy)
end