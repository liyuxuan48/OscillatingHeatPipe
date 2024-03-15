export merging_affect!,merging_condition,nucleateboiling,merging

function merging_affect!(integrator)

    p = getcurrentsys(integrator.u,integrator.p);
    δv = 0.5*p.wall.L_newbubble

    merge_flags = getmerge_flags(δv,p)
    indexmergingsite = sort(findall(x->x == true, merge_flags),rev = true)

    for i in indexmergingsite
        p = merging(p,i)
    end

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
    @unpack PtoD = p.propconvert
    M = PtoD.(p.vapor.P) .* volume_vapor

    unew=[XMδLtovec(p.liquid.Xp,p.liquid.dXdt,M,δstart,δend,Lfilm_start,Lfilm_end); liquidθtovec(p.liquid.θarrays)];

    resize!(integrator.u,size(unew,1)::Int)
    integrator.u = deepcopy(unew)
end

function merging_condition(u,t,integrator)     # only for closed loop tube

    p = getcurrentsys!(integrator.u,integrator.p);
    # δv = p.tube.d > (integrator.dt*maximum(p.liquid.dXdt)[1]) ? p.tube.d : (integrator.dt*maximum(p.liquid.dXdt)[1])
    δv = 0.5*p.wall.L_newbubble

    sys = p

    merge_flags = getmerge_flags(δv,sys)

    return sum(merge_flags) != 0
    # return true
end

function merging(p,i)

    closedornot = p.tube.closedornot

    # get the liquid interface velocities and lengthes for merging
    Lliquidslug = XptoLliquidslug(p.liquid.Xp,p.tube.L)
    Lvaporplug  = XptoLvaporplug(p.liquid.Xp,p.tube.L,p.tube.closedornot)

    numofliquidslug = length(Lliquidslug)
    numofvaporbubble = length(Lvaporplug)

     
    L = p.tube.L
    Ac = p.tube.Ac
    ρₗ = p.liquid.ρ

    Mvapor = getMvapor(p)
    Mfilm = getMfilm(p);

    @unpack PtoD = p.propconvert

    # get compensated L of merged liquid slug for mass conservation
    if closedornot == true
        left_index = i > 1 ? i-1 : length(Lvaporplug)
        right_index = i < length(Lvaporplug) ? i+1 : 1
        liquid_right_i = i #  operating liquid i
        liquid_left_i = left_index

        splice_liquid_i_1 = i > 1 ? left_index : numofliquidslug
        splice_liquid_i_2 = i > 1 ? left_index : 1
        insert_i_1 = i > 1 ? left_index : numofliquidslug-1

    else
        left_index = i-1
        right_index = i+1
        liquid_right_i = i+1 #  operating liquid i
        liquid_left_i = i

        splice_liquid_i_1 = i
        splice_liquid_i_2 = i
        insert_i_1 = i

   
    end
 

  
    if closedornot == false && i == 1
        Linsert = (Mvapor[i] + Mfilm[1][i] + Mfilm[2][i] - 0.5 .* Ac .* Lvaporplug[i] .* (PtoD(p.vapor.P[right_index]))) ./ (ρₗ .* Ac .- 0.5 .* Ac .* (PtoD(p.vapor.P[right_index])))
        Xpnewone = (p.liquid.Xp[liquid_left_i][1], p.liquid.Xp[liquid_right_i][end] - Lvaporplug[i] + Linsert)
        dXdtnewonevalue = 0.0
    elseif i == numofvaporbubble
        Linsert = (Mvapor[i] + Mfilm[1][i] + Mfilm[2][i] - 0.5 .* Ac .* Lvaporplug[i] .* (PtoD(p.vapor.P[left_index]))) ./ (ρₗ .* Ac .- 0.5 .* Ac .* (PtoD(p.vapor.P[left_index])))
        Xpnewone = (p.liquid.Xp[liquid_left_i][1]+Lvaporplug[i] - Linsert, p.liquid.Xp[liquid_right_i][end])
        dXdtnewonevalue = 0.0
    else
        Linsert = (Mvapor[i] + Mfilm[1][i] + Mfilm[2][i] - 0.5 .* Ac .* Lvaporplug[i] .* (PtoD(p.vapor.P[left_index]) .+ PtoD(p.vapor.P[right_index]))) ./ (ρₗ .* Ac .- 0.5 .* Ac .* (PtoD(p.vapor.P[left_index]) .+ PtoD(p.vapor.P[right_index])))
        Xpnewone = (mod(p.liquid.Xp[liquid_left_i][1]+Lvaporplug[i]/2 - Linsert/2,L), mod(p.liquid.Xp[liquid_right_i][end] - Lvaporplug[i]/2 + Linsert/2,L))
        dXdtnewonevalue = (p.liquid.dXdt[liquid_left_i][1]*Lliquidslug[liquid_left_i] + p.liquid.dXdt[liquid_right_i][end]*Lliquidslug[liquid_right_i])/(Lliquidslug[liquid_left_i]+Lliquidslug[liquid_right_i]) 
    end

    # dXdtnewonevalue = (i != 1) ? (p.liquid.dXdt[i-1][1]*Lliquidslug[i-1] + p.liquid.dXdt[i][end]*Lliquidslug[i])/(Lliquidslug[i-1]+Lliquidslug[i]) : (p.liquid.dXdt[end][1]*Lliquidslug[end] + p.liquid.dXdt[i][end]*Lliquidslug[i])/(Lliquidslug[end]+Lliquidslug[i])
    # dXdtnewonevalue = (p.liquid.dXdt[liquid_left_i][1]*Lliquidslug[liquid_left_i] + p.liquid.dXdt[liquid_right_i][end]*Lliquidslug[liquid_right_i])/(Lliquidslug[liquid_left_i]+Lliquidslug[liquid_right_i]) 

    systemp = deepcopy(p)

    # up to here

    splice!(systemp.liquid.Xp,splice_liquid_i_1)
    splice!(systemp.liquid.Xp,splice_liquid_i_2)
    insert!(systemp.liquid.Xp,insert_i_1,Xpnewone)

    splice!(systemp.liquid.dXdt,splice_liquid_i_1)
    splice!(systemp.liquid.dXdt,splice_liquid_i_2)
    insert!(systemp.liquid.dXdt,insert_i_1,(dXdtnewonevalue,dXdtnewonevalue))

    # if i != 1
    #     splice!(systemp.liquid.Xp,i-1:i,[Xpnewone])
    # else
    #     splice!(systemp.liquid.Xp,length(systemp.liquid.Xp));
    #     splice!(systemp.liquid.Xp,1);
    #     insert!(systemp.liquid.Xp,length(systemp.liquid.Xp)+1,Xpnewone)
    # end

    # if i != 1
    #     splice!(systemp.liquid.dXdt,i-1:i,[(dXdtnewonevalue,dXdtnewonevalue)])
    # else
    #     splice!(systemp.liquid.dXdt,length(systemp.liquid.dXdt));
    #     splice!(systemp.liquid.dXdt,1);
    #     insert!(systemp.liquid.dXdt,length(systemp.liquid.dXdt)+1,(dXdtnewonevalue,dXdtnewonevalue))
    # end

    splice!(systemp.vapor.δstart,i)
    splice!(systemp.vapor.δend,i)
    splice!(systemp.vapor.Lfilm_start,i)
    splice!(systemp.vapor.Lfilm_end,i)
    splice!(systemp.vapor.P,i)

    # Nliquids = (i != 1) ? length([systemp.liquid.Xarrays[i-1]; systemp.liquid.Xarrays[i]]) : length([systemp.liquid.Xarrays[end]; systemp.liquid.Xarrays[i]])
    Nliquids = length([systemp.liquid.Xarrays[liquid_left_i]; systemp.liquid.Xarrays[liquid_right_i]]) 

    # Xarraysnewone = constructoneXarray((i != 1) ? systemp.liquid.Xp[i-1] : systemp.liquid.Xp[end],Nliquids-1,p.tube.L)
    Xarraysnewone = constructoneXarray(systemp.liquid.Xp[insert_i_1],Nliquids-1,p.tube.L)
    splice!(systemp.liquid.Xarrays,splice_liquid_i_1);
    splice!(systemp.liquid.Xarrays,splice_liquid_i_2);
    insert!(systemp.liquid.Xarrays,insert_i_1,Xarraysnewone);

    # splice!(systemp.liquid.Xarrays,i);
    # (i != 1) ? splice!(systemp.liquid.Xarrays,i-1) : splice!(systemp.liquid.Xarrays,length(systemp.liquid.Xarrays));
    # (i != 1) ? insert!(systemp.liquid.Xarrays,i-1,Xarraysnewone) : insert!(systemp.liquid.Xarrays,length(systemp.liquid.Xarrays)+1,Xarraysnewone);

    # θarraysnewone = (i != 1) ? [p.liquid.θarrays[i-1][1:end-1]; (p.liquid.θarrays[i-1][end-1]+p.liquid.θarrays[i][2])/2 ;p.liquid.θarrays[i][2:end]] : [p.liquid.θarrays[end][1:end-1]; (p.liquid.θarrays[end][end-1]+p.liquid.θarrays[i][2]) / 2 ;p.liquid.θarrays[i][2:end]]
    θarraysnewone = [p.liquid.θarrays[liquid_left_i][1:end-1]; (p.liquid.θarrays[liquid_left_i][end-1]+p.liquid.θarrays[liquid_right_i][2])/2 ;p.liquid.θarrays[liquid_right_i][2:end]]
    splice!(systemp.liquid.θarrays,splice_liquid_i_1);
    splice!(systemp.liquid.θarrays,splice_liquid_i_2);
    insert!(systemp.liquid.θarrays,insert_i_1,θarraysnewone);
    
    # splice!(systemp.liquid.θarrays,i);
    # (i != 1) ? splice!(systemp.liquid.θarrays,i-1) : splice!(systemp.liquid.θarrays,length(systemp.liquid.θarrays));
    # (i != 1) ? insert!(systemp.liquid.θarrays,i-1,θarraysnewone) : insert!(systemp.liquid.θarrays,length(systemp.liquid.θarrays)+1,θarraysnewone);

    return deepcopy(systemp)
end

function getmerge_flags(δv,sys)

    Xpvapor = getXpvapor(sys.liquid.Xp,sys.tube.L,sys.tube.closedornot)
    numofmergingsite = length(Xpvapor)
    merge_flags = Array{Bool,1}(undef, numofmergingsite)

    for i in 1:numofmergingsite
        merge_flags[i] = sys.tube.closedornot ? (mod(Xpvapor[i][2] - Xpvapor[i][1],sys.tube.L) < δv) : ((Xpvapor[i][2] - Xpvapor[i][1]) < δv)
    end

    return merge_flags
end
