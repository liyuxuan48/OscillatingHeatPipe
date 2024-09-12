export constructmapping,sys_interpolation,slug_interp

function sys_interpolation(sys)
    if sys.tube.closedornot == true
            return sys_interpolation_closedloop(sys)
        else
            return sys_interpolation_openloop(sys)
    end
end

function sys_interpolation_openloop(sys)
    X_inner = Array{Float64}(undef, 0)
    θ_inner = Array{Float64}(undef, 0)
    H_inner = Array{Float64}(undef, 0)
    X_inner_pres  = Array{Float64}(undef, 0)
    P_inner = Array{Float64}(undef, 0)

    Xp  = sys.liquid.Xp

    @unpack PtoT = sys.propconvert
    θ = PtoT.(sys.vapor.P)
    P = sys.vapor.P
    δstart = sys.vapor.δstart
    δend = sys.vapor.δend
    Lfilm_start = sys.vapor.Lfilm_start
    Lfilm_end = sys.vapor.Lfilm_end
    Xpvapor = getXpvapor(Xp,sys.tube.closedornot)

    H_film_start = Hfilm.(δstart,[sys])
    H_film_end = Hfilm.(δend,[sys])

    Tavg = median(θ)

    H_vapor = sys.vapor.Hᵥ
    H_liquid = sys.liquid.Hₗ(Tavg)

    Nvapor = length(P)


    append!(X_inner,sys.liquid.Xarrays[1])
    append!(θ_inner,sys.liquid.θarrays[1])
    append!(H_inner,H_liquid .* ones(length(sys.liquid.Xarrays[1])))
    append!(X_inner_pres,[Xp[1][1],Xp[1][end]])
    append!(P_inner,[P[1], P[1]])

    for i = 1:Nvapor

        # vapors
                append!(X_inner,[Xpvapor[i][1],Xpvapor[i][1]+Lfilm_start[i],
                Xpvapor[i][1]+Lfilm_start[i],Xpvapor[i][end]-Lfilm_end[i],
                Xpvapor[i][end]-Lfilm_end[i],Xpvapor[i][end]])

                append!(θ_inner,[θ[i],θ[i],θ[i],θ[i],θ[i],θ[i]])
                append!(H_inner,[H_film_start[i], H_film_start[i],
                H_vapor,H_vapor,
                H_film_end[i],H_film_end[i]])

                append!(X_inner_pres,[Xpvapor[i][1],Xpvapor[i][end]])
                append!(P_inner,[P[i], P[i]])
    
        # liquids
            append!(X_inner,sys.liquid.Xarrays[i+1])
            append!(θ_inner,sys.liquid.θarrays[i+1])
            append!(H_inner,H_liquid .* ones(length(sys.liquid.Xarrays[i+1])))
            append!(X_inner_pres,[Xp[i+1][1],Xp[i+1][end]])
            append!(P_inner,[P[i], P[i]])
        end

    # ``` extend wall Xarray by adding its 0.0 point (plan to write this as a seperate function)```
    extend_wall_Xarray = deepcopy(sys.wall.Xarray)
    extend_wall_θarray = deepcopy(sys.wall.θarray)

    extend_wall_Xarray = [0.0;sys.wall.Xarray;sys.tube.L]
    extend_wall_θarray = [(sys.wall.θarray[1]+sys.wall.θarray[end])/2;sys.wall.θarray;(sys.wall.θarray[1]+sys.wall.θarray[end])/2]

    Interpolations.deduplicate_knots!(X_inner,move_knots = true)
    Interpolations.deduplicate_knots!(extend_wall_Xarray,move_knots = true)
    Interpolations.deduplicate_knots!(X_inner_pres,move_knots = true)

    θ_interp_liquidtowall = LinearInterpolation(X_inner, θ_inner);

    H_interp_liquidtowall = LinearInterpolation(X_inner, H_inner);

    θ_interp_walltoliquid = LinearInterpolation(extend_wall_Xarray, extend_wall_θarray);

    P_interp_liquidtowall = LinearInterpolation(X_inner_pres, P_inner);


    return θ_interp_walltoliquid, θ_interp_liquidtowall, H_interp_liquidtowall, P_interp_liquidtowall

end


function sys_interpolation_closedloop(sys)
    X_inner = Array{Float64}(undef, 0)
    θ_inner = Array{Float64}(undef, 0)
    H_inner = Array{Float64}(undef, 0)
    X_inner_pres  = Array{Float64}(undef, 0)
    P_inner = Array{Float64}(undef, 0)

    Xp  = sys.liquid.Xp

    @unpack PtoT = sys.propconvert
    θ = PtoT.(sys.vapor.P)
    P = sys.vapor.P
    δstart = sys.vapor.δstart
    δend = sys.vapor.δend
    Lfilm_start = sys.vapor.Lfilm_start
    Lfilm_end = sys.vapor.Lfilm_end
    Xpvapor = getXpvapor(Xp,sys.tube.closedornot)

    H_film_start = Hfilm.(δstart,[sys])
    H_film_end = Hfilm.(δend,[sys])

    Tavg = median(θ)

    H_vapor = sys.vapor.Hᵥ
    H_liquid = sys.liquid.Hₗ(Tavg)
    # H_vapor = sys.vapor.Hᵥ
    # H_liquid = sys.liquid.Hₗ

    Nvapor = length(P)

    loop_plus_index = [2:Nvapor;1]

    max_i = argmax(Xpvapor)
    max_j = argmax(Xpvapor[max_i])
  

    for i = 1:Nvapor

            if i == max_i && max_j != 2

                X_inner_loop_append,H_inner_loop_append = XHloop_append(i,H_film_start[i],H_film_end[i],sys)

                append!(X_inner,X_inner_loop_append)
                append!(θ_inner,[θ[i], θ[i], θ[i], θ[i], θ[i], θ[i], θ[i], θ[i]])
                append!(H_inner,H_inner_loop_append)

                append!(X_inner_pres,[Xpvapor[i][1], sys.tube.L, 0.0, Xpvapor[i][end]])
                append!(P_inner,[P[i], P[i], P[i], P[i]])
            else
                append!(X_inner,[Xpvapor[i][1],Xpvapor[i][1]+Lfilm_start[i],
                Xpvapor[i][1]+Lfilm_start[i],Xpvapor[i][end]-Lfilm_end[i],
                Xpvapor[i][end]-Lfilm_end[i],Xpvapor[i][end]])

                append!(θ_inner,[θ[i],θ[i],θ[i],θ[i],θ[i],θ[i]])
                append!(H_inner,[H_film_start[i], H_film_start[i],
                H_vapor,H_vapor,
                H_film_end[i],H_film_end[i]])

                append!(X_inner_pres,[Xpvapor[i][1],Xpvapor[i][end]])
                append!(P_inner,[P[i], P[i]])
    
        end
        

    if (sys.liquid.Xarrays[i][1] > sys.liquid.Xarrays[i][end])
        period_index = findmin(sys.liquid.Xarrays[i])[2]
        Xarrays_temp = deepcopy(sys.liquid.Xarrays[i])
        θarrays_temp = deepcopy(sys.liquid.θarrays[i])
        
        
        insert!(Xarrays_temp,period_index, sys.tube.L)
        insert!(Xarrays_temp,period_index+1, 0.0)
        insert!(θarrays_temp,period_index, sys.liquid.θarrays[i][period_index])
        insert!(θarrays_temp,period_index+1, sys.liquid.θarrays[i][period_index])
        append!(X_inner,Xarrays_temp)
        append!(θ_inner,θarrays_temp)
        
        append!(H_inner,H_liquid .* ones(length(sys.liquid.Xarrays[i]) + 2))
        

        append!(X_inner_pres,[sys.liquid.Xarrays[i][1], sys.tube.L, 0.0, sys.liquid.Xarrays[i][end]])
        P_inner_end = (sys.tube.L-sys.liquid.Xarrays[i][1])/mod(sys.liquid.Xarrays[i][end]-sys.liquid.Xarrays[i][1],sys.tube.L) * (P[loop_plus_index[i]]-P[i]) + P[i]
        append!(P_inner,[P[i], P_inner_end, P_inner_end, P[loop_plus_index[i]]])

        else
            append!(X_inner,sys.liquid.Xarrays[i])
            append!(θ_inner,sys.liquid.θarrays[i])
            append!(H_inner,H_liquid .* ones(length(sys.liquid.Xarrays[i])))

    end
end

        imin = argmin(X_inner)
        X_inner_final =Array{Float64}(undef, 0)
        θ_inner_final =Array{Float64}(undef, 0)
        H_inner_final =Array{Float64}(undef, 0)
        if imin != 1
            X_inner_final = [X_inner[imin:end];X_inner[1:imin-1]]
            θ_inner_final = [θ_inner[imin:end];θ_inner[1:imin-1]]
            H_inner_final = [H_inner[imin:end];H_inner[1:imin-1]]
        end

        imin = argmin(X_inner_pres)
        X_inner_pres_final =Array{Float64}(undef, 0)
        P_inner_final =Array{Float64}(undef, 0)
        if imin != 1
            X_inner_pres_final = [X_inner_pres[imin:end];X_inner_pres[1:imin-1]]
            P_inner_final = [P_inner[imin:end];P_inner[1:imin-1]]
        end
    
    ``` extend wall Xarray by adding its 0.0 point```
    extend_wall_Xarray = deepcopy(sys.wall.Xarray)
    extend_wall_θarray = deepcopy(sys.wall.θarray)

    extend_wall_Xarray = [0.0;sys.wall.Xarray;sys.tube.L]
    extend_wall_θarray = [(sys.wall.θarray[1]+sys.wall.θarray[end])/2;sys.wall.θarray;(sys.wall.θarray[1]+sys.wall.θarray[end])/2]

    Interpolations.deduplicate_knots!(X_inner_final,move_knots = true)
    Interpolations.deduplicate_knots!(extend_wall_Xarray,move_knots = true)
    Interpolations.deduplicate_knots!(X_inner_pres_final,move_knots = true)

    θ_interp_liquidtowall = LinearInterpolation(X_inner_final, θ_inner_final);

    H_interp_liquidtowall = LinearInterpolation(X_inner_final, H_inner_final);

    θ_interp_walltoliquid = LinearInterpolation(extend_wall_Xarray, extend_wall_θarray);

    P_interp_liquidtowall = LinearInterpolation(X_inner_pres_final, P_inner_final);


    return θ_interp_walltoliquid, θ_interp_liquidtowall, H_interp_liquidtowall, P_interp_liquidtowall

end

function XHloop_append(i,H_film_start,H_film_end,sys)


    Xp = sys.liquid.Xp[i]
    Xpvapor = getXpvapor(sys.liquid.Xp,sys.tube.closedornot)[i]
    Lfilm_start = sys.vapor.Lfilm_start[i]
    Lfilm_end = sys.vapor.Lfilm_end[i]

    H_vapor = sys.vapor.Hᵥ

    L = sys.tube.L

    if Xpvapor[2] > Xpvapor[1]
        # println("Xp error 1")
        return error("Xp error 1")
        
    elseif mod(Xpvapor[1],L) + Lfilm_start > L
        
        X_inner_loop_append = mod.([Xpvapor[1],L,
                0.0,Xpvapor[1]+Lfilm_start,
                Xpvapor[1]+Lfilm_start,Xpvapor[end]-Lfilm_end,
                Xpvapor[end]-Lfilm_end,Xpvapor[end]],L)
        X_inner_loop_append[2] = L
        
        H_inner_loop_append = [H_film_start, H_film_start,
                H_film_start, H_film_start,
                H_vapor,H_vapor,
                H_film_end,H_film_end]
        
    elseif Xpvapor[2] - Lfilm_end > 0.0
        
        X_inner_loop_append = mod.([Xpvapor[1],Xpvapor[1]+Lfilm_start,
                Xpvapor[1]+Lfilm_start, L,
                0.0, Xpvapor[end]-Lfilm_end,
                Xpvapor[end]-Lfilm_end,Xpvapor[end]],L)
        X_inner_loop_append[4] = L
        
        H_inner_loop_append = [H_film_start, H_film_start,
                H_vapor, H_vapor,
                H_vapor, H_vapor,
                H_film_end,H_film_end]
        
    elseif Xpvapor[2] > 0
        
        X_inner_loop_append = mod.([Xpvapor[1],Xpvapor[1]+Lfilm_start,
                Xpvapor[1]+Lfilm_start,  Xpvapor[end]-Lfilm_end,
                Xpvapor[end]-Lfilm_end,L,
                0.0,Xpvapor[end]],L)
        X_inner_loop_append[6] = L
        
        H_inner_loop_append = [H_film_start, H_film_start,
                H_vapor, H_vapor,
                H_film_end,H_film_end,
                H_film_end,H_film_end]
        
    else 
        println("Xp error 2")
        return "Xp error 2"
    end

return X_inner_loop_append,H_inner_loop_append
end

function slug_interp(sys::PHPSystem)
        filmLend = sys.vapor.Lfilm_end
    filmLstart = sys.vapor.Lfilm_start
    # Xpstart = [elem[1] for elem in sys.liquid.Xp]
    # Xpend = [elem[2] for elem in sys.liquid.Xp]
    Xpvapor = getXpvapor(sys.liquid.Xp,sys.tube.closedornot)
    Xvaporend = [elem[2] for elem in Xpvapor]
    Xvaporstart = [elem[1] for elem in Xpvapor]
    Xfilmstart = Xvaporstart .+ filmLstart
    Xfilmend = Xvaporend .- filmLend
    Xslugs = zeros(8*size(Xvaporstart)[1])
    slug_flags = zeros(8*size(Xvaporstart)[1])
    Xfinalslugs = zeros(8*size(Xvaporstart)[1])
    slug_flags = zeros(8*size(Xvaporstart)[1])

    for i in eachindex(Xvaporstart)
        Xslugs[8i-7:8i-6] = [Xvaporstart[i];Xvaporstart[i]]
        Xslugs[8i-5:8i-4] = [Xfilmstart[i];Xfilmstart[i]]
        Xslugs[8i-3:8i-2] = [Xfilmend[i];Xfilmend[i]]
        Xslugs[8i-1:8i] = [Xvaporend[i];Xvaporend[i]]

        slug_flags[8i-7:8i-6] = [1.0; 2.0]
        slug_flags[8i-5:8i-4] = [2.0; 0.0]
        slug_flags[8i-3:8i-2] = [0.0; 2.0]
        slug_flags[8i-1:8i] =   [2.0; 1.0]

        Xslugs = mod.(Xslugs,sys.tube.L)
    end

    Xmax, imax = findmax(Xslugs)
    #     Xmin, imin = findmax(Xslugs)

    edge_status = slug_flags[imax+1]

    xtemp = splice!(Xslugs,1:imax+1,0.0)
    slug_temp = splice!(slug_flags,1:imax+1,edge_status)
    append!(Xslugs,[xtemp;sys.tube.L])
    append!(slug_flags,[slug_temp;edge_status])

    Interpolations.deduplicate_knots!(Xslugs,move_knots = true)

    sluginterp = LinearInterpolation(Xslugs, slug_flags);
end