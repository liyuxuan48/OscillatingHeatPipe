export dMdtdynamicsmodel,dynamicsmodel_steadyfilm,wallmodel,liquidmodel,dynamicsmodel,sys_to_heatflux,sys_to_Harray,integrator_to_heatflux,integrator_to_Harray

# this is a equation takes the current u and returns du/dt, p was already updated from the current u.
function dynamicsmodel(u::Array{Float64,1},p::PHPSystem)
    
    du = zeros(size(u))

    # extracts essenstial information stored at sys::PHPSystem
    sys = p
    @unpack d,peri,Ac,g,L,closedornot = sys.tube
    @unpack σ,μₗ,ρ,Xp,dXdt = sys.liquid
    ρₗ = ρ
    @unpack P,Eratio_plus,Eratio_minus,δstart,δend,Lfilm_start,Lfilm_end,ad_fac = sys.vapor
    @unpack L_newbubble = sys.wall

    # println(Xp)

    # number of liquid slugs
    numofliquidslug = length(Xp)

    # characteristic interface velocities for each liquid slug
    V = [mean(elem) for elem in dXdt]

    # get a characteristic Capilarry number based on the average velocities
    Vavg = mean(abs.(V))
    Ca = getCa.(μₗ,σ,Vavg)

    # get liquid film deposition area
    δdeposit = Catoδ(d,Ca,adjust_factor=ad_fac)
    Adeposit = getAdeposit(sys,δdeposit)
    Adeposit_left = [elem[1] for elem in Adeposit]
    Adeposit_right = [elem[2] for elem in Adeposit]

    # get actrual liquid film area
    Astart = getδarea(Ac,d,δstart)
    Aend = getδarea(Ac,d,δend)

    Lvaporplug = XptoLvaporplug(Xp,L,closedornot)
    Lliquidslug = XptoLliquidslug(Xp,L)

    # get gravity potential
    heightg_interp = sys.mapping.heightg_interp
    Xp1 = [elem[1] for elem in Xp]
    Xp2 = [elem[2] for elem in Xp]
    heightg = map(tuple,heightg_interp(Xp1),heightg_interp(Xp2))

    Xpvapor = getXpvapor(Xp,L,closedornot)

    # get differential equation factors
    lhs = ρₗ*Ac .* Lliquidslug
    rhs_press = Ac ./ lhs
    Re_list = ρₗ .* abs.(V) .* d ./ μₗ
    f_coefficient = f_churchill.(Re_list .+ 1e-4)
    dXdt_to_stress = -0.125 .* f_coefficient .* ρₗ .* V .* abs.(V)
    rhs_dXdt = peri .* Lliquidslug .* dXdt_to_stress ./ lhs
    rhs_g = Ac*ρₗ ./ lhs

    # parameters for boolean flags
    L0threshold_film = 0.15*L_newbubble
    L0threshold_pure_vapor = 0.5*L_newbubble

    if closedornot == false
        numofvaporbubble = numofliquidslug - 1

        # initialize liquid and vapor velocities for left(smaller ξ) and right(larger ξ) side
        # v_liquid_left_normal=zeros(numofliquidslug)
        # v_liquid_right_normal=zeros(numofliquidslug)
        v_vapor_left_normal=zeros(numofvaporbubble)
        v_vapor_right_normal=zeros(numofvaporbubble)

        # get average velocities of the liquid slugs
        v_momentum = @view u[2*numofliquidslug+1:2:4*numofliquidslug]
        v_momentum_vapor_end = v_momentum[2:end]
        v_momentum_vapor_start = v_momentum[1:end-1]
        v_liquid_left_normal  = v_momentum .+ v_momentum .* Adeposit_left ./ (Ac .- Adeposit_left)
        v_liquid_right_normal = v_momentum .+ v_momentum .* Adeposit_right ./(Ac .- Adeposit_right)
        
        rhs_dLdt = -v_momentum .* (v_liquid_right_normal .- v_liquid_left_normal) ./ Lliquidslug

        # vapor δ
        v_vapor_left_normal = v_liquid_right_normal[1:end-1]
        v_vapor_right_normal = v_liquid_left_normal[2:end]

        A_dδdt_right_vapor = Adeposit_left[2:end]
        A_dδdt_left_vapor = Adeposit_right[1:end-1]
  

        dMdt_latent_start,dMdt_latent_end,dMdt_latent_start_positive,dMdt_latent_end_positive = dMdtdynamicsmodel(Xpvapor,sys)
        dMdt_latent_start_negative = dMdt_latent_start .- dMdt_latent_start_positive
        dMdt_latent_end_negative =   dMdt_latent_end   .- dMdt_latent_end_positive

        # println(dMdt_latent_end)
        # println(dMdt_latent_end_positive)
        # println(dMdt_latent_end_negative)

        # println(dMdt_latent_start)
        # println(dMdt_latent_start_positive)
        # println(dMdt_latent_start_negative)


        dLdt_start,dLdt_end,dδdt_start,dδdt_end,v_vapor_start_final,v_vapor_end_final = film_dynamics(ρₗ, Ac,d,δstart,δend,dMdt_latent_start,dMdt_latent_end,dMdt_latent_start_positive,dMdt_latent_start_negative,Eratio_plus,Eratio_minus,
        v_vapor_left_normal,dMdt_latent_end_positive,dMdt_latent_end_negative,v_vapor_right_normal,A_dδdt_left_vapor,A_dδdt_right_vapor,
        Lvaporplug,Lfilm_start,Lfilm_end,L0threshold_film,L0threshold_pure_vapor,v_momentum_vapor_start,v_momentum_vapor_end,Astart,Aend)

        v_liquid_left_final = [0.0;vcat(v_vapor_end_final...)]
        v_liquid_right_final = [vcat(v_vapor_start_final...);0.0]

        # up tp here
        for i = 1:numofliquidslug
                du[2*i-1] = v_liquid_left_final[i]
                du[2*i] = v_liquid_right_final[i]

                if (i != 1) && (i != numofliquidslug)
                    du[2*numofliquidslug + 2*i-1] = rhs_dXdt[i] + rhs_g[i]*(heightg[i][1]-heightg[i][end]) + rhs_press[i] * (P[i-1]-P[i]) + rhs_dLdt[i]
                else
                    du[2*numofliquidslug + 2*i-1] = 0.0
                end                

                du[2*numofliquidslug + 2*i] = du[2*numofliquidslug + 2*i-1]

            end


            du[4*numofliquidslug+1:4*numofliquidslug+numofvaporbubble] = dMdt_latent_start+dMdt_latent_end
            du[4*numofliquidslug+numofvaporbubble+1:4*numofliquidslug+2*numofvaporbubble] .= 0.0 # equals to 0 for now
            du[4*numofliquidslug+2*numofvaporbubble+1:4*numofliquidslug+3*numofvaporbubble] .= 0.0 # equals to 0 for now
            # du[4*numofliquidslug+numofvaporbubble+1:4*numofliquidslug+2*numofvaporbubble] = dδdt_start # equals to 0 for now
            # du[4*numofliquidslug+2*numofvaporbubble+1:4*numofliquidslug+3*numofvaporbubble] = dδdt_end # equals to 0 for now
            du[4*numofliquidslug+3*numofvaporbubble+1:4*numofliquidslug+4*numofvaporbubble] = dLdt_start # equals to 0 for now
            du[4*numofliquidslug+4*numofvaporbubble+1:4*numofliquidslug+5*numofvaporbubble] = dLdt_end # equals to 0 for now


            # println(dδdt_end)
            # println(δdeposit)
            # println(Vavg)
            # println(v_momentum)
            # println(Adeposit)
            # println(Adeposit_left)
            # println(v_liquid_left_normal)
            # println(δstart)
            # println(δend)
            # println(Lfilm_end)
            # println(dLdt_start)
            # println(dLdt_end)
            # println(v_liquid_left_final)
            # println(v_liquid_right_final)

            return du

    end

    if closedornot == true

        numofvaporbubble = numofliquidslug

        # initialize liquid and vapor velocities for left(smaller ξ) and right(larger ξ) side
        # v_liquid_left_normal=zeros(numofliquidslug)
        # v_liquid_right_normal=zeros(numofliquidslug)
        v_vapor_left_normal=zeros(numofvaporbubble)
        v_vapor_right_normal=zeros(numofvaporbubble)

        # get average velocities of the liquid slugs
        v_momentum = @view u[2*numofliquidslug+1:2:4*numofliquidslug]
        v_momentum_vapor_end = v_momentum
        v_momentum_vapor_start = [v_momentum[end];v_momentum[1:end-1]]
        v_liquid_left_normal  = v_momentum .+ v_momentum .* Adeposit_left ./ (Ac .- Adeposit_left)
        v_liquid_right_normal = v_momentum .+ v_momentum .* Adeposit_right ./(Ac .- Adeposit_right)
        
        rhs_dLdt = -v_momentum .* (v_liquid_right_normal .- v_liquid_left_normal) ./ Lliquidslug

        # vapor δ
        v_vapor_left_normal[2:end] = @view v_liquid_right_normal[1:end-1]
        v_vapor_left_normal[1] = v_liquid_right_normal[end]
        v_vapor_right_normal = v_liquid_left_normal

        A_dδdt_right_vapor = Adeposit_left
        A_dδdt_right_liquid = Adeposit_right
        A_dδdt_left_vapor = zeros(size(A_dδdt_right_vapor))
        A_dδdt_left_vapor[2:end] = @view A_dδdt_right_liquid[1:end-1]
        A_dδdt_left_vapor[1] = A_dδdt_right_liquid[end]

        dMdt_latent_start,dMdt_latent_end,dMdt_latent_start_positive,dMdt_latent_end_positive = dMdtdynamicsmodel(Xpvapor,sys)
        dMdt_latent_start_negative = dMdt_latent_start .- dMdt_latent_start_positive
        dMdt_latent_end_negative =   dMdt_latent_end   .- dMdt_latent_end_positive

        dLdt_start,dLdt_end,dδdt_start,dδdt_end,v_vapor_start_final,v_vapor_end_final = film_dynamics(ρₗ, Ac,d,δstart,δend,dMdt_latent_start,dMdt_latent_end,dMdt_latent_start_positive,dMdt_latent_start_negative,Eratio_plus,Eratio_minus,
        v_vapor_left_normal,dMdt_latent_end_positive,dMdt_latent_end_negative,v_vapor_right_normal,A_dδdt_left_vapor,A_dδdt_right_vapor,
        Lvaporplug,Lfilm_start,Lfilm_end,L0threshold_film,L0threshold_pure_vapor,v_momentum_vapor_start,v_momentum_vapor_end,Astart,Aend)

        v_liquid_left_final = v_vapor_end_final
        v_liquid_right_final = [@view v_vapor_start_final[2:end];v_vapor_start_final[1]]

        for i = 1:numofliquidslug
                du[2*i-1] = v_liquid_left_final[i]
                du[2*i] = v_liquid_right_final[i]

                if i != numofliquidslug
                    du[2*numofliquidslug + 2*i-1] = rhs_dXdt[i] + rhs_g[i]*(heightg[i][1]-heightg[i][end]) + rhs_press[i] * (P[i]-P[i+1]) + rhs_dLdt[i]
                else
                    du[2*numofliquidslug + 2*i-1] = rhs_dXdt[i] + rhs_g[i]*(heightg[i][1]-heightg[i][end]) + rhs_press[i] * (P[i]-P[1]) + rhs_dLdt[i]
                end                

                du[2*numofliquidslug + 2*i] = du[2*numofliquidslug + 2*i-1]

            end

            du[4*numofliquidslug+1:5*numofliquidslug] = dMdt_latent_start+dMdt_latent_end
            du[5*numofliquidslug+1:6*numofliquidslug] = dδdt_start # equals to 0 for now
            du[6*numofliquidslug+1:7*numofliquidslug] = dδdt_end # equals to 0 for now
            du[7*numofliquidslug+1:8*numofliquidslug] = dLdt_start # equals to 0 for now
            du[8*numofliquidslug+1:9*numofliquidslug] = dLdt_end # equals to 0 for now

            return du
    end


end

function film_dynamics(ρₗ, Ac,d,δstart,δend,dMdt_latent_start,dMdt_latent_end,dMdt_latent_start_positive,dMdt_latent_start_negative,Eratio_plus,Eratio_minus,
    v_vapor_left_normal,dMdt_latent_end_positive,dMdt_latent_end_negative,v_vapor_right_normal,A_dδdt_left_vapor,A_dδdt_right_vapor,
    Lvaporplug,Lfilm_start,Lfilm_end,L0threshold_film,L0threshold_pure_vapor,v_momentum_vapor_start,v_momentum_vapor_end,Astart,Aend)

    F_start = ρₗ .* Ac .* 4 .* δstart .* (d .- δstart) ./ (d^2)
    C_start = ρₗ .* Ac .* 4 .* (d .- 2δstart) ./ (d^2)

    F_end = ρₗ .* Ac .* 4 .* δend .* (d .- δend) ./ (d^2)
    C_end = ρₗ .* Ac .* 4 .* (d .- 2δend) ./ (d^2)

    dLdt_start_normal = -(dMdt_latent_start_positive .* Eratio_plus .+ dMdt_latent_start_negative .* Eratio_minus) ./ F_start .- v_vapor_left_normal
    dLdt_end_normal = -(dMdt_latent_end_positive .* Eratio_plus .+ dMdt_latent_end_negative .* Eratio_minus) ./ F_end .+ v_vapor_right_normal

    # get boolean flags for five cases
    he_matrix_start,he_matrix_end = case_flags(Lvaporplug,Lfilm_start,Lfilm_end,L0threshold_film,L0threshold_pure_vapor,dLdt_start_normal,dLdt_end_normal)

    v_vapor_left_case5 = v_momentum_vapor_start .+ v_momentum_vapor_start .* Aend ./ (Ac .- Aend)
    v_vapor_right_case5 = v_momentum_vapor_end .+ v_momentum_vapor_end .* Astart ./ (Ac .- Astart)

    V_vapor_matrix_start = hcat(v_vapor_left_normal,v_momentum_vapor_start,v_vapor_left_normal,v_vapor_left_normal,v_vapor_left_case5)'
    V_vapor_matrix_end   = hcat(v_vapor_right_normal,v_momentum_vapor_end,v_vapor_right_normal,v_vapor_right_normal,v_vapor_right_case5)'

    v_vapor_start_final = sum(he_matrix_start .* V_vapor_matrix_start,dims=1)
    v_vapor_end_final = sum(he_matrix_end .* V_vapor_matrix_end,dims=1)

    dLdt_matrix_start = hcat(dLdt_start_normal,0 .* dLdt_start_normal,-v_vapor_left_normal,(v_vapor_right_case5 .- v_vapor_left_normal),0 .* dLdt_start_normal)'
    dLdt_matrix_end   = hcat(dLdt_end_normal,0 .* dLdt_end_normal,v_vapor_right_normal,(v_vapor_right_normal .- v_vapor_left_case5),0 .* dLdt_end_normal)'

    dLdt_start = sum(he_matrix_start .* dLdt_matrix_start,dims=1)
    dLdt_end = sum(he_matrix_end .* dLdt_matrix_end,dims=1)

    dδdt_start_normal = (-dMdt_latent_start .- F_start .* dLdt_start' .- ρₗ .* A_dδdt_left_vapor  .* v_vapor_left_normal) ./ (C_start .* Lfilm_start) 
    dδdt_end_normal = (-dMdt_latent_end     .- F_end   .* dLdt_end'   .+ ρₗ .* A_dδdt_right_vapor .* v_vapor_right_normal) ./ (C_end .* Lfilm_end)

    dδdt_start_case4 = (-dMdt_latent_start .- F_start .* dLdt_start' .- ρₗ .* A_dδdt_left_vapor  .* v_vapor_left_normal .+ ρₗ .* Astart  .* v_vapor_right_case5) ./ (C_start .* Lfilm_start) 
    dδdt_end_case4 = (-dMdt_latent_end     .- F_end   .* dLdt_end'   .+ ρₗ .* A_dδdt_right_vapor .* v_vapor_right_normal .- ρₗ .* Aend  .* v_vapor_left_case5) ./ (C_end .* Lfilm_end)

    dδdt_matrix_start = hcat(dδdt_start_normal,0 .* dδdt_start_normal,dδdt_start_normal,dδdt_start_case4,0 .* dδdt_start_normal)'
    dδdt_matrix_end   = hcat(dδdt_end_normal,0 .* dδdt_end_normal,dδdt_end_normal,dδdt_end_case4,0 .* dLdt_end_normal)'

    dδdt_start = sum(he_matrix_start .* dδdt_matrix_start,dims=1)
    dδdt_end = sum(he_matrix_end .* dδdt_matrix_end,dims=1)

    # println((-dMdt_latent_end     .- F_end   .* dLdt_end'   .+ ρₗ .* A_dδdt_right_vapor .* v_vapor_right_normal))
    # println(-dMdt_latent_end)
    # println( F_end   .* dLdt_end')
    # println(C_end .* Lfilm_end)
    # println(dδdt_end_normal)


    return dLdt_start,dLdt_end,dδdt_start,dδdt_end,v_vapor_start_final,v_vapor_end_final
end
    
function case_flags(Lvaporplug,Lfilm_start,Lfilm_end,L0threshold_film,L0threshold_pure_vapor,dLdt_start_normal,dLdt_end_normal)
    # get boolean flags for five cases
    he_start_short = Bool.(heaviside.(-Lfilm_start .+ L0threshold_film))
    he_end_short = Bool.(heaviside.(-Lfilm_end .+ L0threshold_film))
    he_start_positive = Bool.(heaviside.(dLdt_start_normal))
    he_end_positive = Bool.(heaviside.(dLdt_end_normal))
    he_meet= Bool.(heaviside.(-Lvaporplug .+ Lfilm_start .+ Lfilm_end .+ L0threshold_pure_vapor))


    # zero Lfilm and negative dLdt case
    case2_flag_start = Bool.((1 .- he_meet) .* he_start_short .* (1 .- he_start_positive))
    # two ends meet and both nonzero case
    case3_flag_start = Bool.(he_meet .* he_start_short .* he_start_positive .+ he_meet .* (1 .- he_start_short) .* (1 .- he_end_short))
    # two ends meet and other side zero case
    case4_flag_start = Bool.(he_meet .* (1 .- he_start_short) .* he_end_short .* (1 .- he_end_positive))
    # two ends meet and this side zero case
    case5_flag_start = Bool.(he_meet .* he_start_short .* (1 .- he_start_positive))
    # normal case
    case1_flag_start = Bool.((1 .- case3_flag_start) .* (1 .- case2_flag_start) .* (1 .- case4_flag_start) .* (1 .- case5_flag_start))

    # zero Lfilm and negative dLdt case
    case2_flag_end = Bool.((1 .- he_meet) .* he_end_short .* (1 .- he_end_positive))
    # two ends meet and both nonzero case
    case3_flag_end = Bool.(he_meet .* he_end_short .* he_end_positive .+ he_meet .* (1 .- he_end_short) .* (1 .- he_start_short))
    # two ends meet and other side zero case
    case4_flag_end = Bool.(he_meet .* (1 .- he_end_short) .* he_start_short .* (1 .- he_start_positive))
    # two ends meet and this side zero case
    case5_flag_end = Bool.(he_meet .* he_end_short .* (1 .- he_end_positive))            
    # normal case
    case1_flag_end = Bool.((1 .- case3_flag_end) .* (1 .- case2_flag_end) .* (1 .- case4_flag_end) .* (1 .- case5_flag_end))

    # he_matrix_start = hcat([case1_flag_start';case2_flag_start';case3_flag_start';case4_flag_start';case5_flag_start'])
    # he_matrix_end   = hcat([case1_flag_end';case2_flag_end';case3_flag_end';case4_flag_end';case5_flag_end'])

    he_matrix_start = hcat(case1_flag_start,case2_flag_start,case3_flag_start,case4_flag_start,case5_flag_start)'
    he_matrix_end   = hcat(case1_flag_end,case2_flag_end,case3_flag_end,case4_flag_end,case5_flag_end)'

    return he_matrix_start,he_matrix_end
end

function dMdtdynamicsmodel(Xpvapor::Array{Tuple{Float64,Float64},1},sys::PHPSystem)

    dMdt_latent_start=zeros(length(Xpvapor))
    dMdt_latent_end=zeros(length(Xpvapor))
    dMdt_latent_start_positive=zeros(length(Xpvapor))
    dMdt_latent_end_positive=zeros(length(Xpvapor))

    L = sys.tube.L

    peri = sys.tube.peri
    Ac = sys.tube.Ac
    k = sys.vapor.k

    Lfilm_start = sys.vapor.Lfilm_start
    Lfilm_end = sys.vapor.Lfilm_end

    δstart = sys.vapor.δstart
    δend = sys.vapor.δend
    δmin = sys.vapor.δmin

    P = sys.vapor.P
    θarrays = sys.liquid.θarrays
    Xarrays = sys.liquid.Xarrays

    @unpack PtoT,PtoHfg = sys.tube
    θ = PtoT.(P)

    Hfg = PtoHfg.(P)

    
    H_interp = sys.mapping.H_interp_liquidtowall
    θ_wall_interp = sys.mapping.θ_interp_walltoliquid

    dx_wall = sys.tube.L/sys.tube.N

    # Lvapor = XptoLvaporplug(sys.liquid.Xp,sys.tube.L,sys.tube.closedornot)
    # Lvapor_pure = max.(Lvapor - Lfilm_start - Lfilm_end,0.0)

    # axialrhs_total = 0

    for i = 1:length(Xpvapor)

        Nstart = Int64(max(2 , div(Lfilm_start[i],dx_wall)))
        
        heatflux_start = quad_trap(H_interp,θ_wall_interp,θ[i],Xpvapor[i][1],mod(Xpvapor[i][1]+Lfilm_start[i],L),L,Nstart)
        heatflux_start_positive = quad_trap_positive(H_interp,θ_wall_interp,θ[i],Xpvapor[i][1],mod(Xpvapor[i][1]+Lfilm_start[i],L),L,Nstart)

        Nend = Int64(max(2 , div(Lfilm_end[i],dx_wall)))

        heatflux_end = quad_trap(H_interp,θ_wall_interp,θ[i],mod(Xpvapor[i][2]-Lfilm_end[i],L),Xpvapor[i][2],L,Nend)
        heatflux_end_positive = quad_trap_positive(H_interp,θ_wall_interp,θ[i],mod(Xpvapor[i][2]-Lfilm_end[i],L),Xpvapor[i][2],L,Nend)

        slope_r = getslope(θarrays[i][2],θarrays[i][1],Xarrays[i][2],Xarrays[i][1])
        slope_l = (i == 1) ? getslope(θarrays[1][end],θarrays[1][end-1],Xarrays[1][end],Xarrays[1][end-1]) : getslope(θarrays[i-1][end],θarrays[i-1][end-1],Xarrays[i-1][end],Xarrays[i-1][end-1])


        axial_rhs_end = Ac*k*slope_r /Hfg[i]
        axial_rhs_start = Ac*k*(-slope_l) /Hfg[i]

        dMdt_latent_start[i] = heatflux_start*peri/Hfg[i] + axial_rhs_start
        dMdt_latent_end[i] = heatflux_end*peri/Hfg[i] + axial_rhs_end
        dMdt_latent_start_positive[i] = heatflux_start_positive*peri/Hfg[i] + axial_rhs_start
        dMdt_latent_end_positive[i] = heatflux_end_positive*peri/Hfg[i]+ axial_rhs_end

        end

        dMdt_latent_start = dMdt_latent_start.*heaviside.(δstart .- δmin)
        dMdt_latent_end = dMdt_latent_end.*heaviside.(δend .- δmin)
        dMdt_latent_start_positive = dMdt_latent_start_positive.*heaviside.(δstart .- δmin)
        dMdt_latent_end_positive = dMdt_latent_end_positive.*heaviside.(δend .- δmin)

            return dMdt_latent_start,dMdt_latent_end,dMdt_latent_start_positive,dMdt_latent_end_positive
end

function liquidmodel(p::PHPSystem)
    sys = p
    θarrays = sys.liquid.θarrays
    # nondihv_tonondihl = 0.0046206704347650325 # temperary variable to fix different nondimensionlaization

    du = 0 .* θarrays

    # γ = sys.vapor.γ
    Ac = sys.tube.Ac

    Hₗ = sys.liquid.Hₗ
    peri = sys.tube.peri
    α = sys.liquid.α
    Cpₗ = sys.liquid.Cp
    ρₗ = sys.liquid.ρ

    H_rhs = peri / (ρₗ*Cpₗ*Ac)

    for i in eachindex(θarrays)
        
        xs = sys.liquid.Xarrays[i];
        dx = mod(xs[2] - xs[1], sys.tube.L)

        θ_wall_inter = sys.mapping.θ_interp_walltoliquid

        fx = map(θ_wall_inter, xs) - θarrays[i]
        du[i] = α .* laplacian(θarrays[i]) ./ dx ./ dx + Hₗ .* fx .* H_rhs
    end
    return du
end


"""
    This is a function to get the laplacian of a vector field u

    For now zero-gradient boundary condition is used.

    u    ::  an array
"""


function laplacian(u)
    unew = zeros(size(u))

    dl = ones(length(u)-1)
    dr = dl
    d  = -2*ones(length(u))
    d[1] = -1
    d[end]= -1

    A = Tridiagonal(dl, d, dr)

    unew = A*u

    return (unew)
end


# q'
function sys_to_heatflux(p::PHPSystem)

    sys = p

    θarray = sys.wall.θarray
    # γ = sys.vapor.γ
    Hₗ = sys.liquid.Hₗ
    # He = sys.evaporator.He
    # k = sys.vapor.k
    # δ = sys.vapor.δ
    # Hvapor = k ./ δ

    peri = sys.tube.peri

    # dx = sys.wall.Xarray[2]-sys.wall.Xarray[1]

    # Xarray = sys.wall.Xarray
    θ_interp_liquidtowall = sys.mapping.θ_interp_liquidtowall
    H_interp_liquidtowall = sys.mapping.H_interp_liquidtowall

    xs =  sys.wall.Xarray

    dθarray = map(θ_interp_liquidtowall, xs) .- θarray
    Harray  = map(H_interp_liquidtowall, xs)

    # qwallarray = -Harray.*dθarray
    qwallarray = -Harray.*dθarray*peri
end


function integrator_to_heatflux(inte)
    sys_to_heatflux(inte.p)
end

function sys_to_Harray(p::PHPSystem)

    sys = p

    H_interp_liquidtowall = sys.mapping.H_interp_liquidtowall

    xs =  sys.wall.Xarray

    Harray  = map(H_interp_liquidtowall, xs)

    Harray
end

function quad_trap(H_interp,θ_interp,θvapor_one, a,b,L,N) 
    h = mod((b-a),L)/N
    int = h * ( H_interp(a)*(θ_interp(a)-θvapor_one) + H_interp(b)*(θ_interp(b)-θvapor_one) ) / 2
    for k=1:N-1
        xk = mod(mod((b-a),L) * k/N + a,L)
        int = int + h*H_interp(xk)*(θ_interp(xk)-θvapor_one)
    end
    return int
end

function quad_trap_positive(H_interp,θ_interp,θvapor_one, a,b,L,N) 
    h = mod((b-a),L)/N
    int = maximum([h * ( H_interp(a)*(θ_interp(a)-θvapor_one) + H_interp(b)*(θ_interp(b)-θvapor_one) ) / 2, 0.0])
    for k=1:N-1
        xk = mod(mod((b-a),L) * k/N + a,L)
        int = int + maximum([h*H_interp(xk)*(θ_interp(xk)-θvapor_one),0.0])
    end
    return int
end

function integrator_to_Harray(inte)
    sys_to_Harray(inte.p)
end

getslope(y2,y1,x2,x1) = (y2-y1)/(x2-x1)

heaviside(x::AbstractFloat) = ifelse(x < 0, zero(x), ifelse(x > 0, one(x), oftype(x,0.0)))