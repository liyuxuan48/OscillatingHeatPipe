import ConstrainedSystems: init, solve

export ODE_innertube,ODE_steadyfilm,timemarching!

function ODE_innertube(u,p,t)

    sys_init = p

    index_dynamics_end = findfirst(x->abs(x+1e10) <= 10^(-1), u)

    newsys = getcurrentsys!(u,sys_init)

    dynamicsdu = dynamicsmodel(u[1:index_dynamics_end-1],newsys)

    liquiddu = duliquidθtovec(liquidmodel(newsys))

    du = [dynamicsdu;liquiddu]
    
    return(du)

end

# function ODE_steadyfilm(u,p,t)

#     index_dynamics_end = findfirst(x->abs(x+1e10) <= 10^(-1), u)

#     newsys = getcurrentsys!(u,p)

#     dynamicsdu = dynamicsmodel_steadyfilm(u[1:index_dynamics_end-1],newsys)

#     liquiddu = duliquidθtovec(liquidmodel(newsys))

#     du = [dynamicsdu;liquiddu]

#     return(du)

# end

function temperature_linesource(integrator_plate)
    T = ComputationalHeatTransfer.temperature(integrator_plate)
    ohp_E = integrator_plate.p.extra_cache.fcache[end].region_cache.cache.E;
    ohp_E*T
end

# weakly coupled alternate time marching
function timemarching!(integrator_tube,integrator_plate,tstep::Float64,Tstart::Float64)

    currentsys = getcurrentsys!(integrator_tube.u,integrator_tube.p)
    currentsys.wall.θarray = temperature_linesource(integrator_plate)

    # println(currentsys.wall.θarray)
    # currentsys.wall.θarray = temperature_linesource(integrator_plate) .+ Tstart

    # println(currentsys.wall.θarray)
    qtmp = sys_to_heatflux(currentsys)
    # println(qtmp)
 
    integrator_plate.p.phys_params["ohp_flux"] = -qtmp # q source per length
    
    step!(integrator_tube,tstep,true);
    step!(integrator_plate,tstep,true);
    # ADI_timemarching!(temperature(integrator_plate),sys_plate,tstep)
    # integrator_plate.t += tstep
    integrator_tube,integrator_plate
end

function init(u_tube::Vector{Float64},tspan::Tuple{Any, Any},sys_tube::PHPSystem,kwargs...)
    
    cb_boiling =  DiscreteCallback(boiling_condition,boiling_affect!)
    cb_liquidmerging =  DiscreteCallback(merging_condition,merging_affect!)
    cb_vapormerging = DiscreteCallback(vaporMergingCondition,vaporMergingAffect!)
    cb_fixdx =  DiscreteCallback(fixdx_condition,fixdx_affect!)
    cb_slugbc = DiscreteCallback(slugbc_condition,slugbc_affect!)
    cbst = CallbackSet(cb_fixdx,cb_boiling,cb_vapormerging,cb_liquidmerging,cb_slugbc);
    
    # cb_liquidmerging =  DiscreteCallback(merging_condition,merging_affect!)
    # cbst = CallbackSet(cb_liquidmerging);

    prob = ODEProblem(ODE_innertube, u_tube, tspan, sys_tube,kwargs...) # construct integrator_tube problem
    return init(prob, alg=RK4(),dt=1e-3,save_on=false, callback=cbst,maxiters=1e10,kwargs...)
    # return init(prob, alg=RK4(),dt=1e-3,save_on=false, maxiters=1e10,kwargs...)
end