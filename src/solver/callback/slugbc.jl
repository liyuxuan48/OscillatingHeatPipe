export slugbc_affect!,slugbc_condition
# boiling_condition,
function slugbc_condition(u,t,integrator)
    return true
end

function slugbc_affect!(integrator)
    sys = getcurrentsys!(integrator.u,integrator.p);
    unew = deepcopy(integrator.u)

    θarrays = sys.liquid.θarrays
    
    adjacentT = getadjacentT.([sys],1:length(θarrays))
    
    for i in eachindex(θarrays)
        θarrays[i][1] = adjacentT[i][1]
        θarrays[i][end] = adjacentT[i][end]
    end
    
    index_dynamics_end = findfirst(x->abs(x+1e10) <= 10^(-1), unew)

    unew[index_dynamics_end:end] = liquidθtovec(θarrays)
    
    resize!(integrator.u,size(unew,1)::Int)
    integrator.u = deepcopy(unew)
end


function getadjacentT(p::PHPSystem,i::Int64)
    @unpack PtoT = p.tube
    Tfirst = PtoT(p.vapor.P[i])

    if p.tube.closedornot == true
        Tlast = i >= length(p.vapor.P) ? PtoT(p.vapor.P[1]) : PtoT(p.vapor.P[i+1])
    else
        Tlast = PtoT(p.vapor.P[i+1])
    end

    (Tfirst,Tlast)
end
