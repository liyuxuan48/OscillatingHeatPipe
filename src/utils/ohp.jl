export construct_oneloop_curve,construct_ohp_curve

function construct_ohp_curve(OHPshape::String,Δx::Real)


    if OHPshape == "ASETS"
        ds = 1.5Δx
        nturn = 16
        width_ohp = 46.25*1e-3
        length_ohp = 133.83*1e-3
        gap = 1e-3
        pitch = width_ohp/(2*nturn+1)
        x0, y0 = -length_ohp/2 -2e-3, -width_ohp/2

        return x, y, xf, yf = construct_ohp_curve(nturn,pitch,length_ohp,gap,ds,x0,y0,false,false,pi/2)
    end
    

    if OHPshape == "openloop"
        ds = 1.5Δx
        nturn = 16
        width_ohp = 46.25*1e-3
        length_ohp = 100*1e-3
        gap = 10*1e-3
        pitch = width_ohp/(2*nturn+1)
        x0, y0 = -length_ohp/2 -2e-3, -width_ohp/2

        return x, y, xf, yf = construct_ohp_curve_open(nturn,pitch,length_ohp,gap,ds,x0,y0,false,false,pi/2)
    end

    return error("Type wrong")
end

function construct_ohp_curve_open(nturn,pitch,height,gap,ds,x0,y0,flipx,flipy,angle)
    rad = 0.5*pitch
    len = height
    x, y = Float64[], Float64[]
    xf, yf = 0.0, -gap-rad
    for i in 1:nturn-1
        xi, yi, xf, yf = OscillatingHeatPipe._channel_pair(pitch,len,ds,xf,yf,false,false,0.0)
        append!(x,xi)
        append!(y,yi)
    end
       
    xi, yi, xf, yf = _closure_open(pitch,len,gap,2*nturn*pitch,ds,xf,yf,false,false,0.0)
    append!(x,xi)
    append!(y,yi)
    push!(x,xf)
    push!(y,yf)
    
    OscillatingHeatPipe._transform!(x,y,x0,y0,flipx,flipy,angle)
    x[1:end-1], y[1:end-1], x[end], y[end]
end

function construct_ohp_curve(nturn,pitch,height,gap,ds,x0,y0,flipx,flipy,angle)
    rad = 0.5*pitch
    len = height
    x, y = Float64[], Float64[]
    xf, yf = 0.0, -gap-rad
    for i in 1:nturn
        xi, yi, xf, yf = _channel_pair(pitch,len,ds,xf,yf,false,false,0.0)
        append!(x,xi)
        append!(y,yi)
    end
    xi, yi, xf, yf = _closure(pitch,len,gap,2*nturn*pitch,ds,xf,yf,false,false,0.0)
    append!(x,xi)
    append!(y,yi)
    push!(x,xf)
    push!(y,yf)
    _transform!(x,y,x0,y0,flipx,flipy,angle)
    x[1:end-1], y[1:end-1], x[end], y[end]
end



function _transform!(x,y,x0,y0,flipx,flipy,angle)
    xtmp, ytmp = copy(x), copy(y)
    if flipx
        ytmp .= -ytmp
    end
    if flipy
        xtmp .= -xtmp
    end
    ca, sa = cos(angle), sin(angle)
    x .= xtmp.*ca .- ytmp.*sa
    y .= xtmp.*sa .+ ytmp.*ca
    x .+= x0
    y .+= y0
    x, y
 end

 function _channel_pair(pitch,len,ds,x0,y0,flipx,flipy,angle)
     rad = 0.5*pitch
     x, y = Float64[], Float64[]
     xi, yi, xf, yf = _line(len,ds,0.0,0.0,false,false,-π/2)
     append!(x,xi)
     append!(y,yi)
     xi, yi, xf, yf = _halfturn(rad,ds,xf,yf,true,true,0.0)
     append!(x,xi)
     append!(y,yi)
     xi, yi, xf, yf = _line(len,ds,xf,yf,false,false,π/2)
     append!(x,xi)
     append!(y,yi)
     xi, yi, xf, yf = _halfturn(rad,ds,xf,yf,false,true,0.0)
     append!(x,xi)
     append!(y,yi)
     push!(x,xf)
     push!(y,yf)
     _transform!(x,y,x0,y0,flipx,flipy,angle)
     x[1:end-1], y[1:end-1], x[end], y[end]
 end

 function _closure(pitch,channellen,closuregap,closurelen,ds,x0,y0,flipx,flipy,angle)
    rad = 0.5*pitch
    x, y = Float64[], Float64[]
    xi, yi, xf, yf = _line(channellen,ds,0.0,0.0,false,false,-π/2)
    append!(x,xi)
    append!(y,yi)
    xi, yi, xf, yf = _halfturn(rad,ds,xf,yf,true,true,0.0)
    append!(x,xi)
    append!(y,yi)
    xi, yi, xf, yf = _line(channellen+closuregap,ds,xf,yf,false,false,π/2)
    append!(x,xi)
    append!(y,yi)
    xi, yi, xf, yf = _quarterturn(rad,ds,xf,yf,false,false,0.0)
    append!(x,xi)
    append!(y,yi)
    xi, yi, xf, yf = _line(closurelen,ds,xf,yf,false,false,π)
    append!(x,xi)
    append!(y,yi)
    xi, yi, xf, yf = _quarterturn(rad,ds,xf,yf,false,false,π/2)
    append!(x,xi)
    append!(y,yi)
    xi, yi, xf, yf = _line(closuregap,ds,xf,yf,false,false,-π/2)
    append!(x,xi)
    append!(y,yi)
    push!(x,xf)
    push!(y,yf)
    _transform!(x,y,x0,y0,flipx,flipy,angle)
    x[1:end-1], y[1:end-1], x[end], y[end]
end


function _closure_open(pitch,channellen,closuregap,closurelen,ds,x0,y0,flipx,flipy,angle)
    rad = 0.5*pitch
    x, y = Float64[], Float64[]
    xi, yi, xf, yf = OscillatingHeatPipe._line(channellen,ds,0.0,0.0,false,false,-π/2)
    append!(x,xi)
    append!(y,yi)
    xi, yi, xf, yf = OscillatingHeatPipe._halfturn(rad,ds,xf,yf,true,true,0.0)
    append!(x,xi)
    append!(y,yi)
    xi, yi, xf, yf = OscillatingHeatPipe._line(channellen,ds,xf,yf,false,false,π/2)
    append!(x,xi)
    append!(y,yi)
    push!(x,xf)
    push!(y,yf)
    OscillatingHeatPipe._transform!(x,y,x0,y0,flipx,flipy,angle)
    x[1:end-1], y[1:end-1], x[end], y[end]
end

####  Basic elements  ####

# in each element, adjust the spacing so that there is an integral number of points, uniformly spaced
# along entire shape. Sets the first point at the origin.

#=
create a set of points in an arc of radius `radius` and turning angle `turnangle`,
with a nominal spacing `ds` at (0,0) and proceeding counterclockwise about the center
at (-radius,0). Then the initial point is moved to (x0,y0). If either `flipx` and
`flipy` are set to true, then it flips the curve about that axis. The curve is
rotated by angle `angle`.
=#
function _turn(radius,turnangle,ds,x0,y0,flipx,flipy,angle)
    np1 = ceil(Int,turnangle*radius/ds)
    θ = range(0.0,turnangle,length=np1)
    x, y = Float64[], Float64[]
    for θi in θ
        push!(x,radius*cos(θi)-radius)
        push!(y,radius*sin(θi))
    end
    _transform!(x,y,x0,y0,flipx,flipy,angle)
    x[1:end-1], y[1:end-1], x[end], y[end]
end

_halfturn(radius,ds,x0,y0,flipx,flipy,angle) = _turn(radius,π,ds,x0,y0,flipx,flipy,angle)

_quarterturn(radius,ds,x0,y0,flipx,flipy,angle) = _turn(radius,π/2,ds,x0,y0,flipx,flipy,angle)


#=
create a set of points in a straight line of length `length` with nominal spacing `ds`
at (0,0) and proceeding in the x direction.
The initial point is moved to (x0,y0). If either `flipx` and `flipy` are set to
true, then it flips the curve about that axis. The curve is rotated by angle `angle`.
=#

function _line(length,ds,x0,y0,flipx,flipy,angle)
    np1 = ceil(Int,length/ds)
    s = range(0.0,length,length=np1)
    x, y = Float64[], Float64[]
    for si in s
        push!(x,si)
        push!(y,0.0)
    end
    _transform!(x,y,x0,y0,flipx,flipy,angle)
    x[1:end-1], y[1:end-1], x[end], y[end]
end
