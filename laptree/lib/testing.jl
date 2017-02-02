function get_acc(particles, x,y,z, eps2)
    ax=0.0
    ay=0.0
    az=0.0
    @inbounds for p in particles
        dx = getx(p)-x
        dy = gety(p)-y
        dz = getz(p)-z
        dr2 = dx*dx+dy*dy+dz*dz+eps2
        dr3 = dr2*sqrt(dr2)
        ax += dx/dr3
        ay += dy/dr3
        az += dz/dr3
    end
    ax,ay,az
end

function get_err(ph,particles, x,y,z, alpha2,eps2)
    tax,tay,taz = get_acc(ph,length(ph),1, x,y,z, alpha2,eps2)
    eax,eay,eaz = get_acc(particles, x,y,z, eps2)
    dax = tax-eax
    day = tay-eay
    daz = taz-eaz
    da = sqrt(dax*dax+day*day+daz*daz)
    da/sqrt(eax*eax+eay*eay+eaz*eaz)*100.0
end

function get_err_arr(ph,particles,N, alpha2,eps2)
    arr = Float64[]
    ixs = Int64[]
    for i in 1:N
        ix = rand(1:length(particles))
        push!(ixs, ix)
        p = particles[ix]
        push!(arr, get_err(ph,particles, p._x,p._y,p._z, alpha2,eps2))
    end
    arr, ixs
end

function make_radial(x,y,z,ax,ay,az)
    r2 = x*x+y*y+z*z
    r = sqrt(r2)
    nx = x/r
    ny = y/r
    nz = z/r
    ax*nx + ay*ny + az*nz
end

function get_radial_err(ph,particles, x,y,z, alpha2,eps2)
    ta = make_radial(x,y,z,get_acc(ph,length(ph),1, x,y,z, alpha2,eps2)...)
    ea = make_radial(x,y,z,get_acc(particles, x,y,z, eps2)...)
    da = ta-ea
    da/ea
end

function get_radial_err_arr(ph,particles,N, alpha2,eps2)
    arr = Float64[]
    ixs = Int64[]
    for i in 1:N
        ix = rand(1:length(particles))
        push!(ixs, ix)
        p = particles[ix]
        push!(arr, get_radial_err(ph,particles, p._x,p._y,p._z, alpha2,eps2))
    end
    arr, ixs
end
