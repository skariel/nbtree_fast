abstract AbstractParticle

immutable Particle <: AbstractParticle
    _x::Float64
    _y::Float64
    _z::Float64
    _m::Float64
    _pix::Int64 # parent index
end

getx(p::Particle) = p._x
gety(p::Particle) = p._y
getz(p::Particle) = p._z
getm(p::Particle) = p._m
getpix(p::Particle) = p._pix
@inline function setpix(v,i,val)
    @inbounds p = v[i]
    @inbounds v[i] = Particle(p._x, p._y, p._z, p._m, val)
end

function myrand()
    r=rand()^(1/2) #abs(randn()/5)^2
    x=0.0
    y=0.0
    z=0.0
    while true
        x1 = 2*rand()-1
        x2 = 2*rand()-1
        dx2 = x1*x1+x2*x2
        dx2 >= 1.0 && continue
        s = 2*sqrt(1-dx2)*r
        x = x1*s
        y = x2*s
        z = (1-2*dx2)*r
        return x,y,z
    end
end

make_particles(N) =
    [Particle(myrand()..., 1.0, -1) for i in 1:N]
