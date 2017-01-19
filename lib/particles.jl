immutable Particle
    x::Float64
    y::Float64
    z::Float64
    m::Float64
    pix::Int64 # parent index
end

function myrand()
    r=abs(randn()/5)^5
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

