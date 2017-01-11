immutable Particle
    x::Float64
    y::Float64
    z::Float64
    m::Float64
    pix::Int64 # parent index
end

function myrand()
    x=-1
    while x<0 || x>1
        x=randn()^3+0.5
    end
    x
end

make_particles(N) =
    [Particle(myrand(),myrand(),myrand(), 1.0, -1) for i in 1:N]

