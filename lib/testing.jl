function get_acc(particles, ix, eps2)
    ax=0.0
    ay=0.0
    az=0.0
    @inbounds for p in particles
        dx = getx(p)-getx(particles[ix])
        dy = gety(p)-gety(particles[ix])
        dz = getz(p)-getz(particles[ix])
        dr2 = dx*dx+dy*dy+dz*dz+eps2
        dr3 = dr2*sqrt(dr2)
        ax += dx/dr3
        ay += dy/dr3
        az += dz/dr3
    end
    ax,ay,az
end

function perf(particles, ax,ay,az, N, eps2)
    ixs = randperm(length(particles))[1:N]
    rax = zeros(N)
    ray = zeros(N)
    raz = zeros(N)
    vax = ax[ixs]
    vay = ay[ixs]
    vaz = az[ixs]
    @threads for i in 1:N
        tax,tay,taz = get_acc(particles, ixs[i], eps2)
        rax[i] = tax;
        ray[i] = tay;
        raz[i] = taz;
    end
    @show mean(abs(rax))
    dax = vax-rax
    day = vay-ray
    daz = vaz-raz
    da = sqrt(dax.^2+day.^2+daz.^2)
    a = sqrt(rax.^2+ray.^2+raz.^2)
    ee = abs(da./a.*100)
    hx = linspace(-2,2,35)
    mx = 0.5*(hx[2:end]+hx[1:(end-1)])
    my = zeros(length(hx)-1)
    for ei in ee
        for i in 2:length(hx)
            if log10(ei)<hx[i] && log10(ei)>hx[1]
                my[i-1] += 1.0
                break
            end
        end
    end

    ee50 = sort(ee)[round(Int64, N*0.5)]
    ee90 = sort(ee)[round(Int64, N*0.9)]
    ee95 = sort(ee)[round(Int64, N*0.9)]
    ee99 = sort(ee)[round(Int64, N*0.99)]
    mx,my, ee50,ee90,ee95,ee99, mean(log10(ee)), std(log10(ee))
end

