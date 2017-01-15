using Base.Test

function get_acc(particles, ix, eps2)
    ax=0.0
    ay=0.0
    az=0.0
    @inbounds for p in particles
        dx = p.x-particles[ix].x
        dy = p.y-particles[ix].y
        dz = p.z-particles[ix].z
        dr2 = dx*dx+dy*dy+dz*dz+eps2
        dr3 = dr2*sqrt(dr2)
        ax += dx/dr3
        ay += dy/dr3
        az += dz/dr3
    end
    ax,ay,az
end

function perf(particles, ax, ay, az, N, eps2)
    ixs = randperm(length(t.particles))[1:N]
    rax = zeros(N)
    vax = ax[ixs]
    ray = zeros(N)
    vay = ay[ixs]
    raz = zeros(N)
    vaz = az[ixs]
    @threads for i in 1:N
        tax,tay,taz = get_acc(particles, ixs[i], eps2)
        rax[i] = tax;
        ray[i] = tay;
        raz[i] = taz;
    end
    eex = abs((vax-rax)./rax*100)
    eey = abs((vay-ray)./ray*100)
    eez = abs((vaz-raz)./raz*100)
    ix = round(Int64, N*0.995)
    pex = sort(eex)[ix]
    pey = sort(eey)[ix]
    pez = sort(eez)[ix]
    @show pex, pey, pez
end

function test_in_cell_mass(t::Tree)
    for p in t.particles
        @assert p.pix > 0
        i=0
        n = t.nodes[p.pix]
        while true
            i += 1
            @assert i<256
            if n.pix == 1
                @assert i>4
                break
            end
            n = t.nodes[n.pix]
        end
    end
    nn = zeros(Int64, length(t.particles))
    six = 1
    t.stack[1] = 1
    while six > 0
        n = t.nodes[t.stack[six]]
        six -= 1
        if n.cix1<0 && n.cix2<0
            for j in n.iix:n.fix
                nn[j] += 1
            end
            continue
        end
        if n.cix1 > 0
            six += 1
            t.stack[six] = n.cix1
        end
        if n.cix2 > 0
            six += 1
            t.stack[six] = n.cix2
        end
    end
    for x in nn
        @assert x == 1
    end

    for nix in 1:t.num_nodes_used
        n = t.nodes[nix]
        @assert n.l == max(n.maxx-n.minx, n.maxy-n.miny, n.maxz-n.minz)
        m = 0.0
        for pix in n.iix:n.fix
            p = t.particles[pix]
            @assert p.x <= n.maxx && p.x >= n.minx
            @assert p.y <= n.maxy && p.y >= n.miny
            @assert p.z <= n.maxz && p.z >= n.minz
            m += p.m
        end
        @test_approx_eq n.m m
        m = 0.0
        if n.cix1 > 0
            n1 = t.nodes[n.cix1]
            if n.dir==0
                @assert n1.minx == n.minx
                @test_approx_eq n1.maxx (n.minx+n.maxx)/2
                @assert n1.miny == n.miny && n1.maxy == n.maxy
                @assert n1.minz == n.minz && n1.maxz == n.maxz
            elseif n.dir==1
                @assert n1.minx == n.minx && n1.maxx == n.maxx
                @assert n1.miny == n.miny
                @test_approx_eq n1.maxy (n.miny+n.maxy)/2
                @assert n1.minz == n.minz && n1.maxz == n.maxz
            elseif n.dir==2
                @assert n1.minx == n.minx && n1.maxx == n.maxx
                @assert n1.miny == n.miny && n1.maxy == n.maxy
                @assert n1.minz == n.minz
                @test_approx_eq n1.maxz (n.minz+n.maxz)/2
            else
                error("bad direction!")
            end
            m += n1.m
        end
        if n.cix2 > 0
            n2 = t.nodes[n.cix2]
            if n.dir==0
                @test_approx_eq n2.minx (n.minx+n.maxx)/2
                @assert n2.maxx == n.maxx
                @assert n2.miny == n.miny && n2.maxy == n.maxy
                @assert n2.minz == n.minz && n2.maxz == n.maxz
            elseif n.dir==1
                @assert n2.minx == n.minx && n2.maxx == n.maxx
                @test_approx_eq n2.miny (n.miny+n.maxy)/2
                @assert n2.maxy == n.maxy
                @assert n2.minz == n.minz && n2.maxz == n.maxz
            elseif n.dir==2
                @assert n2.minx == n.minx && n2.maxx == n.maxx
                @assert n2.miny == n.miny && n2.maxy == n.maxy
                @test_approx_eq n2.minz (n.minz+n.maxz)/2
                @assert n2.maxz == n.maxz
            else
                error("bad direction!")
            end
            m += n2.m
        end
        if n.cix1>0 || n.cix2>0
            @test_approx_eq n.m m
        end
    end
end