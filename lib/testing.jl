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

function perf(particles, ax, N, eps2)
    ixs = randperm(length(t.particles))[1:N]
    rax = zeros(N)
    vax = ax[ixs]
    @inbounds for i in 1:N
        tax,tay,taz = get_acc(particles, ixs[i], eps2)
        rax[i] = tax;
    end
    ee = abs((vax-rax)./rax*100)
    sort(ee)[round(Int64, N*0.99)]
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
    t.stack1[1] = 1
    while six > 0
        n = t.nodes[t.stack1[six]]
        six -= 1
        if n.cix1<0 && n.cix2<0
            for j in n.iix:n.fix
                nn[j] += 1
            end
            continue
        end
        if n.cix1 > 0
            six += 1
            t.stack1[six] = n.cix1
        end
        if n.cix2 > 0
            six += 1
            t.stack1[six] = n.cix2
        end
    end
    for x in nn
        @assert x == 1
    end

    for nix in 1:t.num_nodes_used
        n = t.nodes[nix]
        m = 0.0
        for pix in n.iix:n.fix
            p = t.particles[pix]
            @assert p.x <= n.maxx && p.x >= n.minx
            @assert p.y <= n.maxy && p.y >= n.miny
            @assert p.z <= n.maxz && p.z >= n.minz
            dx = p.x-n.x
            dy = p.y-n.y
            dz = p.z-n.z
            dr2 = dx*dx+dy*dy+dz*dz
            try
                @assert (n.l*n.l - dr2) >= -1.0e-13
            catch
                @show n.l*n.l
                @show dr2
                @show nix, n.m
                return
            end
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

function get_random_node()
    Node(
        randn(), # x::Float64
        randn(), # y::Float64
        randn(), # z::Float64
        1.0, # m::Float64
        abs(randn()+1.0), # l::Float64 # maximal geometrical dimension
        rand()+1, # maxx::Float64
        rand()-1, # minx::Float64
        rand()+1, # maxy::Float64
        rand()-1, # miny::Float64
        rand()+1, # maxz::Float64
        rand()-1, # minz::Float64
        0 ,  # dir::Int64 # direction node should be splitted (mod 3, 0==x, 1==y, 2==z)
        -1,  # pix::Int64 # parent index
        -1,  # iix::Int64 # first particle index
        -1,  # fix::Int64 # final particle index
        -1,  # cix1::Int64 # first child index
        -1,  # cix2::Int64 # second child index

        randn(), # px::Float64
        randn(), # pxx::Float64
        randn(), # pxxx::Float64
        randn(), # pxxy::Float64
        randn(), # pxxz::Float64
        randn(), # pxy::Float64
        randn(), # pxyy::Float64
        randn(), # pxyz::Float64
        randn(), # pxz::Float64
        randn(), # pxzz::Float64
        randn(), # py::Float64
        randn(), # pyy::Float64
        randn(), # pyyy::Float64
        randn(), # pyyz::Float64
        randn(), # pyz::Float64
        randn(), # pyzz::Float64
        randn(), # pz::Float64
        randn(), # pzz::Float64
        randn(), # pzzz::Float64
    )
end

function calc_pot(n, x,y,z)
    dx = x-n.x
    dy = y-n.y
    dz = z-n.z
    dx2 = dx*dx
    dy2 = dy*dy
    dz2 = dz*dz

    return n.px*dx + 
        n.pxx*dx2/2 +
        n.pxxx*dx2*dx/6 +
        n.pxxy*dx2/2*dy +
        n.pxxz*dx2/2*dz +
        n.pxy*dx*dy +
        n.pxyy*dx*dy2/2 +
        n.pxyz*dx*dy*dz +
        n.pxz*dx*dz +
        n.pxzz*dx*dz2/2 +
        n.py*dy +
        n.pyy*dy2/2 +
        n.pyyy*dy2*dy/6 +
        n.pyyz*dy2/2*dz +
        n.pyz*dy*dz +
        n.pyzz*dy*dz2/2 +
        n.pz*dz +
        n.pzz*dz2/2 +
        n.pzzz*dz2*dz/6    
end

function get_numerical_ax(n, x,y,z)
    dx = 1.0e-8*(x-n.x)
    (calc_pot(n,x+dx,y,z) - calc_pot(n,x,y,z))/dx
end

function get_numerical_ay(n, x,y,z)
    dy = 1.0e-8*(y-n.y)
    (calc_pot(n,x,y+dy,z) - calc_pot(n,x,y,z))/dy
end

function get_numerical_az(n, x,y,z)
    dz = 1.0e-8*(z-n.z)
    (calc_pot(n,x,y,z+dz) - calc_pot(n,x,y,z))/dz
end

function test_acc()
    tot = 0
    for i in 1:1000000
        n = get_random_node()
        x=randn(); y=randn(); z=randn();
        dax,day,daz = get_accel_from_node(n, x,y,z)
        abs(x-n.x)<1e-1&&continue
        abs(y-n.y)<1e-1&&continue
        abs(z-n.z)<1e-1&&continue
        abs(dax)<1e-3 && continue
        abs(day)<1e-3 && continue
        abs(daz)<1e-3 && continue
        tot+=1
        nax = get_numerical_ax(n, x,y,z)
        nay = get_numerical_ay(n, x,y,z)
        naz = get_numerical_az(n, x,y,z)
        dx = abs((nax-dax)/dax)
        dy = abs((nay-day)/day)
        dz = abs((naz-daz)/daz)
        try
            @assert dx<1.0e-2
        catch
            @show dx
            @show n.x,n.y,n.z
            @show x,y,z
            @show dax, nax
            break
        end
        try
            @assert dy<1.0e-2
        catch
            @show dy
            @show n.x,n.y,n.z
            @show x,y,z
            @show day, nay
            break
        end
        try
            @assert dz<1.0e-2
        catch
            @show dz
            @show n.x,n.y,n.z
            @show x,y,z
            @show daz, naz
            break
        end
    end
    tot
end



