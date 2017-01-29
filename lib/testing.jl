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

get_random_node() = Node(
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
    )

get_random_exp() = NodeExp(
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

function calc_pot(n,e, x,y,z)
    dx = x-n.x
    dy = y-n.y
    dz = z-n.z
    dx2 = dx*dx
    dy2 = dy*dy
    dz2 = dz*dz

    return e.px*dx + 
        e.pxx*dx2/2 +
        e.pxxx*dx2*dx/6 +
        e.pxxy*dx2/2*dy +
        e.pxxz*dx2/2*dz +
        e.pxy*dx*dy +
        e.pxyy*dx*dy2/2 +
        e.pxyz*dx*dy*dz +
        e.pxz*dx*dz +
        e.pxzz*dx*dz2/2 +
        e.py*dy +
        e.pyy*dy2/2 +
        e.pyyy*dy2*dy/6 +
        e.pyyz*dy2/2*dz +
        e.pyz*dy*dz +
        e.pyzz*dy*dz2/2 +
        e.pz*dz +
        e.pzz*dz2/2 +
        e.pzzz*dz2*dz/6    
end

function get_numerical_ax(n,e, x,y,z)
    dx = 1.0e-8*(x-n.x)
    (calc_pot(n,e, x+dx,y,z) - calc_pot(n,e, x,y,z))/dx
end

function get_numerical_ay(n,e, x,y,z)
    dy = 1.0e-8*(y-n.y)
    (calc_pot(n,e, x,y+dy,z) - calc_pot(n,e, x,y,z))/dy
end

function get_numerical_az(n,e, x,y,z)
    dz = 1.0e-8*(z-n.z)
    (calc_pot(n,e, x,y,z+dz) - calc_pot(n,e, x,y,z))/dz
end

function test_acc()
    tot = 0
    for i in 1:1000000
        n = get_random_node()
        e = get_random_exp()
        x=randn(); y=randn(); z=randn();
        dax,day,daz = get_accel_from_node(n,e, x,y,z)
        abs(x-n.x)<1e-1&&continue
        abs(y-n.y)<1e-1&&continue
        abs(z-n.z)<1e-1&&continue
        abs(dax)<1e-3 && continue
        abs(day)<1e-3 && continue
        abs(daz)<1e-3 && continue
        tot+=1
        nax = get_numerical_ax(n,e, x,y,z)
        nay = get_numerical_ay(n,e, x,y,z)
        naz = get_numerical_az(n,e, x,y,z)
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


function test_mov()
    tot = 0
    for i in 1:1000000
        n1 = get_random_node()
        e1 = get_random_exp()
        x=randn(); y=randn(); z=randn();
        dax1,day1,daz1 = get_accel_from_node(n1, e1, x,y,z)
        n2 = get_random_node()
        e2 = get_random_exp()
        dax2,day2,daz2 = get_accel_from_node(n2, e2, x,y,z)

        abs(x-n1.x)<1e-1&&continue
        abs(y-n1.y)<1e-1&&continue
        abs(z-n1.z)<1e-1&&continue
        abs(dax1)<1e-3 && continue
        abs(day1)<1e-3 && continue
        abs(daz1)<1e-3 && continue

        abs(x-n2.x)<1e-1&&continue
        abs(y-n2.y)<1e-1&&continue
        abs(z-n2.z)<1e-1&&continue
        abs(dax2)<1e-3 && continue
        abs(day2)<1e-3 && continue
        abs(daz2)<1e-3 && continue

        tot+=1

        dax= dax1 + dax2
        day= day1 + day2
        daz= daz1 + daz2
         
        e1 = add_expansion_to_n1(n1,e1, n2,e2)
        nax, nay, naz = get_accel_from_node(n1,e1, x,y,z)

        dx = abs((nax-dax)/dax)
        dy = abs((nay-day)/day)
        dz = abs((naz-daz)/daz)
        try
            @assert dx<1.0e-5
        catch
            @show dx
            @show n1.x,n1.y,n1.z
            @show x,y,z
            @show dax, nax
            @show dax1, dax2
            break
        end
        try
            @assert dy<1.0e-5
        catch
            @show dy
            @show n1.x,n1.y,n1.z
            @show x,y,z
            @show day, nay
            break
        end
        try
            @assert dz<1.0e-5
        catch
            @show dz
            @show n1.x,n1.y,n1.z
            @show x,y,z
            @show daz, naz
            break
        end
    end
    tot
end


function test_exp()
    n1 = Node(
        -2.0, 2.0, 1.0, 1.0, 1.0,
        1.0, -1.0, 1.0, -1.0, 1.0, -1.0,
        0, -1, -1, -1, -1, -1, 
    )   
    e1 = NodeExp()

    p = Particle(5.0,0.0,0.0,1.0,-1)

    # make them interact...
    x=p.x
    y=p.y
    z=p.z
    m=p.m

    dx = x-n1.x
    dy = y-n1.y
    dz = z-n1.z
    dx2 = dx*dx
    dy2 = dy*dy
    dz2 = dz*dz
    dr2 = dx2 + dy2 + dz2
    dr = sqrt(dr2)
    M = n1.m+m

    fac = 1.0
    alpha = 0.5
    l = 0.0
    @assert (n1.l + l)/dr < alpha/fac # MAC test

    # MAC succesful, execute interaction
    dr3 = dr2*dr
    dr5 = dr3*dr2
    dr7 = dr2*dr5
    dx3 = dx2*dx
    dy3 = dy2*dy
    dz3 = dz2*dz

    @show dx
    @show dx3
    px = dx/dr3
    py = dy/dr3
    pz = dz/dr3

    pxx = (3*dx2-dr2)/dr5
    pyy = (3*dy2-dr2)/dr5
    pzz = (3*dz2-dr2)/dr5

    pxy = 3*dx*dy/dr5
    pxz = 3*dx*dz/dr5
    pyz = 3*dy*dz/dr5

    pxxx = 3*dx*(5*dx2-3*dr2)/dr7
    pyyy = 3*dy*(5*dy2-3*dr2)/dr7
    pzzz = 3*dz*(5*dz2-3*dr2)/dr7

    pxxy = -3*dy*(dr2-5*dx2)/dr7
    pxxz = -3*dz*(dr2-5*dx2)/dr7
    pyyz = -3*dz*(dr2-5*dy2)/dr7
    pyzz = -3*dy*(dr2-5*dz2)/dr7
    pxyy = -3*dx*(dr2-5*dy2)/dr7
    pxzz = -3*dx*(dr2-5*dz2)/dr7

    pxyz = 15*dx*dy*dz/dr7

    e1 = NodeExp(
        m*px,
        m*pxx,
        m*pxxx,
        m*pxxy,
        m*pxxz,
        m*pxy,
        m*pxyy,
        m*pxyz,
        m*pxz,
        m*pxzz,
        m*py,
        m*pyy,
        m*pyyy,
        m*pyyz,
        m*pyz,
        m*pyzz,
        m*pz,
        m*pzz,
        m*pzzz,
    )   
    tx = n1.minx
    ty = n1.miny
    tz = n1.maxz
    @show get_accel_from_node(n1,e1, n1.x,n1.y,n1.z)[1] 
    @show get_accel_from_node(n1,e1, tx,ty,tz)[1]

    dx = x-tx
    dy = y-ty
    dz = z-tz
    dr2 = dx*dx+dy*dy+dz*dz
    @show dx*m/dr2/sqrt(dr2)

    @show px
    @show pxx
    @show pxxx

    @show pxy
end

