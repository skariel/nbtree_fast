
function find_nearest_neighbors(particles)
    data = zeros(3, length(particles))
    data[1,:] = [p._x for p in particles]
    data[2,:] = [p._y for p in particles]
    data[3,:] = [p._z for p in particles]
    kdtree = KDTree(data; leafsize = 10)
    ixs, d = knn(kdtree, data, 2, true)
    [ix[2] for ix in ixs]
end

function get_parents!(particles)
    nixs = find_nearest_neighbors(particles)
    parents = Particle[] 
    grouped = falses(length(particles))
    for ix in eachindex(particles)
        nix = nixs[ix]
        
        p = particles[ix]
        np = particles[nix]

        grouped[ix] && continue
        grouped[nix] && continue
        nixs[nix] != ix && continue

        # add a new parent for these two particles
        grouped[ix] = true
        grouped[nix] = true

        newm = p._m+np._m
        newx = (p._x*p._m + np._x*np._m)/newm
        newy = (p._y*p._m + np._y*np._m)/newm
        newz = (p._z*p._m + np._z*np._m)/newm

        dx = p._x-newx
        dy = p._y-newy
        dz = p._z-newz
        d = sqrt(dx*dx+dy*dy+dz*dz)
        dx = np._x-newx
        dy = np._y-newy
        dz = np._z-newz
        nd = sqrt(dx*dx+dy*dy+dz*dz)
        newr = max(d+p._r, nd+np._r)

        push!(parents, Particle(
            newx,newy,newz,
            -1,
            newm,
            ix,nix,newr,
        ))

        # update the two particles
        particles[ix] = Particle(
            p._x,p._y,p._z,
            length(parents),
            p._m,p._c1,p._c2,p._r,
        )
        particles[nix] = Particle(
            np._x,np._y,np._z,
            length(parents),
            np._m,np._c1,np._c2,np._r,
        )
    end

    # update ungrouped particles
    for ix in eachindex(particles)
        grouped[ix] && continue

        p = particles[ix]

        # add a new parent for this particle
        push!(parents, Particle(
            p._x,p._y,p._z,
            -1,
            p._m,ix,-1,p._r,
        ))

        # update the particle
        particles[ix] = Particle(
            p._x,p._y,p._z,
            length(parents),
            p._m,p._c1,p._c2,p._r,
        )
    end
    parents
end

function get_parent_hierarchy(particles)
    v = Vector{Particle}[]
    push!(v, particles)
    while true
        length(v[end]) < 2 && break
        push!(v, get_parents!(v[end]))
    end
    v
end

function get_leafs(v, i)
    leafs = Particle[]
    stack_l = Int64[length(v)]
    stack_ix = Int64[i]
    level=0
    ix=0
    while length(stack_l) > 0
        try
            level = pop!(stack_l)
            ix = pop!(stack_ix)

            p = v[level][ix]
            if p._c1<0 && p._c2<0
                push!(leafs, p)
                continue
            end
            if p._c1>0
                push!(stack_l, level-1)
                push!(stack_ix, p._c1)
            end
            if p._c2>0
                push!(stack_l, level-1)
                push!(stack_ix, p._c2)
            end
        catch e
            @show level
            @show ix
            @show e
            return
        end
    end
    leafs
end

function get_acc(v,l,i, x,y,z, alpha2, eps2)
    ax=0.0
    ay=0.0
    az=0.0

    p = v[l][i]

    # MAC
    dx = p._x-x
    dy = p._y-y
    dz = p._z-z
    r2 = dx*dx+dy*dy+dz*dz
    
    if p._r*p._r/r2 < alpha2
        # Success!
        r2 += eps2
        r = sqrt(r2)
        r3 = r2*r

        fac = p._m/r3
        ax += dx*fac
        ay += dy*fac
        az += dz*fac

        return ax,ay,az
    end

    # Failure!
    if p._c1>0
        _ax, _ay, _az = get_acc(v,l-1,p._c1, x,y,z, alpha2, eps2)
        ax += _ax
        ay += _ay
        az += _az
    end
    if p._c2>0
        _ax, _ay, _az = get_acc(v,l-1,p._c2, x,y,z, alpha2, eps2)
        ax += _ax
        ay += _ay
        az += _az
    end
        
    ax,ay,az
end