
immutable Node
    x::Float64
    y::Float64
    z::Float64
    m::Float64
    l::Float64 # maximal geometrical dimension
    maxx::Float64
    minx::Float64
    maxy::Float64
    miny::Float64
    maxz::Float64
    minz::Float64
    dir::Int64 # direction node should be splitted (mod 3, 0==x, 1==y, 2==z)
    pix::Int64 # parent index
    iix::Int64 # first particle index
    fix::Int64 # final particle index
    cix1::Int64 # first child index
    cix2::Int64 # second child index

    qxx::Float64
    qyy::Float64
    qxy::Float64
    qxz::Float64
    qyz::Float64
    
end

type Tree
    nodes::Vector{Node}
    particles::Vector{Particle}
    stack::Vector{Int64}
    num_nodes_used::Int64
end

function Tree(particles)
    nodes = Array{Node}(round(Int64, 4.8*length(particles)))
    Tree(nodes, particles, Array{Int64}(10000), 0)
end

function getminmax(t::Tree)
    minx = 1.0e30 # infinity, ha!
    maxx = -1.0e30 # minus infinity, ha!
    miny = 1.0e30 # infinity, ha!
    maxy = -1.0e30 # minus infinity, ha!
    minz = 1.0e30 # infinity, ha!
    maxz = -1.0e30 # minus infinity, ha!
    for p in t.particles
        if p.x<minx
            minx=p.x
        elseif p.x>maxx
            maxx=p.x
        end
        if p.y<miny
            miny=p.y
        elseif p.y>maxy
            maxy=p.y
        end
        if p.z<minz
            minz=p.z
        elseif p.z>maxz
            maxz=p.z
        end
    end
    dx = (maxx-minx)*0.01
    dy = (maxy-miny)*0.01
    dz = (maxz-minz)*0.01
    minx-dx,maxx+dx, miny-dy,maxy+dy, minz-dz,maxz+dz
end

function get_minmax_low(n::Node)
    if n.dir==0
        n.minx,n.maxx-(n.maxx-n.minx)/2, n.miny,n.maxy, n.minz,n.maxz
    elseif n.dir==1
        n.minx,n.maxx, n.miny,n.maxy-(n.maxy-n.miny)/2, n.minz,n.maxz
    elseif n.dir==2
        n.minx,n.maxx, n.miny,n.maxy, n.minz,n.maxz-(n.maxz-n.minz)/2
    else
        error("bad direction!")
    end
end

function get_minmax_high(n::Node)
    if n.dir==0
        n.minx+(n.maxx-n.minx)/2,n.maxx, n.miny,n.maxy, n.minz,n.maxz
    elseif n.dir==1
        n.minx,n.maxx, n.miny+(n.maxy-n.miny)/2,n.maxy, n.minz,n.maxz
    elseif n.dir==2
        n.minx,n.maxx, n.miny,n.maxy, n.minz+(n.maxz-n.minz)/2,n.maxz
    else
        error("bad direction!")
    end
end

function group!(t::Tree, S::Int64)
    minx,maxx, miny,maxy, minz,maxz = getminmax(t)
    l = max(maxx-minx, maxy-miny, maxz-minz)

    # create root node
    node_ix = 1
    t.nodes[node_ix] = Node(
        0.0,                 # x::Float64
        0.0,                 # y::Float64
        0.0,                 # z::Float64
        0.0,                 # m::Float64
        l,
        maxx,
        minx,
        maxy,
        miny,
        maxz,
        minz,
        0,                   # direction
        -1,                  # pix::Int64 # parent index
        1,                   # iix::Int64 # first particle index
        length(t.particles), # fix
        -1,                  # cix1::Int64 # first child index
        -1,                  # cix2::Int64 # second child index        
        0,0,0,0,0,
    )
    # stack contains parents to be splitted
    # add root to stack
    stack_ix = 1
    t.stack[stack_ix] = 1

    @inbounds while stack_ix > 0
        # pop node to split
        pix = t.stack[stack_ix]  # parent index to be splitted
        pn = t.nodes[pix]        # parent node to be splitted
        stack_ix -= 1

        split = splitdir!(t.particles, pn)

        if split >= pn.iix
            # we have a left child!
            # lets create it...
            node_ix += 1
            minx,maxx, miny,maxy, minz,maxz = get_minmax_low(pn)
            l = max(maxx-minx, maxy-miny, maxz-minz)
            t.nodes[node_ix] = Node(
                0.0,               # x::Float64
                0.0,               # y::Float64
                0.0,               # z::Float64
                0.0,               # m::Float64
                l,
                maxx,
                minx,
                maxy,
                miny,
                maxz,
                minz,
                (pn.dir+1)%3,
                pix,               # pix::Int64 # parent index
                pn.iix,            # iix::Int64 # first particle index
                split,             # pnum::Int64 # number of particles
                -1,                # cix1::Int64 # first child index
                -1,                # cix2::Int64 # second child index        
                0,0,0,0,0,        
            )
            if split-pn.iix+1 > S # comparing number of particles
                # we have enough particles to split this node
                # push it to the stack!
                stack_ix += 1
                t.stack[stack_ix] = node_ix
            else
                # node has not enough particles to be splitted
                # tell the particles of their new parent!
                for i in pn.iix:split
                    p = t.particles[i]
                    t.particles[i] = Particle(p.x,p.y,p.z,p.m,node_ix)
                end
            end
            # update parent node that it has a new child
            t.nodes[pix] = Node(
                pn.x,
                pn.y,
                pn.z,
                pn.m,
                pn.l,
                pn.maxx,
                pn.minx,
                pn.maxy,
                pn.miny,
                pn.maxz,
                pn.minz,
                pn.dir,
                pn.pix,
                pn.iix,
                pn.fix,
                node_ix, ### <<<--- this is the update!!!
                pn.cix2, 
                0,0,0,0,0,                
            )
        end

        if split < pn.fix
            # we have a right child!
            # lets create it...
            node_ix += 1            
            minx,maxx, miny,maxy, minz,maxz = get_minmax_high(pn)
            l = max(maxx-minx, maxy-miny, maxz-minz)
            t.nodes[node_ix] = Node(
                0.0,                 # x::Float64
                0.0,                 # y::Float64
                0.0,                 # z::Float64
                0.0,                 # m::Float64
                l,
                maxx,
                minx,
                maxy,
                miny,
                maxz,
                minz,
                (pn.dir+1)%3,
                pix,                 # pix::Int64 # parent index
                split+1,             # iix::Int64 # first particle index
                pn.fix,                # pnum::Int64 # number of particles
                -1,                  # cix1::Int64 # first child index
                -1,                  # cix2::Int64 # second child index        
                0,0,0,0,0,        
            )
            if pn.fix-split > S # comparing number of particles
                # we have enough particles to pslit this node
                # push it to the stack!
                stack_ix += 1
                t.stack[stack_ix] = node_ix
            else
                # node has not enough particles to be splitted
                # tell the particles of their new parent!
                for i in (split+1):pn.fix
                    p = t.particles[i]
                    t.particles[i] = Particle(p.x,p.y,p.z,p.m,node_ix)
                end
            end
            # update parent node that it has a new child
            pn = t.nodes[pix]
            t.nodes[pix] = Node(
                pn.x,
                pn.y,
                pn.z,
                pn.m,
                pn.l,
                pn.maxx,
                pn.minx,
                pn.maxy,
                pn.miny,
                pn.maxz,
                pn.minz,
                pn.dir,
                pn.pix,
                pn.iix,
                pn.fix,
                pn.cix1,
                node_ix, ### <<<--- this is the update!!!
                0,0,0,0,0,        
            )
        end
    end
    t.num_nodes_used = node_ix
    nothing
end

@inline function split!(from_ix, to_ix, x, at, cmp)
    lix = from_ix
    rix = to_ix
    @inbounds while lix < rix
        if cmp(x[lix], at)
            lix += 1
            continue
        end
        x[rix], x[lix] = x[lix], x[rix]
        rix -= 1
    end
    !cmp(x[lix], at) ? (lix-1) : lix
end
@inline splitx!(from_ix, to_ix, x, at) = split!(from_ix, to_ix, x, at, (u,v)->u.x<v)
@inline splity!(from_ix, to_ix, x, at) = split!(from_ix, to_ix, x, at, (u,v)->u.y<v)
@inline splitz!(from_ix, to_ix, x, at) = split!(from_ix, to_ix, x, at, (u,v)->u.z<v)

@inline function splitdir!(x, node)
    if node.dir==0
        splitx!(node.iix, node.fix, x, (node.minx+node.maxx)/2)
    elseif node.dir==1
        splity!(node.iix, node.fix, x, (node.miny+node.maxy)/2)
    elseif node.dir==2
        splitz!(node.iix, node.fix, x, (node.minz+node.maxz)/2)
    else
        error("bad direction!")
    end
end













