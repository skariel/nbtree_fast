
immutable NodeExp
    px::Float64
    pxx::Float64
    pxxx::Float64
    pxxy::Float64
    pxxz::Float64
    pxy::Float64
    pxyy::Float64
    pxyz::Float64
    pxz::Float64
    pxzz::Float64
    py::Float64
    pyy::Float64
    pyyy::Float64
    pyyz::Float64
    pyz::Float64
    pyzz::Float64
    pz::Float64
    pzz::Float64
    pzzz::Float64
end
NodeExp() = NodeExp(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

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
end

Node() = Node(
        0.0,                 # x::Float64
        0.0,                 # y::Float64
        0.0,                 # z::Float64
        0.0,                 # m::Float64
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0,                   # direction
        -1,                  # pix::Int64 # parent index
        -1,                  # iix::Int64 # first particle index
        -1,                  # fix
        -1,                  # cix1::Int64 # first child index
        -1,                  # cix2::Int64 # second child index        
    )

type Tree{T, U, V}
    total_mass::Float64
    nodes::U
    exps::V
    particles::T
    stack1::Vector{Int64}
    stack2::Vector{Int64}
    stack3::Vector{Int64}
    num_nodes_used::Int64
    S::Int64                 # maximum number of particles per node
end

function Tree(particles, nodes, exps, S)
    stack1 = zeros(Int64,10000)
    stack2 = zeros(Int64,10000)
    stack3 = zeros(Int64,10000)
    num_nodes_used = 0
    Tree(
        0.0,
        nodes,
        exps,
        particles,
        stack1,
        stack2,
        stack3,
        num_nodes_used,
        S)    
end

function Tree(particles, S)
    nodec = round(Int64, 0.8*length(particles))
    nodes = Node[Node() for i in 1:nodec]
    exps = NodeExp[NodeExp() for i in 1:nodec]
    Tree(particles, nodes, exps, S)
end

function getminmax(t::Tree)
    minx = 1.0e30 # infinity, ha!
    maxx = -1.0e30 # minus infinity, ha!
    miny = 1.0e30 # infinity, ha!
    maxy = -1.0e30 # minus infinity, ha!
    minz = 1.0e30 # infinity, ha!
    maxz = -1.0e30 # minus infinity, ha!
    for p in t.particles
        if getx(p)<minx
            minx=getx(p)
        elseif getx(p)>maxx
            maxx=getx(p)
        end
        if gety(p)<miny
            miny=gety(p)
        elseif gety(p)>maxy
            maxy=gety(p)
        end
        if getz(p)<minz
            minz=getz(p)
        elseif getz(p)>maxz
            maxz=getz(p)
        end
    end
    dx = (maxx-minx)*0.001
    dy = (maxy-miny)*0.001
    dz = (maxz-minz)*0.001
    minx-dx,maxx+dx, miny-dy,maxy+dy, minz-dz,maxz+dz
end

function get_minmax_low(n::Node)
    if n.dir==0
        n.minx,n.maxx-(n.maxx-n.minx)/2, n.miny,n.maxy, n.minz,n.maxz
    elseif n.dir==1
        n.minx,n.maxx, n.miny,n.maxy-(n.maxy-n.miny)/2, n.minz,n.maxz
    else
        n.minx,n.maxx, n.miny,n.maxy, n.minz,n.maxz-(n.maxz-n.minz)/2
    end
end

function get_minmax_high(n::Node)
    if n.dir==0
        n.minx+(n.maxx-n.minx)/2,n.maxx, n.miny,n.maxy, n.minz,n.maxz
    elseif n.dir==1
        n.minx,n.maxx, n.miny+(n.maxy-n.miny)/2,n.maxy, n.minz,n.maxz
    else
        n.minx,n.maxx, n.miny,n.maxy, n.minz+(n.maxz-n.minz)/2,n.maxz
    end
end

function group!(t::Tree)
    group!(t,t)
end

function group!(t::Tree, tminmax::Tree)
    minx,maxx, miny,maxy, minz,maxz = getminmax(tminmax)

    # create root node
    node_ix = 1
    t.nodes[node_ix] = Node(
        0.0,                 # x::Float64
        0.0,                 # y::Float64
        0.0,                 # z::Float64
        0.0,                 # m::Float64
        0.0,                 # l
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
    )
    # stack contains parents to be splitted
    # add root to stack
    stack_ix = 1
    t.stack1[stack_ix] = 1

    @fastmath @inbounds while stack_ix > 0
        # pop node to split
        pix = t.stack1[stack_ix]  # parent index to be splitted
        pn = t.nodes[pix]        # parent node to be splitted
        stack_ix -= 1

        split = splitdir!(t.particles, pn)
        cix1 = -1
        cix2 = -1

        if split >= pn.iix
            # we have a left child!
            # lets create it...
            node_ix += 1
            minx,maxx, miny,maxy, minz,maxz = get_minmax_low(pn)
            t.nodes[node_ix] = Node(
                0.0,               # x::Float64
                0.0,               # y::Float64
                0.0,               # z::Float64
                0.0,               # m::Float64
                0.0,               # l
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
            )
            if split-pn.iix+1 > t.S # comparing number of particles
                # we have enough particles to split this node
                # push it to the stack!
                stack_ix += 1
                t.stack1[stack_ix] = node_ix
            else
                # node has not enough particles to be splitted
                # tell the particles of their new parent!
                for i in pn.iix:split
                    setpix(t.particles, i, node_ix)
                end
            end
            cix1 = node_ix
        end

        if split < pn.fix
            # we have a right child!
            # lets create it...
            node_ix += 1            
            minx,maxx, miny,maxy, minz,maxz = get_minmax_high(pn)
            t.nodes[node_ix] = Node(
                0.0,                 # x::Float64
                0.0,                 # y::Float64
                0.0,                 # z::Float64
                0.0,                 # m::Float64
                0.0,
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
            )
            if pn.fix-split > t.S # comparing number of particles
                # we have enough particles to pslit this node
                # push it to the stack!
                stack_ix += 1
                t.stack1[stack_ix] = node_ix
            else
                # node has not enough particles to be splitted
                # tell the particles of their new parent!
                for i in (split+1):pn.fix
                    setpix(t.particles, i, node_ix)
                end
            end
            cix2 = node_ix
        end
        if cix1>0 || cix2>0
            # update parent node that it has a new child/s
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
                cix1, ### <<<--- this is the update!!!
                cix2, ### <<<--- this is the update!!!
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
@inline splitx!(from_ix, to_ix, x, at) = split!(from_ix, to_ix, x, at, (p,_x)->getx(p)<_x)
@inline splity!(from_ix, to_ix, y, at) = split!(from_ix, to_ix, y, at, (p,_y)->gety(p)<_y)
@inline splitz!(from_ix, to_ix, z, at) = split!(from_ix, to_ix, z, at, (p,_z)->getz(p)<_z)

@inline function splitdir!(x, node)
    if node.dir==0
        splitx!(node.iix, node.fix, x, (node.minx+node.maxx)/2)
    elseif node.dir==1
        splity!(node.iix, node.fix, x, (node.miny+node.maxy)/2)
    else
        splitz!(node.iix, node.fix, x, (node.minz+node.maxz)/2)
    end
end













