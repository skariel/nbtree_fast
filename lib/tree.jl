
immutable Node
    x::Float64
    y::Float64
    z::Float64
    m::Float64
    l::Float64 # maximal cell side length
    pix::Int64 # parent index
    iix::Int64 # first particle index
    pnum::Int64 # number of particles
    cix1::Int64 # first child index
    cix2::Int64 # second child index
end

immutable _GroupingStackData
    iix::Int64 # initial particle index
    fix::Int64 # final particle index
    dir::Int64 # direction, mod 3 (0==x, 1==y, 2==z)
    midx::Float64 # split coordinate
    midy::Float64
    midz::Float64
    pix::Int64 # parent node index
    slx::Float64 # cell side length (nothing to do with actual particle inside!)
    sly::Float64
    slz::Float64
end


type Tree
    nodes::Vector{Node}
    particles::Vector{Particle}
    grouping_stack::Vector{_GroupingStackData}
    nodes_stack::Vector{Int64}
    num_nodes_used::Int64
end

function Tree(particles)
    nodes = Array{Node}(round(Int64, 4.8*length(particles)))
    Tree(nodes, particles,
        Array{_GroupingStackData}(10000), Array{Int64}(10000), 0)
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
    minx,maxx, miny,maxy, minz,maxz
end

function get_sl_mid_low(data, dir)
    if dir==0
        data.slx/2, data.midx-data.slx/4, data.sly, data.midy, data.slz, data.midz
    elseif dir==1
        data.slx, data.midx, data.sly/2, data.midy-data.sly/4, data.slz, data.midz        
    elseif dir==2
        data.slx, data.midx, data.sly, data.midy, data.slz/2, data.midz-data.slz/4
    else
        error("bad direction!")
    end
end

function get_sl_mid_high(data, dir)
    if dir==0
        data.slx/2, data.midx+data.slx/4, data.sly, data.midy, data.slz, data.midz
    elseif dir==1
        data.slx, data.midx, data.sly/2, data.midy+data.sly/4, data.slz, data.midz        
    elseif dir==2
        data.slx, data.midx, data.sly, data.midy, data.slz/2, data.midz+data.slz/4
    else
        error("bad direction!")
    end
end

function group!(t::Tree, S::Int64)
    minx,maxx, miny,maxy, minz,maxz = getminmax(t)
    newsl = max(maxx-minx, maxy-miny, maxz-minz)

    # create root node
    node_ix = 1
    @inbounds t.nodes[node_ix] = Node(
        0.0,                 # x::Float64
        0.0,                 # y::Float64
        0.0,                 # z::Float64
        0.0,                 # m::Float64
        newsl,               # l::Float64 # maximal cell side length
        -1,                  # pix::Int64 # parent index
        1,                   # fpix::Int64 # first particle index
        length(t.particles), # pnum::Int64 # number of particles
        -1,                  # cix1::Int64 # first child index
        -1,                  # cix2::Int64 # second child index        
    )

    # node in the stack are splitted, add root to stack
    stack_ix = 1
    @inbounds t.grouping_stack[stack_ix] = _GroupingStackData(
       1,                    # iix::Int64 # initial particle index
       length(t.particles),  # fix::Int64 # final particle index
       0,                    # dir::Int64 # direction, mod 3 (0==x, 1==y, 2==z)
       0.5*(minx+maxx),      # mid::Float64 # split coordinate
       0.5*(miny+maxy),      # mid::Float64 # split coordinate
       0.5*(minz+maxz),      # mid::Float64 # split coordinate
       node_ix,              # pix::Int64 # parent node index
       maxx-minx,            # sl::Float64 # cell side length (nothing to do with actual particle inside!)
       maxy-miny,            # sl::Float64 # cell side length (nothing to do with actual particle inside!)
       maxz-minz,            # sl::Float64 # cell side length (nothing to do with actual particle inside!)
    )

    @inbounds while stack_ix > 0
        stack = t.grouping_stack[stack_ix]
        stack_ix -= 1

        split = splitdir!(stack.iix, stack.fix, t.particles, stack, stack.dir)

        if split >= stack.iix
            # we have a left child!
            # lets create it...
            node_ix += 1
            pnum = split-stack.iix+1 # particle number
            slx,midx, sly,midy, slz,midz = get_sl_mid_low(stack, stack.dir)
            newsl = max(slz,sly,slz)
            t.nodes[node_ix] = Node(
                0.0,               # x::Float64
                0.0,               # y::Float64
                0.0,               # z::Float64
                0.0,               # m::Float64
                newsl,             # l::Float64 # maximal cell side length
                stack.pix,         # pix::Int64 # parent index
                stack.iix,         # iix::Int64 # first particle index
                pnum,              # pnum::Int64 # number of particles
                -1,                # cix1::Int64 # first child index
                -1,                # cix2::Int64 # second child index        
            )
            if pnum > S
                # we have enough particles to pslit this node
                # push it to the stack!
                stack_ix += 1
                t.grouping_stack[stack_ix] = _GroupingStackData(
                    stack.iix,           # iix::Int64 # initial particle index
                    split,               # fix::Int64 # final particle index
                    (stack.dir+1)%3,     # dir::Int64 # direction (1==x, 2==y, 3==z)
                    midx,                # mid::Float64 # split coordinate
                    midy,                # mid::Float64 # split coordinate
                    midz,                # mid::Float64 # split coordinate
                    node_ix,             # pix::Int64 # parent node index
                    slx,                 # sl::Float64 # cell side length (nothing to do with actual particle inside!)
                    sly,                 # sl::Float64 # cell side length (nothing to do with actual particle inside!)
                    slz,                 # sl::Float64 # cell side length (nothing to do with actual particle inside!)
                )
            else
                # node has not enough particles to be splitted
                # tell the particles of their new parent!
                for i in stack.iix:split
                    p = t.particles[i]
                    t.particles[i] = Particle(p.x,p.y,p.z,p.m,node_ix)
                end
            end
            # update parent node that it has a new child
            pn = t.nodes[stack.pix]
            t.nodes[stack.pix] = Node(
                pn.x,
                pn.y,
                pn.z,
                pn.m,
                pn.l,
                pn.pix,
                pn.iix,
                pn.pnum,
                node_ix, ### <<<--- this is the update!!!
                pn.cix2, 
            )
        end

        if split < stack.fix
            # we have a right child!
            # lets create it...
            node_ix += 1
            pnum = stack.fix-split # particle number
            slx,midx, sly,midy, slz,midz = get_sl_mid_high(stack, stack.dir)
            newsl = max(slx,sly,slz)
            t.nodes[node_ix] = Node(
                0.0,               # x::Float64
                0.0,               # y::Float64
                0.0,               # z::Float64
                0.0,               # m::Float64
                newsl,             # l::Float64 # maximal 1d separation between particles (not the node side length!)
                stack.pix,         # pix::Int64 # parent index
                split+1,           # iix::Int64 # first particle index
                pnum,              # pnum::Int64 # number of particles
                -1,                # cix1::Int64 # first child index
                -1,                # cix2::Int64 # second child index        
            )
            if pnum > S
                # we have enough particles to pslit this node
                # push it to the stack!
                stack_ix += 1
                t.grouping_stack[stack_ix] = _GroupingStackData(
                    split+1,             # iix::Int64 # initial particle index
                    stack.fix,           # fix::Int64 # final particle index
                    (stack.dir+1)%3,     # dir::Int64 # direction (1==x, 2==y, 3==z)
                    midx,                # mid::Float64 # split coordinate
                    midy,                # mid::Float64 # split coordinate
                    midz,                # mid::Float64 # split coordinate
                    node_ix,             # pix::Int64 # parent node index
                    slx,                 # sl::Float64 # cell side length (nothing to do with actual particle inside!)
                    sly,                 # sl::Float64 # cell side length (nothing to do with actual particle inside!)
                    slz,                 # sl::Float64 # cell side length (nothing to do with actual particle inside!)
                )
            else
                # node has not enough particles to be splitted
                # tell the particles of their new parent!
                for i in (split+1):stack.fix
                    p = t.particles[i]
                    t.particles[i] = Particle(p.x,p.y,p.z,p.m,node_ix)
                end
            end
            # update parent node that it has a new child
            pn = t.nodes[stack.pix]
            t.nodes[stack.pix] = Node(
                pn.x,
                pn.y,
                pn.z,
                pn.m,
                pn.l,
                pn.pix,
                pn.iix,
                pn.pnum,
                pn.cix1,
                node_ix, ### <<<--- this is the update!!!
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

@inline function splitdir!(from_ix, to_ix, x, stack_at, dir)
    if dir==0
        splitx!(from_ix, to_ix, x, stack_at.midx)
    elseif dir==1
        splity!(from_ix, to_ix, x, stack_at.midy)
    elseif dir==2
        splitz!(from_ix, to_ix, x, stack_at.midz)
    else
        error("bad direction!")
    end
end













