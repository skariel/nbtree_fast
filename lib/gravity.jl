function inform!(t::Tree)
    @inbounds for i in t.num_nodes_used:-1:1
        x = 0.0
        y = 0.0
        z = 0.0
        m = 0.0

        n = t.nodes[i]
        if n.cix1<0 && n.cix2<0
            # leaf node, just use particles
            for j in n.iix:n.fix
                p = t.particles[j]
                x += p.x*p.m
                y += p.y*p.m
                z += p.z*p.m
                m += p.m
            end
            x /= m
            y /= m
            z /= m
        elseif n.cix1<0
            # single child, just transfer properties
            n2 = t.nodes[n.cix2]
            x = n2.x
            y = n2.y
            z = n2.z
            m = n2.m
        elseif n.cix2<0
            # single child, just transfer properties
            n1 = t.nodes[n.cix1]
            x = n1.x
            y = n1.y
            z = n1.z
            m = n1.m
        else
            # two children, merge properties
            n1 = t.nodes[n.cix1]
            n2 = t.nodes[n.cix2]
            m = n1.m + n2.m
            x = (n1.x*n1.m + n2.x*n2.m)/m
            y = (n1.y*n1.m + n2.y*n2.m)/m
            z = (n1.z*n1.m + n2.z*n2.m)/m
        end
        t.nodes[i] = Node(
            x, # x::Float64
            y, # y::Float64
            z, # z::Float64
            m, # m::Float64
            n.l,
            n.maxx,
            n.minx,
            n.maxy,
            n.miny,
            n.maxz,
            n.minz,
            n.dir,
            n.pix, # pix::Int64 # parent index
            n.iix, # iix::Int64 # first particle index
            n.fix, # 
            n.cix1, # cix1::Int64 # first child index
            n.cix2, # cix2::Int64 # second child index            
        )
    end
end

@inline function get_accel(t::Tree, pix::Int64, alpha2::Float64, eps2::Float64)
    stack_ix = 1
    const tid = threadid()
    t.stack[stack_ix, tid] = 1
    ax = 0.0
    ay = 0.0
    az = 0.0
    p = t.particles[pix]
    @inbounds while stack_ix > 0
        nix = t.stack[stack_ix, tid]
        stack_ix -= 1
        n = t.nodes[nix]
        dx = n.x - p.x
        dy = n.y - p.y
        dz = n.z - p.z
        dr2 = dx*dx + dy*dy + dz*dz + eps2
        if n.l*n.l/dr2 > alpha2
            # open criterion failed, we should try to open this node
            if n.cix1 > 0
                stack_ix += 1
                t.stack[stack_ix, tid] = n.cix1
            end
            if n.cix2 > 0
                stack_ix += 1
                t.stack[stack_ix, tid] = n.cix2
            end
            if n.cix1 < 0 && n.cix2 < 0
                # try direct summation
                for j in n.iix:n.fix
                    p2 = t.particles[j]
                    dx = p2.x - p.x
                    dy = p2.y - p.y
                    dz = p2.z - p.z
                    dr2 = dx*dx + dy*dy + dz*dz + eps2
                    fac = p2.m/dr2/sqrt(dr2)
                    ax += dx*fac
                    ay += dy*fac
                    az += dz*fac                    
                end
            end
            continue
        end
        # open criterion succeeded
        fac = n.m/dr2/sqrt(dr2)
        ax += dx*fac
        ay += dy*fac
        az += dz*fac
    end
    ax, ay, az
end

@inline function get_accel_rel(t::Tree, pix::Int64, alpha2::Float64, eps2::Float64, oax,oay,oaz)
    stack_ix = 1
    const tid = threadid()    
    t.stack[stack_ix, tid] = 1
    ax = 0.0
    ay = 0.0
    az = 0.0
    p = t.particles[pix]
    oa = sqrt(oax*oax+oay*oay+oaz*oaz)
    @inbounds while stack_ix > 0
        nix = t.stack[stack_ix, tid]
        stack_ix -= 1
        n = t.nodes[nix]
        dx = n.x - p.x
        dy = n.y - p.y
        dz = n.z - p.z
        dr2 = dx*dx + dy*dy + dz*dz + eps2
        na = n.m/dr2
        if na*n.l*n.l/dr2 > oa*alpha2
            # open criterion failed, we should try to open this node
            if n.cix1 > 0
                stack_ix += 1
                t.stack[stack_ix, tid] = n.cix1
            end
            if n.cix2 > 0
                stack_ix += 1
                t.stack[stack_ix, tid] = n.cix2
            end
            if n.cix1 < 0 && n.cix2 < 0
                # try direct summation
                for j in n.iix:n.fix
                    p2 = t.particles[j]
                    dx = p2.x - p.x
                    dy = p2.y - p.y
                    dz = p2.z - p.z
                    dr2 = dx*dx + dy*dy + dz*dz + eps2
                    fac = p2.m/dr2/sqrt(dr2)
                    ax += dx*fac
                    ay += dy*fac
                    az += dz*fac                    
                end
            end
            continue
        end
        # open criterion succeeded
        fac = n.m/dr2/sqrt(dr2)
        ax += dx*fac
        ay += dy*fac
        az += dz*fac
    end
    ax, ay, az
end

function get_all_accel!(t::Tree, alpha2::Float64, eps2::Float64, ax::Vector{Float64}, ay::Vector{Float64}, az::Vector{Float64})
    @threads for i in 1:length(t.particles)
        tax, tay, taz = get_accel(t, i, alpha2, eps2)
        ax[i] = tax
        ay[i] = tay
        az[i] = taz
    end
    nothing
end    

function get_all_accel_rel!(t::Tree, alpha2::Float64, eps2::Float64, ax::Vector{Float64}, ay::Vector{Float64}, az::Vector{Float64})
    @threads for i in 1:length(t.particles)
        tax, tay, taz = get_accel_rel(t, i, alpha2, eps2, ax[i],ay[i],az[i])
        ax[i] = tax
        ay[i] = tay
        az[i] = taz
    end
    nothing
end    
