function inform!(t::Tree)
    @inbounds for i in t.num_nodes_used:-1:1
        x = 0.0
        y = 0.0
        z = 0.0
        m = 0.0

        qxx = 0.0
        qyy = 0.0
        qxy = 0.0
        qxz = 0.0
        qyz = 0.0        

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
            for j in n.iix:n.fix
                p = t.particles[j]
                dx= p.x-x; dy=p.y-y; dz=p.z-z;
                dx2 = dx*dx
                dy2 = dy*dy
                dr2=dx2+dy2+dz*dz;
                qxx += p.m*(3*dx2 - dr2)
                qyy += p.m*(3*dy2 - dr2)
                qxy += p.m*dx*dy
                qxz += p.m*dx*dz
                qyz += p.m*dy*dz                
            end
            qxy *= 3
            qxz *= 3
            qyz *= 3
            qzz = -qxx-qyy
        elseif n.cix1<0
            # single child, just transfer properties
            n2 = t.nodes[n.cix2]
            x = n2.x
            y = n2.y
            z = n2.z
            m = n2.m
            qxx = n2.qxx
            qyy = n2.qyy
            qxy = n2.qxy
            qxz = n2.qxz
            qyz = n2.qyz
        elseif n.cix2<0
            # single child, just transfer properties
            n1 = t.nodes[n.cix1]
            x = n1.x
            y = n1.y
            z = n1.z
            m = n1.m
            qxx = n1.qxx
            qyy = n1.qyy
            qxy = n1.qxy
            qxz = n1.qxz
            qyz = n1.qyz
        else
            # two children, merge properties
            n1 = t.nodes[n.cix1]
            n2 = t.nodes[n.cix2]
            m = n1.m + n2.m
            x = (n1.x*n1.m + n2.x*n2.m)/m
            y = (n1.y*n1.m + n2.y*n2.m)/m
            z = (n1.z*n1.m + n2.z*n2.m)/m

            dx= n1.x-x; dy=n1.y-y; dz=n1.z-z;
            dx2 = dx*dx
            dy2 = dy*dy
            dr2=dx2+dy2+dz*dz;
            qxx += n1.m*(3*dx2 - dr2)
            qyy += n1.m*(3*dy2 - dr2)
            qxy += n1.m*dx*dy
            qxz += n1.m*dx*dz
            qyz += n1.m*dy*dz                

            dx= n2.x-x; dy=n2.y-y; dz=n2.z-z;
            dx2 = dx*dx
            dy2 = dy*dy
            dr2=dx2+dy2+dz*dz;
            qxx += n2.m*(3*dx2 - dr2)
            qyy += n2.m*(3*dy2 - dr2)
            qxy += n2.m*dx*dy
            qxz += n2.m*dx*dz
            qyz += n2.m*dy*dz                

            qxy *= 3
            qxz *= 3
            qyz *= 3
            qzz = -qxx-qyy

            qxx += n1.qxx
            qyy += n1.qyy
            qxy += n1.qxy
            qxz += n1.qxz
            qyz += n1.qyz
            qxx += n2.qxx
            qyy += n2.qyy
            qxy += n2.qxy
            qxz += n2.qxz
            qyz += n2.qyz
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
            qxx,
            qyy,
            qxy,
            qxz,
            qyz,
        )
    end
end

@inline function get_accel(t::Tree, pix::Int64, alpha2::Float64, eps2::Float64)
    stack_ix = 1
    t.stack[stack_ix] = 1
    ax = 0.0
    ay = 0.0
    az = 0.0
    p = t.particles[pix]
    @inbounds while stack_ix > 0
        nix = t.stack[stack_ix]
        stack_ix -= 1
        n = t.nodes[nix]
        dx = n.x - p.x
        dy = n.y - p.y
        dz = n.z - p.z
        dx2 = dx*dx
        dy2 = dy*dy
        dz2 = dz*dz
        dr2 = dx2+dy2+dz2 + eps2
        if n.l*n.l/dr2 > alpha2
            # open criterion failed, we should try to open this node
            if n.cix1 > 0
                stack_ix += 1
                t.stack[stack_ix] = n.cix1
            end
            if n.cix2 > 0
                stack_ix += 1
                t.stack[stack_ix] = n.cix2
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
        dr = sqrt(dr2)
        dr3 = dr2*dr
        fac = n.m/dr3
        ax += dx*fac
        ay += dy*fac
        az += dz*fac     
        dr5 = dr3*dr2
        dr7 = dr5*dr2
        qzz = -n.qyy-n.qxx
        fac1 = 5.0*(0.5*dx2*n.qxx + dx*dy*n.qxy + dx*dz*n.qxz + 0.5*dy2*n.qyy + dy*dz*n.qyz + 0.5*dz2*qzz)/dr7
        ax += (-dx*n.qxx-dy*n.qxy-dz*n.qxz)/dr5 + dx*fac1
        ay += (-dx*n.qxy-dy*n.qyy-dz*n.qyz)/dr5 + dy*fac1
        az += (-dx*n.qxz-dy*n.qyz-dz*qzz)/dr5 + dz*fac1
    end
    ax, ay, az
end

@inline function get_accel_rel(t::Tree, pix::Int64, alpha2::Float64, eps2::Float64, oax,oay,oaz)
    stack_ix = 1
    t.stack[stack_ix] = 1
    ax = 0.0
    ay = 0.0
    az = 0.0
    p = t.particles[pix]
    oa = sqrt(oax*oax+oay*oay+oaz*oaz)
    @inbounds while stack_ix > 0
        nix = t.stack[stack_ix]
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
                t.stack[stack_ix] = n.cix1
            end
            if n.cix2 > 0
                stack_ix += 1
                t.stack[stack_ix] = n.cix2
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
    @inbounds for i in 1:length(t.particles)
        tax, tay, taz = get_accel(t, i, alpha2, eps2)
        ax[i] = tax
        ay[i] = tay
        az[i] = taz
    end
    nothing
end    

function get_all_accel_rel!(t::Tree, alpha2::Float64, eps2::Float64, ax::Vector{Float64}, ay::Vector{Float64}, az::Vector{Float64})
    @inbounds for i in 1:length(t.particles)
        tax, tay, taz = get_accel_rel(t, i, alpha2, eps2, ax[i],ay[i],az[i])
        ax[i] = tax
        ay[i] = tay
        az[i] = taz
    end
    nothing
end    
