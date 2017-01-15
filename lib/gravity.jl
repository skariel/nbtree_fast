function inform!(t::Tree)
    @inbounds for i in t.num_nodes_used:-1:1
        x = 0.0
        y = 0.0
        z = 0.0
        m = 0.0
        l = 0.0

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
                dx = x-p.x
                dy = y-p.y
                dz = z-p.z
                dr2 = dx*dx+dy*dy+dz*dz
                if dr2>l
                    l=dr2
                end
            end
            l = sqrt(l)
        elseif n.cix1<0
            # single child, just transfer properties
            n2 = t.nodes[n.cix2]
            x = n2.x
            y = n2.y
            z = n2.z
            m = n2.m
            l = n2.l
        elseif n.cix2<0
            # single child, just transfer properties
            n1 = t.nodes[n.cix1]
            x = n1.x
            y = n1.y
            z = n1.z
            m = n1.m
            l = n1.l
        else
            # two children, merge properties
            n1 = t.nodes[n.cix1]
            n2 = t.nodes[n.cix2]
            m = n1.m + n2.m
            x = (n1.x*n1.m + n2.x*n2.m)/m
            y = (n1.y*n1.m + n2.y*n2.m)/m
            z = (n1.z*n1.m + n2.z*n2.m)/m

            dx1 = x-n1.x
            dy1 = y-n1.y
            dz1 = z-n1.z
            dr1 = sqrt(dx1*dx1+dy1*dy1+dz1*dz1)
            l1 = dr1+n1.l

            dx2 = x-n2.x
            dy2 = y-n2.y
            dz2 = z-n2.z
            dr2 = sqrt(dx2*dx2+dy2*dy2+dz2*dz2)
            l2 = dr2+n2.l

            lch = max(l1,l2)

            dx1 = x-n.maxx
            dy1 = y-n.maxy
            dz1 = z-n.maxz
            l1 = dx1*dx1+dy1*dy1+dz1*dz1            
            dx2 = x-n.maxx
            dy2 = y-n.maxy
            dz2 = z-n.maxz
            l2 = dx2*dx2+dy2*dy2+dz2*dz2
            dx3 = x-n.maxx
            dy3 = y-n.maxy
            dz3 = z-n.maxz
            l3 = dx3*dx3+dy3*dy3+dz3*dz3
            dx4 = x-n.maxx
            dy4 = y-n.maxy
            dz4 = z-n.maxz
            l4 = dx4*dx4+dy4*dy4+dz4*dz4
            dx5 = x-n.maxx
            dy5 = y-n.maxy
            dz5 = z-n.maxz
            l5 = dx5*dx5+dy5*dy5+dz5*dz5
            dx6 = x-n.maxx
            dy6 = y-n.maxy
            dz6 = z-n.maxz
            l6 = dx6*dx6+dy6*dy6+dz6*dz6
            dx7 = x-n.maxx
            dy7 = y-n.maxy
            dz7 = z-n.maxz
            l7 = dx7*dx7+dy7*dy7+dz7*dz7
            dx8 = x-n.maxx
            dy8 = y-n.maxy
            dz8 = z-n.maxz
            l8 = dx8*dx8+dy8*dy8+dz8*dz8

            lco = sqrt(max(l1,l2,l3,l4,l5,l6,l7,l8))
            l = min(lco,lch)
        end
        alpha = (m/t.total_mass)^(-0.1)
        t.nodes[i] = Node(
            x, # x::Float64
            y, # y::Float64
            z, # z::Float64
            m, # m::Float64
            l,
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
            0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, alpha,  
        )
    end
end

function interact!(t::Tree, alpha2::Float64)
    # pushing root into the stack
    six = 1
    t.stack1[six]=1
    t.stack2[six]=1
    while six > 0
        ix1 = t.stack1[six]
        ix2 = t.stack2[six]
        six -= 1
        if ix1==ix2
            # expand self interactions
            n = t.nodes[ix1]
            if n.cix1>0 && n.cix2>0
                six += 1
                t.stack1[six] = n.cix1
                t.stack2[six] = n.cix2
                six += 1
                t.stack1[six] = n.cix1
                t.stack2[six] = n.cix1
                six += 1
                t.stack1[six] = n.cix2
                t.stack2[six] = n.cix2
                continue
            end
            if n.cix1>0
                six += 1
                t.stack1[six] = n.cix1
                t.stack2[six] = n.cix1
                continue
            elseif n.cix2>0
                six += 1
                t.stack1[six] = n.cix2
                t.stack2[six] = n.cix2
                continue
            end
            continue
        end

        n1 = t.nodes[ix1]
        n2 = t.nodes[ix2]
        
        # check for direct summation, including leafs!
        if n1.cix1<0 && n1.cix2<0 && n2.cix1<0 && n2.cix2<0
            
        end

        dx = n2.x-n1.x
        dy = n2.y-n1.y
        dz = n2.z-n1.z
        dr2 = dx*dx+dy*dy+dz*dz
        if n.l*n.l/dr2 > alpha2*n.alpha*n.alpha
            # failed MAC, opening node

            if (n1.l>n2.l && (n1.cix1>0 || n1.cix2>0)) ||
                (n2.cix1<0 && n2.cix2<0)
                # opnening n1
                if n1.cix1>0
                    six += 1
                    t.stack1[six] = n1.cix1
                    t.stack2[six] = ix2
                end
                if n1.cix2>0
                    six += 1
                    t.stack1[six] = n1.cix2
                    t.stack2[six] = ix2
                end
                continue
            else
                # opening n2
                if n2.cix1>0
                    six += 1
                    t.stack1[six] = n2.cix1
                    t.stack2[six] = ix1
                end
                if n2.cix2>0
                    six += 1
                    t.stack1[six] = n2.cix2
                    t.stack2[six] = ix1
                end
                continue
            end
        end

        # nodes are well separated...

        dr = sqrt(dr2)
        dr3 = dr2*dr
        dr5 = dr3*dr2
        dr7 = dr3*dr5
        dx2 = dx*dx
        dy2 = dy*dy
        dz2 = dz*dz
        dx3 = dx2*dx
        dy3 = dy2*dy
        dz3 = dz2*dz

        px  = -dx/dr3
        py  = -dy/dr3
        pz  = -dz/dr3

        pxx = (3*dx2-dr2)/dr5
        pyy = (3*dy2-dr2)/dr5
        pzz = (3*dz2-dr2)/dr5

        pxxx = (-6*dx3 + 9*dx*(dr2-dx2))/dr7
        pyyy = (-6*dy3 + 9*dy*(dr2-dy2))/dr7
        pzzz = (-6*dz3 + 9*dz*(dr2-dz2))/dr7

        pxxy = 3*dy*(dr2-5*dx2)/dr7
        pxxz = 3*dz*(dr2-5*dx2)/dr7
        pxzz = 3*dx*(dr2-5*dz2)/dr7
        pyyz = 3*dz*(dr2-5*dy2)/dr7
        pyzz = 3*dy*(dr2-5*dz2)/dr7
        pxyy = 3*dx*(dr2-5*dy2)/dr7

        pxyz = 3*dy*(dr2-5*dx2)/dr7

        pxy = (3*dx*dy)/dr5
        pxz = (3*dx*dz)/dr5
        pyz = (3*dz*dy)/dr5
        
        t.nodes[ix1] = Node(
            n1.x,
            n1.y,
            n1.z,
            n1.m,
            n1.l,
            n1.maxx,
            n1.minx,
            n1.maxy,
            n1.miny,
            n1.maxz,
            n1.minz,
            n1.dir,
            n1.pix,
            n1.iix,
            n1.fix,
            n1.cix1,
            n1.cix2,
            -n2.m*px,
            -n2.m*pxx,
            -n2.m*pxxx,
            -n2.m*pxxy,
            -n2.m*pxxz,
            -n2.m*pxy,
            -n2.m*pxyy,
            -n2.m*pxyz,
            -n2.m*pxz,
            -n2.m*pxzz,
            -n2.m*py,
            -n2.m*pyy,
            -n2.m*pyyy,
            -n2.m*pyyz,
            -n2.m*pyz,
            -n2.m*pyzz,
            -n2.m*pz,
            -n2.m*pzz,
            -n2.m*pzzz,
            n1.alpha,
        )
        t.nodes[ix2] = Node(
            n2.x,
            n2.y,
            n2.z,
            n2.m,
            n2.l,
            n2.maxx,
            n2.minx,
            n2.maxy,
            n2.miny,
            n2.maxz,
            n2.minz,
            n2.dir,
            n2.pix,
            n2.iix,
            n2.fix,
            n2.cix1,
            n2.cix2,
            n1.m*px,
            n1.m*pxx,
            n1.m*pxxx,
            n1.m*pxxy,
            n1.m*pxxz,
            n1.m*pxy,
            n1.m*pxyy,
            n1.m*pxyz,
            n1.m*pxz,
            n1.m*pxzz,
            n1.m*py,
            n1.m*pyy,
            n1.m*pyyy,
            n1.m*pyyz,
            n1.m*pyz,
            n1.m*pyzz,
            n1.m*pz,
            n1.m*pzz,
            n1.m*pzzz,
            n2.alpha,
        )
    end
end

@inline function get_accel(t::Tree, pix::Int64, alpha2::Float64, eps2::Float64)
    stack_ix = 1
    t.stack1[stack_ix] = 1
    ax = 0.0
    ay = 0.0
    az = 0.0
    p = t.particles[pix]
    @inbounds while stack_ix > 0
        nix = t.stack1[stack_ix]
        stack_ix -= 1
        n = t.nodes[nix]
        dx = n.x - p.x
        dy = n.y - p.y
        dz = n.z - p.z
        dr2 = dx*dx + dy*dy + dz*dz + eps2
        if n.l*n.l/dr2 > alpha2*n.alpha*n.alpha
            # open criterion failed, we should try to open this node
            if n.cix1 > 0
                stack_ix += 1
                t.stack1[stack_ix] = n.cix1
            end
            if n.cix2 > 0
                stack_ix += 1
                t.stack1[stack_ix] = n.cix2
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
    t.stack1[stack_ix] = 1
    ax = 0.0
    ay = 0.0
    az = 0.0
    p = t.particles[pix]
    oa = sqrt(oax*oax+oay*oay+oaz*oaz)
    @inbounds while stack_ix > 0
        nix = t.stack1[stack_ix]
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
                t.stack1[stack_ix] = n.cix1
            end
            if n.cix2 > 0
                stack_ix += 1
                t.stack1[stack_ix] = n.cix2
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
