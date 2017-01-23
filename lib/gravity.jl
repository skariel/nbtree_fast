function inform!(t::Tree)
    p = t.particles[1]
    n = t.nodes[1]
    n1=n
    n2=n
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
            dx2 = x-n.minx
            dy2 = y-n.maxy
            dz2 = z-n.maxz
            l2 = dx2*dx2+dy2*dy2+dz2*dz2
            dx3 = x-n.maxx
            dy3 = y-n.miny
            dz3 = z-n.maxz
            l3 = dx3*dx3+dy3*dy3+dz3*dz3
            dx4 = x-n.minx
            dy4 = y-n.miny
            dz4 = z-n.maxz
            l4 = dx4*dx4+dy4*dy4+dz4*dz4
            dx5 = x-n.maxx
            dy5 = y-n.maxy
            dz5 = z-n.minz
            l5 = dx5*dx5+dy5*dy5+dz5*dz5
            dx6 = x-n.minx
            dy6 = y-n.maxy
            dz6 = z-n.minz
            l6 = dx6*dx6+dy6*dy6+dz6*dz6
            dx7 = x-n.maxx
            dy7 = y-n.miny
            dz7 = z-n.minz
            l7 = dx7*dx7+dy7*dy7+dz7*dz7
            dx8 = x-n.minx
            dy8 = y-n.miny
            dz8 = z-n.minz
            l8 = dx8*dx8+dy8*dy8+dz8*dz8

            lco = sqrt(max(l1,l2,l3,l4,l5,l6,l7,l8))
            l = min(lco,lch)
        end

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
        )
    end
end

@inline function add_expansion_to_n1(n1, e1, n2, e2)
    dx = n1.x-n2.x
    dy = n1.y-n2.y
    dz = n1.z-n2.z
    dx2 = dx*dx
    dy2 = dy*dy
    dz2 = dz*dz

    new_px = e2.px+e1.px +
        e2.pxx  * dx    +
        e2.pxxx * dx2/2 +
        e2.pxxy * dx*dy +
        e2.pxxz * dx*dz +
        e2.pxy  * dy    +
        e2.pxyy * dy2/2 +
        e2.pxyz * dy*dz + 
        e2.pxz  * dz    +
        e2.pxzz * dz2/2
    new_pxx = e2.pxx+e1.pxx    +
        e2.pxxx * dx    +
        e2.pxxy * dy    +
        e2.pxxz * dz    
    new_pxxx = e2.pxxx+e1.pxxx
    new_pxxy = e2.pxxy+e1.pxxy
    new_pxxz = e2.pxxz+e1.pxxz
    new_pxy = e2.pxy+e1.pxy    +
        e2.pxxy * dx    +
        e2.pxyy * dy    +
        e2.pxyz * dz
    new_pxyy = e2.pxyy+e1.pxyy 
    new_pxyz = e2.pxyz+e1.pxyz
    new_pxz = e2.pxz+e1.pxz +
        e2.pxxz * dx    +
        e2.pxyz * dy    + 
        e2.pxzz * dz
    new_pxzz = e2.pxzz+e1.pxzz
    new_py = e2.py+e1.py      +
        e2.pxxy * dx2/2 +
        e2.pxy  * dx    +
        e2.pxyy * dx*dy +
        e2.pxyz * dx*dz +
        e2.pyy  * dy    + 
        e2.pyyy * dy2/2 +
        e2.pyyz * dy*dz +
        e2.pyz  * dz    +
        e2.pyzz * dz2/2
    new_pyy = e2.pyy+e1.pyy    +
        e2.pxyy * dx    +
        e2.pyyy * dy    +
        e2.pyyz * dz 
    new_pyyy = e2.pyyy+e1.pyyy
    new_pyyz = e2.pyyz+e1.pyyz
    new_pyz = e2.pyz+e1.pyz    +
        e2.pxyz * dx    +
        e2.pyyz * dy    +
        e2.pyzz * dz
    new_pyzz = e2.pyzz+e1.pyzz
    new_pz = e2.pz+e1.pz      +
        e2.pxxz * dx2/2 +
        e2.pxyz * dx*dy +
        e2.pxz  * dx    +
        e2.pxzz * dx*dz +
        e2.pyyz * dy2/2 +
        e2.pyz  * dy    +
        e2.pyzz * dy*dz +
        e2.pzz  * dz    +
        e2.pzzz * dz2/2
    new_pzz = e2.pzz+e1.pzz    +
        e2.pxzz * dx    +
        e2.pyzz * dy    +
        e2.pzzz * dz
    new_pzzz = e2.pzzz+e1.pzzz
    return NodeExp(
        new_px,
        new_pxx,
        new_pxxx,
        new_pxxy,
        new_pxxz,
        new_pxy,
        new_pxyy,
        new_pxyz,
        new_pxz,
        new_pxzz,
        new_py,
        new_pyy,
        new_pyyy,
        new_pyyz,
        new_pyz,
        new_pyzz,
        new_pz,
        new_pzz,
        new_pzzz,
    )
end

@inline function get_accel_from_node(n, e, x,y,z)
    dx = x-n.x
    dy = y-n.y
    dz = z-n.z
    dx2 = dx*dx
    dy2 = dy*dy
    dz2 = dz*dz

    ax = e.px +
        e.pxx*dx +
        e.pxxx*dx2/2 +
        e.pxxy*dx*dy +
        e.pxxz*dx*dz +
        e.pxy*dy +
        e.pxyy*dy2/2 +
        e.pxyz*dy*dz +
        e.pxz*dz +
        e.pxzz*dz2/2

    ay = e.py +
        e.pxxy*dx2/2 +
        e.pxy*dx +
        e.pxyy*dx*dy +
        e.pxyz*dx*dz +
        e.pyy*dy +
        e.pyyy*dy2/2 +
        e.pyyz*dy*dz +
        e.pyz*dz +
        e.pyzz*dz2/2
    
    az = e.pz +
        e.pxxz*dx2/2 +
        e.pxyz*dx*dy +
        e.pxz*dx +
        e.pxzz*dx*dz +
        e.pyyz*dy2/2 +
        e.pyz*dy +
        e.pyzz*dy*dz +
        e.pzz*dz +
        e.pzzz*dz2/2

    ax, ay, az
end

const I_CC = 0
const I_CS = 1
const I_CB = 2
function interact!(t::Tree, alpha::Float64, ax,ay,az, eps2)
    # stack1 -> 1st index
    # stack2 -> 2nd index (maybe particle)
    # stack3 -> interaction type (CC, CS, CB)

    # pushing root into the stack_ix
    six = 1
    t.stack1[six]=1
    t.stack2[six]=1
    t.stack3[six]=I_CS # root is a self interaction
    n = t.nodes[1]
    n1 = n
    p1 = Particle(0.0,0.0,0.0,0.0,-1)
    p2 = p1
    e = t.exps[1]
    @fastmath @inbounds while six > 0
        ix1 = t.stack1[six]
        ix2 = t.stack2[six]
        itype = t.stack3[six]
        six -= 1
        if itype==I_CS
            # we have a self interaction
            n = t.nodes[ix1]
            pnum = n.fix-n.iix+1
            if (pnum <= 32) || (n.cix1<0 && n.cix2<0)
                # its either a leaf or a small particle count
                # perform direct summation
                for i1 in n.iix:(n.fix-1)
                    p1 = t.particles[i1]
                    @fastmath @inbounds @simd for i2 in (i1+1):n.fix
                        p2 = t.particles[i2]
                        dx = p2.x - p1.x
                        dy = p2.y - p1.y
                        dz = p2.z - p1.z
                        dr2 = dx*dx + dy*dy + dz*dz + eps2
                        dr3 = dr2*sqrt(dr2)
                        fac1 = p2.m/dr3
                        fac2 = p1.m/dr3
                        ax[i1] += dx*fac1
                        ay[i1] += dy*fac1
                        az[i1] += dz*fac1
                        ax[i2] -= dx*fac2
                        ay[i2] -= dy*fac2
                        az[i2] -= dz*fac2
                    end            
                end                
                continue
            end
            # split self interactions
            if n.cix1>0 && n.cix2>0
                six += 1
                t.stack1[six] = n.cix1
                t.stack2[six] = n.cix2
                t.stack3[six] = I_CC
            end
            if n.cix1>0
                six += 1
                t.stack1[six] = n.cix1
                t.stack2[six] = n.cix1
                t.stack3[six] = I_CS
            end
            if n.cix2>0
                six += 1
                t.stack1[six] = n.cix2
                t.stack2[six] = n.cix2
                t.stack3[six] = I_CS
            end
            continue
        end

        n1 = t.nodes[ix1]
        nbody1 = n1.fix-n1.iix+1
        if itype==I_CB && nbody1<16
            # just do direct summation
            p2 = t.particles[ix2]
            @fastmath @inbounds @simd for i1 in n1.iix:n1.fix
                p1 = t.particles[i1]
                dx = p2.x - p1.x
                dy = p2.y - p1.y
                dz = p2.z - p1.z
                dr2 = dx*dx + dy*dy + dz*dz + eps2
                dr3 = dr2*sqrt(dr2)
                fac1 = p2.m/dr3
                fac2 = p1.m/dr3
                ax[i1] += dx*fac1
                ay[i1] += dy*fac1
                az[i1] += dz*fac1
                ax[ix2] -= dx*fac2
                ay[ix2] -= dy*fac2
                az[ix2] -= dz*fac2
            end            
            continue
        end

        # Do a MAC test
        x=0.0; y=0.0; z=0.0; m=0.0; l=0.0;
        if itype==I_CB
            p2 = t.particles[ix2]
            x=p2.x
            y=p2.y
            z=p2.z
            m=p2.m
        else
            n = t.nodes[ix2]
            x=n.x
            y=n.y
            z=n.z
            m=n.m
            l=n.l
        end                    
        dx = x-n1.x
        dy = y-n1.y
        dz = z-n1.z
        dx2 = dx*dx
        dy2 = dy*dy
        dz2 = dz*dz
        dr2 = dx2 + dy2 + dz2
        dr = sqrt(dr2)
        M = n1.m+m
        @fastmath fac = (M/t.total_mass)^0.15
        @fastmath if (n1.l + l)/dr < alpha/fac
            # MAC succesful, execute interaction
            dr2 += eps2
            dr = sqrt(dr2)

            dr3 = dr2*dr
            dr5 = dr3*dr2
            dr7 = dr2*dr5
            dr73 = dr7/3
            dx3 = dx2*dx
            dy3 = dy2*dy
            dz3 = dz2*dz

            px = dx/dr3
            py = dy/dr3
            pz = dz/dr3

            pxx = (3*dx2-dr2)/dr5
            pyy = (3*dy2-dr2)/dr5
            pzz = (3*dz2-dr2)/dr5

            pxy = 3*dx*dy/dr5
            pxz = 3*dx*dz/dr5
            pyz = 3*dy*dz/dr5

            pxxx = dx*(5*dx2-3*dr2)/dr73
            pyyy = dy*(5*dy2-3*dr2)/dr73
            pzzz = dz*(5*dz2-3*dr2)/dr73

            pxxy = dy*(5*dx2-dr2)/dr73
            pxxz = dz*(5*dx2-dr2)/dr73
            pyyz = dz*(5*dy2-dr2)/dr73
            pyzz = dy*(5*dz2-dr2)/dr73
            pxyy = dx*(5*dy2-dr2)/dr73
            pxzz = dx*(5*dz2-dr2)/dr73

            pxyz = 15*dx*dy*dz/dr7

            if itype==I_CB
                fac = -n1.m/dr3
                ax[ix2] += dx*fac
                ay[ix2] += dy*fac
                az[ix2] += dz*fac                
            else
                e = t.exps[ix2]
                t.exps[ix2] = NodeExp(
                    e.px-n1.m*px,
                    e.pxx+n1.m*pxx,
                    e.pxxx-n1.m*pxxx,
                    e.pxxy-n1.m*pxxy,
                    e.pxxz-n1.m*pxxz,
                    e.pxy+n1.m*pxy,
                    e.pxyy-n1.m*pxyy,
                    e.pxyz-n1.m*pxyz,
                    e.pxz+n1.m*pxz,
                    e.pxzz-n1.m*pxzz,
                    e.py-n1.m*py,
                    e.pyy+n1.m*pyy,
                    e.pyyy-n1.m*pyyy,
                    e.pyyz-n1.m*pyyz,
                    e.pyz+n1.m*pyz,
                    e.pyzz-n1.m*pyzz,
                    e.pz-n1.m*pz,
                    e.pzz+n1.m*pzz,
                    e.pzzz-n1.m*pzzz,
                )
            end
            e = t.exps[ix1]
            t.exps[ix1] = NodeExp(
                e.px+m*px,
                e.pxx+m*pxx,
                e.pxxx+m*pxxx,
                e.pxxy+m*pxxy,
                e.pxxz+m*pxxz,
                e.pxy+m*pxy,
                e.pxyy+m*pxyy,
                e.pxyz+m*pxyz,
                e.pxz+m*pxz,
                e.pxzz+m*pxzz,
                e.py+m*py,
                e.pyy+m*pyy,
                e.pyyy+m*pyyy,
                e.pyyz+m*pyyz,
                e.pyz+m*pyz,
                e.pyzz+m*pyzz,
                e.pz+m*pz,
                e.pzz+m*pzz,
                e.pzzz+m*pzzz,
            ) 
            continue
        end

        # failed MAC

        if itype==I_CC
            nbody2 = n.fix-n.iix+1
            # postconditions to direct summation
            if nbody1*nbody2 < 64
                @inbounds for i1 in n1.iix:n1.fix
                    p1 = t.particles[i1]
                    @fastmath @inbounds @simd for i2 in n.iix:n.fix
                        p2 = t.particles[i2]
                        dx = p2.x - p1.x
                        dy = p2.y - p1.y
                        dz = p2.z - p1.z
                        dr2 = dx*dx + dy*dy + dz*dz + eps2
                        dr3 = dr2*sqrt(dr2)
                        fac1 = p2.m/dr3
                        fac2 = p1.m/dr3
                        ax[i1] += dx*fac1
                        ay[i1] += dy*fac1
                        az[i1] += dz*fac1
                        ax[i2] -= dx*fac2
                        ay[i2] -= dy*fac2
                        az[i2] -= dz*fac2
                    end            
                end
                continue
            end
            # the interaction cannot be executed, splitting the bigger node
            if n.l > n1.l
                n1 = n
                ix1,ix2 = ix2,ix1
            end
            if n1.cix1>0
                six += 1
                t.stack1[six] = n1.cix1
                t.stack2[six] = ix2
                t.stack3[six] = I_CC
            end
            if n1.cix2>0
                six += 1
                t.stack1[six] = n1.cix2
                t.stack2[six] = ix2
                t.stack3[six] = I_CC
            end
            if n1.cix1<0 && n1.cix2<0
                # splitting into particles
                @inbounds for i1 in n1.iix:n1.fix
                    six += 1
                    t.stack1[six] = ix2
                    t.stack2[six] = i1
                    t.stack3[six] = I_CB
                end
            end
            continue
        end
        
        # we have a c-b interaction
        # test for postconditions for direct summation
        if nbody1<64 || (n1.cix1<0 && n1.cix2<0)
            p2 = t.particles[ix2]
            @fastmath @inbounds @simd for i1 in n1.iix:n1.fix
                p1 = t.particles[i1]
                dx = p2.x - p1.x
                dy = p2.y - p1.y
                dz = p2.z - p1.z
                dr2 = dx*dx + dy*dy + dz*dz + eps2
                dr3 = dr2*sqrt(dr2)
                fac1 = p2.m/dr3
                fac2 = p1.m/dr3
                ax[i1] += dx*fac1
                ay[i1] += dy*fac1
                az[i1] += dz*fac1
                ax[ix2] -= dx*fac2
                ay[ix2] -= dy*fac2
                az[ix2] -= dz*fac2
            end            
            continue
        end

        # cannot execute interaction, split cell
        if n1.cix1>0
            six += 1
            t.stack1[six] = n1.cix1
            t.stack2[six] = ix2
            t.stack3[six] = I_CB
        end
        if n1.cix2>0
            six += 1
            t.stack1[six] = n1.cix2
            t.stack2[six] = ix2
            t.stack3[six] = I_CB
        end

        # thats it, loop for next interaction!
    end
end

function collect!(t::Tree, ax::Vector{Float64},ay::Vector{Float64},az::Vector{Float64})
    n = t.nodes[1]
    _n = n
    p = t.particles[1]
    e = t.exps[1]
    _e = e
    @fastmath @inbounds for i in 1:t.num_nodes_used
        e = t.exps[i]
        n = t.nodes[i]
        if n.cix1 > 0
            _e = t.exps[n.cix1]
            _n = t.nodes[n.cix1]
            t.exps[n.cix1] = add_expansion_to_n1(_n,_e, n,e)
        end
        if n.cix2 > 0
            _e = t.exps[n.cix2]
            _n = t.nodes[n.cix2]
            t.exps[n.cix2] = add_expansion_to_n1(_n,_e, n,e)
        end
        if n.cix1<0 && n.cix2<0
            for i in n.iix:n.fix
                p = t.particles[i]
                dax,day,daz = get_accel_from_node(n, e, p.x,p.y,p.z)
                ax[i] += dax
                ay[i] += day
                az[i] += daz
            end
        end
    end
end


