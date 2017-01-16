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
            0.0, 0.0, 0.0, 0.0,  
        )
    end
end

function add_expansion_to_n1(n1, n2)
    dx = n1.x-n2.x
    dy = n1.y-n2.y
    dz = n1.z-n2.z
    dx2 = dx*dx
    dy2 = dy*dy
    dz2 = dz*dz

    new_px = n1.px +
        n2.pxx  * dx    +
        n2.pxxx * dx2/2 +
        n2.pxxy * dx*dy +
        n2.pxxz * dx*dz +
        n2.pxy  * dy    +
        n2.pxyy * dy2/2 +
        n2.pxyz * dy*dz + 
        n2.pxz  * dz    +
        n2.pxzz * dz2/2
    new_pxx = n1.pxx    +
        n2.pxxx * dx    +
        n2.pxxy * dy    +
        n2.pxxz * dz    
    new_pxxx = n1.pxxx
    new_pxxy = n1.pxxy
    new_pxxz = n1.pxxz
    new_pxy = n1.pxy    +
        n2.pxxy * dx    +
        n2.pxyy * dy    +
        n2.pxyz * dz
    new_pxyy = n1.pxyy 
    new_pxyz = n1.pxyz
    new_pxz = n1.pxz +
        n2.pxxz * dx    +
        n2.pxyz * dy    + 
        n2.pxzz * dz
    new_pxzz = n1.pxzz
    new_py = n1.py      +
        n2.pxxy * dx2/2 +
        n2.pxy  * dx    +
        n2.pxyy * dx*dy +
        n2.pxyz * dx*dz +
        n2.pyy  * dy    + 
        n2.pyyy * dy2/2 +
        n2.pyyz * dy*dz +
        n2.pyz  * dz    +
        n2.pyzz * dz2/2
    new_pyy = n1.pyy    +
        n2.pxyy * dx    +
        n2.pyyy * dy    +
        n2.pyyz * dz 
    new_pyyy = n1.pyyy
    new_pyyz = n1.pyyz
    new_pyz = pyz    +
        n2.pxyz * dx    +
        n2.pyyz * dy    +
        n2.pyzz * dz
    new_pyzz = n1.pyzz
    new_pz = n1.pz      +
        n2.pxxz * dx2/2 +
        n2.pxyz * dx*dy +
        n2.pxz  * dx    +
        n2.pxzz * dx*dz +
        n2.pyyz * dy2/2 +
        n2.pyz  * dy    +
        n2.pyzz * dy*dz +
        n2.pzz  * dz    +
        n2.pzzz * dz2/2
    new_pzz = n1.pzz    +
        n2.pxzz * dx    +
        n2.pyzz * dy    +
        n2.pzzz * dz
    new_pzzz = n1.pzzz
    Node(
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

function get_accel_from_node(n, x,y,z)
    dx = x-n.x
    dy = y-n.y
    dz = z-n.z
    dx2 = dx*dx
    dy2 = dy*dy
    dz2 = dz*dz

    ax = n.px +
        n.pxx*dx +
        n.pxxx*dx2/2 +
        n.pxxy*dx*dy +
        n.pxxz*dx*dz +
        n.pxy*dy +
        n.pxyy*dy2/2 +
        n.pxyz*dy*dz +
        n.pxz*dz +
        n.pxzz*dz2/2


    ay = n.py +
        n.pxxy*dx2/2 +
        n.pxy*dx +
        n.pxyy*dx*dy +
        n.pxyz*dx*dz +
        n.pyy*dy +
        n.pyyy*dy2/2 +
        n.pyyz*dy*dz +
        n.pyz*dz +
        n.pyzz*dz2/2
    
    az = n.pz +
        n.pxxz*dx2/2 +
        n.pxyz*dx*dy +
        n.pxz*dx +
        n.pxzz*dx*dz +
        n.pyyz*dy2/2 +
        n.pyz*dy +
        n.pyzz*dy*dz +
        n.pzz*dz +
        n.pzzz*dz2/2

    ax, ay, az
end

const I_CC = 0
const I_CS = 1
const I_CB = 2
function interact!(t::Tree, alpha::Float64, ax,ay,az)
    # stack1 -> 1st index
    # stack2 -> 2nd index (maybe particle)
    # stack3 -> interaction type (CC, CS, CB)

    # pushing root into the stack_ix
    six = 1
    t.stack1[six]=1
    t.stack2[six]=1
    t.stack3[six]=I_CS # root is a self interaction

    @inbounds while six > 0
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
                for i1 in n.iix:n.fix
                    p1 = t.particles[i1]
                    for i2 in i1:n.fix
                        i1 == i2 && continue
                        p2 = t.particles[i2]
                        dx = p2.x - p1.x
                        dy = p2.y - p1.y
                        dz = p2.z - p1.z
                        dr2 = dx*dx + dy*dy + dz*dz
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
        if itype==I_CB && nbody1<512
            # just do self interactions
            p2 = t.particles[ix2]
            @inbounds for i1 in n1.iix:n1.fix
                p1 = t.particles[i1]
                dx = p2.x - p1.x
                dy = p2.y - p1.y
                dz = p2.z - p1.z
                dr2 = dx*dx + dy*dy + dz*dz
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
            p = t.particles[ix2]
            x=p.x
            y=p.y
            z=p.z
            m=p.m
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
        fac = (M/t.total_mass)^0.1
        if (n1.l + l)/dr < alpha/fac
            # MAC succesful, execute interaction
            dr3 = dr2*dr
            dr5 = dr3*dr2
            dr7 = dr3*dr5
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

            if itype==I_CB
                fac = -n.m/dr3
                ax[ix2] += dx/fac
                ay[ix2] += dy/fac
                az[ix2] += dz/fac                
            else
                n2 = t.nodes[ix2]
                t.nodes[ix2] = Node(
                    n2.x, n2.y, n2.z, n2.m, n2.l,
                    n2.maxx, n2.minx, n2.maxy, n2.miny, n2.maxz, n2.minz,
                    n2.dir, n2.pix, n2.iix, n2.fix,
                    n2.cix1, n2.cix2,
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
                )                
            end
            t.nodes[ix1] = Node(
                n1.x, n1.y, n1.z, n1.m, n1.l,
                n1.maxx, n1.minx, n1.maxy, n1.miny, n1.maxz, n1.minz,
                n1.dir, n1.pix, n1.iix, n1.fix, n1.cix1, n1.cix2,
                -m*px,
                -m*pxx,
                -m*pxxx,
                -m*pxxy,
                -m*pxxz,
                -m*pxy,
                -m*pxyy,
                -m*pxyz,
                -m*pxz,
                -m*pxzz,
                -m*py,
                -m*pyy,
                -m*pyyy,
                -m*pyyz,
                -m*pyz,
                -m*pyzz,
                -m*pz,
                -m*pzz,
                -m*pzzz,
            )
            continue
        end

        # failed MAC

        if itype==I_CC
            n2 = t.nodes[ix2]
            nbody2 = n2.fix-n2.iix+1
            # postconditions to direct summation
            if nbody1*nbody2 < 64
                @inbounds for i1 in n1.iix:n1.fix
                    p1 = t.particles[i1]
                    @inbounds for i2 in n2.iix:n2.fix
                        p2 = t.particles[i2]
                        dx = p2.x - p1.x
                        dy = p2.y - p1.y
                        dz = p2.z - p1.z
                        dr2 = dx*dx + dy*dy + dz*dz
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
            if n2.l > n1.l
                n1 = n2
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
            @inbounds for i1 in n1.iix:n1.fix
                p1 = t.particles[i1]
                dx = p2.x - p1.x
                dy = p2.y - p1.y
                dz = p2.z - p1.z
                dr2 = dx*dx + dy*dy + dz*dz
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

        # thats it, loop for next interaciton!
    end
end







