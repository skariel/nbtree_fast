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

        maxx = -1.0e30 # infinity, ha!
        maxy = -1.0e30 # infinity, ha!
        maxz = -1.0e30 # infinity, ha!
        minx =  1.0e30 # infinity, ha!
        miny =  1.0e30 # infinity, ha!
        minz =  1.0e30 # infinity, ha!

        n = t.nodes[i]
        if n.cix1<0 && n.cix2<0
            # leaf node, just use particles
            for j in n.iix:n.fix
                p = t.particles[j]
                x += getx(p)*getm(p)
                y += gety(p)*getm(p)
                z += getz(p)*getm(p)
                m += getm(p)

                if getx(p) < minx
                    minx = getx(p)
                end
                if getx(p) > maxx
                    maxx = getx(p)
                end
                if gety(p) < miny
                    miny = gety(p)
                end
                if gety(p) > maxy
                    maxy = gety(p)
                end
                if getz(p) < minz
                    minz = getz(p)
                end
                if getz(p) > maxz
                    maxz = getz(p)
                end
            end
            x /= m
            y /= m
            z /= m
            for j in n.iix:n.fix
                p = t.particles[j]
                dx = x-getx(p)
                dy = y-gety(p)
                dz = z-getz(p)
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
            maxx=n2.maxx
            maxy=n2.maxy
            maxz=n2.maxz
            minx=n2.minx
            miny=n2.miny
            minz=n2.minz
        elseif n.cix2<0
            # single child, just transfer properties
            n1 = t.nodes[n.cix1]
            x = n1.x
            y = n1.y
            z = n1.z
            m = n1.m
            l = n1.l
            maxx=n1.maxx
            maxy=n1.maxy
            maxz=n1.maxz
            minx=n1.minx
            miny=n1.miny
            minz=n1.minz
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

            maxx = max(n1.maxx, n2.maxx)
            maxy = max(n1.maxy, n2.maxy)
            maxz = max(n1.maxz, n2.maxz)
            minx = min(n1.minx, n2.minx)
            miny = min(n1.miny, n2.miny)
            minz = min(n1.minz, n2.minz)            

            dx1 = x-maxx
            dy1 = y-maxy
            dz1 = z-maxz
            l1 = dx1*dx1+dy1*dy1+dz1*dz1            
            dx2 = x-minx
            dy2 = y-maxy
            dz2 = z-maxz
            l2 = dx2*dx2+dy2*dy2+dz2*dz2
            dx3 = x-maxx
            dy3 = y-miny
            dz3 = z-maxz
            l3 = dx3*dx3+dy3*dy3+dz3*dz3
            dx4 = x-minx
            dy4 = y-miny
            dz4 = z-maxz
            l4 = dx4*dx4+dy4*dy4+dz4*dz4
            dx5 = x-maxx
            dy5 = y-maxy
            dz5 = z-minz
            l5 = dx5*dx5+dy5*dy5+dz5*dz5
            dx6 = x-minx
            dy6 = y-maxy
            dz6 = z-minz
            l6 = dx6*dx6+dy6*dy6+dz6*dz6
            dx7 = x-maxx
            dy7 = y-miny
            dz7 = z-minz
            l7 = dx7*dx7+dy7*dy7+dz7*dz7
            dx8 = x-minx
            dy8 = y-miny
            dz8 = z-minz
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
            maxx,
            minx,
            maxy,
            miny,
            maxz,
            minz,
            n.dir,
            n.pix, # pix::Int64 # parent index
            n.iix, # iix::Int64 # first particle index
            n.fix, # 
            n.cix1, # cix1::Int64 # first child index
            n.cix2, # cix2::Int64 # second child index
        )
    end
    t.total_mass = t.nodes[1].m
end

@inline function add_expansion_to_n1(n1, e1, n2, e2)
    dx = n1.x-n2.x
    dy = n1.y-n2.y
    dz = n1.z-n2.z

    new_px = e2.px+e1.px +
        e2.pxx  * dx    +
        e2.pxy  * dy    +
        e2.pxz  * dz    
    new_pxx = e2.pxx+e1.pxx    
    new_pxy = e2.pxy+e1.pxy    
    new_pxz = e2.pxz+e1.pxz 
    new_py = e2.py+e1.py      +
        e2.pxy  * dx    +
        e2.pyy  * dy    + 
        e2.pyz  * dz    
    new_pyy = e2.pyy+e1.pyy    
    new_pyz = e2.pyz+e1.pyz    
    new_pz = e2.pz+e1.pz      +
        e2.pxz  * dx    +
        e2.pyz  * dy    +
        e2.pzz  * dz    
    new_pzz = e2.pzz+e1.pzz    
    return NodeExp(
        new_px,
        new_pxx,
        new_pxy,
        new_pxz,
        new_py,
        new_pyy,
        new_pyz,
        new_pz,
        new_pzz,
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
        #e.pxxx*dx2/2 +
        #e.pxxy*dx*dy +
        #e.pxxz*dx*dz +
        e.pxy*dy +
        #e.pxyy*dy2/2 +
        #e.pxyz*dy*dz +
        e.pxz*dz #+
        #e.pxzz*dz2/2

    ay = e.py +
        #e.pxxy*dx2/2 +
        e.pxy*dx +
        #e.pxyy*dx*dy +
        #e.pxyz*dx*dz +
        e.pyy*dy +
        #e.pyyy*dy2/2 +
        #e.pyyz*dy*dz +
        e.pyz*dz #+
        #e.pyzz*dz2/2
    
    az = e.pz +
        #e.pxxz*dx2/2 +
        #e.pxyz*dx*dy +
        e.pxz*dx +
        #e.pxzz*dx*dz +
        #e.pyyz*dy2/2 +
        e.pyz*dy +
        #e.pyzz*dy*dz +
        e.pzz*dz #+
        #e.pzzz*dz2/2

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
    p1 = t.particles[1]
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
                        dx = getx(p2) - getx(p1)
                        dy = gety(p2) - gety(p1)
                        dz = getz(p2) - getz(p1)
                        dr2 = dx*dx + dy*dy + dz*dz + eps2
                        dr3 = dr2*sqrt(dr2)
                        fac1 = getm(p2)/dr3
                        fac2 = getm(p1)/dr3
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
                dx = getx(p2) - getx(p1)
                dy = gety(p2) - gety(p1)
                dz = getz(p2) - getz(p1)
                dr2 = dx*dx + dy*dy + dz*dz + eps2
                dr3 = dr2*sqrt(dr2)
                fac1 = getm(p2)/dr3
                fac2 = getm(p1)/dr3
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
            x=getx(p2)
            y=gety(p2)
            z=getz(p2)
            m=getm(p2)
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
        @fastmath fac = (M/t.total_mass)^MASS_POWER_FAC
        @fastmath if (n1.l + l)/dr < alpha/fac
            # MAC succesful, execute interaction
            dr2 += eps2
            dr = sqrt(dr2)

            dr3 = dr2*dr
            dr5 = dr3*dr2

            px = dx/dr3
            py = dy/dr3
            pz = dz/dr3

            pxx = (3*dx2-dr2)/dr5
            pyy = (3*dy2-dr2)/dr5
            pzz = (3*dz2-dr2)/dr5

            pxy = 3*dx*dy/dr5
            pxz = 3*dx*dz/dr5
            pyz = 3*dy*dz/dr5

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
                    e.pxy+n1.m*pxy,
                    e.pxz+n1.m*pxz,
                    e.py-n1.m*py,
                    e.pyy+n1.m*pyy,
                    e.pyz+n1.m*pyz,
                    e.pz-n1.m*pz,
                    e.pzz+n1.m*pzz,
                )
            end
            e = t.exps[ix1]
            t.exps[ix1] = NodeExp(
                e.px+m*px,
                e.pxx+m*pxx,
                e.pxy+m*pxy,
                e.pxz+m*pxz,
                e.py+m*py,
                e.pyy+m*pyy,
                e.pyz+m*pyz,
                e.pz+m*pz,
                e.pzz+m*pzz,
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
                        dx = getx(p2) - getx(p1)
                        dy = gety(p2) - gety(p1)
                        dz = getz(p2) - getz(p1)
                        dr2 = dx*dx + dy*dy + dz*dz + eps2
                        dr3 = dr2*sqrt(dr2)
                        fac1 = getm(p2)/dr3
                        fac2 = getm(p1)/dr3
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
                dx = getx(p2) - getx(p1)
                dy = gety(p2) - gety(p1)
                dz = getz(p2) - getz(p1)
                dr2 = dx*dx + dy*dy + dz*dz + eps2
                dr3 = dr2*sqrt(dr2)
                fac1 = getm(p2)/dr3
                fac2 = getm(p1)/dr3
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

function collect!(t::Tree, ax,ay,az)
    collect!(t)
    accel!(t, ax,ay,az)
end

function collect!(t::Tree)
    n = t.nodes[1]
    _n = n
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
    end
end

function accel!(t::Tree, ax,ay,az)
    n = t.nodes[1]
    p = t.particles[1]
    e = t.exps[1]
    old_pix = -1
    @fastmath @inbounds for i in eachindex(t.particles)
        p = t.particles[i]
        if old_pix != getpix(p)
            old_pix = getpix(p)
            e = t.exps[getpix(p)]
            n = t.nodes[getpix(p)]
        end
        dax,day,daz = get_accel_from_node(n, e, getx(p),gety(p),getz(p))
        ax[i] += dax
        ay[i] += day
        az[i] += daz
    end

end
