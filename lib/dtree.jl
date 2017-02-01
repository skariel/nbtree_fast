immutable DTree
    trees::Vector{Tree}
    tree::Tree
end

function _rng(vec, i, N)
    rat = length(vec)/N
    left = round(Int64, ceil(rat*(i-1)+1))
    right = round(Int64, ceil(rat*i))
    left:right
end

function DTree(particles, S::Int64)
    tree = Tree(view(particles,1:length(particles)), S)
    group!(tree)
    trees = Tree[
        Tree(
            view(tree.particles, _rng(tree.particles, i,nthreads())),
            view(tree.nodes    , _rng(tree.nodes    , i,nthreads())),
            view(tree.exps     , _rng(tree.exps     , i,nthreads())),
            S,
        )
        for i in 1:nthreads()
    ]
    DTree(trees,tree)
end

function group!(t::DTree, glbl=false)
    if glbl
        group!(t.tree)
    end
    @threads for st in t.trees
        group!(st, t.tree)
    end
    nothing
end

function inform!(t::DTree)
    total_mass = sum(getm(p) for p in t.tree.particles)
    @threads for st in t.trees
        inform!(st)
        st.total_mass = total_mass
    end
    nothing    
end

function collect!(t::DTree)
    for i in eachindex(t.trees)
        collect!(t.trees[i])
    end
end

function accel!(t::DTree, ax,ay,az)
    @threads for i in eachindex(t.trees)
        vax = view(ax,_rng(ax,i,length(t.trees)))
        vay = view(ay,_rng(ay,i,length(t.trees)))
        vaz = view(az,_rng(az,i,length(t.trees)))
        accel!(t.trees[i], vax,vay,vaz)
    end
    nothing    
end

function collect!(t::DTree, ax,ay,az)
    collect!(t)
    accel!(t, ax,ay,az)
end

function interact!(t::DTree, alpha::Float64, ax,ay,az, eps2)
    @threads for i in eachindex(t.trees)
        t1 = t.trees[i]
        vax = view(ax,_rng(ax,i,length(t.trees)))
        vay = view(ay,_rng(ay,i,length(t.trees)))
        vaz = view(az,_rng(az,i,length(t.trees)))
        interact!(t1, alpha, vax,vay,vaz, eps2)
        for j in eachindex(t.trees)
            i==j && continue
            t2 = t.trees[j]
            interact!(t1, t2, alpha, vax,vay,vaz, eps2)
        end
    end
    nothing        
end



const I_BC = 3
function interact!(t1::Tree, t2::Tree, alpha::Float64, ax,ay,az, eps2)
    # stack1 -> 1st index (maybe particle)
    # stack2 -> 2nd index (maybe particle)
    # stack3 -> interaction type (CC, BC, CB) -- no CS type here!

    # pushing root into the stack_ix
    six = 1
    t1.stack1[six]=1
    t1.stack2[six]=1
    t1.stack3[six]=I_CC
    n1 = t1.nodes[1]
    n2 = t2.nodes[1]
    p1 = t1.particles[1]
    p2 = p1
    e = t1.exps[1]
    @fastmath @inbounds while six > 0
        ix1   = t1.stack1[six]
        ix2   = t1.stack2[six]
        itype = t1.stack3[six]
        six  -= 1

        if itype==I_BC
            n2 = t2.nodes[ix2]
            nbody2 = n2.fix-n2.iix+1
            if nbody2<16
                # just do direct summation
                p1 = t1.particles[ix1]
                @simd for i2 in n2.iix:n2.fix
                    p2 = t2.particles[i2]
                    dx = getx(p2) - getx(p1)
                    dy = gety(p2) - gety(p1)
                    dz = getz(p2) - getz(p1)
                    dr2 = dx*dx + dy*dy + dz*dz + eps2
                    dr3 = dr2*sqrt(dr2)
                    fac = getm(p2)/dr3
                    ax[ix1] += dx*fac
                    ay[ix1] += dy*fac
                    az[ix1] += dz*fac
                end            
                continue
            end
        end

        # Do a MAC test
        x1=0.0; y1=0.0; z1=0.0; m1=0.0; l1=0.0;
        x2=0.0; y2=0.0; z2=0.0; m2=0.0; l2=0.0;
        if itype==I_CB
            n1 = t1.nodes[ix1]
            nbody1 = n1.fix-n1.iix+1
            x1=n1.x
            y1=n1.y
            z1=n1.z
            m1=n1.m
            l1=n1.l
            p2 = t2.particles[ix2]
            x2=getx(p2)
            y2=gety(p2)
            z2=getz(p2)
            m2=getm(p2)
        elseif itype==I_BC
            p1 = t1.particles[ix1]
            x1=getx(p1)
            y1=gety(p1)
            z1=getz(p1)
            m1=getm(p1)

            x2=n2.x
            y2=n2.y
            z2=n2.z
            m2=n2.m
            l2=n2.l
        else
            # CC
            n1 = t1.nodes[ix1]
            nbody1 = n1.fix-n1.iix+1
            x1=n1.x
            y1=n1.y
            z1=n1.z
            m1=n1.m
            l1=n1.l
            n2 = t2.nodes[ix2]
            nbody2 = n2.fix-n2.iix+1
            x2=n2.x
            y2=n2.y
            z2=n2.z
            m2=n2.m
            l2=n2.l
        end             
        dx = x2-x1
        dy = y2-y1
        dz = z2-z1
        dx2 = dx*dx
        dy2 = dy*dy
        dz2 = dz*dz
        dr2 = dx2 + dy2 + dz2
        dr = sqrt(dr2)
        M = m1+m2
        @fastmath fac = (M/t1.total_mass)^0.15
        @fastmath if (l1 + l2)/dr < alpha/fac
            # MAC succesful, execute interaction
            dr2 += eps2
            dr = sqrt(dr2)
            dr3 = dr2*dr
            if itype!=I_BC
                # expansion needed
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

                e = t1.exps[ix1]
                t1.exps[ix1] = NodeExp(
                    e.px+m2*px,
                    e.pxx+m2*pxx,
                    e.pxxx+m2*pxxx,
                    e.pxxy+m2*pxxy,
                    e.pxxz+m2*pxxz,
                    e.pxy+m2*pxy,
                    e.pxyy+m2*pxyy,
                    e.pxyz+m2*pxyz,
                    e.pxz+m2*pxz,
                    e.pxzz+m2*pxzz,
                    e.py+m2*py,
                    e.pyy+m2*pyy,
                    e.pyyy+m2*pyyy,
                    e.pyyz+m2*pyyz,
                    e.pyz+m2*pyz,
                    e.pyzz+m2*pyzz,
                    e.pz+m2*pz,
                    e.pzz+m2*pzz,
                    e.pzzz+m2*pzzz,
                ) 
                continue
            end
            fac = m2/dr3
            ax[ix1] += dx*fac
            ay[ix1] += dy*fac
            az[ix1] += dz*fac                
            continue
        end

        # failed MAC

        if itype==I_CC
            # postconditions to direct summation
            if nbody1*nbody2 < 64
                @inbounds for i1 in n1.iix:n1.fix
                    p1 = t1.particles[i1]
                    @fastmath @inbounds @simd for i2 in n2.iix:n2.fix
                        p2 = t2.particles[i2]
                        dx = getx(p2) - getx(p1)
                        dy = gety(p2) - gety(p1)
                        dz = getz(p2) - getz(p1)
                        dr2 = dx*dx + dy*dy + dz*dz + eps2
                        dr3 = dr2*sqrt(dr2)
                        fac1 = getm(p2)/dr3
                        ax[i1] += dx*fac1
                        ay[i1] += dy*fac1
                        az[i1] += dz*fac1
                    end            
                end
                continue
            end
            # the interaction cannot be executed, splitting the bigger node
            if n1.l > n2.l
                if n1.cix1>0
                    six += 1
                    t1.stack1[six] = n1.cix1
                    t1.stack2[six] = ix2
                    t1.stack3[six] = I_CC
                end
                if n1.cix2>0
                    six += 1
                    t1.stack1[six] = n1.cix2
                    t1.stack2[six] = ix2
                    t1.stack3[six] = I_CC
                end
                if n1.cix1<0 && n1.cix2<0
                    # splitting into particles
                    @inbounds for i1 in n1.iix:n1.fix
                        six += 1
                        t1.stack1[six] = i1
                        t1.stack2[six] = ix2
                        t1.stack3[six] = I_BC
                    end
                end
            else
                if n2.cix1>0
                    six += 1
                    t1.stack1[six] = ix1
                    t1.stack2[six] = n2.cix1
                    t1.stack3[six] = I_CC
                end
                if n2.cix2>0
                    six += 1
                    t1.stack1[six] = ix1
                    t1.stack2[six] = n2.cix2
                    t1.stack3[six] = I_CC
                end
                if n2.cix1<0 && n2.cix2<0
                    # splitting into particles
                    @inbounds for i2 in n2.iix:n2.fix
                        six += 1
                        t1.stack1[six] = ix1
                        t1.stack2[six] = i2
                        t1.stack3[six] = I_CB
                    end
                end
            end
            continue
        end
                
        # we have a c-b (or b-c) interaction
        # test for postconditions for direct summation
        if itype==I_CB
            if nbody1<64 || (n1.cix1<0 && n1.cix2<0)
                p2 = t2.particles[ix2]
                @fastmath @inbounds @simd for i1 in n1.iix:n1.fix
                    p1 = t1.particles[i1]
                    dx = getx(p2) - getx(p1)
                    dy = gety(p2) - gety(p1)
                    dz = getz(p2) - getz(p1)
                    dr2 = dx*dx + dy*dy + dz*dz + eps2
                    dr3 = dr2*sqrt(dr2)
                    fac = getm(p2)/dr3
                    ax[i1] += dx*fac
                    ay[i1] += dy*fac
                    az[i1] += dz*fac
                end            
                continue
            end
            # cannot execute interaction, split cell
            if n1.cix1>0
                six += 1
                t1.stack1[six] = n1.cix1
                t1.stack2[six] = ix2
                t1.stack3[six] = I_CB
            end
            if n1.cix2>0
                six += 1
                t1.stack1[six] = n1.cix2
                t1.stack2[six] = ix2
                t1.stack3[six] = I_CB
            end
            continue
        end

        # interaction is BC
        if nbody2<64 || (n2.cix1<0 && n2.cix2<0)
            p1 = t1.particles[ix1]
            @fastmath @inbounds @simd for i2 in n2.iix:n2.fix
                p2 = t2.particles[i2]
                dx = getx(p2) - getx(p1)
                dy = gety(p2) - gety(p1)
                dz = getz(p2) - getz(p1)
                dr2 = dx*dx + dy*dy + dz*dz + eps2
                dr3 = dr2*sqrt(dr2)
                fac1 = getm(p2)/dr3
                ax[ix1] += dx*fac1
                ay[ix1] += dy*fac1
                az[ix1] += dz*fac1
            end            
            continue
        end
        # cannot execute interaction, split cell
        if n2.cix1>0
            six += 1
            t1.stack1[six] = ix1
            t1.stack2[six] = n2.cix1
            t1.stack3[six] = I_BC
        end
        if n2.cix2>0
            six += 1
            t1.stack1[six] = ix1
            t1.stack2[six] = n2.cix2
            t1.stack3[six] = I_BC
        end
        # thats it, loop for next interaction!
    end
end
