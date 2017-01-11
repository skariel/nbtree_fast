function get_acc(particles, ix, eps2)
    ax=0.0
    ay=0.0
    az=0.0
    @inbounds for p in particles
        dx = p.x-particles[ix].x
        dy = p.y-particles[ix].y
        dz = p.z-particles[ix].z
        dr2 = dx*dx+dy*dy+dz*dz+eps2
        dr3 = dr2*sqrt(dr2)
        ax += dx/dr3
        ay += dy/dr3
        az += dz/dr3
    end
    ax,ay,az
end

function perf(particles, ax, N, eps2)
    ixs = randperm(length(t.particles))[1:N]
    rax = zeros(N)
    vax = ax[ixs]
    @inbounds for i in 1:N
        tax,tay,taz = get_acc(particles, ixs[i], eps2)
        rax[i] = tax;
    end
    ee = abs((vax-rax)./rax*100)
    sort(ee)[round(Int64, N*0.99)]
end

# function test_particles_in_cell_mass_and_cm(t::Tree)
#     @inbounds for nix in 1:t.num_nodes_used
#         n = t.nodes[nix]
#         for pix in n.iix:(n.iix+n.pnum-1)
#             p = t.particles[pix]
#             @assert 
#         end
#     end
# end