import Base:push!

immutable LockVec
    locks::Vector{Atomic{Int64}}
    LockVec(N) = new(Atomic{Int64}[Atomic{Int64}(0) for i in 1:N])
end

@inline function is_locked(l::LockVec, i::Int64)
    @inbounds p = l.locks[i][]==1
    p
end

@inline function try_lock!(l::LockVec, i::Int64)
    @inbounds if is_locked(l,i)
        return false
    end
    @inbounds p = atomic_xchg!(l.locks[i], 1)
    p == 0
end

@inline function lock!(l::LockVec, i::Int64)
    while !try_lock(l,i)
        # I'm busy!
    end
    nothing
end

@inline function unlock!(l::LockVec, i::Int64)
    @inbounds l.locks[i][] = 0
    nothing
end

@inline function lock_something!(l::LockVec)
    i = 1
    while !try_lock!(l, i)
        i = i%length(l.locks)+1
    end
    i
end



# immutable AKStack
#     six_i::Atomic{Int64}
#     six_f::Atomic{Int64}
#     count::Atomic{Int64}
#     stack1::Vector{Int64}
#     stack2::Vector{Int64}
#     stack3::Vector{Int64}
#     function AKStack(N)
#         new(Atomic{Int64}(0), Atomic{Int64}(0), Atomic{Int64}(0), zeros(Int64,N), zeros(Int64,N), zeros(Int64,N))
#     end
# end

# # TODO: batch push
# @inline function push!(s::AKStack, a::Int64, b::Int64, c::Int64)
#     s.count[] > 900000 && error("cannot push to stack!")
#     ix = atomic_add!(s.six_f,1) % length(s.stack1) + 1
#     @inbounds s.stack1[ix] = a
#     @inbounds s.stack2[ix] = b
#     @inbounds s.stack3[ix] = c
#     atomic_add!(s.count,1)
#     nothing
# end

# @inline function try_pop!(s::AKStack)
#     s.count[]==0 && return -1,-1,-1,false
#     atomic_add!(s.count,-1)
#     ix = atomic_add!(s.six_i,1) % length(s.stack1) + 1
#     a = s.stack1[ix]
#     b = s.stack2[ix]
#     c = s.stack3[ix]
#     return a,b,c,true
# end

# @inline function empty!(s::AKStack)
#     s.count[] = 0
#     s.six_i[] = 0
#     s.six_f[] = 0
#     nothing
# end




immutable AKStack
    six::Vector{Int64}
    locks::LockVec
    stack1::Array{Int64,2}
    stack2::Array{Int64,2}
    stack3::Array{Int64,2}
    function AKStack(N)
        sz = (nthreads(),N)
        new(zeros(Int64,nthreads()), LockVec(nthreads()), zeros(Int64,sz), zeros(Int64,sz), zeros(Int64,sz))
    end
end

# TODO: batch push
@inline function push!(s::AKStack, a::Int64, b::Int64, c::Int64)
    i = lock_something!(s.locks)
    @inbounds six = s.six[i]+1
    @inbounds s.six[i] = six
    @inbounds s.stack1[i,six] = a
    @inbounds s.stack2[i,six] = b
    @inbounds s.stack3[i,six] = c
    unlock!(s.locks, i)
    nothing
end

@inline function try_pop!(s::AKStack)
    i=0
    @inbounds while true
        i = i%nthreads()+1
        if !try_lock!(s.locks,i)
            continue
        end
        six = s.six[i]
        if six==0
            unlock!(s.locks,i)
            continue
        end
        a = s.stack1[i,six]
        b = s.stack2[i,six]
        c = s.stack3[i,six]
        s.six[i] -= 1
        unlock!(s.locks,i)
        return a,b,c,true
    end
    return -1,-1,-1,false
end

@inline function empty!(s::AKStack)
    s.six[:] = 0
    nothing
end


