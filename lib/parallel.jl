immutable LockVec
    locks::Vector{Atomic{Int64}}
    LockVec(N) = new(Atomic{Int64}[Atomic{Int64}(0) for i in 1:N])
end

@inline function try_lock!(l::LockVec, i::Int64)
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
    const rng = 1:nthreads()
    i = 1
    while !try_lock!(l, i)
        i = i%NTH+1
    end
    i
end



const NTH = nthreads()

immutable AKStack
    count::Atomic{Int64}
    six::Vector{Int64}
    locks::LockVec
    stack1::Array{Int64,2}
    stack2::Array{Int64,2}
    stack3::Array{Int64,2}
    function AKStack(N)
        sz = (nthreads(),N)
        new(Atomic{Int64}(0), zeros(Int64,nthreads()), LockVec(nthreads()), zeros(Int64,sz), zeros(Int64,sz), zeros(Int64,sz))
    end
end

@inline function push!(s::AKStack, a::Int64, b::Int64, c::Int64)
    i = lock_something!(s.locks)
    @inbounds six = s.six[i]+1
    @inbounds s.six[i] = six
    @inbounds s.stack1[i,six] = a
    @inbounds s.stack2[i,six] = b
    @inbounds s.stack3[i,six] = c
    atomic_add!(s.count, 1)
    unlock!(s.locks, i)
    nothing
end

@inline function try_pop!(s::AKStack)
    i = 0
    @inbounds while true
        i = i%NTH+1
        while !try_lock!(s.locks,i)
            i = i%NTH+1
        end
        if s.count[]==0
            unlock!(s.locks, i)
            return -1,-1,-1,false
        end
        six = s.six[i]
        if six==0
            unlock!(s.locks, i)
            continue
        end
        a = s.stack1[i,six]
        b = s.stack2[i,six]
        c = s.stack3[i,six]
        s.six[i] -= 1
        atomic_add!(s.count, -1)
        unlock!(s.locks, i)
        return a,b,c,true
    end
end

@inline function empty!(s::AKStack)
    @inbounds for i in 1:nthreads()
        s.six[i] = 0
    end
    s.count[] = 0
    nothing
end


