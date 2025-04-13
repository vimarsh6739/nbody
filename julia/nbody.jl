using Random

const SOFTENING = 1e-9
const POSMAX = 100

mutable struct Body
    x::Int
    y::Int
    z::Int
    vx::Float64
    vy::Float64
    vz::Float64
    m::Float64
    key::Int
    
    function Body(x::Int, y::Int, z::Int, vx::Float64, vy::Float64, vz::Float64, m::Float64)
        new(x, y, z, vx, vy, vz, m, 0)
    end
    
    function Body(x::Int, y::Int, z::Int, vx::Float64, vy::Float64, vz::Float64, m::Float64, key::Int)
        new(x, y, z, vx, vy, vz, m, key)
    end
end

mutable struct Node
    key::Int
    isLeaf::Bool
    parent::Union{Node, Nothing}
    children::Vector{Union{Node, Nothing}}
    whichChildren::UInt8
    
    function Node(key::Int, isLeaf::Bool)
        new(key, isLeaf, nothing, fill(nothing, 8), 0)
    end
end

mutable struct Octree
    root::Node
    leafLength::Int
    nLevels::Int
    
    function Octree(maxKeyLength::Int)
        root = Node(1, false)
        leafLength = maxKeyLength
        nLevels = div(maxKeyLength, 3)
        println("Octree with $(nLevels) levels")
        new(root, leafLength, nLevels)
    end
end

mutable struct LCGRandomEngine
    state::UInt32
    
    function LCGRandomEngine(seed::Int)
        new(UInt32(seed))
    end
end

function Base.rand(rng::LCGRandomEngine)
    a = UInt32(1664525)
    c = UInt32(1013904223)
    
    rng.state = a * rng.state + c
    
    return rng.state / (1.0 * typemax(UInt32))
end

function Base.rand(rng::LCGRandomEngine, range::UnitRange{Int})
    width = range.stop - range.start + 1
    
    r = rand(rng) * width
    
    return range.start + floor(Int, r)
end

function getKeyNoPrepend(body::Body)
    key = 0
    
    x, y, z = body.x, body.y, body.z
    
    i = 0
    while x != 0 || y != 0 || z != 0
        if (body.x & (1 << i)) != 0
            key |= (1 << (3*i))
        end
        if (body.y & (1 << i)) != 0
            key |= (1 << (3*i + 1))
        end
        if (body.z & (1 << i)) != 0
            key |= (1 << (3*i + 2))
        end
        
        x >>= 1
        y >>= 1
        z >>= 1
        i += 1
    end
    
    return key
end

function binaryLength(n::Int)
    if n == 0
        return 0
    end
    
    length = floor(Int, log2(n)) + 1
    if length % 3 != 0
        length += 3 - length % 3
    end
    
    @assert length % 3 == 0
    return length
end

function binaryString(k::Int)
    if k <= 1
        return string(k % 2)
    else
        return binaryString(div(k, 2)) * string(k % 2)
    end
end

function addChild!(octree::Octree, parent::Node, index::Int, isLeaf::Bool, mykey::Int)
    @assert parent.isLeaf == false
    @assert mykey > 0
    
    parent.children[index+1] = Node(mykey, isLeaf)
    parent.children[index+1].parent = parent
    parent.whichChildren |= 1 << index
end

function insert!(octree::Octree, body::Body)
    current = octree.root
    key = body.key
    
    @assert key >> (octree.leafLength - 1) == 1
    @assert octree.nLevels * 3 == octree.leafLength - 1
    
    level = 0
    while level < octree.nLevels - 1
        shifted = key >> (3 * (octree.nLevels - level - 1))
        index = shifted & 0x7
        
        if current.children[index+1] === nothing
            addChild!(octree, current, index, false, shifted)
        end
        
        current = current.children[index+1]
        level += 1
    end
    
    @assert level == octree.nLevels - 1
    index = key & 0x7
    
    @assert current.children[index+1] === nothing
    addChild!(octree, current, index, true, key)
end

function printTree(node::Union{Node, Nothing}, level::Int)
    if node === nothing
        return
    end
    
    print("  " ^ level)
    println("0b$(binaryString(node.key)) ($(node.key))")
    
    if !node.isLeaf
        for i in 0:7
            if node.children[i+1] !== nothing
                printTree(node.children[i+1], level + 1)
            end
        end
    end
end

function bodyForceRange!(bodies::Vector{Body}, dt::Float64, start::Int, finish::Int, n::Int)
    for i in start:finish
        Fx, Fy, Fz = 0.0, 0.0, 0.0
        
        for j in 1:n
            dx = Float64(bodies[j].x - bodies[i].x)
            dy = Float64(bodies[j].y - bodies[i].y)
            dz = Float64(bodies[j].z - bodies[i].z)
            distSqr = dx^2 + dy^2 + dz^2 + SOFTENING
            invDist = 1.0 / sqrt(distSqr)
            invDist3 = invDist^3
            
            Fx += dx * bodies[j].m * invDist3
            Fy += dy * bodies[j].m * invDist3
            Fz += dz * bodies[j].m * invDist3
        end
        
        bodies[i] = Body(
            bodies[i].x, bodies[i].y, bodies[i].z,
            bodies[i].vx + dt * Fx, bodies[i].vy + dt * Fy, bodies[i].vz + dt * Fz,
            bodies[i].m, bodies[i].key
        )
    end
end

function bodyForce!(bodies::Vector{Body}, dt::Float64, n::Int)
    mid = div(n, 2)
    bodyForceRange!(bodies, dt, 1, mid, n)
    bodyForceRange!(bodies, dt, mid+1, n, n)
end

function integratePositionsRange!(bodies::Vector{Body}, dt::Float64, start::Int, finish::Int)
    for i in start:finish
        fx = bodies[i].x + bodies[i].vx * dt
        fy = bodies[i].y + bodies[i].vy * dt
        fz = bodies[i].z + bodies[i].vz * dt
        
        newX = Int(round(fx))
        newY = Int(round(fy))
        newZ = Int(round(fz))
        
        newX = clamp(newX, 0, POSMAX-1)
        newY = clamp(newY, 0, POSMAX-1)
        newZ = clamp(newZ, 0, POSMAX-1)
        
        bodies[i] = Body(
            newX, newY, newZ,
            bodies[i].vx, bodies[i].vy, bodies[i].vz,
            bodies[i].m, bodies[i].key
        )
    end
end

function integratePositions!(bodies::Vector{Body}, dt::Float64, n::Int)
    mid = div(n, 2)
    integratePositionsRange!(bodies, dt, 1, mid)
    integratePositionsRange!(bodies, dt, mid+1, n)
end

function nbody_simulation(nBodies::Int = 2)
    println("NBody using Octree")
    
    rng = LCGRandomEngine(2025)
    
    dt = 0.01
    nIters = 10
    
    println("Simulating $nBodies bodies")
    
    bodies = Vector{Body}(undef, nBodies)
    
    for i in 1:nBodies
        bodies[i] = Body(0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0)
    end
    
    maxKeyLength = randomizeBodies(rng, bodies, nBodies)
    println("Randomized bodies with maxKeyLength = $maxKeyLength")
    
    octree = Octree(maxKeyLength)
    
    for i in 1:nBodies
        insert!(octree, bodies[i])
    end
    
    printTree(octree.root, 0)

    GC.gc()
    
    totalTime = 0.0
    
    for iter in 1:nIters
        startTime = time()
        
        bodyForce!(bodies, dt, nBodies)
        
        integratePositions!(bodies, dt, nBodies)
        
        newMaxKeyLength = 0
        for i in 1:nBodies
            body = bodies[i]
            key = getKeyNoPrepend(body)
            keyLength = binaryLength(key)
            
            if keyLength > newMaxKeyLength
                newMaxKeyLength = keyLength
            end
            
            bodies[i] = Body(
                body.x, body.y, body.z,
                body.vx, body.vy, body.vz,
                body.m, key
            )
        end
        
        prepend = 1 << newMaxKeyLength
        octree = Octree(newMaxKeyLength + 1)
        
        for i in 1:nBodies
            body = bodies[i]
            bodies[i] = Body(
                body.x, body.y, body.z,
                body.vx, body.vy, body.vz,
                body.m, body.key + prepend
            )
            
            insert!(octree, bodies[i])
        end
        
        endTime = time()
        tElapsed = endTime - startTime
        
        if iter > 1
            totalTime += tElapsed
        end
        
        println("Iteration $iter: $tElapsed seconds")
    end
    
    if nIters > 1
        avgTime = totalTime / (nIters - 1)
        println("Average time for iterations 2 through $nIters: $avgTime seconds.")
        println("$nBodies Bodies: average $(1e-9 * nBodies * nBodies / avgTime) Billion Interactions / second")
    end
end

function randomizeBodies(rng::LCGRandomEngine, bodies::Vector{Body}, n::Int)
    maxKeyLength = 0
    
    for i in 1:n
        x = rand(rng, 0:POSMAX-1)
        y = rand(rng, 0:POSMAX-1)
        z = rand(rng, 0:POSMAX-1)
        
        tempBody = Body(x, y, z, 0.0, 0.0, 0.0, 0.0)
        
        key = getKeyNoPrepend(tempBody)
        keyLength = binaryLength(key)
        
        if keyLength > maxKeyLength
            maxKeyLength = keyLength
        end
        
        vx = rand(rng) * 2 - 1
        vy = rand(rng) * 2 - 1
        vz = rand(rng) * 2 - 1
        m = rand(rng) + 0.1
        
        bodies[i] = Body(x, y, z, vx, vy, vz, m, key)
    end
    
    prepend = 1 << maxKeyLength
    println("Using prepend value: $prepend (2^$maxKeyLength)")
    
    for i in 1:n
        bodies[i] = Body(
            bodies[i].x, bodies[i].y, bodies[i].z,
            bodies[i].vx, bodies[i].vy, bodies[i].vz,
            bodies[i].m, bodies[i].key + prepend
        )
    end
    
    return maxKeyLength + 1
end

function main()
    nBodies = 2
    if length(ARGS) > 0
        nBodies = parse(Int, ARGS[1])
    end
    try
        nbody_simulation(nBodies)
    catch e
        println("Error in simulation: $e")
        println(catch_backtrace())
    end
end

main()