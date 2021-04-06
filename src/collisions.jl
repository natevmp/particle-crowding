# Improved algorithm for finding and performing collisions

struct CollQ
    coll_t_Cid::Vector{SVector{2,Int}}
    _t::Vector{Float64}

    function CollQ()
        coll_t_Cid = SVector{2,Int}[]
        _t = Float64[]
        new(coll_t_Cid, _t)
    end
end


"""
For a collision queue sorted by times _t, find the index at which to insert a collision at time T.
"""
function findIndexInQueue(_t::AbstractVector, T::Real)
    if length(_t)==0
        return 1
    end
    for (i,tt) in enumerate(_t)
        if tt>T
            return i
        end
    end
    return length(_t)+1
end

function findIndexInQueue(_t::AbstractVector, T::Real, ind::Int)
    if length(_t)==0
        return 1
    end
    for (i,tt) in enumerate(_t[ind+1:end])
        if tt>T
            return ind + i
        end
    end
    return length(_t)+1
end

"""
Add a collision with cells l_cid at time t to the sorted queue Q.
"""
function addToCollQueue!(Q::CollQ, t::Float64, l_cid::AbstractVector{Int})
    addToCollQueue!(Q, t, SVector(l_cid...))
end

function addToCollQueue!(Q::CollQ, t::Float64, l_cid::SVector{2,Int})
    index = findIndexInQueue(Q._t, t)
    insert!(Q._t, index, t)
    insert!(Q.coll_t_Cid, index, l_cid)
    return Q
end

function addToCollQueue!(Q::CollQ, t::Float64, l_cid::AbstractVector{Int}, indIn::Int)
    addToCollQueue!(Q, t, SVector(l_cid...), indIn)
end

function addToCollQueue!(Q::CollQ, t::Float64, l_cid::SVector{2,Int}, indIn::Int)
    index = findIndexInQueue(Q._t, t, indIn)
    insert!(Q._t, index, t)
    insert!(Q.coll_t_Cid, index, l_cid)
    return Q
end


"""
Remove the collision from queue Q found at index ind.
"""
function removeFromCollQueue!(Q::CollQ, ind)
    deleteat!(Q.coll_t_Cid, ind)
    deleteat!(Q._t, ind)
    return Q
end

# struct Collision
#     cid::SVector{2, Int}
#     t::Float64
# end
# struct collQ
#     coll_t::Vector{Collision}
#     _t::Vector{Float64}
# end

function collisionCrossSection(arena::Arena)
    2arena.cellsList[1].radius
end


function collTimeCalc(cellA::Cell, cellB::Cell, bounds::Bounds)
    d = cellA.radius + cellB.radius
    if cellDistance(cellA, cellB, bounds) >= d
        return -1.
    end

    r = connectingVector(cellA.pos, cellB.pos, bounds)
    v = cellA.vel .- cellB.vel
    tColl = ( r⋅v + sqrt( d^2*(v⋅v) - (r[2]*v[1]-r[1]*v[2])^2 ) ) / (v⋅v)

    return tColl
end

function collCheck(cellA::Cell, cellB::Cell, bounds::Bounds)
    cellDistance(cellA, cellB, bounds) < (cellA.radius + cellB.radius) ? true : false
end


function getPairs(elems::Vector{Int})
    pairs_p = Vector{Vector{Int}}(undef,0)
    for i in 1:length(elems)
        for el2 in elems[i+1:end]
            push!(pairs_p, [elems[i], el2])
        end
    end
    return pairs_p
end

"""Perform collision of two cells."""
function collideCells!(bounds::Bounds, cellA::Cell, cellB::Cell, tRetc::Real; verbose=true)
    # reverse overlap movement
    moveCell!(cellA, -tRetc)
    moveCell!(cellB, -tRetc)
    #calculate new velocities
    newVelA = collVelCalc(cellA.vel, cellA.pos, cellB.vel, cellB.pos, bounds)
    newVelB = collVelCalc(cellB.vel, cellB.pos, cellA.vel, cellA.pos, bounds)
    cellA.vel .= newVelA
    cellB.vel .= newVelB
    #finish movement with remaining time
    moveCell!(cellA, tRetc)
    moveCell!(cellB, tRetc)
    return nothing
end

# function reverseCellTime!(t::Float64, cells_id::Vararg{Cell})
#     for cell in cells_id
#         moveCell!(cell, -t)
#     end
#     return nothing
# end

"""Calculate velocities of cell pair after collision."""
function collVelCalc(velA::MVector{2, Float64}, posA::MVector{2, Float64},
                        velB::MVector{2, Float64}, posB::MVector{2, Float64},
                            bounds::Bounds)
    posR = connectingVector(posA, posB, bounds)
    return velA .- (velA .- velB)⋅(posR)/norm(posR)^2 * posR
end

function removeFutureCollisions!(collQ::CollQ, qInd::Int, coll_cid)
    qIndF = qInd+1
    while qIndF <= length(collQ._t)
        if !isempty(intersect(coll_cid, collQ.coll_t_Cid[qIndF]))
            removeFromCollQueue!(collQ, qIndF)
            continue
        end
        qIndF += 1
    end
end

function collider!(arena::Arena, tStep::Real; verbose::Bool=false)

    # Initiation: create collision queue and list of collided cells
    collQ = CollQ()
    collidedCells_cid = Int[]

    # Create initial NNTree: treeZero
    pZero_dim_id = cellPositionsPeriodic_DIM_ID(arena)
    treeZero = BallTree(pZero_dim_id, PeriodicEuclidean([arena.bounds.xLen, arena.bounds.yLen]))
    overlapsTZ_nbh_Cid = inrange(treeZero, pZero_dim_id, collisionCrossSection(arena))

    # add collisions from treeZero to queue
    for nbh_cid in overlapsTZ_nbh_Cid
        if length(nbh_cid) < 2
            continue
        elseif length(nbh_cid)==2
            tRetc = collTimeCalc((@view arena.cellsList[nbh_cid])..., arena.bounds)
            addToCollQueue!(collQ, tStep-tRetc, nbh_cid)
        elseif length(nbh_cid)>2
            collPairs_nbh_Cid = getPairs(nbh_cid)
            for collPair_cid in collPairs_nbh_Cid
                tRetc = collTimeCalc((@view arena.cellsList[collPair_cid])..., arena.bounds)
                # ! ===== debug =====
                tRetc > tStep && println("Error: collision time larger than timestep ", tStep, " encountered in treeZero: ", tRetc)
                # ! =================
                tRetc > 0 && addToCollQueue!(collQ, tStep-tRetc, collPair_cid)
            end
        end
    end

    # == move through collisions in queue ==
    for (qInd, coll_cid) in enumerate(collQ.coll_t_Cid)
        
        # resolve collision
        collideCells!(arena.bounds, arena.cellsList[coll_cid[1]], arena.cellsList[coll_cid[2]], tStep - collQ._t[qInd])
        # add colliding cells to collidedCells_cid
        append!(collidedCells_cid, coll_cid)
        
        # === rm future collisions involving currently colliding cells from queue ===
        removeFutureCollisions!(collQ, qInd, coll_cid)

        # === add new collisions involving currently colliding cells to queue ===
        # find overlaps with treeZero -- disregarding cells in collidedCells_cid -- and add to queue
        posCur_dim_collCid = cellPositionsPeriodic_DIM_ID(arena, coll_cid)
        overlaps_pCur_Cid = inrange(treeZero, posCur_dim_collCid, collisionCrossSection(arena)) # get list of new overlaps per cell in current collision
        for (i, ovlNew_cid) in enumerate(overlaps_pCur_Cid) # perform this for both cells in coll_cid
            for cidNew in ovlNew_cid
                cidNew==coll_cid[i] && continue
                if isempty(intersect(cidNew, collidedCells_cid))
                    collNew_cid = @SVector[coll_cid[i], cidNew]
                    tRetc = collTimeCalc((@view arena.cellsList[collNew_cid])..., arena.bounds)
                    # ! ===== debug =====
                    tRetc > tStep && println("Error: collision time larger than timestep ", tStep, " encountered in treeMoved: ", tRetc)
                    tStep-tRetc < collQ._t[qInd] && println("Error: collision occurring in past detected")
                    # ! =================
                    addToCollQueue!(collQ, tStep-tRetc, collNew_cid, qInd)
                end
            end
        end
        # find overlaps with collided cells, calculate collision times, and add to correct position in queue
        for cidCur in coll_cid # do this for both cells in current collision
            # println(length(collidedCells_cid))
            distance_cid = cellDistance(arena.cellsList[cidCur], arena.cellsList[collidedCells_cid], arena.bounds)
            # println("distance_cid: ", distance_cid)
            overlaps_cid = findall(distance_cid .< collisionCrossSection(arena))
            # println("overlaps_cid", overlaps_cid)
            for cid in collidedCells_cid[overlaps_cid]
                cid==cidCur && continue
                tRetc = collTimeCalc(arena.cellsList[cidCur], arena.cellsList[cid], arena.bounds)
                # ! ===== debug =====
                tRetc > tStep && println("Error: collision time larger than timestep ", tStep, " encountered in treeMoved: ", tRetc)
                tStep-tRetc < collQ._t[qInd] && println("Error: collision occurring in past detected")
                # ! =================
                addToCollQueue!(collQ, tStep-tRetc, @SVector[cidCur, cid])
            end
        end


    end

    # ! ===== debug =====
    # println(length(collidedCells_cid), " collided in this timestep.")
    # ! =================

    return collidedCells_cid
end
