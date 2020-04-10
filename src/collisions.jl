
### ======== Collision Functions ========

"""Perform collisions in a given arena. Return index of collided cells."""
function collider!(arena::Arena; nbhcut=50)
    collidedCellsList = Int[]

    counter = 1
    while true
        collidedCellsPass = passThroughNeighborhood(arena)
        if length(collidedCellsPass)==0
            break
        end
        union!(collidedCellsList, collidedCellsPass)
        counter += 1
        if counter > nbhcut
            println(">>> Error: Could not resolve neighborhood collisions in ", nbhcut, " steps ! <<<")
            break
        end
    end
    return collidedCellsList
end

"""Find collisions and perform collisions in arena. """
function passThroughNeighborhood(arena::Arena)
    # create neighborhoods
    nbhoods_id = collisionFinder(arena)
    collidedCellsOut = Int[]
    for nbh_c in nbhoods_id
        # check for collision in neighborhood
        collision = collCheck!(nbh_c, collidedCellsOut)
        if collision === false
            continue
        end
        # compose neighborhoods with shared elements
        nbhComposer!(nbh_c, nbhoods_id)
        if length(nbh_c) == 2
            # get collision time
            tRetc = collTimeCalc(arena.cellsList[nbh_c]..., arena.bounds)
        else
            # get earliest collision an time in composite neighborhood
            tRetc, nbh_c = collTimeCompare(arena, nbh_c...)
        end
        #perform collision
        collideCells!(arena.bounds, arena.cellsList[nbh_c]..., tRetc)
        push!(collidedCellsOut, nbh_c...)
    end
    return collidedCellsOut
end

"""Find earliest collision in higher order (3 or more cells) collisions"""
function collTimeCompare(arena::Arena, nbh_c::Vararg{Int})
    cellPairsInd_p = getPairs(nbh_c...)
    timepair_tp = [(collTimeCalc(arenaIndsToCells(arena, cpi)..., arena.bounds), cpi)
                        for cpi in cellPairsInd_p]
    return sort!(filter!(tp -> tp[1]>0, timepair_tp), rev=true)[1]
end

"""Search for collisions."""
function collisionFinder(positions_dim_id::AbstractArray{T, 2} where T,
                            locations_dim_id::AbstractArray{T, 2} where T,
                            collisionRadius::Float64, bounds::Bounds)
    tree = BallTree(positions_dim_id, PeriodicEuclidean([bounds.xLen, bounds.yLen]))
    return inrange(tree, locations_dim_id, collisionRadius)
end

function collisionFinder(arena::Arena)
    positions_dim_id = cellPositions_DIM_ID(arena)
    collisionRadius = collisionCrossSection(arena)
    return collisionFinder(positions_dim_id, positions_dim_id, collisionRadius, arena.bounds)
end

function collisionFinder(arena::Arena, inds::Vector{Int64})
    positions_dim_id = cellPositions_DIM_ID(arena)
    collisionRadius = collisionCrossSection(arena)
    return collisionFinder(positions_dim_id, positions_dim_id[:,inds], collisionRadius, arena.bounds)
end

function collCheck!(nbh::Vector, collidedCellsList::Array{Int})
    if length(nbh) < 2
        return false
    end
    removeDoubleCounts!(nbh, collidedCellsList)
    if length(nbh) < 2
        return false
    else
        return true
    end
end

# function nbhComposer!(nbh::Vector, collidedCellsList::Array{Int},
#                         nbhoods_id::AbstractVector{V} where V)
#     if length(nbh) < 2
#         return false
#     else
#         removeDoubleCounts!(nbh, collidedCellsList)
#         if length(nbh) < 2
#             return false
#         else
#             # take union of all neighborhoods of cells in current nbh
#             union!(nbh, (nbhoods_id[i] for i in nbh)...)
#             return true
#         end
#     end
# end

function nbhComposer!(nbh::Vector, nbhoods_id::AbstractVector{V} where V)
    # take union of all neighborhoods of cells in current nbh
    union!(nbh, (nbhoods_id[i] for i in nbh)...)
end

function removeDoubleCounts!(nbh::Vector{Int}, pastColls::Vector{Int})
    deleteat!(nbh, [(nbh[i] in pastColls) for i in 1:length(nbh)])
end

function reverseCellTime!(bounds::Bounds, t::Float64, cells_id::Vararg{Cell})
    for cell in cells_id
        moveCell!(cell, bounds, -t)
    end
    return nothing
end

"""Calculate exact time within timestep when collision took place."""
function collTimeCalc(cellA::Cell, cellB::Cell, bounds::Bounds)
    # r = cellA.pos .- cellB.pos
    r = connectingVector(cellA.pos, cellB.pos, [bounds.xLen, bounds.yLen])
    v = cellA.vel .- cellB.vel
    d = cellA.radius + cellB.radius
    if cellDistance(cellA, cellB, bounds) < d   # Only calculate if cells actually overlap. If this is not the case, complex solutions may occur.
        return ( r⋅v + sqrt( d^2*(v⋅v) - (r[2]*v[1]-r[1]*v[2])^2 ) ) / (v⋅v)
    else
        return -1. # any negative value will automatically be discarded later at collTimeCompare
    end
end

"""Perform collision of two cells."""
function collideCells!(bounds::Bounds, cellA::Cell, cellB::Cell, tRetc::Float64)
    #undo cell overlap
    reverseCellTime!(bounds, tRetc, cellA, cellB)

    #calculate new velocities
    newVelA = collVelCalc(cellA.vel, cellA.pos, cellB.vel, cellB.pos, [bounds.xLen, bounds.yLen])
    newVelB = collVelCalc(cellB.vel, cellB.pos, cellA.vel, cellA.pos, [bounds.xLen, bounds.yLen])
    cellA.vel .= newVelA
    cellB.vel .= newVelB

    #finish movement with remaining time
    if tRetc < 1
        moveCell!(cellA, bounds, tRetc)
        moveCell!(cellB, bounds, tRetc)
    else
        println("anomalous collision time ", tRetc, " found in step time: ", 1)
        # println("A position: ", cellA.pos)
        # println("B position: ", cellB.pos)
    end
    return nothing
end

"""Calculate velocities of cell pair after collision."""
function collVelCalc(velA::MVector{2, Float64}, posA::MVector{2, Float64},
                        velB::MVector{2, Float64}, posB::MVector{2, Float64},
                            period::AbstractVector)
    posR = connectingVector(posA, posB, period)
    return velA .- (velA .- velB)⋅(posR)/norm(posR)^2 * posR
end

function unoverlapPair!(cellA::Cell, cellB::Cell, bounds::Bounds)
    rAxis = cellA.pos .- cellB.pos
    ovFrac = ((cellA.radius + cellB.radius)-norm(rAxis))/norm(rAxis)
    if ovFrac < 0
        println("Error unoverlapping cells: Cells do not overlap.")
        return nothing
    end
    # println("initial distance between cells: ", cellDistance(cellA, cellB))
    # println("moving...")
    moveCell!(cellA, bounds, (ovFrac/2 + cellA.radius/100)*rAxis)
    moveCell!(cellB, bounds, -(ovFrac/2 + cellB.radius/100)*rAxis)
    # println("new distance between cells: ", cellDistance(cellA, cellB))
    return nothing
end

# ==== Utilities ====

function collisionCrossSection(arena::Arena)
    2arena.cellsList[1].radius
end

function collisionCrossSection(cell::Cell)
    2cell.radius
end

function collisionCrossSection(cellA::Cell, cellB::Cell)
    cellA.radius + cellB.radius
end

function getPairs(elems::Vararg{Int})
    pairs_p = Vector(undef,0)
    for i in 1:length(elems)
        for el2 in elems[i+1:end]
            push!(pairs_p, [elems[i], el2])
        end
    end
    return pairs_p
end
