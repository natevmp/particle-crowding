
### ======== Collision Functions ========

"""Perform collisions in a given arena. Return index of collided cells."""
function collider!(arena::Arena, tStep::Real; nbhcut=50, verbose=true)
    collidedCellsList = Int[]

    counter = 1
    while true
        collidedCellsPass = passThroughNeighborhood(arena, tStep; verbose=verbose)
        if length(collidedCellsPass)==0
            break
        end
        union!(collidedCellsList, collidedCellsPass)
        counter += 1
        if counter > nbhcut
            if verbose
                println(">>> Error: Could not resolve neighborhood collisions in ", nbhcut, " steps ! <<<")
            end
            break
        end
    end
    return collidedCellsList
end

"""Find and perform collisions in arena."""
function passThroughNeighborhood(arena::Arena, tStep::Real; verbose::Bool=true)
    # ! ======== DEBUG ========
    # println("\n ----- passing through neighborhoods -----\n")
    # ! =======================
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
            tRetc = collTimeCalc(arena.cellsList[nbh_c]..., arena.bounds, tStep; verbose=verbose)
            # ! ======== DEBUG ========
            # println("cell IDs: ", nbh_c)
            # ! =======================
        else
            # get earliest collision and time in composite neighborhood
            tRetc, nbh_c = collTimeCompare(arena, tStep, nbh_c...; verbose=verbose)
            # ! ======== DEBUG ========
            # println("cell IDs: ", nbh_c)
            # ! =======================
        end

        if tRetc>tStep && verbose   #debug
            println("anomalous collision time ", tRetc, " found in time step ", tStep)
            throw(ErrorException("the thing happened"))
        end
        #perform collision

        push!(collidedCellsOut, nbh_c...)
        
        # ! ======== DEBUG ========
        # plotArena(arena, nbh_c, 0; title=string(nbh_c), xlims=(3,5), ylims=(3,5))
        # ! =======================

        collideCells!(arena.bounds, arena.cellsList[nbh_c]..., tRetc, tStep; verbose=verbose)
        
        # ! ======== DEBUG ========
        # plotArena(arena, nbh_c, 0; title=string(nbh_c), xlims=(3,5), ylims=(3,5))
        # ! =======================

    end
    return collidedCellsOut
end

"""Search for collisions."""
function collisionFinder(arena::Arena)
    positionsP_dim_id = cellPositionsPeriodic_DIM_ID(arena)
    collisionRadius = collisionCrossSection(arena)
    return collisionFinder(positionsP_dim_id, positionsP_dim_id, collisionRadius, arena.bounds)
end

function collisionFinder(arena::Arena, inds::Vector{Int64})
    positionsP_dim_id = cellPositionsPeriodic_DIM_ID(arena)
    collisionRadius = collisionCrossSection(arena)
    return collisionFinder(positionsP_dim_id, positionsP_dim_id[:,inds], collisionRadius, arena.bounds)
end

function collisionFinder(positions_dim_id::AbstractArray{T, 2} where T,
                            locations_dim_id::AbstractArray{T, 2} where T,
                            collisionRadius::Float64, bounds::Bounds)
    tree = BallTree(positions_dim_id, PeriodicEuclidean([bounds.xLen, bounds.yLen]))
    return inrange(tree, locations_dim_id, collisionRadius)
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

function nbhComposer!(nbh::Vector, nbhoods_id::AbstractVector{V} where V)
    # take union of all neighborhoods of cells in current nbh
    union!(nbh, (nbhoods_id[i] for i in nbh)...)
end

function removeDoubleCounts!(nbh::Vector{Int}, pastColls::Vector{Int})
    deleteat!(nbh, [(nbh[i] in pastColls) for i in 1:length(nbh)])
end

"""Calculate exact time within timestep when collision took place."""
function collTimeCalc(cellA::Cell, cellB::Cell, bounds::Bounds, tStep::Real; verbose::Bool=true)
    # r = cellA.pos .- cellB.pos
    r = connectingVector(cellA.pos, cellB.pos, bounds)
    v = cellA.vel .- cellB.vel
    d = cellA.radius + cellB.radius
    if cellDistance(cellA, cellB, bounds) <= d   # Only calculate if cells actually overlap. If this is not the case, complex solutions may occur.
        tColl = ( r⋅v + sqrt( d^2*(v⋅v) - (r[2]*v[1]-r[1]*v[2])^2 ) ) / (v⋅v)
        # if tColl > tStep && verbose
        #     println("anomalous collision time ", tColl, " found in time step ", tStep)
        # end
        return tColl
    else
        # this can occur if a higher order collision neighborhood (>2 cells) is detected, but two of the cells are not actually overlapping. This function is still called, so seems to be easily dealt with by simply returning a "skip" value (e.g. -1)
        return -1. # any negative value will automatically be discarded later at collTimeCompare
    end
end

"""Find earliest collision in higher order (3 or more cells) collisions"""
function collTimeCompare(arena::Arena, tStep::Real, nbh_c::Vararg{Int}; verbose=true)
    cellPairsInd_p = getPairs(nbh_c...)

    
    timepair_tp = [
        (collTimeCalc(arenaIndsToCells(arena, cpi)..., arena.bounds, tStep; verbose=verbose), cpi) 
        for cpi in cellPairsInd_p 
    ]
    
    ## ! ==== debug ====
    # println("\nMult collision detected:")
    # for tp in sort!(timepair_tp, rev=true)
    #     println(tp, "; positions ", arena.cellsList[tp[2][1]].pos, " ", arena.cellsList[tp[2][2]].pos)
    # end
    ## ! ===============

    return sort!( filter!(tp -> tp[1]>0, timepair_tp), rev=true )[1]
end

"""Perform collision of two cells."""
function collideCells!(bounds::Bounds, cellA::Cell, cellB::Cell, tRetc::Real, tStep::Real; verbose=true)
    #undo cell overlap
    # reverseCellTime!(bounds, tRetc, cellA, cellB)
    
    ## ! ==== debug ====
    # println( "overlap distance: ", norm(connectingVector(cellA.pos, cellB.pos, bounds)) )
    # println("positions ", cellA.pos, "; ", cellB.pos)
    # println("time to retcon: ", tRetc)
    ## ! ===============

    reverseCellTime!(tRetc, cellA, cellB)

    ## ! ==== debug ====
    # println( "overlap distance: ", norm(connectingVector(cellA.pos, cellB.pos, bounds)))
    # println("positions ", cellA.pos, "; ", cellB.pos)
    ## ! ===============

    #calculate new velocities
    newVelA = collVelCalc(cellA.vel, cellA.pos, cellB.vel, cellB.pos, bounds)
    newVelB = collVelCalc(cellB.vel, cellB.pos, cellA.vel, cellA.pos, bounds)
    cellA.vel .= newVelA
    cellB.vel .= newVelB

    #finish movement with remaining time
    if 0<tRetc<tStep
        # moveCell!(cellA, bounds, tRetc)
        # moveCell!(cellB, bounds, tRetc)
        moveCell!(cellA, tRetc)
        moveCell!(cellB, tRetc)
    else
        # something went wrong if we're here
        nothing
    end

    ## ! ==== debug ====
    # println( "overlap distance: ", norm(connectingVector(cellA.pos, cellB.pos, bounds)) )
    # println("positions ", cellA.pos, "; ", cellB.pos, "\n")
    ## ! ===============

    return nothing
end

function reverseCellTime!(t::Float64, cells_id::Vararg{Cell})
    for cell in cells_id
        moveCell!(cell, -t)
    end
    return nothing
end

"""Calculate velocities of cell pair after collision."""
function collVelCalc(velA::MVector{2, Float64}, posA::MVector{2, Float64},
                        velB::MVector{2, Float64}, posB::MVector{2, Float64},
                            bounds::Bounds)
    posR = connectingVector(posA, posB, bounds)
    return velA .- (velA .- velB)⋅(posR)/norm(posR)^2 * posR
end

function unoverlapPair!(cellA::Cell, cellB::Cell, bounds::Bounds)
    # rAxis = cellA.pos .- cellB.pos
    rAxis = connectingVector(cellA.pos, cellB.pos, bounds)
    ovFrac = ((cellA.radius + cellB.radius)-norm(rAxis))/norm(rAxis)
    if ovFrac < 0
        # println("Error unoverlapping cells: Cells do not overlap.")
        return nothing
    end
    # println("initial distance between cells: ", cellDistance(cellA, cellB))
    # println("moving...")
    # moveCell!(cellA, bounds, (ovFrac/2 + cellA.radius/100)*rAxis)
    # moveCell!(cellB, bounds, -(ovFrac/2 + cellB.radius/100)*rAxis)
    moveCell!(cellA, (ovFrac/2 + cellA.radius/100)*rAxis)
    moveCell!(cellB, -(ovFrac/2 + cellB.radius/100)*rAxis)
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
    pairs_p = Vector{Vector{Int}}(undef,0)
    for i in 1:length(elems)
        for el2 in elems[i+1:end]
            push!(pairs_p, [elems[i], el2])
        end
    end
    return pairs_p
end
