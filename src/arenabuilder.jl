
### ======== Initializing functions ========

"""Build a new arena with a specified number of randomly placed cells."""
function buildRandArena(bounds::Bounds, nCells::Int, cellRadius::Real, s0::Real;
    fixSpeed=true, verbose=false, overlapScans=40, attempts=10)

    verbose && println("Building arena...")
    arena = Arena([randCell(bounds, cellRadius, s0, fixSpeed=fixSpeed) for i=1:nCells], bounds)
    overlaps = overlapScan(arena; verbose=verbose)
    if length(overlaps) > 0
        unique!(sort!(overlaps))
        deleteat!(arena.cellsList, overlaps)
    end
    nRemoved = length(overlaps)
    verbose && println(string(nRemoved)*" overlaps found.")

    if nRemoved == 0
        return arena
    end

    for i in 1:nRemoved
        # make ballTree of positions
        position_dim_id = cellPositionsPeriodic_DIM_ID(arena)
        posTree = BallTree(position_dim_id, PeriodicEuclidean([bounds.xLen, bounds.yLen]))

        noOverlap = false
        while !noOverlap
            # make a new cell
            newCell = randCell(bounds, cellRadius, s0; fixSpeed=fixSpeed)
            # check whether it overlaps any existing cells OR any newly added cells:
            overlaps = inrange(posTree, newCell.pos, collisionCrossSection(arena))
            if length(overlaps) == 0
                push!(arena.cellsList, newCell)
                noOverlap=true
            end
        end
    end

    verbose && println("Arena with "*string(length(arena.cellsList))*" cells built.")
    # test arena overlaps one final time

    return arena
end

"""Build a new arena with a specified number of randomly placed cells."""
function buildRandArena(bounds::Tuple{Tuple{Real, Real}, Tuple{Real, Real}},
     nCells::Int, cellRadius::Real, s0::Real; fixSpeed=false, verbose=false, overlapScans=40, attempts=10)
     bounds = Bounds(bounds[1], bounds[2])
     buildRandArena(bounds, nCells, cellRadius, s0;
        fixSpeed=fixSpeed, attempts=attempts, verbose=verbose, overlapScans=overlapScans)
end

"""Perform a scan to find and delete cells which overlap."""
function deleteArenaOverlaps!(arena::Arena)
    overlaps = overlapScan(arena; verbose=verbose)
    if length(overlaps) > 0
        unique!(sort!(overlaps))
        deleteat!(arena.cellsList, overlaps)
    end
    return length(overlaps)
end

# """Find and unoverlap cells in arena."""
# function overlapScan!(arena::Arena; verbose=false)
#     nbhoods_id = collisionFinder(arena)
#     unoverlappedCellsList = Int[]
#     for nbh in nbhoods_id
#         # overlap = nbhCollCheck!(nbh, unoverlappedCellsList, nbhoods_id)
#         overlap = collCheck!(nbh, unoverlappedCellsList)
#         if overlap === false
#             continue
#         end
#         nbhComposer!(nbh, nbhoods_id)
#         # println(">> Overlap found:")
#         # if length(nbh)==2
#         #     push!(unoverlappedCellsList, nbh[1], nbh[2])
#         #     unoverlapPair!(arena.cellsList[nbh[1]], arena.cellsList[nbh[2]], arena.bounds)
#         # else
#         #     nbh = nbh[1:2]
#         #     push!(unoverlappedCellsList, nbh...)
#         #     unoverlapPair!(arena.cellsList[nbh[1]], arena.cellsList[nbh[2]], arena.bounds)
#         # end
#         for cID in nbh[2:end]
#             randPos!(arena.cellsList[cID], arena.bounds)
#             push!(unoverlappedCellsList, nbh[2:end]...)
#             # if verbose
#             #     println("moved cell ID's: ", nbh[2:end])
#             # end
#         end
#     end
#     return unoverlappedCellsList
# end

"""Find and unoverlap cells in arena."""
function overlapScan(arena::Arena; verbose=false)
    nbhoods_id = collisionFinder(arena)
    removeCellsList = Int[]
    for nbh in nbhoods_id
        overlap = collCheck!(nbh, removeCellsList)
        if overlap === false
            continue
        end
        nbhComposer!(nbh, nbhoods_id)
        for cID in nbh[2:end]
            push!(removeCellsList, nbh[2:end]...)
        end
    end
    return removeCellsList
end



"""Create a random cell within given boundaries"""
function randCell(bounds::Bounds, radius::Real, speed::Real; fixSpeed=false)
    xPos = bounds.x[1] + (bounds.x[2]-bounds.x[1])*rand()
    yPos = bounds.y[1] + (bounds.y[2]-bounds.y[1])*rand()

    # note:
    #if fixSpeed=true speed will not equal the average particle speed at equilibrium
    if fixSpeed
        α = rand() * 2*π
        xVel = speed * cos(α)
        yVel = speed * sin(α)
    else
        s = rand(Rayleigh(speed*√(2/π)))
        α = rand() * 2*π
        xVel = s * cos(α)
        yVel = s * sin(α)
    end
    return Cell(@MVector[xPos,yPos], @MVector[xVel, yVel], radius)
end

"""Change a given cell's position to a new random position within specified bounds."""
function randPos!(cell::Cell, bounds::Bounds)
    cell.pos[1] = bounds.x[1] + (bounds.x[2]-bounds.x[1])*rand()
    cell.pos[2] = bounds.y[1] + (bounds.y[2]-bounds.y[1])*rand()
    return nothing
end


### ======== visualization ========

"""Plot cell positions in a given arena"""
function plotArena(arena::Arena, collidedCellsList_id::Union{Vector{Int}, Nothing}, nDaughters::Int; displayPlot=true, animator::Union{Animation, Nothing}=nothing)
    arenaCellPositions_dim_id = cellPositionsPeriodic_DIM_ID(arena)
    s = scatter(arenaCellPositions_dim_id[1,:], arenaCellPositions_dim_id[2,:], xlims = (0,10), ylims = (0,10), legend=false)
    if !(collidedCellsList_id===nothing)
        collisions_dim_id = arenaCellPositions_dim_id[:, collidedCellsList_id]
        scatter!(collisions_dim_id[1,:], collisions_dim_id[2,:])
    end
    if nDaughters > 0
        daughters_dim_id = arenaCellPositions_dim_id[:, end-nDaughters+1:end]
        scatter!(daughters_dim_id[1,:], daughters_dim_id[2,:])
    end
    if displayPlot
        display(s)
    end
    if animator!==nothing
        frame(animator, s)
    end
end


### ======== Utilities ========

"""Get cells in an arena from their index"""
function arenaIndsToCells(arena::Arena, inds_i::Union{Tuple{N, Int}, Array{Int, N}} where N)
    [arena.cellsList[ind] for ind in inds_i]
end
