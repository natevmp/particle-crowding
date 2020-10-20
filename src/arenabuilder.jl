
### ======== Initializing functions ========

"""Build a new arena with a specified number of randomly placed cells."""
function buildRandArena(bounds::Bounds, nCells::Int, cellRadius::Real, s0::Real;
    fixSpeed=false, overlapScans=40, attempts=10, verbose=false)

    verbose && println("Building arena...")
    for attempt in 1:attempts
        arena = Arena([randCell(bounds, cellRadius, s0, fixSpeed=fixSpeed) for i=1:nCells], bounds)
        success = fixArenaOverlaps!(arena; scanLimit=overlapScans, verbose=verbose)
        if success
            verbose && println("Arena succesfully built")
            return arena
        end
    end
    println("could not resolve overlaps in ", attempts, " attempts. Try a smaller number of cells or increase the number of attempts.")
    return nothing
end

"""Build a new arena with a specified number of randomly placed cells."""
function buildRandArena(bounds::Tuple{Tuple{Real, Real}, Tuple{Real, Real}},
     nCells::Int, cellRadius::Real, s0::Real; fixSpeed=false, overlapScans=40, attempts=10, verbose=false)
     bounds = Bounds(bounds[1], bounds[2])
     buildRandArena(bounds, nCells, cellRadius, s0;
        fixSpeed=fixSpeed, attempts=attempts, verbose=verbose, overlapScans=overlapScans)
end

"""Perform multiple scans to find and fix arena overlaps."""
function fixArenaOverlaps!(arena::Arena; scanLimit::Int=40, verbose=false)
    succeeded = false
    limitReached = false
    counter = 0
    while (!succeeded) && (!limitReached)
        overlaps = overlapScan!(arena; verbose=verbose)
        if length(overlaps)==0
            if verbose
                println("No overlaps found.")
            end
            succeeded = true
        else
            if verbose
                println(Int(ceil(length(overlaps)/2)), " overlaps fixed. Making another pass...")
            end
        end
        counter += 1
        # println(counter)
        # println(scanLimit)
        if counter > scanLimit
            # println("could not resolve overlaps in ", scanLimit, " passes. Try a different cell configuration or a smaller number of cells.")
            limitReached = true
        end
    end
    return succeeded
end

"""Find and unoverlap cells in arena."""
function overlapScan!(arena::Arena; verbose=false)
    nbhoods_id = collisionFinder(arena)
    unoverlappedCellsList = Int[]
    for nbh in nbhoods_id
        # overlap = nbhCollCheck!(nbh, unoverlappedCellsList, nbhoods_id)
        overlap = collCheck!(nbh, unoverlappedCellsList)
        if overlap === false
            continue
        end
        nbhComposer!(nbh, nbhoods_id)
        # println(">> Overlap found:")
        # if length(nbh)==2
        #     push!(unoverlappedCellsList, nbh[1], nbh[2])
        #     unoverlapPair!(arena.cellsList[nbh[1]], arena.cellsList[nbh[2]], arena.bounds)
        # else
        #     nbh = nbh[1:2]
        #     push!(unoverlappedCellsList, nbh...)
        #     unoverlapPair!(arena.cellsList[nbh[1]], arena.cellsList[nbh[2]], arena.bounds)
        # end
        for cID in nbh[2:end]
            randPos!(arena.cellsList[cID], arena.bounds)
            push!(unoverlappedCellsList, nbh[2:end]...)
            # if verbose
            #     println("moved cell ID's: ", nbh[2:end])
            # end
        end
    end
    return unoverlappedCellsList
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
    newCell = Cell(@MVector[xPos,yPos], @MVector[xVel, yVel], radius)
    return newCell
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
