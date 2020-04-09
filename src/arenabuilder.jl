
### ======== Initializing functions ========

function buildRandArena(bounds::Bounds, nCells::Int, cellRadius::Float64, vStd::Float64; fixSpeed=false)
    arena = Arena([randCell(bounds, cellRadius, vStd, fixSpeed=fixSpeed) for i=1:nCells], bounds)
    # fixArenaOverlaps!(arena)
    return arena
end

function buildArena(bounds::Bounds, cells::Array{Cell, 1}=Array{Cell}(undef, 0))
    Arena(bounds, cells)
end

function fixArenaOverlaps!(arena::Arena, scanLimit::Int=20)
    println("Cleaning arena overlaps...")
    counter = 0
    while true
        overlaps = overlapScan!(arena)
        if length(overlaps)==0
            println("No overlaps found.")
            break
        else
            println(Int(ceil(length(overlaps)/2)), " overlaps fixed. Making another pass...")
        end
        counter += 1
        if counter > scanLimit
            error("could not resolve overlaps in ", scanLimit, " passes")
        end
    end

    return nothing
end

function buildArena(bounds::Bounds, cells::Vararg{Cell})
    return Arena([c for c in cells], bounds)
end

"""Create a random cell within given boundaries"""
function randCell(bounds::Bounds, radius::Real, vStd::Float64; fixSpeed=false)
    xPos = bounds.x[1] + (bounds.x[2]-bounds.x[1])*rand()
    yPos = bounds.y[1] + (bounds.y[2]-bounds.y[1])*rand()
    if fixSpeed
        α = rand() * 2*π
        xVel = vStd * cos(α)
        yVel = vStd * sin(α)
    else
        xVel = randn() * vStd
        yVel = randn() * vStd
    end
    newCell = Cell(@MVector[xPos,yPos], @MVector[xVel, yVel], radius)
    return newCell
end

function randPos!(cell::Cell, bounds::Bounds)
    cell.pos[1] = bounds.x[1] + (bounds.x[2]-bounds.x[1])*rand()
    cell.pos[2] = bounds.y[1] + (bounds.y[2]-bounds.y[1])*rand()
    return nothing
end

function overlapScan!(arena::Arena)
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
        if length(nbh)==2
            push!(unoverlappedCellsList, nbh[1], nbh[2])
            unoverlapPair!(arena.cellsList[nbh[1]], arena.cellsList[nbh[2]], arena.bounds)
        else
            nbh = nbh[1:2]
            push!(unoverlappedCellsList, nbh...)
            unoverlapPair!(arena.cellsList[nbh[1]], arena.cellsList[nbh[2]], arena.bounds)
        end
    end
    return unoverlappedCellsList
end

### ======== dynamics functions ========

# function evolveArena(arena::Arena, steps::Int; plotsteps=false, animator::Union{Animation, Nothing}=nothing)

#     posTime_id_t_dim = fill(NaN, length(arena.cellsList), steps, 2)
#     velTime_id_t_dim = fill(NaN, length(arena.cellsList), steps, 2)
#     cells_T_ID = Vector{Vector{Cell}}(undef, steps)

#     for t in 1:steps
#         for cell in arena.cellsList
#             moveCell!(cell, arena.bounds)
#         end

#         collidedCellsList_id = collider!(arena)
#         if plotsteps||(animator!==nothing)
#             plotArena(arena, collidedCellsList_id, displayPlot=plotsteps, animator=animator)
#         end

#         snapshotCells!(posTime_id_t_dim, velTime_id_t_dim, arena, t)
#         snapshotCells!(cells_T_ID, arena, t)
#     end
#     # gif(anim, "anim_2.gif", fps=15)
#     return posTime_id_t_dim, velTime_id_t_dim, cells_T_ID
# end

function evolveArena(arena::Arena, steps::Int, growthParams::Union{Tuple{Real,Real,Real}, Nothing}=nothing; plotsteps=false,
                     animator::Union{Animation, Nothing}=nothing)

    # posTime_id_t_dim = fill(NaN, length(arena.cellsList), steps, 2)
    # velTime_id_t_dim = fill(NaN, length(arena.cellsList), steps, 2)

    posTime_t_dim_id = ElasticArray(fill(NaN, steps, 2, length(arena.cellsList)))
    velTime_t_dim_id = ElasticArray(fill(NaN, steps, 2, length(arena.cellsList)))

    cells_T_ID = Vector{Vector{Cell}}(undef, steps)


    for t in 1:steps

        # == Movement ==
        for cell in arena.cellsList
            moveCell!(cell, arena.bounds)
        end

        # == Collisions ==
        collidedCellsList_id = collider!(arena)
        
        # == Growth ==
        if growthParams != nothing
            nDaughters = cultivateArena!(arena, 1., growthParams...)
        end

        # == Plot data ==
        if plotsteps||(animator!==nothing)
            plotArena(arena, collidedCellsList_id, nDaughters, displayPlot=plotsteps, animator=animator)
        end

        # == record data ==
        snapshotCells!(posTime_t_dim_id, velTime_t_dim_id, arena, t)
        snapshotCells!(cells_T_ID, arena, t)
    end
    # gif(anim, "anim_2.gif", fps=15)
    return posTime_t_dim_id, velTime_t_dim_id, cells_T_ID
end


### ======== visualization ========

function plotArena(arena::Arena, collidedCellsList_id::Union{Vector{Int}, Nothing}, nDaughters::Int; displayPlot=true, animator::Union{Animation, Nothing}=nothing)
    arenaCellPositions_dim_id = cellPositions_DIM_ID(arena)
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

function arenaIndsToCells(arena::Arena, inds_i::Union{Tuple{N, Int}, Array{Int, N}} where N)
    [arena.cellsList[ind] for ind in inds_i]
end
