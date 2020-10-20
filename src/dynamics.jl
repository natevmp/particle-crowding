### ======== dynamics functions ========

function evolveArena!(arena::Arena, time::Real, stepSize::Real, growthParams::Union{Dict, Nothing}=nothing;
    plotsteps=false, animator::Union{Animation, Nothing}=nothing, progress=true, verbose=true, saveTimeStep=1)

    # timepoints to evolve
    timeSteps = 0:stepSize:time

    # timepoints to save
    timeSaves = zeros(Float64, length(0:saveTimeStep:time))

    # preallocate position, velocity, and cell arrays
    posTime_t_dim_id = ElasticArray{Union{Float64, Missing}}(fill(missing, length(timeSaves), 2, length(arena.cellsList)))
    velTime_t_dim_id = ElasticArray{Union{Float64, Missing}}(fill(missing, length(timeSaves), 2, length(arena.cellsList)))
    cells_T_ID = Vector{Vector{Cell}}(undef, length(timeSaves))
    # --> record time 0
    snapshotCells!(posTime_t_dim_id, velTime_t_dim_id, arena, 1)
    snapshotCells!(cells_T_ID, arena, 1)

    # progresMeter
    if progress
        prog = Progress(length(timeSteps), 1)
    end

    nextSaveTime = 0+saveTimeStep
    nextSaveInd = 2

    for t in timeSteps[2:end]

        # == Movement ==
        for cell in arena.cellsList
            # moveCell!(cell, arena.bounds)
            moveCell!(cell, stepSize)
        end

        # == Collisions ==
        collidedCellsList_id = collider!(arena, stepSize; verbose=verbose)

        # == Growth ==
        if growthParams != nothing
            if t>growthParams["waitTime"]
                if growthParams["randGrowth"]
                    nDaughters = cultivateArena!(
                        arena, 1., growthParams["rateFunc"],
                        growthParams["radius"],
                        growthParams["speed"], randGrowth=true)
                else
                    # popSizeNew = growthParams["growthFunc"](t)
                    rateFunc(nPop) = growthParams["growthFunc"](t-growthParams["waitTime"]) - nPop
                    nDaughters = cultivateArena!(
                        arena, 1., rateFunc,
                        growthParams["radius"],
                        growthParams["speed"], randGrowth=false)
                end
            end
        else
            nDaughters = 0
        end

        # savepoint
        if t >= nextSaveTime
            # == Plot data ==
            if plotsteps||(animator!==nothing)
                plotArena(arena, collidedCellsList_id, nDaughters, displayPlot=plotsteps, animator=animator)
            end

            # == record data ==
            snapshotCells!(posTime_t_dim_id, velTime_t_dim_id, arena, nextSaveInd)
            snapshotCells!(cells_T_ID, arena, nextSaveInd)

            # == update save times ==
            timeSaves[nextSaveInd] = t
            nextSaveTime += saveTimeStep
            nextSaveInd += 1
        end

        # progressMeter
        if progress
            next!(prog)
        end
    end
    # gif(anim, "anim_2.gif", fps=15)
    return posTime_t_dim_id, velTime_t_dim_id, cells_T_ID, timeSaves
end

function evolveArena!(arena::Arena, steps::Int, growthParams::Union{Dict, Nothing}=nothing;
    plotsteps=false, animator::Union{Animation, Nothing}=nothing, progress=true, verbose=true)

    evolveArena!(arena, steps, 1, growthParams;
        plotsteps=plotsteps, animator=animator, progress=progress, verbose=verbose)

end
