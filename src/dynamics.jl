### ======== dynamics functions ========

function evolveArena!(
    arena::Arena, 
    time::Real, 
    stepSize::Real, 
    growthParams::Union{Dict, Nothing}=nothing;
    coldGrowth=false, plotsteps=false, animator::Union{Animation, Nothing}=nothing, progress=true, verbose=true, saveTimeStep=1)

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

    # average particle energy
    Eav = avParticleEnergy(arena.cellsList)

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
        if growthParams !== nothing
            if t>growthParams["waitTime"]
                if growthParams["randGrowth"]
                    rateFunc = growthParams["rateFunc"]
                else
                    rateFunc(nPop) = growthParams["growthFunc"](t-growthParams["waitTime"]) - nPop
                end
                nDaughters = cultivateArena!(
                    arena, stepSize, rateFunc,
                    growthParams["radius"],
                    growthParams["speed"];
                    randGrowth=growthParams["randGrowth"])

                if nDaughters > 0
                    rmComDriftArena!(arena.cellsList)
                    # rescaleEnergy!(arena.cellsList, Eav)
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

# ====== COM normalization =====
function avParticleEnergy(cellsList::AbstractVector{C}) where C<:Cell
    ETot = 0.
    for cell in cellsList
        ETot += norm(cell.vel)^2/2
    end
    return ETot / length(cellsList)
end

function comVel(cellsList::AbstractVector{C}) where C<:Cell
    vTot = MVector{2,Float64}([0.,0.])
    for cell in cellsList
        vTot .+= cell.vel
    end
    return vTot / length(cellsList)
end

function rmComDriftArena!(cellsList::AbstractVector{C}) where C<:Cell
    # COM (drift) velocity
    vCom = comVel(cellsList)
    # COM (drift) energy
    ECom = norm(vCom)^2/2
    # average particle energy
    E = avParticleEnergy(cellsList)

    χ = sqrt(E/(E-ECom))

    # subtract drift and rescale each cell velocity
    for cell in cellsList
        cell.vel .= χ*(cell.vel - vCom)
    end
    return nothing
end

function rescaleEnergy!(cellsList::AbstractVector{C}, Etarget::Real) where C<:Cell
    # current average energy
    Ecur = avParticleEnergy(cellsList)
    χ = sqrt(Etarget/Ecur)
    # rescale each cell velocity to obtain target average energy
    for cell in cellsList
        cell.vel .*= χ
    end
    return nothing
end
