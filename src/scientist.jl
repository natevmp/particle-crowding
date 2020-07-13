using Distributions
using PyCall
using Statistics

# ===== Running experiment =====


"""Construct an arena with randomly distributed cells and evolve it for a specified time."""
function randArenaEvolve(nCells::Int, steps::Int, arenaParams::Dict, growthParams::Union{Dict, Nothing}=nothing;
    plotting=false, animating=false, progress=true, verbose=true)

    arena = buildRandArena(arenaParams["bounds"], nCells, arenaParams["radius"], arenaParams["speed"], fixSpeed=true)

    arenaCellPositions_dim_id = BParts.cellPositions_DIM_ID(arena)

    eKin = BParts.kineticEnergy(arena)
    # println("::::: Initial total kinetic energy: ", eKin)

    if animating
        anim = Animation()
    else anim = nothing
    end

    if isnothing(growthParams)
        evolveGrowthParams = nothing
    else
        logisticRate(n::Real, ρ::Real, k::Real) = n*ρ*(1-n/k)

        evolveGrowthParams =
            Dict(
                "rateFunc"=> n->logisticRate(n, growthParams["ρ"], growthParams["k"]),
                "growthFunc"=> t->logisticGrowth(t, growthParams["ρ"], growthParams["k"],
                                                    arenaParams["n0"]),
                "radius"=> arenaParams["radius"],
                "speed"=> arenaParams["speed"],
                "randGrowth"=> growthParams["randGrowth"],
                "waitTime"=>growthParams["waitTime"]
            )
    end

    posTime_t_dim_id, velTime_t_dim_id, cells_T_ID =
        evolveArena!(arena, steps, evolveGrowthParams, plotsteps=plotting, animator=anim, progress=progress, verbose=verbose)

    if animating
        gif(anim, "figures/animation.gif", fps=10)
    end

    eKin = BParts.kineticEnergy(arena)
    # println("::::: Final total kinetic energy: ", eKin)

    return arena, posTime_t_dim_id, velTime_t_dim_id, cells_T_ID
end

# ===== Analysing experiments =====

function kineticEnergy(arena::Arena)
    sum(v -> norm(v)^2/2, [c.vel for c in arena.cellsList])
end

function snapshotPos(arena::Arena, t::Int, posTime_t_dim_id::AbstractArray)
    for (id, cell) in enumerate(arena.cellsList)
        posTime_t_dim_id[t, :, id] = cell.pos
    end
end

function snapshotVel(arena::Arena, t::Int, velTime_t_dim_id::AbstractArray)
    for (id, cell) in enumerate(arena.cellsList)
        velTime_t_dim_id[t, :, id] = cell.vel
    end
end

function snapshotCells!(posTime_t_dim_id::AbstractArray, velTime_t_dim_id::AbstractArray,
                        arena::Arena, t::Int)

    if length(arena.cellsList) > size(posTime_t_dim_id)[3]
        nBirths = length(arena.cellsList) - size(posTime_t_dim_id)[3]
        tDim = size(posTime_t_dim_id)[1]
        pDim = size(posTime_t_dim_id)[2]
        append!(posTime_t_dim_id, fill(missing, tDim, pDim, nBirths))
        append!(velTime_t_dim_id, fill(missing, tDim, pDim, nBirths))
    end
    for (id, cell) in enumerate(arena.cellsList)
        posTime_t_dim_id[t, :, id] = cell.pos
        velTime_t_dim_id[t, :, id] = cell.vel
    end
    return nothing
end

function snapshotCells!(cells_Time_ID::Vector{Vector{Cell}}, arena::Arena, t::Int)
    cells_Time_ID[t] = arena.cellsList
    return nothing
end

function speedCalc(velTime_t_dim_id::AbstractArray)
    # does not yet work with growing population (Missing values cause error)
    speedTime_t_id = Array{Float64}(undef, size(velTime_t_dim_id)[1], size(velTime_t_dim_id)[3])
    for t in 1:size(speedTime_t_id, 1)
        for id in 1:size(speedTime_t_id, 2)
            speedTime_t_id[t, id] = norm(velTime_t_dim_id[t, :, id])
        end
    end
    return speedTime_t_id
end

function speedCalc_Vardims(velTime_vardims_xy::AbstractArray{Float64, N} where N)
    vardims = size(velTime_vardims_xy)[1:end-1]
    speedTime_vardims = Array{Float64}(undef, vardims...)
    function assign!(state_vdims, var_vdims, var_vdims_xy, func)
        var_vdims[state_vdims...] = func(var_vdims_xy[state_vdims..., :])
    end
    nestfor(vardims, assign!, speedTime_vardims, velTime_vardims_xy, norm)
    return speedTime_vardims
end

function rayleighDistFit(velTime_t_dim_id::AbstractArray)
    speed_t_id = speedCalc(velTime_t_dim_id)
    stats = pyimport("scipy.stats")
    _, rshape = stats.rayleigh.fit(speed_t_id, floc=0)
    rDist = Rayleigh(rshape)
    return rDist
end


function velocityAutocorrelation(cells_T_ID::Vector{Vector{Cell}})

end

function velocityAutocorrelation(vel_t_dim_id::AbstractArray)
    # vCorr_id_t = zeros(Float64, size(vel_id_t_dim)[1], size(vel_id_t_dim)[2])
    vCorr_t = zeros(Float64, size(vel_t_dim_id)[1])
    for t in 1:length(vCorr_t)
        vCorr_t[t] = mean([ vel_t_dim_id[1, :, id] ⋅ vel_t_dim_id[t, :, id] for id in 1:(size(vel_t_dim_id)[3]) ])
    end
    return vCorr_t
end

function meanFreePath(vel_t_dim_id::AbstractArray, dt::Float64)
    nCells, nTime = size(vel_t_dim_id)[[3,1]]
    pathlengths_l = Vector{Float64}(undef, 0)
    for cellInd in 1:nCells
        speedCur::Float64 = norm(vel_t_dim_id[1, :, cellInd])
        pathsteps::Integer = 0
        for tInd in 2:nTime
            # if the x velocity does not change, no collision has occurred
            if vel_t_dim_id[tInd, 1, cellInd] == vel_t_dim_id[tInd-1, 1, cellInd]
                pathsteps += 1
                continue
            else
                append!(pathlengths_l, dt * pathsteps * speedCur)
                speedCur = norm(vel_t_dim_id[tInd, :, cellInd])
                pathsteps = 0
            end
        end
    end
    return mean(pathlengths_l)
end

function meanSquaredDisplacement(pos_t_dim_id::AbstractArray, bperiod_xy, tspan::Tuple{Real, Real})
    nCells = size(pos_t_dim_id)[3]
    steps = tspan[2]-tspan[1]+1
    msd_t = Vector{Float64}(undef, steps)
    for (i,t) in enumerate(range(tspan[1], tspan[2], step=1))
        sd_CID_t =
            [peuclidean(pos_t_dim_id[t, :, cellInd], pos_t_dim_id[tspan[1], :, cellInd], bperiod_xy)^2
            for cellInd in 1:nCells]
        msd_t[i] = mean(skipmissing(sd_CID_t))
    end
    return msd_t
end

function meanSquaredDisplacement(pos_t_dim_id::AbstractArray, bounds::Bounds, tspan::Tuple{Real, Real})
    meanSquaredDisplacement(pos_t_dim_id, [bounds.xLen, bounds.yLen], tspan)
end

function logisticGrowth(t, ρ, k, n0)
    return k/( 1 + (k-n0)/n0 * exp(-ρ*t) )
end

function nCellsTime(cells_T_ID)
    nCells_t = map(length, cells_T_ID)
end
