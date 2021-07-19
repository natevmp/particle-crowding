# using Distributions
# # using PyCall
# using Statistics

# ===== Running experiment =====


"""Construct an arena with randomly distributed cells and evolve it for a specified time."""
function randArenaEvolve(
    nCells::Int,
    time::Real,
    stepSize::Real,
    arenaParams::Dict,
    growthParams::Union{Dict, Nothing}=nothing;
    coldGrowth=false, plotting=false, animating=false, progress=true, verbose=false, attempts=1, overlapScans=0)

    arena = buildRandArena(arenaParams["bounds"], nCells, arenaParams["radius"], arenaParams["speed"];
                fixSpeed=true, verbose=verbose)

    arenaCellPositions_dim_id = BParts.cellPositions_DIM_ID(arena)

    eKin = BParts.kineticEnergy(arena)
    verbose && println("::::: Initial total kinetic energy: ", eKin)

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
                "growthFunc"=> t->logisticGrowth(t, growthParams["ρ"], growthParams["k"], arenaParams["n0"]),
                "radius"=> arenaParams["radius"],
                "randGrowth"=> growthParams["randGrowth"],
                "waitTime"=>growthParams["waitTime"]
            )
            if coldGrowth
                evolveGrowthParams["speed"] = 0.
            else
                evolveGrowthParams["speed"] = arenaParams["speed"]
            end
    end

    posTime_t_dim_id, velTime_t_dim_id, cells_T_ID, times_t =
        evolveArena!(arena, time, stepSize, evolveGrowthParams;
            coldGrowth=coldGrowth, plotsteps=plotting, animator=anim, progress=progress, verbose=verbose)

    if animating
        gif(anim, "figures/animation.gif", fps=10)
    end

    eKin = BParts.kineticEnergy(arena)
    # println("::::: Final total kinetic energy: ", eKin)

    return arena, posTime_t_dim_id, velTime_t_dim_id, cells_T_ID, times_t
end

function randArenaEvolve(
    arena::Arena,
    time::Real,
    stepSize::Real,
    arenaParams::Dict,
    growthParams::Union{Dict, Nothing}=nothing;
    coldGrowth=false, plotting=false, animating=false, progress=true, verbose=false, attempts=1, overlapScans=0, saveTimeStep=1)

    arenaCellPositions_dim_id = BParts.cellPositions_DIM_ID(arena)

    eKin = BParts.kineticEnergy(arena)
    verbose && println("::::: Initial total kinetic energy: ", eKin)

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
                "randGrowth"=> growthParams["randGrowth"],
                "waitTime"=>growthParams["waitTime"]
            )
            if coldGrowth
                evolveGrowthParams["speed"] = 0.
            else
                evolveGrowthParams["speed"] = arenaParams["speed"]
            end
    end

    posTime_t_dim_id, velTime_t_dim_id, cells_T_ID, times_t =
        evolveArena!(arena, time, stepSize, evolveGrowthParams;
            coldGrowth=coldGrowth, plotsteps=plotting, animator=anim, progress=progress, verbose=verbose, saveTimeStep=saveTimeStep)

    if animating
        gif(anim, "figures/animation.gif", fps=10)
    end

    eKin = BParts.kineticEnergy(arena)
    # println("::::: Final total kinetic energy: ", eKin)

    return arena, posTime_t_dim_id, velTime_t_dim_id, cells_T_ID, times_t
end

function randArenaEvolve(nCells::Int, steps::Int, arenaParams::Dict, growthParams::Union{Dict, Nothing}=nothing;
    plotting=false, animating=false, progress=true, verbose=true)

    randArenaEvolve(nCells, steps, 1, arenaParams, growthParams;
        plotting=plotting, animating=animating, progress=progress, verbose=verbose)

end

# ===== Analysing experiments =====
function extendParams!(arenaParams::Dict)
    bounds = arenaParams["bounds"]
    arenaParams["volume"] = abs(bounds[1][2]-bounds[1][1])*abs(bounds[2][2]-bounds[2][1])
    arenaParams["bperiod"] = [abs(bounds[1][2]-bounds[1][1]), abs(bounds[2][2]-bounds[2][1])]
    arenaParams["E"] = arenaParams["speed"]^2/2
end

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
                        arena::Arena, tInd::Int)
    # increase pos_t_dim_id and vel_t_dim_id sizes to accomodate new cells
    if length(arena.cellsList) > size(posTime_t_dim_id)[3]
        nBirths = length(arena.cellsList) - size(posTime_t_dim_id)[3]
        tDim = size(posTime_t_dim_id)[1]
        pDim = size(posTime_t_dim_id)[2]
        append!(posTime_t_dim_id, fill(missing, tDim, pDim, nBirths))
        append!(velTime_t_dim_id, fill(missing, tDim, pDim, nBirths))
    end
    # add positions and velocities to arrays at timestep
    for (id, cell) in enumerate(arena.cellsList)
        posTime_t_dim_id[tInd, :, id] = cell.pos
        velTime_t_dim_id[tInd, :, id] = cell.vel
    end
    return nothing
end

function snapshotCells!(cells_Time_ID::Vector{Vector{Cell}}, arena::Arena, tInd::Int)
    cells_Time_ID[tInd] = arena.cellsList
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
    vCorr_t = zeros(Float64, size(vel_t_dim_id,1))
    for t in 1:length(vCorr_t)
        vCorr_t[t] = mean([ vel_t_dim_id[1, :, id] ⋅ vel_t_dim_id[t, :, id] for id in 1:(size(vel_t_dim_id,3)) ])
    end
    return vCorr_t
end

function velocityAutocorrelation(vel_t_dim_id::AbstractArray, times::Tuple{Real, Real}, tStep::Real=1)
    # vCorr_id_t = zeros(Float64, size(vel_id_t_dim)[1], size(vel_id_t_dim)[2])
    # tInds = length(tStep:tStep:time)
    tInds = Int(floor(times[1]/tStep)):Int(floor(times[2]/tStep))
    vCorr_t = velocityAutocorrelation(vel_t_dim_id[tInds, :, :])
    return 0:tStep:times[2]-times[1], vCorr_t
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

function meanSquaredDisplacement(pos_t_dim_id::AbstractArray, tspan::Tuple{Real, Real}, tStep::Real=1)
    nCells = size(pos_t_dim_id)[3]
    tInds = Int( ceil(tspan[1]/tStep) ) : Int( floor(tspan[2]/tStep) )
    msd_t = Vector{Float64}(undef, length(tInds))
    sd_cid = Vector{Union{Float64, Missing}}(undef, size(pos_t_dim_id)[3])
    for (i,tInd) in enumerate(tInds)
        sd_cid .= [euclidean(pos_t_dim_id[tInd, :, cellInd], pos_t_dim_id[tInds[1], :, cellInd])^2
            for cellInd in 1:nCells
        ]
        msd_t[i] = mean(skipmissing(sd_cid))
    end
    return msd_t
end

function meanSquaredDisplacement(pos_t_dim_id::AbstractArray)
    nCells = size(pos_t_dim_id,3)
    msd_t = Vector{Float64}(undef, size(pos_t_dim_id,1))
    sd_cid = Vector{Union{Float64, Missing}}(undef, size(pos_t_dim_id,3))
    for t in 1:size(pos_t_dim_id,1)
        sd_cid .= [euclidean(pos_t_dim_id[t, :, cellInd], pos_t_dim_id[1, :, cellInd])^2
            for cellInd in 1:nCells
        ]
        msd_t[t] = mean(skipmissing(sd_cid))
    end
    return msd_t
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

function photographImagePositions(pos_dim_id::AbstractArray, xBounds::Tuple{Real, Real}, yBounds::Tuple{Real, Real})
    posImage_dim_id = Array{Float64}(undef, 2, size(pos_dim_id)[2])
    println(size(pos_dim_id))
    println(size(posImage_dim_id))

    posImage_dim_id[1,:] = (x -> toBoundsPeriodic(x, xBounds)).(pos_dim_id[1,:])
    posImage_dim_id[2,:] = (x -> toBoundsPeriodic(x, xBounds)).(pos_dim_id[2,:])
    return posImage_dim_id
end



# ========== analysis ==========

function getAngle(v2::AbstractVector, v1::AbstractVector)
    θ1 = atan(v1[2], v1[1])
    θ2 = atan(v2[2], v2[1])
    atan( sin(θ2-θ1), cos(θ2-θ1) )
end

function particleAngleChanges(vel_t_dim)
    dθ_ = Float64[]
    for tInd in 2:size(vel_t_dim, 1)
        if vel_t_dim[tInd, 1] != vel_t_dim[tInd-1, 1]
            dθ = getAngle(vel_t_dim[tInd, :], vel_t_dim[tInd-1, :])
            push!(dθ_, dθ)
        end
    end
    return dθ_
end

function ensembleAngleChanges(vel_t_dim_id)
    dθ_ = Float64[]
    for cid in 1:size(vel_t_dim_id, 3)
        dθ_cid = particleAngleChanges(@view vel_t_dim_id[:,:,cid])
        append!(dθ_, dθ_cid)
    end
    return dθ_
end

function particleFreePaths(vel_t_dim, dt=1)
    freePaths_ = Float64[]
    pathCur = 0.
    for tInd in 2:size(vel_t_dim, 1)
        if vel_t_dim[tInd, 1] == vel_t_dim[tInd-1, 1]
            pathCur += norm(vel_t_dim[tInd, :]) * dt
        else
            pathCur += norm(vel_t_dim[tInd, :]) * dt/2
            push!(freePaths_, pathCur)
            pathCur = 0.
        end
    end
    return freePaths_
end

function particleFreePaths(vel_t_dim, times_t::AbstractVector)
    freePaths_ = Float64[]
    pathCur = 0.
    for tInd in 2:length(times_t)
        dt = (times_t[tInd]-times_t[tInd-1])
        if vel_t_dim[tInd, 1] == vel_t_dim[tInd-1, 1]
            pathCur += norm(vel_t_dim[tInd, :]) * dt
        else
            pathCur += norm(vel_t_dim[tInd, :]) * dt/2
            push!(freePaths_, pathCur)
            pathCur = 0.
        end
    end
    return freePaths_
end

function ensembleFreePaths(vel_t_dim_id, dt::Real=1)
    freePaths_ = Float64[]
    for cid in 1:size(vel_t_dim_id, 3)
        append!(freePaths_, particleFreePaths((@view vel_t_dim_id[:,:,cid]), dt))
    end
    return freePaths_
end

function ensembleFreePaths(vel_t_dim_id, times_t::AbstractVector)
    freePaths_ = Float64[]
    for cid in 1:size(vel_t_dim_id, 3)
        append!(freePaths_, particleFreePaths((@view vel_t_dim_id[:,:,cid]), times_t))
    end
    return freePaths_
end