using Distributions
using PyCall
using Statistics

function kineticEnergy(arena::Arena)
    sum(v -> norm(v)^2/2, [c.vel for c in arena.cellsList])
end

function snapshotPos(arena::Arena, t::Int, posTime_t_dim_id::AbstractArray{Float64, 3})
    for (id, cell) in enumerate(arena.cellsList)
        posTime_t_dim_id[t, :, id] = cell.pos
    end
end

function snapshotVel(arena::Arena, t::Int, velTime_t_dim_id::AbstractArray{Float64, 3})
    for (id, cell) in enumerate(arena.cellsList)
        velTime_t_dim_id[t, :, id] = cell.vel
    end
end

function snapshotCells!(posTime_t_dim_id::AbstractArray{Float64, 3}, velTime_t_dim_id::AbstractArray{Float64, 3}, 
                        arena::Arena, t::Int)
    
    if length(arena.cellsList) > size(posTime_t_dim_id)[3]
        nBirths = length(arena.cellsList) - size(posTime_t_dim_id)[3]
        tDim = size(posTime_t_dim_id)[1]
        pDim = size(posTime_t_dim_id)[2]
        append!(posTime_t_dim_id, fill(NaN, tDim, pDim, nBirths))
        append!(velTime_t_dim_id, fill(NaN, tDim, pDim, nBirths))
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

function speedCalc(velTime_t_dim_id::AbstractArray{Float64, 3})
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

function rayleighDistCompare(velTime_t_dim_id)
    speed_t_id = speedCalc(velTime_t_dim_id)
    # using scipy
    stats = pyimport("scipy.stats")
    _, rshape = stats.rayleigh.fit(speed_t_id, floc=0)
    rDist = Rayleigh(rshape)
    # using Distributions -- problem: doesn't ignore NaNs?
    # rDist = Distributions.fit(Rayleigh, reshape(speed_t_id, :))
    return rDist
end

function velocityAutocorrelation(cells_T_ID::Vector{Vector{Cell}})

end


function velocityAutocorrelation(vel_t_dim_id::Array{Float64, 3})
    # vCorr_id_t = zeros(Float64, size(vel_id_t_dim)[1], size(vel_id_t_dim)[2])
    vCorr_t = zeros(Float64, size(vel_t_dim_id)[1])
    for t in 1:length(vCorr_t)
        vCorr_t[t] = mean([ vel_t_dim_id[1, :, id] ⋅ vel_t_dim_id[t, :, id] for id in 1:(size(vel_t_dim_id)[3]) ])
    end
    return vCorr_t
end

function meanFreePath(vel_t_dim_id::AbstractArray{Float64, 3}, dt::Float64)
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


function meanSquaredDisplacement(pos_t_dim_id::AbstractArray{Float64, 3}, bperiod::Vector{Float64})
    nCells, nTime = size(pos_t_dim_id)[[3,1]]
    msd_t = Vector{Float64}(undef, nTime)
    for t in 1:nTime
        msd_t[t] = mean([peuclidean(pos_t_dim_id[t, :, cellInd], pos_t_dim_id[1, :, cellInd], bperiod)^2
                    for cellInd in 1:nCells])
    end
    return msd_t
end
