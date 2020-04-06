using Distributions
using PyCall
using Statistics

function kineticEnergy(arena::Arena)
    sum(v -> norm(v)^2/2, [c.vel for c in arena.cellsList])
end

function snapshotPos(arena::Arena, t::Int, posTime_id_t_dim::Array{Float64, 3})
    for (id, cell) in enumerate(arena.cellsList)
        posTime_id_t_dim[id, t, :] = cell.pos
    end
end

function snapshotVel(arena::Arena, t::Int, velTime_id_t_dim::Array{Float64, 3})
    for (id, cell) in enumerate(arena.cellsList)
        velTime_id_t_dim[id, t, :] = cell.vel
    end
end

function snapshotCells(arena::Arena, t::Int, posTime_id_t_dim::Array{Float64, 3}, velTime_id_t_dim::Array{Float64, 3})
    for (id, cell) in enumerate(arena.cellsList)
        posTime_id_t_dim[id, t, :] = cell.pos
        velTime_id_t_dim[id, t, :] = cell.vel
    end
end

function snapshotCells!(arena::Arena, t::Int, cells_Time_ID::Vector{Vector{Cell}})
    cells_Time_ID[t] = arena.cellsList
end

function speedCalc(velTime_id_t_dim::Array{Float64, 3})
    speedTime_id_t = Array{Float64}(undef, size(velTime_id_t_dim)[1], size(velTime_id_t_dim)[2])
    for id in 1:size(speedTime_id_t, 1)
        for t in 1:size(speedTime_id_t, 2)
            speedTime_id_t[id, t] = norm(velTime_id_t_dim[id, t, :])
        end
    end
    return speedTime_id_t
end

function speedCalc_Vardims(velTime_vardims_xy::Array{Float64, N} where N)
    vardims = size(velTime_vardims_xy)[1:end-1]
    speedTime_vardims = Array{Float64}(undef, vardims...)
    function assign!(state_vdims, var_vdims, var_vdims_xy, func)
        var_vdims[state_vdims...] = func(var_vdims_xy[state_vdims..., :])
    end
    nestfor(vardims, assign!, speedTime_vardims, velTime_vardims_xy, norm)
    return speedTime_vardims
end

function rayleighDistCompare(velTime_id_t_dim)
    speed_id_t = BParts.speedCalc(velTime_id_t_dim::Array{Float64, 3})
    sMean = [sum(speed_id_t[:, t])/size(speed_id_t, 1) for t in 1:size(speed_id_t, 2)]
    # d = Distributions.fit(Normal, velTime_id_t_dim)
    stats = pyimport("scipy.stats")
    _, rshape = stats.rayleigh.fit(speed_id_t, floc=0)
    rDist = Rayleigh(rshape)
end


function velocityAutocorrelation(cells_T_ID::Vector{Vector{Cell}})

end


function velocityAutocorrelation(vel_id_t_dim::Array{Float64, 3})
    # vCorr_id_t = zeros(Float64, size(vel_id_t_dim)[1], size(vel_id_t_dim)[2])
    vCorr_t = zeros(Float64, size(vel_id_t_dim)[2])
    for t in 1:size(vel_id_t_dim)[2]
        # vCorr_id_t[:, t] = [ Vel_id_t_dim[id, 1, :] ⋅ vel_id_t_dim[id, t, :] for id in range(size(vel_id_t_dim)[1]) ]
        vCorr_t[t] = mean([ vel_id_t_dim[id, 1, :] ⋅ vel_id_t_dim[id, t, :] for id in 1:(size(vel_id_t_dim)[1]) ])
    end
    return vCorr_t
end

function meanFreePath(vel_id_t_dim::Array{Float64, 3}, dt::Float64)
    nCells, nTime = size(vel_id_t_dim)[[1,2]]
    pathlengths_l = Vector{Float64}(undef, 0)
    for cellInd in 1:nCells
        speedCur::Float64 = norm(vel_id_t_dim[cellInd, 1, :])
        pathsteps::Integer = 0
        for tInd in 2:nTime
            if vel_id_t_dim[cellInd, tInd, 1] == vel_id_t_dim[cellInd, tInd-1, 1]
                pathsteps += 1
                # println("vxCur: ", vel_id_t_dim[cellInd, tInd, 1])
                # println("pathsteps: ", pathsteps)
                # println("speedCur: ", speedCur)
                continue
            else
                append!(pathlengths_l, dt * pathsteps * speedCur)
                speedCur = norm(vel_id_t_dim[cellInd, tInd, :])
                pathsteps = 0
                # println("========== Collision! ==========")
                # println("vxCur: ", vel_id_t_dim[cellInd, tInd, 1])
                # println("pathsteps: ", pathsteps)
                # println("speedCur: ", speedCur)
            end
        end
    end
    return mean(pathlengths_l)
end


function meanSquaredDisplacement(pos_id_t_dim::Array{Float64, 3}, bperiod::Vector{Float64})
    nCells, nTime = size(pos_id_t_dim)[[1,2]]
    msd_t = Vector{Float64}(undef, nTime)
    for t in 1:nTime
        msd_t[t] = mean([distance(pos_id_t_dim[cellInd, t, :], pos_id_t_dim[cellInd, 1, :], bperiod)^2
                    for cellInd in 1:nCells])
    end
    return msd_t
end
