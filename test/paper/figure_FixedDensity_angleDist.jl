include("../../src/bmtheory.jl")
using .Theorist
using Distributions
using JLD2, FileIO


## ===== Load data =====

files_f = [
    "./data/FixedDensity/simResultFixedDensity_0.01.jld2",
    "./data/FixedDensity/simResultFixedDensity_0.05.jld2",
    "./data/FixedDensity/simResultFixedDensity_0.1.jld2",
    "./data/FixedDensity/simResultFixedDensity_0.2.jld2",
    "./data/FixedDensity/simResultFixedDensity_0.3.jld2",
    "./data/FixedDensity/simResultFixedDensity_0.5.jld2",
]

volumeDens(n, r, v) = n*π*r^2/v
ρ_sim = Float64[]
thermVals_sim = Dict[]
pos_sim_T_dim_id = Array{Float64, 3}[]
vel_sim_T_dim_id = Array{Float64, 3}[]
units_sim = Dict[]

for (i,file) in enumerate(files_f)
    arenaParams, growthParams, pos_t_dim_id, vel_t_dim_id = load(file, "arenaParams", "growthParams", "pos_t_dim_id", "vel_t_dim_id")

    thermVals = Theorist.thermalValues(arenaParams)
    colUnits = Dict(
        :t => 1 / thermVals["γ"],
        :x => arenaParams["speed"]/thermVals["γ"],
    )
    push!( ρ_sim, volumeDens(arenaParams["n0"], arenaParams["radius"], arenaParams["volume"]) )
    push!(thermVals_sim, thermVals)
    push!(pos_sim_T_dim_id, pos_t_dim_id)
    push!(vel_sim_T_dim_id, vel_t_dim_id)
    push!(units_sim, colUnits)
end

## ===== Collision angle distribution =====

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

function getAngleChanges(vel_t_dim_id)
    dθ_ = Float64[]
    for cid in 1:size(vel_t_dim_id, 3)
        dθ_cid = particleAngleChanges(@view vel_t_dim_id[:,:,cid])
        append!(dθ_, dθ_cid)
    end
    return dθ_
end

dθ_ = getAngleChanges(vel_t_dim_id)

figθ = histogram(dθ_)
display(figθ)

