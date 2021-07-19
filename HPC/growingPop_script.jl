
using JLD2, FileIO
# using JLD
using Statistics
include("../src/bmparticles.jl")
using .BParts

##
simNumber = parse(Int, ARGS[1])
# simNumber = 1

_rho = [0.01, 0.02, 0.05, 0.1]

rho = _rho[simNumber]
display(simNumber)
display(rho)
# time = parse(Int64, ARGS[3])
# simNumber = 1
# rho = 0.002


arenaParams = Dict(
    "n0"=>1000,
    "evolveTime"=>4500,
    # "evolveTime"=>50,
    "bounds"=>((0.,32),(0.,32)), 
    "radius"=>0.08, 
    "speed"=>0.02,
    "timeStep"=> 0.1
)

growthParams =
    Dict(
        "ρ"=> rho,
        "k"=> 10000,
        "randGrowth"=> false,
        "coldGrowth"=> false,
        "waitTime"=> 500
    )

BParts.extendParams!(arenaParams)

println(arenaParams)
println(growthParams)

## ==== Run simulations ====
println("running simulation...")

@time _, pos_t_dim_id, vel_t_dim_id, cells_t_Id, times_t = 
    BParts.randArenaEvolve(
            arenaParams["n0"],
            arenaParams["evolveTime"], 
            arenaParams["timeStep"],
            arenaParams,
            growthParams;
            coldGrowth=growthParams["coldGrowth"],
            progress=false,
            verbose=false);


println("Simulation finished. Calculating msd...")
# ==== Get mean squared displacements ====
msdTimes = (growthParams["waitTime"]+1, arenaParams["evolveTime"])
msdPart_t = BParts.meanSquaredDisplacement(pos_t_dim_id, msdTimes)

## ==== save data ====
println("MSD calculation finished. Saving data...")
filename = "growingPop_rho"*string(growthParams["ρ"])*".jld2"

save(filename, "arenaParams", arenaParams, "growthParams", growthParams, "pos_t_dim_id", pos_t_dim_id, "msdPart_t", msdPart_t, "times_t", times_t)

println("saving complete")
