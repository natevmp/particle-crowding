
# using JLD2
using JLD
using Statistics
include("../src/bmparticles.jl")
using .BParts

##
simNumber = parse(Int, ARGS[1])
rho = parse(Float64, ARGS[2])
time = parse(Int64, ARGS[3])
# simNumber = 1
# rho = 0.002
display(simNumber)
display(rho)

nSims = 1

# arenaParams = 
#     Dict{String, Any}(
#         "n0"=>500,
#         "evolveTime"=>200,
#         "bounds"=>((0.,22.4),(0.,22.4)), 
#         "radius"=>0.08, 
#         "speed"=>0.02,
#         "timeStep"=> 0.1
#     )
# growthParams =
#     Dict{String, Any}(
#         "ρ"=> rho,
#         "k"=> 5000,
#         "randGrowth"=> false,
#         "coldGrowth"=> false,
#         "waitTime"=> 100
#     )
arenaParams = 
    Dict(
        "n0"=>1000,
        # "evolveTime"=>5500,
        "evolveTime"=>time,
        "bounds"=>((0.,32),(0.,32)), 
        "radius"=>0.08, 
        "speed"=>0.02,
        "timeStep"=> 0.05
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

pos_Sim = Array{Array{Union{Float64, Missing},3}}(undef, nSims)
for i in 1:nSims
    succes = false
    while !succes 
        # try/catch construction in case boundserror occurs
        try
            _, posSim_t_dim_id, __, ___ = 
                BParts.randArenaEvolve(
                        arenaParams["n0"],
                        arenaParams["evolveTime"], 
                        arenaParams["timeStep"],
                        arenaParams,
                        growthParams;
                        coldGrowth=growthParams["coldGrowth"],
                        progress=false,
                        verbose=false);
            pos_Sim[i] = posSim_t_dim_id
        catch e
            if isa(e, BoundsError)
                continue
            else
                println(e)
                throw(ErrorException("unidentified error occurred"))
            end
        end
        succes = true
    end
end

println("Simulation finished. Calculating msd...")
# ==== Get mean squared displacements ====
msdTimes = (growthParams["waitTime"]+1, arenaParams["evolveTime"])
msd_Sim_t = (p_t_dim_id -> BParts.meanSquaredDisplacement(p_t_dim_id, msdTimes)).(pos_Sim);
msdPart_t = Array{Float64, 1}(undef, msdTimes[2]-msdTimes[1]+1)
for i in 1:(msdTimes[2]-msdTimes[1]+1)
    msdPart_t[i] = mean([ msd_tt[i] for msd_tt in msd_Sim_t ])
end

## ==== save data ====
println("MSD calculation finished. Saving data...")
saveName = "growingPop_multSims_rho"*string(growthParams["ρ"])*"_sim"*string(simNumber, pad=2)
# saveName = "growingPop_multSims"*"_sim"*string(simNumber, pad=2)

# @JLD2.save saveName*".jld2" arenaParams growthParams msdPart_t
JLD.save(saveName*".jld", "msdPart_t", msdPart_t)
JLD.save(saveName*"params.jld",
    "n0", arenaParams["n0"],
    "evolveTime", arenaParams["evolveTime"],
    # "bounds", arenaParams["bounds"],
    "boundsX", arenaParams["bounds"][1],
    "boundsY", arenaParams["bounds"][2],
    "radius", arenaParams["radius"],
    "speed", arenaParams["speed"],
    "ρ", growthParams["ρ"],
    "k", growthParams["k"],
    "randGrowth", growthParams["randGrowth"],
    "coldGrowth", growthParams["coldGrowth"],
    "waitTime", growthParams["waitTime"],
    )


println("saving complete")
