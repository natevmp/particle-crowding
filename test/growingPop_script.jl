using JLD2, Statistics
include("../src/bmparticles.jl")
include("../src/bmtheory.jl")
using .BParts

nSims = 50
arenaParams = 
    Dict(
        "n0"=>100,
        "evolveTime"=>1000,
        "bounds"=>((0.,10.),(0.,10.)), 
        "radius"=>0.08, 
        "speed"=>0.02,
        "timeStep"=> 0.1
    )
growthParams =
    Dict(
        "ρ"=> 0.002,
        "k"=> 1000,
        "randGrowth"=> false,
        "waitTime"=> 100
    )


BParts.extendParams!(arenaParams)
# timeSteps = arenaParams["timeStep"]:arenaParams["timeStep"]:arenaParams["evolveTime"]

# ==== Run simulations ====
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
                        coldGrowth=false,
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

# ==== Get mean squared displacements ====
msdTimes = (growthParams["waitTime"]+1, arenaParams["evolveTime"])
msd_Sim_t = (p_t_dim_id -> BParts.meanSquaredDisplacement(p_t_dim_id, msdTimes)).(pos_Sim);
msdPart_t = Array{Float64, 1}(undef, msdTimes[2]-msdTimes[1]+1)
for i in 1:(msdTimes[2]-msdTimes[1]+1)
    msdPart_t[i] = mean([ msd_tt[i] for msd_tt in msd_Sim_t ])
end

# ==== save data ====
saveName = "growingPop_multSims_rho"*string(growthParams["ρ"])*".jld2"
@save saveName arenaParams growthParams msd_Sim_t msdPart_t


