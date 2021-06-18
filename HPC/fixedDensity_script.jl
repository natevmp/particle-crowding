include("../src/bmparticles.jl")
include("../src/bmtheory.jl")
using .BParts
using .Theorist
using Distributions
using JLD2, FileIO

##

## ====== Fixed Density ======
# First we compare the particle simulation with the simulations of the Langevin equation (and some analytical results from it) for a fixed population size.

# ---- Parameters for simulation ----


ρIn_ = [0.01, 0.05, 0.1, 0.2, 0.3, 0.5]


simNumber = parse(Int, ARGS[1])
ρIn = ρIn_[simNumber]

timeWaitC = 4
timeEvolveC = 100

arenaParams = Dict(
    # "evolveTime"=>2500,
    "n0"=> 2000,
    "bounds"=>((0.,1),(0.,1)), 
    "radius"=>0.08, 
    "speed"=>0.02,
    "timeStep"=> 0.1
)

arenaSizeFromDens(ρ, n0, r) = √(n0*π*r^2 / ρ)
arenaSize = round(
    arenaSizeFromDens(ρIn, arenaParams["n0"], arenaParams["radius"]),
    digits=3
)
# println(arenaSize)

arenaParams["bounds"] = ((0.,arenaSize), (0., arenaSize))

growthParams = Dict{String, Real}(
    "ρ"=> 0,
    "k"=> arenaParams["n0"],
    "randGrowth"=> false,
    "coldGrowth"=> false
)

Theorist.extendParams!(arenaParams)
for p in arenaParams
    println(p)
end

thermVals = Theorist.thermalValues(arenaParams)
volumeDens(n, r, v) = n*π*r^2/v
println("mean free time: ", 1/thermVals["γ"])
println("initial volume density: ", volumeDens(arenaParams["n0"], arenaParams["radius"], arenaParams["volume"]))
# println("final volume density: ", volumeDens(growthParams["k"], arenaParams["radius"], arenaParams["volume"]))

# save time points at 1/10 of the expected mfp
saveTimeStep = floor(1/thermVals["γ"])/10
growthParams["waitTime"] = timeWaitC / thermVals["γ"]
arenaParams["evolveTime"] = growthParams["waitTime"] + ceil(timeEvolveC / thermVals["γ"])

# ======== Simulations ========
## ---- Build arena ----
arena = BParts.buildRandArena(arenaParams["bounds"], arenaParams["n0"], arenaParams["radius"], arenaParams["speed"], verbose=true)

println("energy per cell: ", BParts.kineticEnergy(arena) / arenaParams["n0"])

## ---- Run simulations ----

println("running simulation...")
@time arena, pos_t_dim_id, vel_t_dim_id, ___, times_t = 
    BParts.randArenaEvolve(
            arena,
            arenaParams["evolveTime"], 
            arenaParams["timeStep"],
            arenaParams,
            growthParams;
            saveTimeStep=saveTimeStep,
            coldGrowth=growthParams["coldGrowth"],
            plotting=false,
            progress=true,
            verbose=true)
            
println("done")

## ===== Save data =====

filename = "simResultFixedDensity_"*string(ρIn)*".jld2"
save(filename, "arenaParams", arenaParams, "growthParams", growthParams, "pos_t_dim_id", pos_t_dim_id, "vel_t_dim_id", vel_t_dim_id, "times_t", times_t)
