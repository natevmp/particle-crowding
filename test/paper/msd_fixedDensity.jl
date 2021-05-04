include("../../src/bmparticles.jl")
include("../../src/bmtheory.jl")
using .BParts
using .Theorist
using Distributions
using JLD2, FileIO

##
using Plots
# plotly()
pyplot()

## ====== Fixed Density ======
# First we compare the particle simulation with the simulations of the Langevin equation (and some analytical results from it) for a fixed population size.

# ---- Parameters for simulation ----


ρIn = 0.2

arenaParams = Dict(
    "evolveTime"=>2500,
    "n0"=> 1000,
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

growthParams = Dict(
    "ρ"=> 0,
    "k"=> arenaParams["n0"],
    "randGrowth"=> false,
    "coldGrowth"=> false,
    "waitTime"=> 500
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

# ======== Simulations ========
## ---- Build arena ----
arena = BParts.buildRandArena(arenaParams["bounds"], arenaParams["n0"], arenaParams["radius"], arenaParams["speed"], verbose=true)

println("energy per cell: ", BParts.kineticEnergy(arena) / arenaParams["n0"])

## ---- Run simulations ----

println("running simulation...")
@time arena, pos_t_dim_id, vel_t_dim_id, ___, ___ = 
    BParts.randArenaEvolve(
            arena,
            arenaParams["evolveTime"], 
            arenaParams["timeStep"],
            arenaParams,
            growthParams;
            coldGrowth=growthParams["coldGrowth"],
            plotting=false,
            progress=true,
            verbose=true)
            
println("done")

## ---- Langevin equation ----
# run an ensemble of Langevin equations
# @time langevinEnsemble = Theorist.runLangevinSims(5000, arenaParams, growthParams)



## ---- Analysis ----

# # -- Velocity autocorrelation
# corrTime = 10* 1/thermVals["γ"]
# # tParIn = growthParams["waitTime"]
# tParIn = 1000
# # _tCorr = growthParams["waitTime"]:growthParams["waitTime"]+corrTime
# _tCorrInd = range(tParIn, Int(floor(tParIn+corrTime-1)), step=1)
# _tCorr = 0:corrTime-1
# # timesCorr_t, vCorrLan_t = Theorist.velCorrelation(langevinEnsemble, (_tCorr[1],_tCorr[end]), steps=corrTime)
# vCorrPar_t = BParts.velocityAutocorrelation(vel_t_dim_id[_tCorrInd, :, :])
# vCorrPred_t = map( t->2*thermVals["E"]*exp(-thermVals["γ"]*t), _tCorr )

# p2 = plot(_tCorr, vCorrPar_t, label="particle simulation")
# # plot!(_tCorr, vCorrLan_t, label="Langevin simulation")
# plot!(_tCorr, vCorrPred_t, label="analytical prediction")
# xlabel!("time")
# ylabel!("v autocorrelation")
# # xlims!(0, )
# display(p2)



## ===== Collision angle distribution =====

# function getAngle(v2::AbstractVector, v1::AbstractVector)
#     θ1 = atan(v1[2], v1[1])
#     θ2 = atan(v2[2], v2[1])
#     atan( sin(θ2-θ1), cos(θ2-θ1) )
# end

# function particleAngleChanges(vel_t_dim)
#     dθ_ = Float64[]
#     for tInd in 2:size(vel_t_dim, 1)
#         if vel_t_dim[tInd, 1] != vel_t_dim[tInd-1, 1]
#             dθ = getAngle(vel_t_dim[tInd, :], vel_t_dim[tInd-1, :])
#             push!(dθ_, dθ)
#         end
#     end
#     return dθ_
# end

# function getAngleChanges(vel_t_dim_id)
#     dθ_ = Float64[]
#     for cid in 1:size(vel_t_dim_id, 3)
#         dθ_cid = particleAngleChanges(@view vel_t_dim_id[:,:,cid])
#         append!(dθ_, dθ_cid)
#     end
#     return dθ_
# end

# # particleAngleChanges(vel_t_dim_id[:,:,1])
# dθ_ = getAngleChanges(vel_t_dim_id)

# ##

# figθ = histogram(dθ_)
# display(figθ)







## ===== MSD =====
# _tMsd = growthParams["waitTime"]+1:arenaParams["evolveTime"]
# msdPar_t = BParts.meanSquaredDisplacement(pos_t_dim_id, (_tMsd[1],_tMsd[end]))
# timesMSD_t, msdLan_t = Theorist.msd(langevinEnsemble, arenaParams, (_tMsd[1],_tMsd[end]))
# _tMsd = _tMsd .- (growthParams["waitTime"]+1)

## ---- Theory ----
# msdTheory_t = Theorist.msdTheory.(_tMsd, arenaParams["n0"]/arenaParams["volume"], thermVals["σc"], thermVals["E"])

##
# p2 = plot(_tMsd, msdLan_t, label="Langevin simulation", legend=:bottomright)
# plot!(_tMsd, msdPar_t, label="particle simulation")
# plot!(_tMsd, msdTheory_t, label="theoretical prediction")
# xlabel!("time")
# ylabel!("msd")
# display(p2)


## ===== Save data =====
# filename = "./data/FixedDensity/simResultFixedDensity_"*string(ρIn)*"_msd.jld2"
# save(filename, "arenaParams", arenaParams, "growthParams", growthParams, "_t", _tMsd, "msdLan_t", msdLan_t, "msdPar_t", msdPar_t)


filename = "./data/FixedDensity/simResultFixedDensity_"*string(ρIn)*".jld2"
save(filename, "arenaParams", arenaParams, "growthParams", growthParams, "pos_t_dim_id", pos_t_dim_id, "vel_t_dim_id", vel_t_dim_id)
