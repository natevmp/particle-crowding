include("../src/bmparticles.jl")
include("../src/bmtheory.jl")
using .BParts
using .Theorist
using Distributions

##
using Plots
# plotly()
pyplot()

## ====== Fixed Density ======
# First we compare the particle simulation with the simulations of the Langevin equation (and some analytical results from it) for a fixed population size.

# ---- Parameters for simulation ----

arenaParams = Dict(
    "n0"=>400,
    "evolveTime"=>2500,
    # "evolveTime"=>1000,
    "bounds"=>((0.,20),(0.,20)), 
    "radius"=>0.08, 
    "speed"=>0.02,
    "timeStep"=> 0.3
)
growthParams = Dict(
    "ρ"=> 0,
    "k"=> 1000,
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
@time arena, posSim_t_dim_id, __, ___ = 
    BParts.randArenaEvolve(
            arena,
            arenaParams["evolveTime"], 
            arenaParams["timeStep"],
            arenaParams,
            growthParams;
            # saveTimeStep=arenaParams["timeStep"],
            coldGrowth=growthParams["coldGrowth"],
            plotting=false,
            progress=true,
            verbose=true)
            
println("done")

## ---- Langevin equation ----
# run an ensemble of Langevin equations
@time langevinEnsemble = Theorist.runLangevinSims(1000, arenaParams)



## ---- Analysis ----
# -- Speed distributions --
# For a well mixed population the particle speeds are predicted to be Rayleigh distributed with shape parameter $\sigma = \sqrt{E}$ where $E$ is the average particle energie.

evolveTime = arenaParams["evolveTime"]
speedPar_t_id = BParts.speedCalc(vel_t_dim_id)
speedLan_t_id = Theorist.speedCalc(langevinEnsemble, 1:evolveTime)
rDistPred = Rayleigh( sqrt(thermVals["E"]) )

p1 = histogram(vec(speedPar_t_id), bins=100, normalize=true, label="particle simulation")
histogram!(vec(speedLan_t_id), bins=100, normalize=true, label="Langevin simulation")
plot!(range(0,0.10, length=100), pdf.(rDistPred, range(0,0.10, length=100)), linestyle=:dash, linewidth=2, label="predicted")
xlabel!("speed")
ylabel!("distribution")
display(p1)

# -- Velocity autocorrelation
corrTime = 200
timesCorr_t, vCorrLan_t = Theorist.velCorrelation(langevinEnsemble, (1,corrTime), steps=corrTime)
vCorrPar_t = BParts.velocityAutocorrelation(vel_t_dim_id[1:corrTime, :, :])
vCorrPred_t = map( t->2*thermVals["E"]*exp(-thermVals["γ"]*t), 0:corrTime-1 )

p2 = plot(timesCorr_t, vCorrPar_t, label="particle simulation")
plot!(timesCorr_t, vCorrLan_t, label="Langevin simulation")
plot!(timesCorr_t, vCorrPred_t, label="analytical prediction")
xlabel!("time")
ylabel!("v autocorrelation")
display(p2)

# -- Mean squared displacement --
msdTime = 400
timesMSD_t, msdLan_t = Theorist.msd(langevinEnsemble, arenaParams,  (0, msdTime-1), msdTime)
msdPar_t = BParts.meanSquaredDisplacement(pos_t_dim_id[1:msdTime,:,:], arenaParams["bperiod"])

p2 = plot(timesMSD_t, msdLan_t, label="Langevin simulation")
plot!(timesMSD_t, msdPar_t, label="particle simulation")
xlabel!("time")
ylabel!("msd")
display(p2)

# The msd may appear sublinear after a particular distance due to the periodic boundaries. The effect can be removed by increasing the arena size.