include("../src/bmparticles.jl")
include("../src/bmtheory.jl")
using .BParts
using .Theorist
using DifferentialEquations, Distributions

using Plots
# plotly()
gr()

# ====== Fixed Density ======
# First we compare the particle simulation with the simulations of the Langevin equation (and some analytical results from it) for a fixed population size.

# ---- Parameters for simulation ----
function extendParams!(arenaParams::Dict)
    bounds = arenaParams["bounds"]
    arenaParams["volume"] = abs(bounds[1][2]-bounds[1][1])*abs(bounds[2][2]-bounds[2][1])
    arenaParams["bperiod"] = [abs(bounds[1][2]-bounds[1][1]), abs(bounds[2][2]-bounds[2][1])]
end
arenaParams = 
    Dict(
        "n0"=>1600,
        "evolveTime"=> 2000,
        "bounds"=>((0.,40.),(0.,40.)), 
        "radius"=>0.08, 
        "speed"=>0.03
    )

extendParams!(arenaParams)
for p in arenaParams
    println(p)
end

thermVals = Theorist.thermalValues(arenaParams)

# ---- Simulations ----

# -- Particle simulation --
arena, pos_t_dim_id, vel_t_dim_id, cells_T_ID = BParts.randArenaEvolve(arenaParams["n0"], arenaParams["evolveTime"], arenaParams)

println("energy per cell: ", BParts.kineticEnergy(arena) / arenaParams["n0"])

# -- Langevin equation --
# run an ensemble of Langevin equations
langevinEnsemble = Theorist.runLangevinSims(1000, arenaParams)

# ---- Analysis ----
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

# While the particle simulation agrees nicely with the predicted Rayleigh distribution, the Langevin simulations appear to deviate slightly. I haven't quite figured out what is causing this. Especially since technically the predicted rayleigh is derived from the Langevin equation...

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