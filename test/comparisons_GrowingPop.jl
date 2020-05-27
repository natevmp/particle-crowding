include("../src/bmparticles.jl")
include("../src/bmtheory.jl")
using .BParts
using .Theorist
using DifferentialEquations, Distributions
using Plots
# plotly()
gr()
# plotlyjs()
# pyplot()

# ====== Growing population ======

# ---- Parameters ----

function extendParams!(arenaParams::Dict)
    bounds = arenaParams["bounds"]
    arenaParams["volume"] = abs(bounds[1][2]-bounds[1][1])*abs(bounds[2][2]-bounds[2][1])
    arenaParams["bperiod"] = [abs(bounds[1][2]-bounds[1][1]), abs(bounds[2][2]-bounds[2][1])]
end

arenaParams = 
    Dict(
        "n0"=>100,
        "evolveTime"=>300,
        "bounds"=>((0.,20.),(0.,20.)), 
        "radius"=>0.08, 
        "speed"=>0.03
    )

growthParams =
    Dict(
        "ρ"=> 0.04,
        "k"=> 2000,
        "randGrowth"=> false
    )

extendParams!(arenaParams)

println(arenaParams)
println(growthParams)

# ---- Logistic Growth ----
# The number of cells in the system increases according to the logistic growth function

# -- Run Langevin simulations: --
langevinEnsemble = Theorist.runLangevinSims(1000, arenaParams, growthParams)

# -- Run a particle simulation: --
arena, pos_t_dim_id, vel_t_dim_id, cells_T_ID = BParts.randArenaEvolve(arenaParams["n0"],arenaParams["evolveTime"], arenaParams, growthParams)

# Plot number of cells over time
f1 = plot(0:250, Theorist.logisticGrowth.(0:250, growthParams["ρ"], growthParams["k"], arenaParams["n0"]), 
    label="expected", linewidth=2)
plot!(0:250, BParts.nCellsTime(cells_T_ID)[1:251], linestyle=:dash, 
    label="particle simulation", linewidth=2)
xlabel!("time")
ylabel!("expected number of cells")
display(f1)

# ---- Mean squared displacement ----
msdTime = 250
timesMSD_t, msdLan_t = Theorist.msd(langevinEnsemble, arenaParams,  (0, msdTime-1), msdTime);
msdPar_t = BParts.meanSquaredDisplacement(pos_t_dim_id[1:msdTime,:,:], arenaParams["bperiod"])

p2 = plot(timesMSD_t, msdLan_t, label="Langevin simulation", legend=:bottomright)
plot!(timesMSD_t, msdPar_t, label="particle simulation")
display(p2)




