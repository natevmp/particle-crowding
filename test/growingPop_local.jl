
using JLD2
include("../src/bmparticles.jl")
using .BParts
include("../src/bmtheory.jl")
using .Theorist

##
# simNumber = parse(Int, ARGS[1])
# rho = parse(Float64, ARGS[2])
# time = parse(Int64, ARGS[3])
# simNumber = 1
# rho = 0.002
# display(simNumber)
# display(rho)

LOADARENA = false

## ===== Load / create Arena =====
if LOADARENA
    # @load "data/arenaInit_n500.jld2" arena arenaParams growthParams
    @load "data/arenaInit_n500preError.jld2" arena arenaParams growthParams
else
    arenaParams = Dict(
            "n0"=>100,
            "evolveTime"=>100,
            # "evolveTime"=>1000,
            "bounds"=>((0.,5),(0.,5)), 
            "radius"=>0.08, 
            "speed"=>0.02,
            "timeStep"=> 0.3
        )
    growthParams = Dict(
            "ρ"=> 0.002,
            "k"=> 100,
            "randGrowth"=> false,
            "coldGrowth"=> false,
            "waitTime"=> 500
        )
    BParts.extendParams!(arenaParams)
    arena = BParts.buildRandArena(arenaParams["bounds"], arenaParams["n0"], arenaParams["radius"], arenaParams["speed"], verbose=true)
    @save "data/arenaInit_n"*string(arenaParams["n0"])*".jld2" arena arenaParams growthParams
end
println(arenaParams)
println(growthParams)

## ==== Run simulations ====

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
#

# @save "data/arenaInit_n"*string(arenaParams["n0"])*"preError.jld2" arena arenaParams growthParams

# ## ==== Get mean squared displacements ====
# msdTimes = (growthParams["waitTime"]+1, arenaParams["evolveTime"])
# msdPart_t = BParts.meanSquaredDisplacement(posSim_t_dim_id, msdTimes)

# ## Run Langevin simulations
# function extendParams!(arenaParams::Dict)
#     bounds = arenaParams["bounds"]
#     arenaParams["volume"] = abs(bounds[1][2]-bounds[1][1])*abs(bounds[2][2]-bounds[2][1])
#     arenaParams["bperiod"] = [abs(bounds[1][2]-bounds[1][1]), abs(bounds[2][2]-bounds[2][1])]
# end
# extendParams!(arenaParams)

# @time langevinEnsemble = Theorist.runLangevinSims(5000, arenaParams, growthParams)

# ## Calculate msd for Langevin simulations
# @time timesLan_, msdLan_t = 
#     Theorist.msd(
#         langevinEnsemble, 
#         arenaParams,
#         (growthParams["waitTime"]+1, arenaParams["evolveTime"])
#     )

# ##
# thermalVals = Theorist.thermalValues(arenaParams)
# fricTstep = 1/thermalVals["γ"]

# # fricTimes_ = 0:fricTstep:timesMSD_[end]
# fricTimes_ = range( 0, step=5*fricTstep, stop=(msdTimes[2]-msdTimes[1]) )

# ##
# using Plots
# pyplot()
# ##
# using LaTeXStrings

# ##
# timesPar_ = 1:msdTimes[2]-msdTimes[1]
# p2 = plot(timesPar_, msdPart_t[2:end], label="particle simulation", legend=:bottomright, linewidth=2,
#     size=(400,300),
#     xticks = (fricTimes_, 5*(0:length(fricTimes_)-1)),
#     dpi=140)
# plot!(timesPar_, msdLan_t[2:end], label="Langevin simulation", linestyle=:dash, linewidth=2)
# xlabel!(L"t \quad [1/\gamma_0]")
# ylabel!(L"\left\langle x(t)^2 \right\rangle")
# title!(latexstring("\\rho="*string(growthParams["ρ"])))
# display(p2)

