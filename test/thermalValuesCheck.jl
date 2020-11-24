using JLD2, Statistics
include("../src/bmparticles.jl")
include("../src/bmtheory.jl")
using .BParts
using .Theorist

using Plots
gr()

arenaParams = 
    Dict(
        "n0"=>500,
        "evolveTime"=>5500,
        "bounds"=>((0.,22.4),(0.,22.4)), 
        "radius"=>0.08, 
        "speed"=>0.02,
        "timeStep"=> 0.1
    )
growthParams =
    Dict(
        "ρ"=> 0.002,
        "k"=> 5000,
        "randGrowth"=> false,
        "waitTime"=> 500
    )

volumeDens(n, r, v) = n*π*r^2/v

BParts.extendParams!(arenaParams)
thermVals = Theorist.thermalValues(arenaParams)

println("mean free time: ", 1/thermVals["γ"])
println("initial volume density: ", volumeDens(arenaParams["n0"], arenaParams["radius"], arenaParams["volume"]))
println("final volume density: ", volumeDens(growthParams["k"], arenaParams["radius"], arenaParams["volume"]))

##
f1 = plot(0:(arenaParams["evolveTime"]-1-growthParams["waitTime"]), 
        volumeDens.(
            Theorist.logisticGrowth.(0:(arenaParams["evolveTime"]-growthParams["waitTime"]-1), growthParams["ρ"], growthParams["k"], arenaParams["n0"]),
            arenaParams["radius"],
            arenaParams["volume"]
        ),
        label="growth function",
        linewidth=2,
        size=(450,300),
        legend=:topleft)

xlabel!("time")
ylabel!("occupied surface density")
xlims!(0,arenaParams["evolveTime"]-growthParams["waitTime"])
display(f1)

# savefig(f1, "../Figures/rho"*string(growthParams["ρ"])*"surfaceDensity_time.pdf")

savefig(f1, "densityGrowth.pdf")