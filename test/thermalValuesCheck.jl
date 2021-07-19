include("../src/bmparticles.jl")
include("../src/bmtheory.jl")
using .BParts
using .Theorist

using Plots
pyplot()

##
arenaParams = 
    Dict(
        "n0"=>10000,
        "evolveTime"=>3500,
        "bounds"=>((0.,32),(0.,32)), 
        "radius"=>0.08, 
        "speed"=>0.02,
        "timeStep"=> 0.4
    )
growthParams =
    Dict(
        "ρ"=> 0,
        "k"=> 500,
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
        legend=:bottomright)

xlabel!("time")
ylabel!("occupied surface density")
xlims!(0,arenaParams["evolveTime"]-growthParams["waitTime"])
display(f1)

# savefig(f1, "../Figures/rho"*string(growthParams["ρ"])*"surfaceDensity_time.pdf")

# savefig(f1, "densityGrowth.pdf")