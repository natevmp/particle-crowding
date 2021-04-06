include("../src/bmparticles.jl")
include("../src/bmtheory.jl")

using .BParts
using .Theorist

using Glob
using JLD
# using JLD2, FileIO
using LaTeXStrings, Statistics
using Plots 
pyplot()

##
# files_f = glob("growingPop_multSims_rho0*.jld2", "./data/many_gp_50613208")
# filesSims_f = glob(r"growingPop_multSims_rho0*/[0-9]/.jld", "./data/manyGP_02-16")
# filesParams_f = glob("growingPop_multSims_*params.jld", "./data/manyGP_04-06")
# filesSims_f = glob(glob"growingPop_multSims_*[0-9].jld", "./data/manyGP_04-06")
filesParams_f = glob("growingPop_multSims_*params.jld", "./data/manyGP_04-06")
filesSims_f = glob(glob"growingPop_multSims_*[0-9].jld", "./data/manyGP_04-06")

##
params = load(filesParams_f[1])
paramsArena = params
paramsGrowth = params
paramsArena["bounds"] = (params["boundsX"], params["boundsY"])

##
msdPart_Sim_t = Vector{Vector{Float64}}(undef, 0)
for filename in filesSims_f
    fileData = load(filename)
    push!(msdPart_Sim_t, fileData["msdPart_t"])
end

## Run Langevin simulations
function extendParams!(arenaParams::Dict)
    bounds = arenaParams["bounds"]
    arenaParams["volume"] = abs(bounds[1][2]-bounds[1][1])*abs(bounds[2][2]-bounds[2][1])
    arenaParams["bperiod"] = [abs(bounds[1][2]-bounds[1][1]), abs(bounds[2][2]-bounds[2][1])]
end
extendParams!(paramsArena)

@time langevinEnsemble = Theorist.runLangevinSims(5000, paramsArena, paramsGrowth)

# Calculate msd for Langevin simulations
@time timesLan_, msdLan_t = 
    Theorist.msd(
        langevinEnsemble, 
        paramsArena,
        (paramsGrowth["waitTime"]+1, paramsArena["evolveTime"])
    )

## Average msd over all simulations

msdTimes = (paramsGrowth["waitTime"]+1, paramsArena["evolveTime"])

msdPart_t = Array{Float64, 1}(undef, msdTimes[2]-msdTimes[1]+1)
for i in 1:(msdTimes[2]-msdTimes[1]+1)
    msdPart_t[i] = mean([ msd_tt[i] for msd_tt in msdPart_Sim_t ])
end

##

thermalVals = Theorist.thermalValues(paramsArena)
fricTstep = 1/thermalVals["γ"]

# fricTimes_ = 0:fricTstep:timesMSD_[end]
fricTimes_ = range( 0, step=5*fricTstep, stop=(msdTimes[2]-msdTimes[1]) )
##

## Save msd data
# saveName = "data/FigData/msdAverageData_rho"*string(paramsGrowth["ρ"])*".jld"
# save(saveName, "timesLan_", timesLan_, "msdLan_t", msdLan_t, "msdPart_t", msdPart_t, "fricTimes_", fricTimes_, "timesPar_", timesPar_, "ρ", paramsGrowth["ρ"])

##

timesPar_ = 1:msdTimes[2]-msdTimes[1]
p2 = plot(timesPar_, msdPart_t[2:end], label="particle simulation", legend=:bottomright, linewidth=2,
    size=(400,300),
    xticks = (fricTimes_, 5*(0:length(fricTimes_)-1)),
    dpi=140)
plot!(timesPar_, msdLan_t[2:end], label="Langevin simulation", linestyle=:dash, linewidth=2)
xlabel!(L"t \quad [1/\gamma_0]")
ylabel!(L"\left\langle x(t)^2 \right\rangle")
title!(latexstring("\\rho="*string(paramsGrowth["ρ"])))
display(p2)

# savefig(p2, "./figures/2021/rho"*string(paramsGrowth["ρ"])*"msd_time.png")
# savefig(p2, "./figures/2021/particles"*string(paramsGrowth["k"])*"_rho"*string(paramsGrowth["ρ"])*"msd_time.pdf")

