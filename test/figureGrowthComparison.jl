include("../src/bmtheory.jl")
using .Theorist

using Glob
using JLD
# using JLD2, FileIO
using LaTeXStrings
using Plots
# using UnicodeFun

##
res1 = load("data/FigData/msdAverageData_rho0.001.jld")
params1 = load("data/manyGP_02-16/Rho0010/growingPop_multSims_rho0.001_sim01params.jld")
res2 = load("data/FigData/msdAverageData_rho0.0015.jld")
params2 = load("data/manyGP_02-16/Rho0015/growingPop_multSims_rho0.0015_sim01params.jld")
res3 = load("data/FigData/msdAverageData_rho0.002.jld")
params3 = load("data/manyGP_02-16/Rho0020/growingPop_multSims_rho0.002_sim01params.jld")

##
paramsArena = params1
paramsGrowth = params1
paramsArena["bounds"] = (paramsArena["boundsX"], paramsArena["boundsY"])
thermalVals = Theorist.thermalValues(paramsArena)
tStepFric = 1/thermalVals["γ"]
xStep = paramsArena["speed"]*tStepFric
##

units = Dict(
    :t => 1/thermalVals["γ"],
    :x => paramsArena["speed"]/thermalVals["γ"],
    :ρ => thermalVals["γ"]
)
##
msdLan_Sim_t = []
msdPart_Sim_t = []
timesPar_Sim_ = []
fricTimes_Sim_ = []
ρ_Sim = []

for res in [res1, res2, res3]
    push!(msdLan_Sim_t, res["msdLan_t"])
    push!(msdPart_Sim_t, res["msdPart_t"])
    push!(timesPar_Sim_, res["timesPar_"])
    push!(fricTimes_Sim_, res["fricTimes_"])
    push!(ρ_Sim, res["ρ"])
end
ρUnits_Sim = ρ_Sim / units[:ρ]

##
wfscale = 0.9
pyplot(tickfontsize=10, guidefontsize=13, legendfontsize=12, size=(fscale*600,fscale*400))

fig = plot(timesPar_Sim_[1]/units[:t], msdPart_Sim_t[1][2:end]/units[:x], 
    label="particle simulation",
    legend=:bottomright,
    color=1,
    linewidth=1.5
    )
plot!(timesPar_Sim_[1]/units[:t], msdLan_Sim_t[1][2:end]/units[:x],
    label="Langevin simulation",
    linestyle=:dash,
    color=2,
    linewidth=1.5
)

for i in 2:3
    plot!(timesPar_Sim_[i]/units[:t], msdPart_Sim_t[i][2:end]/units[:x], 
        label="",
        legend=:bottomright,
        linewidth=1.5,
        color=1
    )
    plot!(timesPar_Sim_[i]/units[:t], msdLan_Sim_t[i][2:end]/units[:x],
        label="",
        linestyle=:dash,
        color=2,
        linewidth=1.5
    )
end
annotate!(
    [
        ( 20, 45, text(latexstring("\\rho = "*string(round(ρUnits_Sim[1],digits=2))), 10) ),
        ( 30, 40, text(latexstring("\\rho = "*string(round(ρUnits_Sim[2],digits=2))), 10) ),
        ( 35, 30, text(latexstring("\\rho = "*string(round(ρUnits_Sim[3],digits=2))), 10) ),
    ]
)
xlabel!(L"t \quad [1/\gamma_0]")
ylabel!(L"\left\langle x(t)^2 \right\rangle \quad [s_0 / \gamma_0]")

# title!(latexstring("\\rho="*string(paramsGrowth["ρ"])))

display(fig)
figname = "particles"*string(paramsGrowth["k"])*"_MultRho_msd.pdf"
savefig(fig, "./figures/2021/"*figname)
