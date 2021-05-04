include("../../src/bmtheory.jl")
using .Theorist

using Glob
using JLD
# using JLD2, FileIO
using LaTeXStrings
using Plots
# using UnicodeFun


## ========= Load data =========

res1 = load("data/FigData/msdAverageData_rho0.001.jld")
params1 = load("data/manyGP_02-16/Rho0010/growingPop_multSims_rho0.001_sim01params.jld")
res2 = load("data/FigData/msdAverageData_rho0.0015.jld")
params2 = load("data/manyGP_02-16/Rho0015/growingPop_multSims_rho0.0015_sim01params.jld")
res3 = load("data/FigData/msdAverageData_rho0.002.jld")
params3 = load("data/manyGP_02-16/Rho0020/growingPop_multSims_rho0.002_sim01params.jld")
res4 = load("data/FigData/msdAverageData_rho0.005.jld")
params4 = load("data/manyGP_04-06/growingPop_multSims_rho0.005_sim01params.jld")

##
function addBounds!(params)
    params["bounds"] = (params["boundsX"], params["boundsY"])
end

addBounds!(params1)
addBounds!(params2)
addBounds!(params3)
addBounds!(params4)
paramsArena = params1
thermalVals = Theorist.thermalValues(paramsArena)
tStepFric = 1/thermalVals["γ"]
xStep = paramsArena["speed"]*tStepFric
Theorist.extendParams!(paramsArena)
Theorist.extendParams!(params2)
Theorist.extendParams!(params3)
Theorist.extendParams!(params4)

##
volumeDens(n, r, v) = n*π*r^2/v

function getParticleDensTime(params)
    Theorist.logisticGrowth.(0:(params["evolveTime"]-params["waitTime"]-1), params["ρ"], params["k"], params["n0"]) / params["volume"]
end

function getVolumeTime(params)
    volumeDens.(
        getParticleDensTime(params)*params["volume"],
        params["radius"],
        params["volume"]
    )
end

volume1_t = getVolumeTime(params1)
volume2_t = getVolumeTime(params2)
volume3_t = getVolumeTime(params3)
volume4_t = getVolumeTime(params4)

density1_t = getParticleDensTime(params1)
density2_t = getParticleDensTime(params2)
density3_t = getParticleDensTime(params3)
density4_t = getParticleDensTime(params4)

##

colUnits = Dict(
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

for res in [res1, res2, res3, res4]
    push!(msdLan_Sim_t, res["msdLan_t"])
    push!(msdPart_Sim_t, res["msdPart_t"])
    push!(timesPar_Sim_, res["timesPar_"])
    push!(fricTimes_Sim_, res["fricTimes_"])
    push!(ρ_Sim, res["ρ"])
end
ρUnits_Sim = ρ_Sim / colUnits[:ρ]


##
pyplot()
default(palette = palette(:seaborn_deep)) 

## ========= Plot MSD Collisional Units ==========
fscale = 0.9
pyplot(tickfontsize=10, guidefontsize=13, legendfontsize=12, size=(fscale*600,fscale*400))

fig1 = plot(1:10, -1*ones(5),
    color=:grey45,
    label="particle simulation"
)
xlims!(0,40)
ylims!(0,55)

plot!(1:10, -1*ones(5),
    color=:black,
    linestyle=:dash,
    label="Langevin simulation"
)

for i in 1:3
    plot!(timesPar_Sim_[i]/colUnits[:t], msdPart_Sim_t[i][2:end]/colUnits[:x], 
        label="",
        legend=:bottomright,
        linewidth=1.5,
        palette=:seaborn_pastel,
        color=i
    )
    plot!(timesPar_Sim_[i]/colUnits[:t], msdLan_Sim_t[i][2:end]/colUnits[:x],
        label="",
        linestyle=:dash,
        palette=:seaborn_deep,
        color=i,
        linewidth=1.5
    )
end
annotate!(
    [
        ( 20, 45, text(latexstring("\\rho = "*string(round(ρUnits_Sim[1],digits=2))), 10) ),
        ( 30, 40, text(latexstring("\\rho = "*string(round(ρUnits_Sim[2],digits=2))), 10) ),
        ( 35, 30, text(latexstring("\\rho = "*string(round(ρUnits_Sim[3],digits=2))), 10) ),
        # ( 36, 25, text(latexstring("\\rho = "*string(round(ρUnits_Sim[4],digits=2))), 10) ),
    ]
)
xlabel!(L"t \quad [1/\gamma_0]")
ylabel!(L"\left\langle x(t)^2 \right\rangle \quad [s_0 / \gamma_0]")

# title!(latexstring("\\rho="*string(paramsGrowth["ρ"])))

display(fig1)
figname = "particles"*string(params1["k"])*"_MultRho_msd_FrictionUnits.pdf"
savefig(fig1, "./figures/2021/"*figname)

## ========= Plot Growth Collisional Units =========


fig2 = plot(legend=:bottomright)
for (i,s_t) in enumerate([volume1_t, volume2_t, volume3_t])
    _t = (0:length(s_t)-1)/colUnits[:t]
    plot!(_t, s_t,
        label=latexstring("\\rho = "*string(round(ρUnits_Sim[i],digits=2)))
    )
end
xlabel!(L"t")
ylabel!("particle surface density")
display(fig2)

##

friction1_t = Theorist.frictionFromParticleDensity.(density1_t, thermalVals["σc"], thermalVals["E"])
friction2_t = Theorist.frictionFromParticleDensity.(density2_t, thermalVals["σc"], thermalVals["E"])
friction3_t = Theorist.frictionFromParticleDensity.(density3_t, thermalVals["σc"], thermalVals["E"])
friction4_t = Theorist.frictionFromParticleDensity.(density4_t, thermalVals["σc"], thermalVals["E"])
fig3 = plot(legend=:topright)
for (i,f_t) in enumerate([friction1_t, friction2_t, friction3_t])
    _t = (0:length(f_t)-1)/colUnits[:t]
    plot!(_t, 1 ./f_t /colUnits[:t],
        label=latexstring("\\rho = "*string(round(ρUnits_Sim[i],digits=2)))
    )
end
xlabel!(L"t")
ylabel!(L"1 / \gamma")
ylims!(0, 1)
display(fig3)

diffCoeff1_t = Theorist.DiffCoeffFromParticleDensity.(density1_t, thermalVals["σc"], thermalVals["E"])
diffCoeff2_t = Theorist.DiffCoeffFromParticleDensity.(density2_t, thermalVals["σc"], thermalVals["E"])
diffCoeff3_t = Theorist.DiffCoeffFromParticleDensity.(density3_t, thermalVals["σc"], thermalVals["E"])
diffCoeff4_t = Theorist.DiffCoeffFromParticleDensity.(density4_t, thermalVals["σc"], thermalVals["E"])
fig4 = plot(legend=:bottomright)
for (i,d_t) in enumerate([diffCoeff1_t, diffCoeff2_t, diffCoeff3_t])
    _t = (0:length(d_t)-1)/colUnits[:t]
    plot!(_t, d_t*colUnits[:t]^2,
        label=latexstring("\\rho = "*string(round(ρUnits_Sim[i],digits=2)))
    )
end
xlabel!(L"t")
ylabel!(L"D")
display(fig4)

##
fullscale = 1.5
fig5 = plot(fig2, fig4, fig3, fig1, layout=4, size=(fullscale*600,fullscale*400))
figname = "particles"*string(params1["k"])*"_MultRho_complete_frictionUnits.pdf"
savefig(fig5, "./figures/2021/"*figname)
display(fig5)

## ===== Plot MSD in simulation units =====

fscale = 0.9

fig1b = plot(1:10, -1*ones(5),
    color=:grey45,
    label="particle simulation"
)
xlims!(0,5000)
ylims!(0,150)

plot!(1:10, -1*ones(5),
    color=:black,
    linestyle=:dash,
    label="Langevin simulation"
)

for i in 1:3
    plot!(timesPar_Sim_[i], msdPart_Sim_t[i][2:end], 
        label="",
        legend=:bottomright,
        linewidth=1.5,
        palette=:seaborn_pastel,
        color=i
    )
    plot!(timesPar_Sim_[i], msdLan_Sim_t[i][2:end],
        label="",
        linestyle=:dash,
        palette=:seaborn_deep,
        color=i,
        linewidth=1.5
    )
end
annotate!(
    [
        ( 2000, 105, text(latexstring("\\rho = "*string(round(ρ_Sim[1],digits=5))), 10) ),
        ( 3000, 93, text(latexstring("\\rho = "*string(round(ρ_Sim[2],digits=5))), 10) ),
        ( 4000, 70, text(latexstring("\\rho = "*string(round(ρ_Sim[3],digits=5))), 10) ),
        # ( 36, 25, text(latexstring("\\rho = "*string(round(ρUnits_Sim[4],digits=2))), 10) ),
    ]
)
xlabel!(L"t \quad")
ylabel!(L"\left\langle x(t)^2 \right\rangle \quad")

display(fig1b)
figname = "particles"*string(params1["k"])*"_MultRho_msd_simulationUnits.pdf"
savefig(fig1b, "./figures/2021/"*figname)

## ========= Plot Growth Simulation Units =========


fig2b = plot(legend=:bottomright)
for (i,s_t) in enumerate([volume1_t, volume2_t, volume3_t])
    _t = (0:length(s_t)-1)
    plot!(_t, s_t,
        label=latexstring("\\rho = "*string(round(ρ_Sim[i],digits=5)))
    )
end
xlabel!(L"t")
ylabel!("particle surface density")
display(fig2b)

##

friction1_t = Theorist.frictionFromParticleDensity.(density1_t, thermalVals["σc"], thermalVals["E"])
friction2_t = Theorist.frictionFromParticleDensity.(density2_t, thermalVals["σc"], thermalVals["E"])
friction3_t = Theorist.frictionFromParticleDensity.(density3_t, thermalVals["σc"], thermalVals["E"])
friction4_t = Theorist.frictionFromParticleDensity.(density4_t, thermalVals["σc"], thermalVals["E"])
fig3b = plot(legend=:topright)
for (i,f_t) in enumerate([friction1_t, friction2_t, friction3_t])
    _t = (0:length(f_t)-1)
    plot!(_t, 1 ./f_t,
        label=latexstring("\\rho = "*string(round(ρ_Sim[i],digits=5)))
    )
end
xlabel!(L"t")
ylabel!(L"1 / \gamma")
ylims!(0, 130)
display(fig3b)

##
diffCoeff1_t = Theorist.DiffCoeffFromParticleDensity.(density1_t, thermalVals["σc"], thermalVals["E"])
diffCoeff2_t = Theorist.DiffCoeffFromParticleDensity.(density2_t, thermalVals["σc"], thermalVals["E"])
diffCoeff3_t = Theorist.DiffCoeffFromParticleDensity.(density3_t, thermalVals["σc"], thermalVals["E"])
diffCoeff4_t = Theorist.DiffCoeffFromParticleDensity.(density4_t, thermalVals["σc"], thermalVals["E"])
fig4b = plot(legend=:bottomright)
for (i,d_t) in enumerate([diffCoeff1_t, diffCoeff2_t, diffCoeff3_t])
    _t = (0:length(d_t)-1)
    plot!(_t, d_t,
        label=latexstring("\\rho = "*string(round(ρ_Sim[i],digits=5)))
    )
end
xlabel!(L"t")
ylabel!(L"D")
display(fig4b)

##
fullscale = 1.5
fig5b = plot(fig2b, fig4b, fig3b, fig1b, layout=4, size=(fullscale*600,fullscale*400))
figname = "particles"*string(params1["k"])*"_MultRho_complete_simulationUnits.pdf"
savefig(fig5b, "./figures/2021/"*figname)
display(fig5b)

