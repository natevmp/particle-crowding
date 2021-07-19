include("../../src/bmtheory.jl")
using .Theorist
using ElasticArrays
using Glob
using JLD2
# using JLD2, FileIO
using LaTeXStrings
using Plots
# using UnicodeFun

SAVEFIGS = true

## ========= Load data =========

filenames_ = glob("growingPop*jld2", "data/GrowingPop/21-06-25")

_gr = Float64[]
msdPart_gr_T = Vector{Float64}[]
_gr_T = Vector{Float64}[]
params_gr = Dict[]
for fname in filenames_
    @load fname arenaParams growthParams msdPart_t times_t
    push!(_gr, growthParams["ρ"])
    push!(msdPart_gr_T, msdPart_t)
    push!(_gr_T, times_t)
    global paramsArena = arenaParams
    push!(params_gr, merge(arenaParams, growthParams))
end

##
thermalVals = Theorist.thermalValues(paramsArena)
tStepFric = 1/thermalVals["γ"]
xStep = paramsArena["speed"]*tStepFric

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

volume_gr_T = getVolumeTime.(params_gr)
density_gr_T = getParticleDensTime.(params_gr)

##

colUnits = Dict(
    :t => 1/thermalVals["γ"],
    :x => paramsArena["speed"]/thermalVals["γ"],
    :ρ => thermalVals["γ"]
)
##
# Calculate msd for Langevin simulations

msdLan_gr_T = Vector{Float64}[]
_gr_TLan = Vector{Float64}[]
for i in 1:length(_gr)
    @time langevinEnsemble = Theorist.runLangevinSims(5000, paramsArena, params_gr[i])
    @time _tLan, msdLan_t = 
        Theorist.msd(
            langevinEnsemble, 
            paramsArena,
            (params_gr[i]["waitTime"]+1, paramsArena["evolveTime"])
        )
    push!(msdLan_gr_T, msdLan_t)
    push!(_gr_TLan, _tLan)
end




##
pyplot()
default(palette = palette(:seaborn_deep)) 

## ========= Plot MSD Collisional Units ==========
fscale = 0.9
# pyplot(tickfontsize=10, guidefontsize=13, legendfontsize=12, size=(fscale*600,fscale*400))
pyplot()
theme(:bright, size=(0.8*600,0.8*400), minorgrid=false, gridstyle=:dash, fontfamily="DejaVu Sans", legendfontsize=10, palette=:seaborn_deep,)

# cPalette = cgrad(:tableau_orange_blue, 4, categorical=true, rev=true)


##
_t = _gr_T[1][params_gr[1]["waitTime"]+1:end] .- (params_gr[1]["waitTime"])

##

plot(msdPart_gr_T[1])

##
fig1 = plot(1:10, -1*ones(5),
    color=:grey45,
    label="particle simulation"
)
# xlims!(0,40)
# ylims!(0,55)

plot!(1:10, -1*ones(5),
    color=:black,
    linestyle=:dash,
    label="Langevin simulation"
)

for i in 1:3
    plot!(_t[2:end]/colUnits[:t], msdPart_gr_T[i]/colUnits[:x], 
        label="",
        legend=:bottomright,
        linewidth=1.5,
        # color=cPalette[i],
        palette=:seaborn_pastel,
        color=i,
    )
    plot!(_t[2:end]/colUnits[:t], msdLan_gr_T[i]/colUnits[:x],
        label="",
        linestyle=:dash,
        # color=cPalette[i],
        palette=:seaborn_deep,
        color=i,
        linewidth=1.5
    )
end
display(fig1)
##

annotate!(
    [
        ( 20, 45, text(latexstring("\\lambda = "*string(round(ρUnits_Sim[1],digits=2))), 10) ),
        ( 30, 40, text(latexstring("\\lambda = "*string(round(ρUnits_Sim[2],digits=2))), 10) ),
        ( 35, 30, text(latexstring("\\lambda = "*string(round(ρUnits_Sim[3],digits=2))), 10) ),
        # ( 36, 25, text(latexstring("\\rho = "*string(round(ρUnits_Sim[4],digits=2))), 10) ),
    ]
)
xlabel!(L"t \quad [1/\gamma_0]")
ylabel!(L"\left\langle x(t)^2 \right\rangle \quad [s_0 / \gamma_0]")

# title!(latexstring("\\rho="*string(paramsGrowth["ρ"])))

display(fig1)
figname = "particles"*string(params1["k"])*"_MultRho_msd_FrictionUnits.pdf"
SAVEFIGS && savefig(fig1, "./figures/2021/"*figname)

## ========= Plot Growth Collisional Units =========


fig2 = plot(legend=:bottomright)
for (i,s_t) in enumerate([volume1_t, volume2_t, volume3_t])
    _t = (0:length(s_t)-1)/colUnits[:t]
    plot!(_t, s_t,
        label=latexstring("\\lambda = "*string(round(ρUnits_Sim[i],digits=2))),
        palette=:seaborn_deep,
        color=i,
    )
end
xlabel!(L"t \quad [1/\gamma_0]")
ylabel!(L"ρ")
display(fig2)

## ----- Friction -----

friction1_t = Theorist.frictionFromParticleDensity.(density1_t, thermalVals["σc"], thermalVals["E"])
friction2_t = Theorist.frictionFromParticleDensity.(density2_t, thermalVals["σc"], thermalVals["E"])
friction3_t = Theorist.frictionFromParticleDensity.(density3_t, thermalVals["σc"], thermalVals["E"])
friction4_t = Theorist.frictionFromParticleDensity.(density4_t, thermalVals["σc"], thermalVals["E"])
fig3 = plot(legend=:topright)
for (i,f_t) in enumerate([friction1_t, friction2_t, friction3_t])
    _t = (0:length(f_t)-1)/colUnits[:t]
    plot!(_t, 1 ./f_t /colUnits[:t],
        label=latexstring("\\lambda = "*string(round(ρUnits_Sim[i],digits=2)))
    )
end
xlabel!(L"t \quad [1/\gamma_0]")
ylabel!(L"1 / \gamma(t) \quad [1/\gamma_0]")
ylims!(0, 1)
display(fig3)

## ----- diffusion coefficient -----
diffCoeff1_t = Theorist.diffCoeffFromParticleDensity.(density1_t, thermalVals["σc"], thermalVals["E"])
diffCoeff2_t = Theorist.diffCoeffFromParticleDensity.(density2_t, thermalVals["σc"], thermalVals["E"])
diffCoeff3_t = Theorist.diffCoeffFromParticleDensity.(density3_t, thermalVals["σc"], thermalVals["E"])
diffCoeff4_t = Theorist.diffCoeffFromParticleDensity.(density4_t, thermalVals["σc"], thermalVals["E"])
fig4 = plot(legend=:bottomright)
for (i,d_t) in enumerate([diffCoeff1_t, diffCoeff2_t, diffCoeff3_t])
    _t = (0:length(d_t)-1)/colUnits[:t]
    plot!(_t, d_t*colUnits[:t]^2,
        label=latexstring("\\lambda = "*string(round(ρUnits_Sim[i],digits=2)))
    )
end
xlabel!(L"t")
ylabel!(L"D")
display(fig4)

## ===== Complete plot =====
# fullscale = 1.5
# fig5 = plot(fig2, fig4, fig3, fig1, layout=4, size=(fullscale*600,fullscale*400))
# figname = "particles"*string(params1["k"])*"_MultRho_complete_frictionUnits.pdf"
# SAVEFIGS && savefig(fig5, "./figures/2021/"*figname)
# display(fig5)
fig5 = plot(fig2, fig3, fig1, layout=(1,3), size=(2.5*600, 400))
figname = "growingPop_complete_frictionUnits.pdf"
SAVEFIGS && savefig(fig5, "./figures/2021/"*figname)
display(fig5)




## ===================== SIMULATION UNITS ===================
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
        ( 2000, 105, text(latexstring("\\lambda = "*string(round(ρ_Sim[1],digits=5))), 10) ),
        ( 3000, 93, text(latexstring("\\lambda = "*string(round(ρ_Sim[2],digits=5))), 10) ),
        ( 4000, 70, text(latexstring("\\lambda = "*string(round(ρ_Sim[3],digits=5))), 10) ),
        # ( 36, 25, text(latexstring("\\rho = "*string(round(ρUnits_Sim[4],digits=2))), 10) ),
    ]
)
xlabel!(L"t \quad")
ylabel!(L"\left\langle x(t)^2 \right\rangle \quad")

display(fig1b)
figname = "particles"*string(params1["k"])*"_MultRho_msd_simulationUnits.pdf"
SAVEFIGS && savefig(fig1b, "./figures/2021/"*figname)

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
diffCoeff1_t = Theorist.diffCoeffFromParticleDensity.(density1_t, thermalVals["σc"], thermalVals["E"])
diffCoeff2_t = Theorist.diffCoeffFromParticleDensity.(density2_t, thermalVals["σc"], thermalVals["E"])
diffCoeff3_t = Theorist.diffCoeffFromParticleDensity.(density3_t, thermalVals["σc"], thermalVals["E"])
diffCoeff4_t = Theorist.diffCoeffFromParticleDensity.(density4_t, thermalVals["σc"], thermalVals["E"])
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
SAVEFIGS &&  savefig(fig5b, "./figures/2021/"*figname)
display(fig5b)

