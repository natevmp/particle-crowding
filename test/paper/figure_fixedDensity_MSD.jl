using JLD2, Glob
using LsqFit
include("../../src/bmtheory.jl")
using .Theorist

## ==================== Load data ====================

# fileNames = glob("simResultFixedDensity_*_msd.jld2", "./data/FixedDensity")
fileNames = [
    "./data/FixedDensity/21-05-10/simResultFixedDensity_0.01_msd.jld2",
    "./data/FixedDensity/21-05-10/simResultFixedDensity_0.05_msd.jld2",
    "./data/FixedDensity/21-05-10/simResultFixedDensity_0.1_msd.jld2",
    "./data/FixedDensity/21-05-10/simResultFixedDensity_0.2_msd.jld2",
    "./data/FixedDensity/21-05-10/simResultFixedDensity_0.3_msd.jld2",
    "./data/FixedDensity/21-05-10/simResultFixedDensity_0.5_msd.jld2",
]
##
ρ_sim = Float64[]
_t_Sim = Vector{Float64}[]
# msdLan_t_Sim = Vector{Float64}[]
msdPar_sim_T = Vector{Float64}[]
msdTheory_sim_T = Vector{Float64}[]
params_sim = Dict[]
units_sim = Dict[]
for fname in fileNames
    @load fname arenaParams growthParams _t msdTheory_t msdPar_t
    push!(msdPar_sim_T, msdPar_t)
    push!(params_sim, arenaParams)
    push!(_t_Sim, _t)
    push!(msdTheory_sim_T, msdTheory_t)

    r = arenaParams["radius"]
    N = arenaParams["n0"]
    V = arenaParams["volume"]
    ρ = N * π*r^2 / V
    push!(ρ_sim, ρ)
    
    thermalVals = Theorist.thermalValues(arenaParams)
    colUnits = Dict(
        :t => 1 / thermalVals["γ"],
        :x => arenaParams["speed"]/thermalVals["γ"],
    )
    push!(units_sim, colUnits)
end


## ==================== Plot figures ====================
using Plots, LaTeXStrings
pyplot()
theme(:bright, size=(0.8*600,0.8*400), minorgrid=false, gridstyle=:dash, fontfamily="DejaVu Sans")
# theme(:default)


## ---------- MSD simulation units ----------

palette = cgrad(:tableau_red_green_gold, 4, categorical=true, rev=true)

fig1 = plot(legend=:topleft)
for (i,simId) in enumerate([1,2,3,6])
    plot!(_t_Sim[simId], msdPar_sim_T[simId],
        # color=i,
        color=palette[i],
        linealpha=0.5,
        linestyle=:solid,
        label="ρ="*string(round(ρ_sim[simId],digits=3))*"\nparticle sim"
    )
    # plot!(_t, msdLan_t_Sim[i],
    plot!(_t_Sim[simId], msdTheory_sim_T[simId],
        color=palette[i],
        linestyle=:dash,
        label="Langevin model"
    )
end
xlabel!(L"t")
ylabel!(L"\left\langle x(t)^2 \right\rangle")
xlims!(0,500)
ylims!(0,60)

display(fig1)
figname = "fixedDensity_msd.pdf"
savefig(fig1, "./figures/2021/"*figname)



## ---------- Friction density dependence ----------
# ----- Calculations -----
DPar_Sim = Float64[]
γPar_Sim = Float64[]
DTheory_Sim = Float64[]
γTheory_Sim = Float64[]
for simId in 1:length(fileNames)
    thermVals = Theorist.thermalValues(params_sim[simId])
    E = thermVals["E"]
    @. msdModel(t, p) = 4E/p[1] * ( t - (1-exp(-p[1]*t))/p[1] )

    γ0 = [thermVals["γ"]]
    fit = curve_fit(msdModel, _t_Sim[simId], msdPar_sim_T[simId], γ0)
    γPar = fit.param[1]
    DPar = fit.param[1] * E
    n = params_sim[simId]["n0"]/params_sim[simId]["volume"]
    σ = thermVals["σc"]
    DTheory = Theorist.diffCoeffFromParticleDensity(n, σ, E)
    γTheory = Theorist.frictionFromParticleDensity(n, σ, E)

    push!(DPar_Sim, DPar)
    push!(γPar_Sim, γPar)
    push!(DTheory_Sim, DTheory)
    push!(γTheory_Sim, γTheory)
end
## ----- Plot figure -----
fig2 = plot(legend=:topleft)
plot!(ρ_sim, γPar_Sim, 
    markershape=:auto,
    label="particle simulation"
)
plot!(ρ_sim, γTheory_Sim, 
    markershape=:auto,
    label="Langevin theory"
)
xlabel!(L"\rho")
ylabel!(L"\gamma")
display(fig2)
figname = "fixedDensity_frictionDensityDependence.pdf"
savefig(fig2, "./figures/2021/"*figname)





## ---------- MSD natural units ----------

palette = cgrad(:tableau_red_green_gold, 6, categorical=true, rev=true)

figMSDNat = plot(legend=:topleft)

for simId in 1:6

    _tCu = _t_Sim[simId] ./ units_sim[simId][:t]
    msd_tCu = msdPar_sim_T[simId] ./ (units_sim[simId][:x])^2

    plot!(
        _tCu, msd_tCu,
        label="ρ="*string(round(ρ_sim[simId], digits=2)),
        color=palette[simId],
    )

end
_tCu = 0:0.1:400
msdColUnitsTheory_t = ( t -> 2*( t - (1-exp(-t)) ) ).(_tCu)
plot!(
    _tCu, msdColUnitsTheory_t, 
    color=:black,
    linealpha=0.7,
    linestyle=:dash,
    linewidth=1.5,
    label="Langevin\nprediction"
)

xlims!(0,100)
ylims!(0,220)
xlabel!(L"t \quad \left[ 1/\gamma \right]")
ylabel!(L"\left\langle \left( x(t) \right)^2 \right\rangle \quad \left[ (s/\gamma)^2 \right]")

display(figMSDNat)

figname = "fixedDensity_msdNaturalUnits.pdf"
savefig(figMSDNat, "./figures/2021/"*figname)
