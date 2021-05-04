using JLD2, Glob
using LsqFit
using LaTeXStrings
using Plots
pyplot()

include("../../src/bmtheory.jl")
using .Theorist


## ========== Load data ==========

fileNames = glob("simResultFixedDensity_*_msd.jld2", "./data/FixedDensity")

##
ρ_Sim = Float64[]
_t = Int[]
msdLan_t_Sim = Vector{Float64}[]
msdPar_t_Sim = Vector{Float64}[]
msdTheory_t_Sim = Vector{Float64}[]
params_Sim = Dict[]
for fname in fileNames
    @load fname arenaParams growthParams _t msdLan_t msdPar_t
    push!(msdLan_t_Sim, msdLan_t)
    push!(msdPar_t_Sim, msdPar_t)
    push!(params_Sim, arenaParams)
    global _t = _t

    r = arenaParams["radius"]
    N = arenaParams["n0"]
    V = arenaParams["volume"]
    ρ = N * π*r^2 / V
    push!(ρ_Sim, ρ)

    thermalVals = Theorist.thermalValues(arenaParams)
    msdTheory_t = Theorist.msdTheory.(_t, arenaParams["n0"]/arenaParams["volume"], thermalVals["σc"], thermalVals["E"])
    push!(msdTheory_t_Sim, msdTheory_t)
end


## ========== Plot figure ==========

fig1 = plot(legend=:topleft)

for (i,simId) in enumerate([1,2,3,4,5,8])
    plot!(_t, msdPar_t_Sim[simId],
        color=i,
        linealpha=0.5,
        linestyle=:solid,
        label="ρ="*string(round(ρ_Sim[simId],digits=3))*" -- Particle sim"
    )
    # plot!(_t, msdLan_t_Sim[i],
    plot!(_t, msdTheory_t_Sim[simId],
        color=i,
        linestyle=:dash,
        label="ρ="*string(round(ρ_Sim[simId],digits=3))*" -- Langevin"
    )
end
xlabel!(L"t")
ylabel!(L"\left\langle x(t)^2 \right\rangle")

display(fig1)
figname = "msdFixedDensity.pdf"
savefig(fig1, "./figures/2021/"*figname)

##

DPar_Sim = Float64[]
DTheory_Sim = Float64[]
for simId in 1:length(fileNames)
    thermVals = Theorist.thermalValues(params_Sim[simId])
    E = thermVals["E"]
    @. msdModel(t, p) = 4E/p[1] * ( t - (1-exp(-p[1]*t))/p[1] )

    γ0 = [thermVals["γ"]]
    fit = curve_fit(msdModel, _t, msdPar_t_Sim[simId], γ0)
    DPar = fit.param[1] * E
    println(DPar)
    n = params_Sim[simId]["n0"]/params_Sim[simId]["volume"]
    σ = thermVals["σc"]
    DTheory = Theorist.diffCoeffFromParticleDensity(n, σ, E)
    push!(DPar_Sim, DPar)
    push!(DTheory_Sim, DTheory)
end

##
fig2 = plot(legend=:topleft)
plot!(ρ_Sim, DPar_Sim, 
    markershape=:auto,
    label="particle simulation"
)
plot!(ρ_Sim, DTheory_Sim, 
    markershape=:auto,
    label="Langevin theory"
)
xlabel!("particle surface density")
ylabel!("D")
display(fig2)
figname = "intensityDensityDependence.pdf"
savefig(fig2, "./figures/2021/"*figname)


