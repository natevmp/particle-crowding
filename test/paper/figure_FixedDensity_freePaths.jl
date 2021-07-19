include("../../src/bmtheory.jl")
include("../../src/bmparticles.jl")
using .Theorist, .BParts
using Distributions
using JLD2, FileIO
using LaTeXStrings

## ===== Load data =====

files_f = [
    "./data/FixedDensity/21-05-10/simResultFixedDensity_0.01.jld2",
    "./data/FixedDensity/21-05-10/simResultFixedDensity_0.05.jld2",
    "./data/FixedDensity/21-05-10/simResultFixedDensity_0.1.jld2",
    "./data/FixedDensity/21-05-10/simResultFixedDensity_0.2.jld2",
    "./data/FixedDensity/21-05-10/simResultFixedDensity_0.3.jld2",
    "./data/FixedDensity/21-05-10/simResultFixedDensity_0.5.jld2",
]

##
volumeDens(n, r, v) = n*π*r^2/v
ρ_sim = Float64[]
n_sim = Float64[]
thermVals_sim = Dict[]
arenaParams_sim = Dict[]
pos_sim_T_dim_id = Array{Float64, 3}[]
vel_sim_T_dim_id = Array{Float64, 3}[]
times_sim_T = Vector{Float64}[]
units_sim = Dict[]

# @load "./data/FixedDensity/21-05-10/simResultFixedDensity_0.5.jld2"


for (i,file) in enumerate(files_f)
    arenaParams, growthParams, pos_t_dim_id, vel_t_dim_id, times_t = load(file, "arenaParams", "growthParams", "pos_t_dim_id", "vel_t_dim_id", "times_t")

    thermVals = Theorist.thermalValues(arenaParams)
    colUnits = Dict(
        :t => 1 / thermVals["γ"],
        :x => arenaParams["speed"]/thermVals["γ"],
    )

    nDens = arenaParams["n0"]/arenaParams["volume"]

    push!(arenaParams_sim, arenaParams)
    push!( ρ_sim, volumeDens(arenaParams["n0"], arenaParams["radius"], arenaParams["volume"]) )
    push!(n_sim, nDens)
    push!(thermVals_sim, thermVals)
    push!(pos_sim_T_dim_id, pos_t_dim_id)
    push!(vel_sim_T_dim_id, vel_t_dim_id)
    push!(units_sim, colUnits)
    push!(times_sim_T, times_t)
end


## ===== mean free paths =====


mfpPart_sim = Float64[]
@time for sim in 1:length(ρ_sim)
    fpPart_ = BParts.ensembleFreePaths(vel_sim_T_dim_id[sim], times_sim_T[sim])
    mfp = mean(fpPart_)
    push!(mfpPart_sim, mfp)
end

mfpTheory_sim = (n -> 1/(sqrt(2)*thermVals_sim[1]["σc"]*n) ).(n_sim)




## ========== Plotting ==========
using Plots
pyplot()
theme(:bright, size=(0.8*600,0.8*400), minorgrid=false, gridstyle=:dash, fontfamily="DejaVu Sans")
palette = cgrad(:tableau_red_green_gold, 6, categorical=true, rev=true)
# theme(:default)

## ===== MFP theory vs particle sim =====
figMFP = plot(
    ρ_sim, mfpPart_sim,
    markershape=:auto,
    label="particle simulation",
)
plot!(
    ρ_sim, mfpTheory_sim,
    markershape=:auto,
    label="model"
)


## ===== Error on MFP =====

mfpRel_sim = mfpTheory_sim ./ [thermVals_sim[sim]["σc"] for sim in 1:length(arenaParams_sim)]
mfpErrorRel_sim = abs.(mfpTheory_sim .- mfpPart_sim) ./ mfpPart_sim

figError = plot(
    mfpRel_sim, mfpErrorRel_sim,
    # series_annotations="ρ=".*string.(round.(ρ_sim, digits=2)),
    # legend=:topleft,
    color=:grey45,
    markershape=:auto,
    label="",
    xscale=:log10,
)
annotate!(
    [ 
        ( mfpRel_sim[i], mfpErrorRel_sim[i]+0.1, Plots.text(L"c="*string(round(ρ_sim[i], digits=2)), 9) ) for i in 1:length(ρ_sim)
    ]
)
ylims!(0,2.4)
xlims!(10^(-0.6), 10^(1.3))
xlabel!(L"l/\sigma")
ylabel!(L"\epsilon_l")
display(figError)

figname = "figures/2021/fixedDensity_mpfError.pdf"
savefig(figError, figname)


## ===== Velocity autocorrelation =====
palette = cgrad(:tableau_red_green_gold, 6, categorical=true, rev=true)

_tCorr = range(0,5, length=100)
vCorrPred_t = map( t->exp(-t), _tCorr )
figVelCor = plot(
    _tCorr, vCorrPred_t,
    label=L"e^{-\tau}",
    color=:black,
    linestyle=:dash,
    linewidth=1.3,
    # yscale=:log10
)
for simId in 1:length(ρ_sim)

    vel_t_dim_id = vel_sim_T_dim_id[simId]
    thermVals = thermVals_sim[simId]
    units = units_sim[simId]

    tMin = 50/thermVals["γ"]
    tMax = tMin + 5/thermVals["γ"]

    tIndMin = findall(t->t>tMin, times_sim_T[simId])[1]
    tIndMax = findall(t->t>tMax, times_sim_T[simId])[1]

    _tCorr = times_sim_T[simId][tIndMin:tIndMax] .- times_sim_T[simId][tIndMin]

    # tParIn = growthParams["waitTime"]
    # tParIn = 800
    # _tCorr = growthParams["waitTime"]:growthParams["waitTime"]+corrTime
    # _tCorrInd = range(tParIn, Int(floor(tParIn+corrTime-1)), step=1)
    # _tCorr = 0:corrTime-1
    # timesCorr_t, vCorrLan_t = Theorist.velCorrelation(langevinEnsemble, (_tCorr[1],_tCorr[end]), steps=corrTime)

    vCorrPar_t = BParts.velocityAutocorrelation(vel_t_dim_id[tIndMin:tIndMax, :, :])
    vCorrPred_t = map( t->2*thermVals["E"]*exp(-thermVals["γ"]*t), _tCorr )

    plot!(
        _tCorr / units[:t], vCorrPar_t / (units[:x]/units[:t])^2, 
        label=L"c="*string(round(ρ_sim[simId], digits=2)),
        color=palette[simId],
        linewidth=0.9,
    )
    # plot!(
    #     _tCorr / units[:t], vCorrPred_t / (units[:x]/units[:t])^2,
    #     label="analytical prediction"
    # )
end

xlabel!(L"\tau \quad \left[ 1/\gamma \right]")
ylabel!(L"\left\langle v(t)v(t+\tau) \right\rangle \quad \left[ s^2 \right]")
xlims!(0, 3)
ylims!(0, 1)
display(figVelCor)

figname = "figures/2021/fixedDensity_velAutoCor.pdf"
savefig(figVelCor, figname)