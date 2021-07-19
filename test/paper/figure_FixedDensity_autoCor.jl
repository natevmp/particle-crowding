include("../../src/bmtheory.jl")
using .Theorist
using Distributions
using JLD2, FileIO
using LaTeXStrings


## ===== Load data =====

files_f = glob("simResultFixedDensity_*.jld2", "./data/FixedDensity/21-05-10")

##
volumeDens(n, r, v) = n*π*r^2/v
ρ_sim = Float64[]
thermVals_sim = Dict[]
pos_sim_T_dim_id = Array{Float64, 3}[]
vel_sim_T_dim_id = Array{Float64, 3}[]
units_sim = Dict[]

for (i,file) in enumerate(files_f)
    arenaParams, growthParams, pos_t_dim_id, vel_t_dim_id = load(file, "arenaParams", "growthParams", "pos_t_dim_id", "vel_t_dim_id")

    thermVals = Theorist.thermalValues(arenaParams)
    colUnits = Dict(
        :t => 1 / thermVals["γ"],
        :x => arenaParams["speed"]/thermVals["γ"],
    )
    push!( ρ_sim, volumeDens(arenaParams["n0"], arenaParams["radius"], arenaParams["volume"]) )
    push!(thermVals_sim, thermVals)
    push!(pos_sim_T_dim_id, pos_t_dim_id)
    push!(vel_sim_T_dim_id, vel_t_dim_id)
    push!(units_sim, colUnits)
end


## ========== Plotting ==========
using Plots
pyplot()
theme(:bright, size=(0.8*600,0.8*400), minorgrid=false, gridstyle=:dash, fontfamily="DejVu Sans")


## ===== Velocity autocorrelation =====
palette = cgrad(:tableau_red_blue, 6, categorical=true, rev=true)

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

    corrTime = 5 * 1 / thermVals["γ"]
    # tParIn = growthParams["waitTime"]
    tParIn = 800
    # _tCorr = growthParams["waitTime"]:growthParams["waitTime"]+corrTime
    _tCorrInd = range(tParIn, Int(floor(tParIn+corrTime-1)), step=1)
    _tCorr = 0:corrTime-1
    # timesCorr_t, vCorrLan_t = Theorist.velCorrelation(langevinEnsemble, (_tCorr[1],_tCorr[end]), steps=corrTime)
    vCorrPar_t = BParts.velocityAutocorrelation(vel_t_dim_id[_tCorrInd, :, :])
    vCorrPred_t = map( t->2*thermVals["E"]*exp(-thermVals["γ"]*t), _tCorr )

    plot!(
        _tCorr / units[:t], vCorrPar_t / (units[:x]/units[:t])^2, 
        label="ρ="*string(round(ρ_sim[simId], digits=2)),
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
# savefig(figVelCor, figname)


