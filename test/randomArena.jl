
include("../src/bmparticles.jl")
include("../src/bmtheory.jl")

module Test

using Revise
using LaTeXStrings
using ..BParts
using ..Theorist
using Plots
theme(:juno)
gr()

function fillDensity(arenaParams)
    xperiod = abs(arenaParams["bounds"][1][2] - arenaParams["bounds"][1][1])
    yperiod = abs(arenaParams["bounds"][2][2] - arenaParams["bounds"][2][1])
    n0 = arenaParams["n0"]
    r = arenaParams["radius"]
    V = xperiod*yperiod
    fillDens = n0*π*r^2 / V
    return fillDens
end

arenaParams =
    Dict(
        "n0"=>1000,
        "evolveTime"=> 100,
        "bounds"=>((0.,10.),(0.,10.)),
        "radius"=>0.08,
        "speed"=>0.02
    )
println("surface covered: ", fillDensity(arenaParams))

# @time myArena = BParts.buildRandArena(arenaParams["bounds"], arenaParams["n0"], arenaParams["radius"], arenaParams["speed"];
#     fixSpeed=true, overlapScans=40, attempts=10, verbose=true)

@time arena, posTime_t_dim_id, velTime_t_dim_id, cells_T_ID, times_t =
    BParts.randArenaEvolve(arenaParams["n0"], arenaParams["evolveTime"], 0.1, arenaParams)


# arena1, posTime_t_dim_id, velTime_t_dim_id, cells_T_ID =
#     BParts.randArenaEvolve(nCells, evolveTime, 0.1, arenaParams, plotting=false)


# speed_t_id = BParts.speedCalc(velTime_t_dim_id)
# sMean = [sum(speed_id_t[:, t])/size(speed_id_t, 1) for t in 1:size(speed_id_t, 2)]

# rDist = BParts.rayleighDistFit(velTime_t_dim_id)

# h = histogram(vec(speed_t_id[100:500, :]), bins=range(0, 0.25, length=50),
#     normalize=true, ylims=(0,14), xlabel="speed", ylabel="distribution",
#     label="simulation", dpi=100)
# plot!(range(0, 0.25, length=100), pdf.(rDist, range(0, 0.25, length=100)), label="fit")
# plot!(range(0, 0.25, length=100), pdf.(Rayleigh(sqrt(eKinAv)), range(0, 0.25, length=100)),
#     line=(:dash), label="theory")
# display(h)
# # savefig(h, "speedDist.png")
# # savefig(h, "speedDist.svg")
# mfpExp = BParts.meanFreePath(velTime_t_dim_id, 1.)
# mfpTheory = Theorist.meanFreePath(nCells/100., 2*2*0.08)
# # println(mfpExp)
# # println(mfpTheory)
# friction = Theorist.friction(mfpExp, eKinAv)
# friction2 = Theorist.friction(mfpTheory, eKinAv)
#
#
# corrTime = 200
# velCorr_t = BParts.velocityAutocorrelation(velTime_t_dim_id[200:(200+corrTime),:,:])
# velCorrTheory_t = [Theorist.velocityAutoCorrelation(t, eKinAv, friction) for t in 1:corrTime]
# velCorrTheory2_t = [Theorist.velocityAutoCorrelation(t, eKinAv, friction2) for t in 1:corrTime]
# # h2 = plot(velCorr_t, label="simulation", xlabel=L"\tau", ylabel=L"\left< v(t)v(t+\tau) \right>", ylims=(-0.0005, 0.005), dpi=100)
# h2 = plot(velCorr_t, label="simulation", xlabel="\tau", ylabel="E[v(t)v(t+\tau)]", ylims=(-0.0005, 0.005), dpi=100)
# plot!(velCorrTheory_t, label="theory (mfp experiment)")
# plot!(velCorrTheory2_t, label="theory (mfp theory)", line=(:dash))
# # savefig(h2, "velautocor250.svg")
# # savefig(h2, "velautocor250.png")
# display(h2)
#
#
# msd_t = BParts.meanSquaredDisplacement(posTime_t_dim_id[1:200,:,:], [10., 10.])
# h3 = plot(msd_t)
# xlabel!("time")
# ylabel!("msd")
# display(h3)
#
#
# h4 = scatter(posTime_t_dim_id[end,1,:], posTime_t_dim_id[end,2,:], xlims = (0,10), ylims = (0,10), legend=false)
# display(h4)


end
