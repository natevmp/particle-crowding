
include("../src/bmparticles.jl")
include("../src/bmtheory.jl")

module Test

using ..BParts
using ..Theorist
using StaticArrays
using LinearAlgebra
using Plots
using Distributions
using LaTeXStrings
using Debugger
# pyplot()

### ======== Random arena animation test ========

function randArenaEvolve(nCells::Int, steps::Int; plotting=true, animating=false)
    arena = buildRandArena(Bounds((0.,10.), (0.,10.)), nCells, 0.08, 0.05, fixSpeed=false)
    arenaCellPositions_dim_id = BParts.cellPositions_DIM_ID(arena)
    # s = scatter(arenaCellPositions_dim_id[1,:], arenaCellPositions_dim_id[2,:],
    #             xlims = (0,10), ylims = (0,10), legend=false)
    # display(s)

    BParts.fixArenaOverlaps!(arena)

    arenaCellPositions_dim_id = BParts.cellPositions_DIM_ID(arena)
    # s = scatter(arenaCellPositions_dim_id[1,:], arenaCellPositions_dim_id[2,:],
    #             xlims = (0,10), ylims = (0,10), legend=false)
    # display(s)

    eKin = BParts.kineticEnergy(arena)
    println("::::: Initial total kinetic energy: ", eKin)

    println("\n")
    println("=== Evolving arena ===")
    println("\n")

    if animating
        anim = Animation()
    else anim = nothing
    end
    posTime_id_t_dim, velTime_id_t_dim, cells_T_ID = BParts.evolveArena(arena, steps, plotsteps=plotting, animator=anim)
    if animating
        gif(anim, "anim_2.gif", fps=10)
    end

    eKin = BParts.kineticEnergy(arena)
    println("::::: Final total kinetic energy: ", eKin)

    return posTime_id_t_dim, velTime_id_t_dim, cells_T_ID, eKin
end

nCells = 150
evolveTime = 500
posTime_id_t_dim, velTime_id_t_dim, cells_T_ID, eKin = randArenaEvolve(nCells, evolveTime, plotting=false)

eKinAv = eKin / nCells

speed_id_t = BParts.speedCalc(velTime_id_t_dim)
# sMean = [sum(speed_id_t[:, t])/size(speed_id_t, 1) for t in 1:size(speed_id_t, 2)]

rDist = BParts.rayleighDistCompare(velTime_id_t_dim)

# h = histogram(vec(speed_id_t[:,100:2000]), bins=range(0, 0.25, length=50), normalize=true, ylims=(0,14), xlabel="speed", ylabel="distribution", label="simulation", dpi=300)
# plot!(range(0, 0.25, length=100), pdf.(rDist, range(0, 0.25, length=100)), label="fit")
# plot!(range(0, 0.25, length=100), pdf.(Rayleigh(sqrt(eKinAv)), range(0, 0.25, length=100)), line=(:dash), label="theory")
# display(h)
# savefig(h, "speedDist.png")
# savefig(h, "speedDist.svg")
mfpExp = BParts.meanFreePath(velTime_id_t_dim, 1.)
mfpTheory = Theorist.meanFreePath(nCells/100., 2*2*0.08)
# println(mfpExp)
# println(mfpTheory)
friction = Theorist.friction(mfpExp, eKinAv)
friction2 = Theorist.friction(mfpTheory, eKinAv)


corrTime = 200
velCorr_t = BParts.velocityAutocorrelation(velTime_id_t_dim[:,200:(200+corrTime),:])
velCorrTheory_t = [Theorist.velocityAutoCorrelation(t, eKinAv, friction) for t in 1:corrTime]
velCorrTheory2_t = [Theorist.velocityAutoCorrelation(t, eKinAv, friction2) for t in 1:corrTime]
h2 = plot(velCorr_t, label="simulation", xlabel=L"\tau", ylabel=L"\left< v(t)v(t+\tau) \right>", ylims=(-0.0005, 0.005), dpi=300)
plot!(velCorrTheory_t, label="theory (mfp experiment)")
plot!(velCorrTheory2_t, label="theory (mfp theory)", line=(:dash))
# savefig(h2, "velautocor250.svg")
# savefig(h2, "velautocor250.png")
display(h2)


msd_t = BParts.meanSquaredDisplacement(posTime_id_t_dim[:,1:200,:], [10., 10.])
h3 = plot(msd_t)
xlabel!("time")
ylabel!("msd")
display(h3)


end