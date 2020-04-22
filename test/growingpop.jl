
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
gr()

function randArenaEvolve(nCells::Int, steps::Int, growthParams::Union{Tuple, Nothing}=nothing; plotting=false, animating=false)
    arena = buildRandArena(Bounds((0.,10.), (0.,10.)), nCells, 0.08, 0.05, fixSpeed=false)

    arenaCellPositions_dim_id = BParts.cellPositions_DIM_ID(arena)

    eKin = BParts.kineticEnergy(arena)
    println("::::: Initial total kinetic energy: ", eKin)

    println("\n")
    println("=== Evolving arena ===")
    println("\n")

    if animating
        anim = Animation()
    else anim = nothing
    end
    posTime_t_dim_id, velTime_t_dim_id, cells_T_ID = 
        evolveArena!(arena, steps, growthParams, plotsteps=plotting, animator=anim)
    if animating
        gif(anim, "figures/animation.gif", fps=10)
    end

    eKin = BParts.kineticEnergy(arena)
    println("::::: Final total kinetic energy: ", eKin)

    return posTime_t_dim_id, velTime_t_dim_id, cells_T_ID, eKin
end






nCells = 50
evolveTime = 100
logisticRate(n::Real, ρ::Real, k::Real) = n*ρ*(1-n/k)
exponentialRate(n::Real, ρ::Real) = n*ρ
growthParams = (n->logisticRate(n, 0.06, 500), 0.08, 0.05)
posTime_t_dim_id, velTime_t_dim_id, cells_T_ID, eKin = 
    randArenaEvolve(nCells, evolveTime, growthParams, animating=false)
# eKinAv = eKin / nCells

speed_t_id = BParts.speedCalc(velTime_t_dim_id)
sMean_t = [sum(speed_t_id[t, :])/size(speed_t_id, 2) for t in 1:size(speed_t_id, 1)]


rDist = BParts.rayleighDistCompare(velTime_t_dim_id)

p2 = plot(map(length, cells_T_ID))
xlabel!("time")
ylabel!("number of cells")
display(p2)


# h = histogram(vec(speed_t_id[:,:]), bins=range(0, 0.25, length=50), 
#     normalize=true, ylims=(0,14), xlabel="speed", 
#     ylabel="distribution", label="simulation")
# plot!(range(0, 0.25, length=100), pdf.(rDist, range(0, 0.25, length=100)), label="fit")
# # plot!(range(0, 0.25, length=100), pdf.(Rayleigh(sqrt(eKinAv)), range(0, 0.25, length=100)), line=(:dash), label="theory")
# display(h)
# savefig(h, "speedDist.png")
# savefig(h, "speedDist.svg")
# mfpExp = BParts.meanFreePath(velTime_id_t_dim, 1.)
# mfpTheory = Theorist.meanFreePath(nCells/100., 2*2*0.08)
# println(mfpExp)
# println(mfpTheory)
# friction = Theorist.friction(mfpExp, eKinAv)
# friction2 = Theorist.friction(mfpTheory, eKinAv)


end