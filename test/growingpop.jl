
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

function randArenaEvolve(nCells::Int, steps::Int, growthParams::Tuple; plotting=true, animating=false)
    arena = buildRandArena(Bounds((0.,10.), (0.,10.)), nCells, 0.08, 0.05, fixSpeed=false)
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
    posTime_t_dim_id, velTime_t_dim_id, cells_T_ID = 
        BParts.evolveArena(arena, steps, growthParams, plotsteps=plotting, animator=anim)
    if animating
        gif(anim, "anim_2.gif", fps=10)
    end

    eKin = BParts.kineticEnergy(arena)
    println("::::: Final total kinetic energy: ", eKin)

    return posTime_t_dim_id, velTime_t_dim_id, cells_T_ID, eKin
end

nCells = 100
evolveTime = 50
growthParams = (1., 0.08, 0.05)
posTime_t_dim_id, velTime_t_dim_id, cells_T_ID, eKin = 
    randArenaEvolve(nCells, evolveTime, growthParams, plotting=false)

# eKinAv = eKin / nCells

speed_t_id = BParts.speedCalc(velTime_t_dim_id)
sMean_t = [sum(speed_t_id[t, :])/size(speed_t_id, 2) 
            for t in 1:size(speed_t_id, 1)]

rDist = BParts.rayleighDistCompare(velTime_t_dim_id)

println(typeof(posTime_t_dim_id))
println(size(posTime_t_dim_id))
println(typeof(velTime_t_dim_id))
println(size(velTime_t_dim_id))
println(typeof(cells_T_ID))
println(size(cells_T_ID[end]))
println(size(speed_t_id))
println(size(sMean_t))


h = histogram(vec(speed_t_id[:,:]), bins=range(0, 0.25, length=50), 
    normalize=true, ylims=(0,14), xlabel="speed", 
    ylabel="distribution", label="simulation")
plot!(range(0, 0.25, length=100), pdf.(rDist, range(0, 0.25, length=100)), label="fit")
# plot!(range(0, 0.25, length=100), pdf.(Rayleigh(sqrt(eKinAv)), range(0, 0.25, length=100)), line=(:dash), label="theory")
display(h)
# savefig(h, "speedDist.png")
# savefig(h, "speedDist.svg")
# mfpExp = BParts.meanFreePath(velTime_id_t_dim, 1.)
# mfpTheory = Theorist.meanFreePath(nCells/100., 2*2*0.08)
# println(mfpExp)
# println(mfpTheory)
# friction = Theorist.friction(mfpExp, eKinAv)
# friction2 = Theorist.friction(mfpTheory, eKinAv)








end