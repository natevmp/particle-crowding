
include("../src/bmparticles.jl")
include("../src/bmtheory.jl")

module Test

using ..BParts
using ..Theorist
using StaticArrays
using LinearAlgebra
using Plots
gr()
using Distributions
using LaTeXStrings

### ======== Two cell collision test ========

function twoCellColl(bounds::Bounds, cells::AbstractArray{C,1} where C<:Cell)
    arena = Arena(cells, bounds)
    println("vel A: ", norm(arena.cellsList[1].vel))
    println("vel B: ", norm(arena.cellsList[2].vel))
    println("total kinetic energy: ",
            (norm(arena.cellsList[1].vel)^2+norm(arena.cellsList[2].vel)^2)/2)
    for i in 1:15
        println("loop: ", i)
        println("distance between cells: ",
                BParts.cellDistance(arena.cellsList[1], arena.cellsList[2], bounds))

        # fix collisions
        collidedCellsList_id = BParts.collider!(arena)
        if length(collidedCellsList_id)>0
               println(">>> Collision detected! <<<")
        end
        #get cell locations
        arenaCellPositions_dim_id = BParts.cellPositions_DIM_ID(arena)
        #get collision locans
        collisions_dim_id = arenaCellPositions_dim_id[:, collidedCellsList_id]
        #plot
        s = scatter((arenaCellPositions_dim_id[1,1], arenaCellPositions_dim_id[2,1]),
                    xlims = (0,10), ylims = (0,10), legend=false, ms=8)
        scatter!((arenaCellPositions_dim_id[1,2], arenaCellPositions_dim_id[2,2]), ms=8)
        scatter!(collisions_dim_id[1,:], collisions_dim_id[2,:], ms=8)

        # evolve normally
        for cell in arena.cellsList
            moveCell!(cell, arena.bounds)
        end
        display(s)
    end
    println("vel A: ", norm(arena.cellsList[1].vel))
    println("vel B: ", norm(arena.cellsList[2].vel))
    println("total kinetic energy: ",
            (norm(arena.cellsList[1].vel)^2+norm(arena.cellsList[2].vel)^2)/2)
end

bounds = Bounds((0.,10.),(0.,10.))
cells_C = [Cell(@MVector([5.,4.8]), @MVector([0.05,0.]), 0.2),
            Cell(@MVector([6.,5.15]), @MVector([-0.09,0.]), 0.2)]
twoCellColl(bounds, cells_C)

end
