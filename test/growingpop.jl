
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



nCells = 50
evolveTime = 100

arenaParams = 
    Dict(
        "n0"=>nCells,
        "evolveTime"=>evolveTime,
        "bounds"=>((0.,20.),(0.,20.)), 
        "radius"=>0.08, 
        "speed"=>0.03
    )
growthParams = 
    Dict(
        "ρ"=> 0.06,
        "k"=> 2000
    )

arena, pos_t_dim_id, vel_t_dim_id, cells_T_ID = 
    BParts.randArenaEvolve(nCells, evolveTime, arenaParams, growthParams)


p2 = plot(length.(cells_T_ID))
xlabel!("time")
ylabel!("number of cells")
display(p2)


end