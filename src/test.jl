
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





# function timetest(arena::Arena, n::Integer)
#     for i in 1:n
#         BParts.cultivateArena!(arena, 1., 1., 0.08, 0.05)
#     end
# end
# function test(n)
#     arena = buildRandArena(Bounds((0.,10.), (0.,10.)), 10, 0.08, 0.05, fixSpeed=false)
#     println(length(arena.cellsList))
    
#     @time timetest(arena, n)

#     println(length(arena.cellsList))
# end
# test(100)


p = plot([1,2,3,4,NaN,8,9])
display(p)



end