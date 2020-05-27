#=
bmparticles:
- Julia version: 1.11
- Author: nathanielmonpere
- Date: 2019-08-20
=#

module BParts
export Cell, Arena, evolveArena!, cultivateArena!, Bounds, moveCell!, buildRandArena

using StaticArrays
using ElasticArrays
using LinearAlgebra
using Distances
using NearestNeighbors
using Random, Distributions, PoissonRandom
using ProgressMeter
using Plots


### ======== Structures ========
struct Cell{T<:MVector{2,Float64}, R<:AbstractFloat}
    pos::T
    vel::T
    radius::R
end

struct Bounds
    x::Tuple{Real, Real}
    y::Tuple{Real, Real}
    xLen::Real
    yLen::Real

    function Bounds(x, y)
        xLen = abs(x[2]-x[1])
        yLen = abs(y[2]-y[1])
        new(x, y, xLen, yLen)
    end
end

struct Arena{T<:AbstractArray{C,1} where C<:Cell, B<:Bounds}
    # Array containing all cells in the population. A cell's index in the array is its ID
    cellsList::T
    bounds::B
end

include("arenabuilder.jl")
include("dynamics.jl")
include("positions.jl")
include("movement.jl")
include("collisions.jl")
include("growth.jl")
include("scientist.jl")
include("util.jl")

end
