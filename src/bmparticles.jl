#=
bmparticles:
- Julia version: 1.11
- Author: nathanielmonpere
- Date: 2019-08-20
=#

module BParts
export Cell, Arena, evolveArena, Bounds, moveCell!, buildArena, buildRandArena, nestfor

using StaticArrays
using LinearAlgebra
using Distances
using NearestNeighbors
using Random
using Plots


### ======== Structures ========
struct Cell{T<:MVector{2,Float64}, R<:AbstractFloat}
    pos::T
    vel::T
    radius::R
end

# struct Bounds{T<:Real}
#     x::Tuple{T, T}
#     y::Tuple{T, T}
#     xLen::T
#     yLen::T
#     function Bounds{T}(x::Tuple{T, T}, y::Tuple{T, T}) where T<:Real
#         xLen = abs(x[2]-x[1])
#         yLen = abs(y[2]-y[1])
#         new(x, y, xLen, yLen)
#     end
# end
struct Bounds
    x::Tuple{Real,Real}
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
    cellsList::T
    bounds::B
end

include("arenabuilder.jl")
include("positions.jl")
include("movement.jl")
include("collisions.jl")
include("scientist.jl")
include("util.jl")

end
