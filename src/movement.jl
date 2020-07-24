### ======== Cell movement functions ========

function moveCell!(cell::Cell)
    movePoint!(cell.pos, cell.vel)
end

function moveCell!(cell::Cell, time::Float64)
    movePoint!(cell.pos, cell.vel*time)
end

function moveCell!(cell::Cell, bounds::Bounds)
    movePoint!(cell.pos, cell.vel, bounds)
end

function moveCell!(cell::Cell, bounds::Bounds, time::Float64)
    movePoint!(cell.pos, cell.vel*time, bounds)
end

function moveCell!(cell::Cell, bounds::Bounds, r::MVector{2,Float64})
    movePoint!(cell.pos, r, bounds)
end

function moveCell!(cell::Cell, r::MVector{2,Float64})
    movePoint!(cell.pos, r)
end

function movePoint!(p::MVector{2,Float64}, d::MVector{2,Float64})
    p[1] += d[1]
    p[2] += d[2]
    nothing
end

function movePoint!(p::MVector{2,Float64}, d::MVector{2,Float64}, bounds::Bounds)
    movePoint!(p, d)
    if !inBounds(p, bounds)
        toBoundsPeriodic!(p, bounds)
    end
    nothing
end

# function movePointAlt!(p::Vector{Float64}, v::Vector{Float64}, bounds::Bounds)
#     if v[1] > 0
#         p[1] = bounds.x[1] + (p[1]-bounds.x[1]+v[1])%xDist(bounds)
#     else
#         p[1] = bounds.x[2] - (-p[1]+bounds.x[2]-v[1])%xDist(bounds)
#     end
#     if v[2] > 0
#         p[2] = bounds.y[1] + (p[2]-bounds.y[1]+v[2])%yDist(bounds)
#     else
#         p[2] = bounds.y[2] - (-p[2]+bounds.y[2]-v[2])%yDist(bounds)
#     end
#     return p
# end
