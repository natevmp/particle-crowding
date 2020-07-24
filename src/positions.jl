
### ======== Positional Functions ========

function cellPositions_ID_DIM(arena::Arena)
    positions_id_dim = Array{Float64,2}(undef, length(arena.cellsList), 2)
    for (i, cell) in enumerate(arena.cellsList)
        positions_id_dim[i,:] = cell.pos
    end
    return positions_id_dim
end

function cellPositions_ID_DIM!(arena::Arena, posList_id_dim::Array{Float64, 2})
    for (i, cell) in enumerate(arena.cellsList)
        posList_id_dim[i,:] = cell.pos
    end
    return nothing
end

function cellPositions_ID_DIM(arena::Arena, cellsList_id::Vector{Int})
    positions_id_dim .= cellPositions_ID_DIM(arena)[cellsList_id, :]
    return positions_id_dim
end

"""
Get _dim_id array of cell positions.
"""
function cellPositions_DIM_ID(arena::Arena)
    positions_dim_id = Array{Float64,2}(undef, 2, length(arena.cellsList))
    for (i, cell) in enumerate(arena.cellsList)
        positions_dim_id[:,i] = cell.pos
    end
    return positions_dim_id
end

"""
Get _dim_id array of cell positions in primary box.
"""
function cellPositionsPeriodic_DIM_ID(arena::Arena)
    positionsP_dim_id = Array{Float64,2}(undef, 2, length(arena.cellsList))
    for (i, cell) in enumerate(arena.cellsList)
        # positions_dim_id[:,i] = cell.pos
        positionsP_dim_id[1,i] = toBoundsPeriodic(cell.pos[1], arena.bounds.x)
        positionsP_dim_id[2,i] = toBoundsPeriodic(cell.pos[2], arena.bounds.y)
    end
    return positionsP_dim_id
end

function cellPositions_DIM_ID(cellsList_c::Vector{C} where C<:Cell)
    positions_dim_id = Array{Float64,2}(undef, 2, length(cellsList_c))
    for (i, cell) in enumerate(cellsList_c)
        positions_dim_id[:,i] = cell.pos
    end
    return positions_dim_id
end

function cellPositions_DIM_ID(arena::Arena, cellsList_id::Vector{Int})
    positions_dim_id = cellPositions_DIM_ID(arena)[:, cellsList_id]
    return positions_dim_id
end

function cellDistance(cellA::Cell, cellB::Cell, bounds::Bounds)
    # norm(cellA.pos - cellB.pos)
    # dx = abs(cellA.pos[1]-cellB.pos[1])
    # dy = abs(cellA.pos[2]-cellB.pos[2])
    #
    # if (dx > bounds.xLen/2)
    #     dx = bounds.xLen - dx
    # end
    # if (dy > bounds.yLen/2)
    #     dy = bounds.yLen - dy
    # end
    #
    # return sqrt(dx^2 + dy^2)
    # evaluate(PeriodicEuclidean([bounds.xLen, bounds.yLen]), cellA.pos, cellB.pos)
    peuclidean(cellA.pos, cellB.pos, [bounds.xLen, bounds.yLen])
end

function cellVelocities_ID_DIM(arena::Arena)
    velocities_id_dim = Array{Float64,2}(undef, length(arena.cellsList), 2)
    for (i, cell) in enumerate(arena.cellsList)
        velocities_id_dim[i,:] = cell.pos
    end
    return velocities_id_dim
end


### ======== Bounds related functions ========
function inBounds(cell::Cell, bounds::Bounds)
    inX = bounds.x[1] <= cell.pos[1] <= bounds.x[2]
    inY = bounds.y[1] <= cell.pos[2] <= bounds.y[2]
    return inX && inY
end

function inBounds(p::MVector{2,Float64}, bounds::Bounds)
    inX = bounds.x[1] <= p[1] <= bounds.x[2]
    inY = bounds.y[1] <= p[2] <= bounds.y[2]
    return inX && inY
end

function inBounds(r::Float64, b::Tuple{Real,Real})
    inB = b[1] <= r <= b[2]
end

"""
Move input vector to equivalent position in bounds.
"""
function toBoundsPeriodic!(p::AbstractVector, bounds)
    if p[1] < bounds.x[1]
        p[1] = bounds.x[2] - (bounds.x[1]-p[1])%bounds.xLen
    elseif p[1] > bounds.x[2]
        p[1] = bounds.x[1] + (p[1]-bounds.x[2])%bounds.xLen
    end
    if p[2] < bounds.y[1]
        p[2] = bounds.y[2] - (bounds.y[1]-p[2])%bounds.yLen
    elseif p[2] > bounds.y[2]
        p[2] = bounds.y[1] + (p[2]-bounds.y[2])%bounds.yLen
    end
    nothing
end

"""
Create vector with in bounds position equivalent to input vector.
"""
function toBoundsPeriodic(p::AbstractVector, bounds)
    if p[1] < bounds.x[1]
        px = bounds.x[2] - (bounds.x[1]-p[1])%bounds.xLen
    elseif p[1] > bounds.x[2]
        px = bounds.x[1] + (p[1]-bounds.x[2])%bounds.xLen
    else
        px = p[1]
    end
    if p[2] < bounds.y[1]
        py = bounds.y[2] - (bounds.y[1]-p[2])%bounds.yLen
    elseif p[2] > bounds.y[2]
        py = bounds.y[1] + (p[2]-bounds.y[2])%bounds.yLen
    else
        py = p[2]
    end
    return @MVector[px, py]
end

"""
Get equivalent point in bounds.
"""
function toBoundsPeriodic(x::Real, xbounds::Tuple{Real, Real})
    if x < xbounds[1]
        x = xbounds[2] - (xbounds[1]-x)%(xbounds[2]-xbounds[1])
    elseif x > xbounds[2]
        x = xbounds[1] + (x-xbounds[2])%(xbounds[2]-xbounds[1])
    end
    return x
end
