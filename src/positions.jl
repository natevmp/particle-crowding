
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

function cellPositions_DIM_ID(arena::Arena)
    positions_dim_id = Array{Float64,2}(undef, 2, length(arena.cellsList))
    for (i, cell) in enumerate(arena.cellsList)
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
    evaluate(PeriodicEuclidean([bounds.xLen, bounds.yLen]), cellA.pos, cellB.pos)
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

function toBoundsPeriodic!(p::MVector{2,Float64}, bounds)
    if p[1] < bounds.x[1]
        p[1] = bounds.x[2] - (bounds.x[1]-p[1])
    elseif p[1] > bounds.x[2]
        p[1] = bounds.x[1] + (p[1]-bounds.x[2])
    end
    if p[2] < bounds.y[1]
        p[2] = bounds.y[2] - (bounds.y[1]-p[2])
    elseif p[2] > bounds.y[2]
        p[2] = bounds.y[1] + (p[2]-bounds.y[2])
    end
    nothing
end
