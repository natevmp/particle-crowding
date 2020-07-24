
function nestfor(looplengths::Tuple{Vararg{Int}}, loopstates::Vector{Int}, loopdepth, func::Function, args...)
    if loopdepth < length(looplengths)
        for i in 1:looplengths[1]
            loopstates[loopdepth] = i
            nestfor(looplengths, loopstates, loopdepth+1, func, args...)
        end
    else
        for i in 1:looplengths[end]
            loopstates[loopdepth] = i
            func(loopstates, args...)
        end
    end
end

function nestfor(looplengths::Tuple{Vararg{Int}}, func::Function, args...)
    nestfor(looplengths, fill(0, length(looplengths)), 1, func, args...)
end

function distance(p1::Vector{Float64}, p2::Vector{Float64}, period::Vector{Float64})
    peuclidean(p1, p2, period)
end

function connectingVector(pA::AbstractVector, pB::AbstractVector, bounds::BParts.Bounds)
    p1 = BParts.toBoundsPeriodic(pA, bounds)
    p2 = BParts.toBoundsPeriodic(pB, bounds)

    dx_InUpDown = [
        p1[1] - p2[1],
        p1[1]+bounds.xLen - p2[1],
        p1[1]-bounds.xLen - p2[1]
        ]
    indx = argmin(abs.(dx_InUpDown))
    dx = dx_InUpDown[indx]

    dy_InUpDown = [
        p1[2] - p2[2],
        p1[2]+bounds.yLen - p2[2],
        p1[2]-bounds.yLen - p2[2]
        ]
    indy = argmin(abs.(dy_InUpDown))
    dy = dy_InUpDown[indy]

    return @MVector[dx, dy]
end
