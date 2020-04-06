
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

function connectingVector(p1::AbstractVector, p2::AbstractVector, period::AbstractVector)
    dx = p1[1]-p2[1]
    dy = p1[2]-p2[2]

    if (abs(dx) > period[1]/2)
        dx = -sign(dx)*(period[1] - abs(dx))
    end
    if (abs(dy) > period[2]/2)
        dy = -sign(dy)*(period[2] - abs(dy))
    end

    return @MVector[dx, dy]
end
