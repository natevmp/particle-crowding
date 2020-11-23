

"""Add new cells to given arena according to specified growth rate."""
function cultivateArena!(arena::Arena, dt::Real, grFunc::Function, r::Real, speed::Real; randGrowth::Bool=false)

    nPop = length(arena.cellsList)
    
    if randGrowth
        # draw number of new cells in timestep
        n = pois_rand(grFunc(nPop)*dt)
    else
        n = Integer(round(grFunc(nPop)))
    end
    if n > 0
        newcells_c = Vector{Cell}(undef, n)
        #create n random cells
        for cid in 1:n
            newcells_c[cid] = randCell(arena.bounds, r, speed, fixSpeed=true)
        end
        
        # === find and fix potential overlaps ===
        # - create nearest neighbor tree of existing cell cellPositions_DIM_ID
        oldCellTree = BallTree(cellPositions_DIM_ID(arena), PeriodicEuclidean([arena.bounds.xLen, arena.bounds.yLen]))
        
        verifiedBool_c = falses(n)

        while sum(verifiedBool_c)<n
            for (i, cell) in enumerate(newcells_c)
                # skip step if cell is already verified
                if verifiedBool_c[i]
                    continue
                end

                # get nearest existing cell
                nearestDistOld = knn(oldCellTree, cell.pos, 1)[2][1]

                # get nearest verified new cell
                newCellsDists_c = map(cellB -> cellDistance(cell, cellB, arena.bounds), newcells_c[verifiedBool_c])
                if length(newCellsDists_c)>0
                    sort!(newCellsDists_c)
                    nearestDistNew = newCellsDists_c[1]
                else nearestDistNew = NaN end

                # check if cell does not overlap existing cells, or verified new cells
                if nearestDistOld < collisionCrossSection(arena) || nearestDistNew < collisionCrossSection(arena)
                    #change cell's location
                    randPos!(newcells_c[i], arena.bounds)
                else
                    verifiedBool_c[i] = true
                end
            end
        end

        # === add cells to arena ===
        append!(arena.cellsList, newcells_c)
    end

    return n
end

logisticRate(n::Real, ρ::Real, k::Real) = n*ρ*(1-n/k)
exponentialRate(n::Real, ρ::Real) = n*ρ
