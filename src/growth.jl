

"""Add new cells to given arena according to specified growth rate."""
function cultivateArena!(arena::Arena, dt::Real, rate::Real, r::Real, vStd::Real)

    nPop = length(arena.cellsList)
    
    # draw number of new cells in timestep
    n = pois_rand(rate*dt)
    if n > 0
        newcells_c = Vector{Cell}(undef, n)
        #create n random cells
        for cid in 1:n
            newcells_c[cid] = randCell(arena.bounds, r, vStd)
        end
        
        # === find and fix potential overlaps ===
        # - create nearest neighbor tree of existing cell cellPositions_DIM_ID
        tree = BallTree(cellPositions_DIM_ID(arena), PeriodicEuclidean([arena.bounds.xLen, arena.bounds.yLen]))
        # - check for overlaps
        verifiedBool_c = falses(n)
        while sum(verifiedBool_c)<n
            unverifiedId_c = findall(.!verifiedBool_c)
            # calling knn for multiple points at the same time is more efficient than calling it for each point seperately
            idxs, dists = knn(tree, cellPositions_DIM_ID(newcells_c[unverifiedId_c]), 1)
            for (dInd, dist) in enumerate(dists)
                if dist[1] > collisionCrossSection(arena)
                    verifiedBool_c[unverifiedId_c[dInd]] = true
                else
                    #change cell's location
                    randPos!(newcells_c[unverifiedId_c[dInd]], arena.bounds)
                end
            end
        end

        # === add cells to arena ===
        append!(arena.cellsList, newcells_c)
    end

    return n
end