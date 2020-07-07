### ======== dynamics functions ========

# function evolveArena!(arena::Arena, steps::Int, growthParams::Union{Tuple, Nothing}=nothing; plotsteps=false,
#     animator::Union{Animation, Nothing}=nothing)
#     # LEGACY

#     posTime_t_dim_id = ElasticArray{Union{Float64, Missing}}(fill(missing, steps, 2, length(arena.cellsList)))
#     velTime_t_dim_id = ElasticArray{Union{Float64, Missing}}(fill(missing, steps, 2, length(arena.cellsList)))

#     cells_T_ID = Vector{Vector{Cell}}(undef, steps)


#     @showprogress for t in 1:steps

#         # == Movement ==
#         for cell in arena.cellsList
#             moveCell!(cell, arena.bounds)
#         end

#         # == Collisions ==
#         collidedCellsList_id = collider!(arena)

#         # == Growth ==
#         if growthParams != nothing
#             nDaughters = cultivateArena!(arena, 1., growthParams...)
#         end

#         # == Plot data ==
#         if plotsteps||(animator!==nothing)
#             plotArena(arena, collidedCellsList_id, nDaughters, displayPlot=plotsteps, animator=animator)
#         end

#         # == record data ==
#         snapshotCells!(posTime_t_dim_id, velTime_t_dim_id, arena, t)
#         snapshotCells!(cells_T_ID, arena, t)
#     end
#     # gif(anim, "anim_2.gif", fps=15)
#     return posTime_t_dim_id, velTime_t_dim_id, cells_T_ID
# end


function evolveArena!(arena::Arena, steps::Int, growthParams::Union{Dict, Nothing}=nothing; plotsteps=false,
    animator::Union{Animation, Nothing}=nothing, progress=true, verbose=true)

    posTime_t_dim_id = ElasticArray{Union{Float64, Missing}}(fill(missing, steps, 2, length(arena.cellsList)))
    velTime_t_dim_id = ElasticArray{Union{Float64, Missing}}(fill(missing, steps, 2, length(arena.cellsList)))

    cells_T_ID = Vector{Vector{Cell}}(undef, steps)

    # progresMeter
    if progress
        prog = Progress(steps, 1)
    end
    for t in 1:steps

        # == Movement ==
        for cell in arena.cellsList
            moveCell!(cell, arena.bounds)
        end

        # == Collisions ==
        collidedCellsList_id = collider!(arena; verbose=verbose)

        # == Growth ==
        if growthParams != nothing
            if growthParams["randGrowth"]
                nDaughters = cultivateArena!(
                    arena, 1., growthParams["rateFunc"],
                    growthParams["radius"],
                    growthParams["speed"], randGrowth=true)
            else
                # popSizeNew = growthParams["growthFunc"](t)
                rateFunc(nPop) = growthParams["growthFunc"](t) - nPop
                nDaughters = cultivateArena!(
                    arena, 1., rateFunc,
                    growthParams["radius"],
                    growthParams["speed"], randGrowth=false)
            end
        end

        # == Plot data ==
        if plotsteps||(animator!==nothing)
            plotArena(arena, collidedCellsList_id, nDaughters, displayPlot=plotsteps, animator=animator)
        end

        # == record data ==
        snapshotCells!(posTime_t_dim_id, velTime_t_dim_id, arena, t)
        snapshotCells!(cells_T_ID, arena, t)

        # progressMeter
        if progress
            next!(prog)
        end
    end
    # gif(anim, "anim_2.gif", fps=15)
    return posTime_t_dim_id, velTime_t_dim_id, cells_T_ID
end
