### ======== dynamics functions ========

# function evolveArena(arena::Arena, steps::Int; plotsteps=false, animator::Union{Animation, Nothing}=nothing)

#     posTime_id_t_dim = fill(NaN, length(arena.cellsList), steps, 2)
#     velTime_id_t_dim = fill(NaN, length(arena.cellsList), steps, 2)
#     cells_T_ID = Vector{Vector{Cell}}(undef, steps)

#     for t in 1:steps
#         for cell in arena.cellsList
#             moveCell!(cell, arena.bounds)
#         end

#         collidedCellsList_id = collider!(arena)
#         if plotsteps||(animator!==nothing)
#             plotArena(arena, collidedCellsList_id, displayPlot=plotsteps, animator=animator)
#         end

#         snapshotCells!(posTime_id_t_dim, velTime_id_t_dim, arena, t)
#         snapshotCells!(cells_T_ID, arena, t)
#     end
#     # gif(anim, "anim_2.gif", fps=15)
#     return posTime_id_t_dim, velTime_id_t_dim, cells_T_ID
# end

function evolveArena!(arena::Arena, steps::Int, growthParams::Union{Tuple{Real,Real,Real}, Nothing}=nothing; plotsteps=false,
    animator::Union{Animation, Nothing}=nothing)

    # posTime_id_t_dim = fill(NaN, length(arena.cellsList), steps, 2)
    # velTime_id_t_dim = fill(NaN, length(arena.cellsList), steps, 2)

    posTime_t_dim_id = ElasticArray(fill(NaN, steps, 2, length(arena.cellsList)))
    velTime_t_dim_id = ElasticArray(fill(NaN, steps, 2, length(arena.cellsList)))

    cells_T_ID = Vector{Vector{Cell}}(undef, steps)


    for t in 1:steps

    # == Movement ==
    for cell in arena.cellsList
    moveCell!(cell, arena.bounds)
    end

    # == Collisions ==
    collidedCellsList_id = collider!(arena)

    # == Growth ==
    if growthParams != nothing
    nDaughters = cultivateArena!(arena, 1., growthParams...)
    end

    # == Plot data ==
    if plotsteps||(animator!==nothing)
    plotArena(arena, collidedCellsList_id, nDaughters, displayPlot=plotsteps, animator=animator)
    end

    # == record data ==
    snapshotCells!(posTime_t_dim_id, velTime_t_dim_id, arena, t)
    snapshotCells!(cells_T_ID, arena, t)
    end
    # gif(anim, "anim_2.gif", fps=15)
    return posTime_t_dim_id, velTime_t_dim_id, cells_T_ID
end