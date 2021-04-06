using NearestNeighbors, Distances

##

struct CollQ
    coll_t_Cid::Vector{SVector{2,Int}}
    _t::Vector{Float64}

    function CollQ()
        coll_t_Cid = SVector{2,Int}[]
        _t = Float64[]
        new(coll_t_Cid, _t)
    end
end
Q = CollQ()


##


function findPosInQueue(_t, T)
    if length(_t)==0
        return 1
    end
    for (i,tt) in enumerate(_t)
        if tt>T
            return i
        end
    end
    return length(_t)+1
end

_t = [1.,2,3,4,5,6,7,8]

tt = 14
ind = findPosInQueue(_t, tt)

insert!(_t, ind, tt)
println(_t)

##
pZero_dim_id = rand(2, 20)
# Create trees
# kdtree = KDTree(data; leafsize = 10)
balltree = BallTree(pZero_dim_id, PeriodicEuclidean([0.1,0.1]); reorder = false)
# brutetree = BruteTree(data)
##
overlaps_l_Cid = inrange(balltree, pZero_dim_id, 0.01)
display(overlaps_l_Cid)

for (l, l_cid) in enumerate(overlaps_l_Cid)
    if length(l_cid) < 2
        continue
    elseif length(l_cid)==2
        tRetc = collTimeCalc(arena.cellsList[l_cid]..., arena.bounds)
    elseif length(l_cid)>2

    end
end


# collQ_lid_Cid
##

@time for i in 1:100 BallTree(data, PeriodicEuclidean([1,1]); reorder = false) end

# @time for i in 1:200 BallTree(data, PeriodicEuclidean([1,1]); reorder = true) end


##

testa = [1,2,3,4,5,6,7,8,9]

for (i, el) in enumerate(testa)
    el % 2 == 0 && deleteat!(testa, i+1)
    println(i)
end

println(testa)


##

testa = [3,8]

testb = [7,9]

intersect(testa, testb)

##
testa = [1,2,3,4,5]
append!(testa, [6,7,8])
println(testa)

##

testp = [1,2]
testa_p = [testp, testp]

testp[1] = 5

display(testa_p)