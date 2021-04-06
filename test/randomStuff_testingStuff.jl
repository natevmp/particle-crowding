include("../src/bmparticles.jl")

using StaticArrays
using .BParts


function testfunc(vArg_Id::Vararg{String, 3})
    for id in vArg_Id
        println(id)
    end
end

testa = (1,2,3,4)
testfunc(testa...)
