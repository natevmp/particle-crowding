using StaticArrays, LinearAlgebra

##
include("../src/bmparticles.jl")
using .BParts
##

myArena = BParts.buildRandArena(((0,10),(0,1)), 50, 0.08, 0.02, fixSpeed=true);

comVel = BParts.comVel(myArena.cellsList)
E = BParts.avParticleEnergy(myArena.cellsList)
println("before energy: ", E)
println("before COM speed: ", norm(comVel))

BParts.rmComDriftArena!(myArena.cellsList)

comVel = BParts.comVel(myArena.cellsList)
E = BParts.avParticleEnergy(myArena.cellsList)
println("after energy: ", E)
println("after COM speed: ", norm(comVel))

println([norm(c.vel) for c in myArena.cellsList])
# function incrementA!(a, b)
#     for i in 1:10^5
#         a += b
#     end
#     return nothing
# end

# function incrementB!(a, b)
#     for i in 1:10^5
#         a .+= b
#     end
#     return nothing
# end

# testa = [0.,0.]
# testinc = [0.1,0.2]
# println("intiate:")
# @time incrementA!(testa, testinc)
# testb = [0.,0.]
# @time incrementB!(testb, testinc)

# println("===== real tests =====")
# testa = [0.,0.]
# testinc = [0.1,0.2]
# @time incrementA!(testa, testinc)

# testb = [0.,0.]
# @time incrementB!(testb, testinc)

# println(testa)
# println(testb)