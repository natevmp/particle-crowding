"""
Calculate mean squared displacement from simulations. Saves results in *_msd.jld2 files.
"""

##
include("../../src/bmparticles.jl")
include("../../src/bmtheory.jl")
using .BParts
using .Theorist
using JLD2, FileIO

## ===== MSD =====

volumeDens(n, r, v) = n*π*r^2/v
files_f = [
    "./data/FixedDensity/21-05-10/simResultFixedDensity_0.01.jld2",
    "./data/FixedDensity/21-05-10/simResultFixedDensity_0.05.jld2",
    "./data/FixedDensity/21-05-10/simResultFixedDensity_0.1.jld2",
    "./data/FixedDensity/21-05-10/simResultFixedDensity_0.2.jld2",
    "./data/FixedDensity/21-05-10/simResultFixedDensity_0.3.jld2",
    "./data/FixedDensity/21-05-10/simResultFixedDensity_0.5.jld2",
]

for (i,file) in enumerate(files_f)
    arenaParams, growthParams, pos_t_dim_id, vel_t_dim_id, times_t = load(file, "arenaParams", "growthParams", "pos_t_dim_id", "vel_t_dim_id", "times_t")

    thermVals = Theorist.thermalValues(arenaParams)
    colUnits = Dict(
        :t => 1 / thermVals["γ"],
        :x => arenaParams["speed"]/thermVals["γ"],
    )

    t0Ind = findall(t->t>growthParams["waitTime"], times_t)[1]
    msdPar_t = BParts.meanSquaredDisplacement(pos_t_dim_id[t0Ind:end,:,:])
    _t = times_t[times_t .> growthParams["waitTime"]] .- times_t[t0Ind]

    # ---- Theory ----
    msdTheory_t = Theorist.msdTheory.(_t, arenaParams["n0"]/arenaParams["volume"], thermVals["σc"], thermVals["E"])
    ρIn = round(volumeDens(arenaParams["n0"], arenaParams["radius"], arenaParams["volume"]), digits=2)

    # ===== Save data =====
    filename = "./data/FixedDensity/21-05-10/simResultFixedDensity_"*string(ρIn)*"_msd.jld2"
    save(filename, "arenaParams", arenaParams, "growthParams", growthParams, "_t", _t, "msdTheory_t", msdTheory_t, "msdPar_t", msdPar_t)

end