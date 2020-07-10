module Theorist
using DifferentialEquations, LinearAlgebra, Distributions, Distances

# ===== Analytical solutions =====
function meanFreePath(n::Number, σ::Number)
    1. / (sqrt(2.)*σ*n)
end

function meanFreePath(t::Number, n::Function, σ::Number)
    1. / ( sqrt(2.)*σ*n(t) )
end

function friction(l::Number, E::Number)
    sqrt(E)/l * sqrt(π/2)
end

function friction(t::Number, l::Function, E::Number)
    sqrt(E)/l(t) * sqrt(π/2)
end

function velocityAutoCorrelation(τ::Number, E::Number, γ::Number)
    2E * exp(-γ*τ)
end

function velocityAutoCorrelation(τ::Number, t::Number, E::Number, γ::Function)
    2E * exp(-γ(t)*τ)
end

function velocityAutoCorrelation(τ::Number, t::Number, E::Function, γ::Function)
    2E(t) * exp(-γ(t)*τ)
end


# ====== System data and values ======

function thermalValues(arenaParams::Dict, growthFunc::Union{Function, Nothing}=nothing)
    thermVals = Dict()
    radius = arenaParams["radius"]
    s0Av = arenaParams["speed"]
    bounds = arenaParams["bounds"]
    volume = arenaParams["volume"]
    n0 = arenaParams["n0"]
    σc = 2*2radius
    E = s0Av^2 / 2
    thermVals["σc"] = σc
    thermVals["E"] = E
    if isnothing(growthFunc)
        n = n0 / volume
        ρ = n*π*radius^2
        l = 1/(√2*σc*n)
        γ = √(E)/l * √(π/2)
        D = E*γ
        thermVals["n"] = n
        thermVals["l"] = l
        thermVals["γ"] = γ
        thermVals["D"] = D
        thermVals["ρ"] = ρ
    end
    return thermVals
end

# ===== Langevin numeric simulation =====

# --- Problem definitions ---
function brownianLangevin(n::Function, E, σ, u0, tspan)

    l(t) = 1/(√2*σ*n(t))
    γ(t) = √(E)/l(t) * √(π/2)
    D(t) = E*γ(t)

    function drift(du, u, p, t)
        du[1] = - γ(t)*u[1]
        du[2] = - γ(t)*u[2]
    end

    function diff(du, u, p, t)
        du[1] = √(2*D(t))
        du[2] = √(2*D(t))
    end

    return SDEProblem(drift, diff, u0, tspan)
end

function logisticGrowth(t, ρ, k, n0)
    return k/( 1 + (k-n0)/n0 * exp(-ρ*t) )
end

# --- Experiments ---
function runLangevinSims(runs, arenaParams::Dict, growthParams::Union{Dict, Nothing}=nothing)
    n0 = arenaParams["n0"]
    evolveTime = arenaParams["evolveTime"]
    radius = arenaParams["radius"]
    s0 = arenaParams["speed"]
    bounds = arenaParams["bounds"]
    volume = arenaParams["volume"]
    σc = 2*2radius
    E = s0^2 / 2
    println("Energy: "*string(E))

    nDensity(t) = n0 / volume
    if isnothing(growthParams)
        ρ = 0
        k = n0
    else
        ρ = growthParams["ρ"]
        k = growthParams["k"]
        nDensity(t) = logisticGrowth(t, ρ, k, n0) / volume
    end

    u0 = [s0*cos(π), s0*sin(π)]
    tspan = (0., evolveTime)
    langProb = brownianLangevin(nDensity, E, σc, u0, tspan);

    function randSpeedProblem(prob,i,repeat)
        # s = rand(  Rayleigh( norm(prob.u0)*√(2/π) )  )
        s = s0
        # println(s)
        α = rand()*2π
        @. prob.u0 = [s*cos(α), s*sin(α)]
        return prob
    end

    ensembleprob = EnsembleProblem(langProb; prob_func=randSpeedProblem)
    ensSol = solve(ensembleprob, SOSRI(), EnsembleThreads(); trajectories=runs)
    return ensSol
end


# --- get positions from velocities ---
function toBounds(val, b1, b2)
    if val < b1
        return b2 + (val - b1)%(b2-b1)
    elseif val > b2
        return b1 + (val - b1)%(b2-b1)
    else
        return val
    end
end

function velToPosition(vSol::RODESolution, pos0, times)
    pos_t_xy = Array{Float64, 2}(undef, length(times), 2)
    pos_t_xy[1, :] = pos0
    for (i, t) in enumerate(times[2:end])
        pos_t_xy[1+i, :] = pos_t_xy[i, :] .+ vSol(t) .*(times[i+1]-times[i])
    end
    return pos_t_xy
end

function velToPosition(vSol::RODESolution, pos0, bounds_Dim_Val, times)
    pos_t_xy = Array{Float64, 2}(undef, length(times), 2)
    pos_t_xy[1, :] = pos0
    # println(bounds_Dim_Val)
    for (i, t) in enumerate(times[2:end])
        pos_t_xy[1+i, :] =
        [
        toBounds(pos_t_xy[i, 1] + vSol(t)[1] *(times[i+1]-times[i]),
            bounds_Dim_Val[1][1], bounds_Dim_Val[1][2]),
        toBounds(pos_t_xy[i, 2] + vSol(t)[2] *(times[i+1]-times[i]),
            bounds_Dim_Val[2][1], bounds_Dim_Val[2][2]),
        ]
    end
    return pos_t_xy
end

function velToPosition(ensSol::EnsembleSolution, p0, times)
    pos_Traj_t_xy = Array{Array{Float64, 2}, 1}(undef, length(ensSol))
    for (i, vTraj) in enumerate(ensSol)
        pos_Traj_t_xy[i] = velToPosition(vTraj, p0, times)
    end
    return pos_Traj_t_xy
end

function velToPosition(ensSol::EnsembleSolution, p0, bounds_Dim_Val, times)
    pos_Traj_t_xy = Array{Array{Float64, 2}, 1}(undef, length(ensSol))
    for (i, vTraj) in enumerate(ensSol)
        pos_Traj_t_xy[i] = velToPosition(vTraj, p0, bounds_Dim_Val, times)
    end
    return pos_Traj_t_xy
end

function speedCalc(vEnsSol::EnsembleSolution, times)
    speed_t_id = Array{Float64, 2}(undef, length(times), length(vEnsSol))
    for (tt, time) in enumerate(times)
        for (j, vTraj) in enumerate(vEnsSol)
            speed_t_id[tt, j] = norm(vTraj(time))
        end
    end
    return speed_t_id
end

"""Calculate mean squared displacement over time of ensemble solutions"""
function msd(ensSol, tspan::Tuple{Real, Real}, steps::Int=100)

    times_t = range(tspan[1], tspan[2]; length=steps)

    p0 = [0., 0.]
    pos_Traj_t_xy = velToPosition(ensSol, p0, times_t)
    msd_t = Array{Float64, 1}(undef, steps)

    for (i,t) in enumerate(times_t)
        msd_t[i] = mean([norm(pos_t_xy[i,:]-pos_t_xy[1,:])^2 for pos_t_xy in pos_Traj_t_xy])
    end

    return times_t, msd_t
end

function msd(ensSol, arenaParams::Dict, tspan::Tuple{Real, Real}, steps::Int=100)

    times_t = range(tspan[1], tspan[2]; length=steps)

    p0 = [0., 0.]

    pos_Traj_t_xy = velToPosition(ensSol, p0, arenaParams["bounds"], times_t)
    msd_t = Array{Float64, 1}(undef, steps)

    for (i,t) in enumerate(times_t)
        msd_t[i] = mean([peuclidean(pos_t_xy[1,:], pos_t_xy[i,:], arenaParams["bperiod"])^2
                        for pos_t_xy in pos_Traj_t_xy])
    end

    return times_t, msd_t
end

function velCorrelation(ensVSol::EnsembleSolution, tspan::Tuple{Real, Real}; steps::Int=100)
    times_t = range(tspan[1], tspan[2]; length=steps)
    vCorr_t = Array{Float64, 1}(undef, steps)
    for (i,t) in enumerate(times_t)
        vCorr_t[i] = mean( [ vSol(t)⋅vSol(tspan[1]) for vSol in ensVSol ] )
    end
    return times_t, vCorr_t
end

end
