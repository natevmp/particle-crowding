module Theorist

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

function velocityAutoCorrelation(τ::Number, E::Number, γ::Number)
    2E * exp(-γ*τ)
end

function velocityAutoCorrelation(τ::Number, t::Number, E::Number, γ::Function)
    2E * exp(-γ(t)*τ)
end

function velocityAutoCorrelation(τ::Number, t::Number, E::Function, γ::Function)
    2E(t) * exp(-γ(t)*τ)
end




end
