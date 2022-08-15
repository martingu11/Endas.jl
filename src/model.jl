
export 
    AbstractModel, LinearizedModel



"""
Abstract forward model for propagating state from current time to time t + δt.

AbstractModel represents the most generic models such that only define a (possibly non-linear)
function `apply(model, x, t, δt)`. 

# Methods

Abstract models must implement the following methods

    apply(ℳ::AbstractModel, x::AbstractStateVector, t::Real, δt::Real)

Applies model ℳ to state vector x, propagating the state from time t to t + δt. 
The state vector x is updated in-place.
"""
abstract type AbstractModel end



"""
    apply!(X::AbstractMatrix{<:Real}, ℳ::AbstractModel, t::Real, δt::Real)

Apply model to an ensemble. This is a generic implementation that calls `apply!()` on each 
ensemble member (matrix column). The ensemble `X` is updated in-place. Returns `nothing` 
even if the model returned linearization data since ensemble methods do not rely on model
linearization.
"""
function apply!(X::AbstractMatrix{<:Real}, ℳ::AbstractModel, t::Real, δt::Real)
    m = size(X)[2]
    @assert size(X)[1] == statesize(ℳ)

    for i in 1:m
        apply!(view(X, :, 1), ℳ, t, δt)
    end
    return nothing
end



"""
Abstract linearized forward model for propagating state from current time to time t + δt.

This is a subtype of the `AbstractModel`` and in addition to the `apply()` function, 
linearized models must also implement the tangent linear (via `apply_tl()`) and adjoint 
(via `apply_adj()`) of the model. 
"""
abstract type LinearizedModel <: AbstractModel end



"""
    apply_tl!(X::AbstractMatrix{<:Real}, ℳ::LinearizedModel, Mlin::AbstractMatrix{<:Real}) 

Default implementation of tangent linear operator for a linearized model.

If the linearization data returned from `apply!()` is a matrix, it is assumed it represents 
the linearized operator as a matrix transform.
"""
function apply_tl!(X::AbstractMatrix{<:Real}, ℳ::AbstractModel, Mlin::AbstractMatrix{<:Real}) 
    X .= Mlin * X
end


"""
    apply_adj!(X::AbstractMatrix{<:Real}, ℳ::LinearizedModel, Mlin::AbstractMatrix{<:Real}) 

Default implementation of adjoint operator for a linearized model.

If the linearization data returned from `apply!()` is a matrix, it is assumed it 
represents the linearized operator as a matrix transform.
"""
function apply_adj!(X::AbstractMatrix{<:Real}, ℳ::AbstractModel, Mlin::AbstractMatrix{<:Real}) 
    X .= X * Mlin'
end




# Model implementations

include("models/lorenz95.jl")

