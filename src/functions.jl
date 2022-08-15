# 
# Common EnDAS functions.
#

export 
    statesize, obssize, getmodel, gettime, getstate, getcovariance,
    forecast!, assimilate!


"""
    statesize(f::AbstractFilter)
    statesize(ℳ::AbstractModel)
    statesize(ℋ::ObservationOperator)
    statesize(C::CovarianceOperator)

Return the size of the state space of the given filter or smoother, model, covariance
or observation operator.
"""
function statesize end

"""
    obssize(ℋ::ObservationOperator)

Return the size of the observation space of the given observation operator.
"""
function obssize end


"""
    getmodel(f::AbstractFilter)

Return the model instance of the given filter or smoother.
"""
function getmodel end


"""
    getstep(f::AbstractFilter)

Return the current time step of a sequential filter or smoother. Returns named tuple 
`(i, t)` where `i` is the index of the current time step and `t` is the corresponding
(physical) time. 
"""
function gettime end


"""
    getstate(f::AbstractFilter)

Return the current state vector of the given filter or smoother. The returned value is a 
one-dimensional array of size `statesize(f)`
"""
function getstate end

"""
    getcovariance(f::AbstractFilter)

Return the current covariance representation of 
"""
function getcovariance end


"""
    forecast!(f::AbstractFilter, δt::Real, 𝒬::Union{CovarianceOperator, Nothing)

Forecast step of sequential data assimilation.
"""
function forecast! end


"""
    assimilate!(f::AbstractFilter, y::AbstractVector{<:Number}, 
                H::ObservationOperator, R::CovarianceOperator)

Assimilation or analysis step in of sequential data assimilation.
"""
function assimilate! end

