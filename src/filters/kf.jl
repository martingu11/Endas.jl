using LinearAlgebra


export KalmanFilter

"""
Internal data for the Kalman filter and smoother.
"""
mutable struct KFData
    x::Vector{Float64}
    P::Matrix{Float64}
    t::Float64
end


"""
Kalman filter and smoother algorithm.
"""
struct KalmanFilter <: AbstractFilter
    ℳ::LinearizedModel
    _data::KFData

    function KalmanFilter(ℳ::LinearizedModel, x₀::Vector{<:Real}, P₀::Matrix{<:Real}, t::Real)
        data = KFData(x₀, P₀, t)
        return new(ℳ, data)
    end
end

getmodel(kf::KalmanFilter) = kf.ℳ
getstate(kf::KalmanFilter) = kf._data.x
putstate!(kf::KalmanFilter, x::Vector{<:Real}) = kf._data.x = x
getcovariance(kf::KalmanFilter) = kf._data.P
gettime(kf::KalmanFilter) = kf._data.t


function forecast!(kf::KalmanFilter, δt::Real, Q::Union{CovarianceOperator, Nothing} = nothing)

    # Apply the model to the state vector (in-place)
    Mlin = apply!(kf._data.x, kf.ℳ, kf._data.t, δt, true)

    # Update the error covariance matrix: P = M P Mᵀ
    apply_tl!(kf._data.P, kf.ℳ, Mlin)
    apply_adj!(kf._data.P, kf.ℳ, Mlin)

    # Add the model error covariance term, if given
    if !isnothing(Q) 
        add!(kf._data.P, Q)
    end

    kf._data.t += δt
end


function assimilate!(kf::KalmanFilter, y::AbstractVector{<:Number}, 
                     ℋ::LinearizedObservationOp, R::CovarianceOperator)

    n = statesize(ℋ)
    k = obssize(ℋ)

    # Nothing to do
    if length(y) == 0 
        return nothing 
    end

    Hx = similar(kf._data.x, k)
    Hlin = apply!(Hx, ℋ, kf._data.x, kf._data.t)
    
    # Innovation covariance F = HPHᵀ + R 
    PHᵀ = similar(Hx, n, k)
    F = similar(Hx, k, k)

    apply_adj!(PHᵀ, ℋ, Hlin, kf._data.P)
    apply_tl!(F, ℋ, Hlin, PHᵀ)
    add!(F, R)
        
    # And its Cholesky factorization
    F = cholesky!(Symmetric(F))

    # State update x = x + PHᵀ F⁻¹ (y - Hx)
    z = y - Hx
    kf._data.x += PHᵀ * (F \ z)

    # Error covariance update P = P - PHᵀ F⁻¹ HP
    HP = similar(Hx, k, n)
    apply_tl!(HP, ℋ, Hlin, kf._data.P)
    kf._data.P -= PHᵀ * (F \ HP)
end





