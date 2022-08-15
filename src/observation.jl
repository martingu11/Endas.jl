

export ObservationOperator, LinearizedObservationOp, LinearObservationOp


"""
Abstract observation operator.
"""
abstract type ObservationOperator end


"""
Abstract linearized observation operator.
"""
abstract type LinearizedObservationOp <: ObservationOperator end



"""
    apply_tl!(out::AbstractVecOrMat{<:Real}, H::LinearizedObservationOp, 
              δH::AbstractMatrix{<:Real}, x::AbstractVecOrMat{<:Real})

Default implementation of tangent linear operator for a linearized observation operator.

If the linearization data returned from `apply(H, x, t)` is a `Matrix`, it is assumed it 
represents the linearized operator as a matrix transform.
"""
function apply_tl!(out::AbstractVecOrMat{<:Real}, ℋ::LinearizedObservationOp, 
                   Hlin::AbstractMatrix{<:Real}, x::AbstractVecOrMat{<:Real})

    mul!(out, Hlin, x)
end


"""
Default implementation of adjoint operator for a linearized observation operator.

If the linearization data returned from `apply!` is a matrix, it is assumed it 
represents the linearized operator as a matrix transform.
"""
function apply_adj!(out::AbstractVecOrMat{<:Real}, ℋ::LinearizedObservationOp, 
                    Hlin::AbstractMatrix{<:Real}, x::AbstractVecOrMat{<:Real})
    mul!(out, x, Hlin')
end


"""
Linear observation operator wrapping a (possibly sparse) matrix that implements the 
transformation.

"""
struct LinearObservationOp <: LinearizedObservationOp
    coeffs::AbstractMatrix{<:Real}
end


# Basic info
function Base.show(io::IO, H::LinearObservationOp)  
    print(io, "LinearObservationOp(nobs:$(size(H)[1]), nstate:$(size(H)[2]))")
end

obssize(H::LinearObservationOp) = size(H.coeffs)[1]
statesize(H::LinearObservationOp) = size(H.coeffs)[2]
asmatrix(H::LinearObservationOp) = H.coeffs

function apply!(out::AbstractVecOrMat{<:Real}, H::LinearObservationOp, 
                x::AbstractVecOrMat{<:Real}, t::Real)
    
    mul!(out, H.coeffs, x)

    # The matrix itself is the linearized version, let the base type handle the tl and 
    # adjoint operators
    return H.coeffs  
end
