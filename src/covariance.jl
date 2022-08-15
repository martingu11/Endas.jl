import LinearAlgebra
import Random


export 
    CovarianceOperator, DenseCovariance, DiagonalCovariance
    

"""
Abstract covariance operator.

The operator represents a covariance matrix.
"""
abstract type CovarianceOperator end 

"""
Return true if the covariance operator is 'diagonal', i.e. it represents uncorrelated
errors.
"""
isdiagonal(C::CovarianceOperator) = false


"""
Return true if the covariance operator is exact and thus supports all operations.
"""
isexact(C::CovarianceOperator) = false


"""
Full-rank covariance represented explicitly by a covariance matrix.

This is a wrapper around a dense matrix that implements the `CovarianceOperator` API.
Full-rank covariance matrices should only be used for sufficiently small state or 
observation spaces.
"""
struct DenseCovariance <: CovarianceOperator
    cov_chol::LinearAlgebra.Cholesky

    DenseCovariance(cov_chol::LinearAlgebra.Cholesky) = new(cov_chol)
    
    function DenseCovariance(cov::AbstractMatrix{<:Real})
        LinearAlgebra.checksquare(cov)
        op = new(LinearAlgebra.cholesky(cov))
    end

end


"""
Create covariance matrix from symmetric covariance function.
    
Assuming one-dimensional space spanning 1…n, the variance function `varfn(d)` is evaluated
for all i, j ∈ 1…n where `d` is the distance between cells `i` and `j`. The distance is 
computed as `d(i, j) = cellsize * abs(i - j)`.

The variance function `varfn` must have a method `varfn(d::Number) -> Number`. The size of 
the returned covariance matrix is n×n.
"""
function DenseCovariance(varfn, n::Integer, cellsize::Number = 1)
    C = Matrix{Float64}(undef, n, n)
    for i = 1:n, j=i:n
        C[i, j] = C[j, i] = varfn(abs(i - j) * cellsize)
    end
    return DenseCovariance(C)
end

"""
Create covariance matrix from symmetric covariance function.
    
This version is similar as the one-dimensional case but works on a two-dimensional space 
spanning 1…n₁ and 1…n₂. The size of the returned covariance matrix is n×n where n = n₁ * n₂.
"""
function DenseCovariance(varfn, n::Tuple{Integer, Integer}, cellsize::Tuple{Number, Number} = (1, 1))
    N = n[1] * n[2]    
    C = Matrix{Float64}(undef, N, N)
    for i = 1:N, j=i:N
        iy, ix = (i-1) % n[1], (i-1) ÷ n[2]
        jy, jx = (j-1) % n[1], (j-1) ÷ n[2]
        dy = cellsize[1]*(iy-jy)
        dx = cellsize[2]*(ix-jx)
        C[i, j] = C[j, i] = varfn(sqrt(dx^2 + dy^2))
    end
    return DenseCovariance(C)
end


# Basic info
Base.show(io::IO, C::DenseCovariance) = (n = size(C)[1]; print(io, "DenseCovariance($(n)x$(n))"))
Base.size(C::DenseCovariance) = size(C.cov_chol)
isdiagonal(C::DenseCovariance) = false
isexact(C::DenseCovariance) = true
asmatrix(C::DenseCovariance) = C.cov_chol.U' * C.cov_chol.U 

# Solve
solve!(b::AbstractMatrix, C::DenseCovariance) = ldiv!(C.cov_chol, b)

# Add
add!(b::AbstractMatrix, C::DenseCovariance) = broadcast!(+, asmatrix(C), b, b)

# Sampling from the covariance
function random_multivariate_normal!(out::AbstractVecOrMat{<:Real}, C::DenseCovariance)
    n = size(out)[1]
    @assert size(C) == (n, n)
    Random.randn!(Random.GLOBAL_RNG, out)
    # Using U' since L needs to be constructed (=expensive)
    return LinearAlgebra.lmul!(C.cov_chol.U', out)
end



"""
Covariance represented by a sparse diagonal matrix.

This is a wrapper around a diagonal matrix that implements the `CovarianceOperator` API.
"""
struct DiagonalCovariance <: CovarianceOperator
    diag::AbstractVector{<:Number}
end

function DiagonalCovariance(var::Number, n::Integer)
    return DiagonalCovariance(fill(var, n))
end


# Basic info
Base.show(io::IO, C::DiagonalCovariance) = (n = size(C.diag)[1]; print(io, "DiagonalCovariance($(n)x$(n))"))
Base.size(C::DiagonalCovariance) = (n = size(C.diag)[1]; (n, n))
isdiagonal(C::DiagonalCovariance) = true
isexact(C::DiagonalCovariance) = true
asmatrix(C::DiagonalCovariance) = LinearAlgebra.Diagonal(C.diag)

# Solve
solve!(b::AbstractMatrix, C::DiagonalCovariance) = ldiv!(LinearAlgebra.Diagonal(C.diag), b)

# Add
function add!(b::AbstractMatrix, C::DiagonalCovariance) 
    for i in 1:length(C.diag)
        b[i, i] += C.diag[i]
    end
end

# Sampling from the covariance
function random_multivariate_normal!(out::AbstractVecOrMat{<:Real}, C::DiagonalCovariance)
    n = size(out)[1]
    #m = ndims(out) == 1 ? 1 : size(out)[2]
    @assert size(C) == (n, n)

    Random.randn!(Random.GLOBAL_RNG, out)
    out .*= C.diag
end
