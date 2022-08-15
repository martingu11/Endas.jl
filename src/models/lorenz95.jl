using LinearAlgebra

export Lorenz95Model


mutable struct Lorenz95LinData
    δt::Float64
    x::Vector{Float64}
    k₁::Vector{Float64}
    k₂::Vector{Float64}
    k₃::Vector{Float64}
    k₄::Vector{Float64}
end
    


mutable struct Lorenz95ModelData

end


@doc raw"""

"""
struct Lorenz95Model <: LinearizedModel
    n::Integer
    F::Real

end

Lorenz95Model(n::Integer) = Lorenz95Model(n, 8)
Lorenz95Model() = Lorenz95Model(40, 8)


# Basic properties
statesize(ℳ::Lorenz95Model) = ℳ.n


l95_index(i::Integer, n::Integer) = i < 1 ? n+i : (i > n ? i-n : i) 

function l95!(out::AbstractVector{<:Real}, x::AbstractVector{<:Real}, δt::Real, F::Real)
    n = length(x)
    for i in 1:n
        i₋₂ = l95_index(i-2, n)
        i₋₁ = l95_index(i-1, n)
        i₊₁ = l95_index(i+1, n)
        out[i] = (-x[i₋₂]*x[i₋₁] + x[i₋₁]*x[i₊₁] - x[i] + F) * δt
    end
end


function l95tl!(out::AbstractVector{<:Real}, x::AbstractVector{<:Real}, 
                δx::AbstractVector{<:Real}, δt::Real)
    n = length(x)
    for i in 1:n
        i₋₂ = l95_index(i-2, n)
        i₋₁ = l95_index(i-1, n)
        i₊₁ = l95_index(i+1, n)
        out[i] = (-x[i₋₁]*δx[i₋₂] + (x[i₊₁]-x[i₋₂])*δx[i₋₁] - δx[i] + x[i₋₁]*δx[i₊₁]) * δt
    end
end

function l95ad!(out::AbstractVector{<:Real}, x::AbstractVector{<:Real}, 
                δx::AbstractVector{<:Real}, δt::Real=1.)
    n = length(x)
    for i in 1:n
        i₋₂ = l95_index(i-2, n)
        i₋₁ = l95_index(i-1, n)
        i₊₁ = l95_index(i+1, n)
        i₊₂ = l95_index(i+2, n)
        out[i] = (x[i₋₂]*δx[i₋₁] + (x[i₊₂]-x[i₋₁])*δx[i₊₁] - δx[i] - x[i₊₁]*δx[i₊₂]) * δt
    end
end



function apply!(x::AbstractVector{<:Real}, M::Lorenz95Model, t::Real, δt::Real, linearize::Bool)
    
    @assert length(x) == M.n

    k₁ = Vector{Float64}(undef, M.n)
    k₂ = Vector{Float64}(undef, M.n)
    k₃ = Vector{Float64}(undef, M.n)
    k₄ = Vector{Float64}(undef, M.n)

    l95!(k₁, x, δt, M.F);
    l95!(k₂, x + k₁/2, δt, M.F);
    l95!(k₃, x + k₂/2, δt, M.F);
    l95!(k₄, x + k₃, δt, M.F);

    @. x += (k₁ + 2k₂ + 2k₃ + k₄) / 6

    return linearize ? Lorenz95LinData(δt, copy(x), k₁, k₂, k₃, k₄) : nothing
end


# Tangent-linear operator
function apply_tl!(X::AbstractMatrix{<:Real}, ℳ::Lorenz95Model, Mlin::Lorenz95LinData)
    n, m = size(X)
    @assert statesize(ℳ) == n

    δk₁ = similar(X, n)
    δk₂ = similar(X, n)
    δk₃ = similar(X, n)
    δk₄ = similar(X, n)

    for i in 1:m
        x = view(X, :, i)

        l95tl!(δk₁, Mlin.x, x, Mlin.δt)
        l95tl!(δk₂, Mlin.x + Mlin.k₁/2, x + δk₁/2, Mlin.δt)
        l95tl!(δk₃, Mlin.x + Mlin.k₂/2, x + δk₂/2, Mlin.δt)
        l95tl!(δk₄, Mlin.x + Mlin.k₃, x + δk₃, Mlin.δt)

        @. x += (δk₁ + 2δk₂ + 2δk₃ + δk₄) / 6
    end
end


# Adjoint operator
function apply_adj!(X::AbstractMatrix{<:Real}, ℳ::Lorenz95Model, Mlin::Lorenz95LinData)
    n, m = size(X)
    @assert statesize(ℳ) == n

    function k₁ad!(out::AbstractVector{<:Real}, δx::AbstractVector{<:Real}) 
        l95ad!(out, Mlin.x, δx, Mlin.δt)
    end

    function k₂ad!(out::AbstractVector{<:Real}, δx::AbstractVector{<:Real}) 
        tmp = similar(out)
        l95ad!(out, Mlin.x + Mlin.k₁/2, δx, 1)
        k₁ad!(tmp, out)
        @. out = (out + tmp / 2) * Mlin.δt
    end

    function k₃ad!(out::AbstractVector{<:Real}, δx::AbstractVector{<:Real}) 
        tmp = similar(out)
        l95ad!(out, Mlin.x + Mlin.k₂/2, δx, 1)
        k₂ad!(tmp, out)
        @. out = (out + tmp / 2) * Mlin.δt
    end

    function k₄ad!(out::AbstractVector{<:Real}, δx::AbstractVector{<:Real}) 
        tmp = similar(out)
        l95ad!(out, Mlin.x + Mlin.k₃/2, δx, 1)
        k₃ad!(tmp, out)
        @. out = (out + tmp / 2) * Mlin.δt
    end

    outk₁ = similar(X, n)
    outk₂ = similar(X, n)
    outk₃ = similar(X, n)
    outk₄ = similar(X, n)

    for i in 1:m
        x = view(X, :, i)

        k₁ad!(outk₁, x)
        k₂ad!(outk₂, x)
        k₃ad!(outk₃, x)
        k₄ad!(outk₄, x)

        @. x += (outk₁ + 2outk₂ + 2outk₃ + outk₄) / 6
    end
end
    
