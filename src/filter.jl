
export AbstractFilter, AbstractSmoother


"""
"""
abstract type AbstractFilter end

"""
"""
abstract type AbstractSmoother <: AbstractFilter end



# Filter and smoother algorithm implementations

include("filters/kf.jl")
