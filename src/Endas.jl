"""
Endas package.
"""
module Endas


#abstract type MatrixType end
#struct AnyMatrix <: MatrixType end
#struct DenseMatrix <: MatrixType end
#asmatrix(x) = asmatrix(AnyMatrix(), x)

include("functions.jl")
include("internal.jl")

include("model.jl")
include("covariance.jl")
include("observation.jl")
include("filter.jl")


include("utils.jl")





end
