# 
# Common EnDAS functions.
#


"""
Return column vector or matrix.

If `x` is a 1-dimensional array, it is reinterpreted as a column vector of size (n, 1).
Otherwise `x` is returned as it is.
"""
function colvec_or_mat(x::AbstractVecOrMat)
   return ndims(x) == 1 ? reshape(x, :, 1) : x
end