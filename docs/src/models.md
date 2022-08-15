# Implementing own process models


## Implementing the process function


At minimum, all generic models must implement

    statesize(ℳ::<modeltype>)

where `<modeltype>` is the concrete model type.


    apply!(X::AbstractVector{<:Real}, ℳ::<modeltype>, t::Real, δt::Real)


With the above definition, ensembles are handled automatically by Endas by calling the 
`apply!` function on each ensemble member individually.


## Implementing linearization


If the linearized model can be represented by a (sparse) matrix, it is sufficient to return 
the matrix if place of the linearization data. Endas implements `apply_tl!` and `apply_adj` 
for this case which simply apply the linear transformation to a state vector or ensemble.
While this may be sufficient for simple toy models, in real use one typically want to 
implement custom tangent-linear and adjoint operators. The methods to implement are

    apply_tl!(x::AbstractMatrix{<:Real}, ℳ::<modeltype>, Mlin::<lindata>) 

for the tangent-linear operator and

    apply_adj!(x::AbstractMatrix{<:Real}, ℳ::<modeltype>, Mlin::<lindata>) 

for the adjoint operator. Here `<lindata>` is the type returned from `apply` that holds 
data needed to compute the model linearization.


## Parallel execution

TODO!
