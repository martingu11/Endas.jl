

export generate_data


"""
    generate_data(ℳ::AbstractModel, 
                  x::AbstractVector{<:Real}, 
                  t::Integer, δt::Number, nsteps::Integer, 
                  H::ObservationOperator,
                  Q::Union{Nothing, CovarianceOperator} = nothing, 
                  R::Union{Nothing, CovarianceOperator} = nothing, 
                  obsskip::Integer = 1)

Generate data for a twin experiment using the given model.

# Arguments:

- `ℳ::AbstractModel` - Process model
- `x::AbstractVector{<:Real}` - Initial state vector
- `t::Integer` - Initial time
- `δt::Number` - Time differenmce between steps
- `nsteps::Integer`: Number of time steps to generate
- `H::ObservationOperator` - Observation operator for generating observations from the 
                             state
- `Q::Union{Nothing, CovarianceOperator}` - Process noise (optional)
- `R::Union{Nothing, CovarianceOperator}` - Observation noise (optional)
- `obsskip::Integer` - If greater than 1, observations are only generated for every 
                       `obsskip`ᵗʰ step

Return named tuple `(xt, times, y)` where `xt` is an n×nsteps matrix of the true states 
(stored in columns), times is an array of the physical times corresponding to the time
steps and y is a vector of observation vectors (`Vector{Union{Nothing, Vector{Float64}}}`)
containing containing observations for steps for which observations are generated, or
`nothing`.
"""
function generate_data(ℳ::AbstractModel, 
                       x::AbstractVector{<:Real}, 
                       t::Integer, δt::Number, nsteps::Integer, 
                       H::ObservationOperator;
                       Q::Union{Nothing, CovarianceOperator}=nothing, 
                       R::Union{Nothing, CovarianceOperator}=nothing, 
                       obsskip::Integer=1,
                       nspin::Integer=0)
 
    k = obssize(H)
    n = statesize(H)
    @assert statesize(ℳ) == n

    times = zeros(nsteps)
    times[1] = t

    obs = Vector{Union{Nothing, Vector{Float64}}}(nothing, nsteps)

    tmpₓ = Vector{Float64}(undef, n)
    tmpₒ = Vector{Float64}(undef, k)

    x = copy(x)

    # Empty model run (spin-up) if requested
    for i in 0:nspin
        apply!(x, ℳ, t, δt, false)
        t += δt
    end

    xₜ = zeros(n, nsteps)
    xₜ[:, 1] .= x
    
    # Sample the true state by running the model and generate appropriate observations
    for i in 2:nsteps
        apply!(x, ℳ, t, δt, false)
        t += δt

        # Add process noise to the true state
        if !isnothing(Q) 
            random_multivariate_normal!(tmpₓ, Q)
            x .+= tmpₓ
        end

        xₜ[:, i] .= x
        times[i] = t

        # Have observation at this step
        if mod(i, obsskip) == 0

            # Generate the observation from true state
            y = Vector{Float64}(undef, k)
            apply!(y, H, x, t)

            # Add observation noise to the observations
            if !isnothing(R) 
                random_multivariate_normal!(tmpₒ, R)
                y .+= tmpₒ
            end

            obs[i] = y
        end
    end

    return (xt = xₜ, times = times, y = obs)
end