using Endas

using Plots
using LinearAlgebra

const σclim = 3.6414723


ℳ = Endas.Lorenz95Model(40, 8.)
n = statesize(ℳ)


# Number of data assimilation steps. The model has a period of 12 time steps.
nsteps = 1000

# Model integration time step.
# Note: Too long integration step may cause difficulties with convergence, especially
# for KalmanSmoother that relies on the linearization of the model. If this happens, 
# decrease the integration step.
δt = 0.025 / 6

# Interval at which observations are assimilated. Set to 1 to assimilate at every k
obsskip = 5


# Observation operator
# We will observe the last 3 variables in every 5 state vector variables, i.e. 24 out of 
# the 40 state variables in this case
k = 24
H = zeros(k, n)
for i in 1:24
    H[i, 5*((i-1) ÷ 3) + ((i-1) % 3) + 1] = 1
end
ℋ = LinearObservationOp(H)

# Covariance matrix of the process noise
Q = DiagonalCovariance((0.05σclim)^2, n)

# Observation error covariance
R = DiagonalCovariance((0.15σclim)^2, k)


x₀ = fill(8., n)
x₀[20] = 8.004        # Perturb the 20-th coordinate

P₀ = diagm(fill((0.5*σclim)^2, n))


print("Generating test data...\n")
data = Endas.generate_data(ℳ, x₀, 0, δt, nsteps, ℋ, Q=Q, R=R, obsskip=obsskip, nspin=300)


kf = KalmanFilter(ℳ, copy(x₀), copy(P₀), 0)
xres = zeros(n, nsteps)
print("Running Kalman filter...\n")


for i in 2:nsteps
    
    forecast!(kf, δt, Q)
    xres[:,i] .= getstate(kf)

    obs = data.y[i]
    if !isnothing(obs)
        assimilate!(kf, obs, ℋ, R)
    end
end

p = plot(data.xt[21,:], lw=1, seriescolor=:black, label="xₜ")

obs = [y[13] for y in data.y if !isnothing(y)]
obssteps = [i for (i, y) in enumerate(data.y) if !isnothing(y)]

plot!(p, obssteps, obs, seriestype=:scatter)
plot!(p, xres[21,:], seriescolor=:green, label="KF")




