# Scenario initializations useful for either testing to for building a
# customized model based on a pre-defined example

# test problem with a step function in 1d
function initialize1d_step(n=32)

    n2 = int(floor(n/2))
    modelparams = ModelParams(1e-5, (1.0/n,), (n,))
    physics = PhysicalParams(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0)

    # set up a simple temperature field
    T = Array(Float64, n)
    T[1:n2] = 1.0
    T[n2+1:end] = -1.0

    # set up a global phi
    x = linspace(0.5*modelparams.dx[1], modelparams.dx[1]*n-0.5*modelparams.dx[1], n)
    φ = x .- modelparams.dx[1]*n2

    return ModelState1d(T, φ, modelparams), physics
end

# test problem from Hill (1987), section 1.3
function initialize1d_hill(n=32)

    modelparams = ModelParams(1e-5, (1.0/n,), (n,))
    physics = PhysicalParams(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0)

    # set up a simple temperature field
    T = Array(Float64, n)
    T[1] = -1.0
    T[2:end] = 0.0

    # set up a global phi
    x = linspace(0.0, modelparams.dx[1]*(modelparams.nx[1]-1), modelparams.nx[1])
    φ = x .- modelparams.dx[1]*2

    return ModelState1d(T, -φ, modelparams), physics
end

# test problem in 1d with multiple fronts
function initialize1d_nfronts(n=64, nfronts=3)

    modelparams = ModelParams(1e-5, (2pi/n,), (n,))
    physics = PhysicalParams(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0)

    # set up a simple temperature field
    x = linspace(0.0, 2pi, modelparams.nx[1])
    T = cos(nfronts*x/2)

    # set up a global phi
    #φ = x .- modelparams.dx[1]*2
    φ = copy(T)
    reinitialize!(φ, 5)

    return ModelState1d(T, -φ, modelparams), physics
end

# test problem in 2d
function initialize2d(n=32)

end
