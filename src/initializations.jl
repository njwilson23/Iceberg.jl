# Scenario initializations useful for either testing to for building a
# customized model based on a pre-defined example

# test problem with a step function in 1d
function initialize1d_step(n=32)

    n2 = int(floor(n/2))
    modelparams = ModelParams(1e-4, (1.0/n,), (n,))
    physics = PhysicalParams(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0)

    # set up a simple temperature field
    T = Array(Float64, n)
    T[1:n2] = 1.0
    T[n2+1:end] = -1.0

    # set up a global phi
    dx = modelparams.dx[1]
    x = linspace(0.5*dx, dx*n - 0.5*dx, n)
    φ = x .- modelparams.dx[1]*n2

    return ModelState1d(T, φ, modelparams), physics
end

# test problem from Hill (1987), section 1.3
function initialize1d_hill(n=32)

    modelparams = ModelParams(1e-4, (1.0/n,), (n,))
    physics = PhysicalParams(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0)

    # set up a simple temperature field
    T = Array(Float64, n)
    T[1] = -1.0
    T[2:end] = 0.0

    # set up a global phi
    dx = modelparams.dx[1]
    nx = modelparams.nx[1]
    x = linspace(0.0, dx*(nx-1), nx)
    φ = x .- modelparams.dx[1]*2

    return ModelState1d(T, -φ, modelparams), physics
end

# test problem in 1d with multiple fronts
function initialize1d_nfronts(n=64, nfronts=3)

    modelparams = ModelParams(1e-4, (2pi/n,), (n,))
    physics = PhysicalParams(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0)

    # set up a simple temperature field
    x = linspace(0.0, 2pi, modelparams.nx[1])
    T = cos(nfronts*x/2)

    # set up a global phi
    state = ModelState1d(T, -copy(T), modelparams)
    reinitialize!(state, 5)

    return state, physics
end

# 2d problem with a vertical freezing front akin to the Hill (1.3) scenario
function initialize2d_frontv(n=64)

    physics = PhysicalParams(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0)
    modelparams = ModelParams(2.5e-3, (1.0/(n-1), 1.0/(n-1)), (n, n))

    x = linspace(0.0, 1.0, modelparams.nx[1])
    y = linspace(0.0, 1.0, modelparams.nx[2])

    T = Array(Float64, (n,n))
    for i=1:n
        for j=1:n
            if j <= n/2
                T[i,j] = 0.0
            else
                T[i,j] = 1.0
            end
        end
    end

    phi = (0.5 .- T)
    state = ModelState2d(T, phi, modelparams)
    reinitialize!(state, 2n)
    return state, physics
end

# 2d problem with a diagonal freezing front akin to the Hill (1.3) scenario
function initialize2d_frontd(n=64)

    physics = PhysicalParams(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0)
    modelparams = ModelParams(2.5e-3, (1.0/(n-1), 1.0/(n-1)), (n, n))

    x = linspace(0.0, 1.0, modelparams.nx[1])
    y = linspace(0.0, 1.0, modelparams.nx[2])

    T = Array(Float64, (n,n))
    for i=1:n
        for j=1:n
            if i + j <= n
                T[i,j] = 0.0
            else
                T[i,j] = 1.0
            end
        end
    end

    phi = (0.5 .- T)
    state = ModelState2d(T, phi, modelparams)
    reinitialize!(state, 2n)
    return state, physics
end


# test problem in 2d
function initialize2d_square(n=32)

    physics = PhysicalParams(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0)
    modelparams = ModelParams(1e-3, (1.0/(n-1), 1.0/(n-1)), (n, n))

    x = linspace(0.0, 1.0, modelparams.nx[1])
    y = linspace(0.0, 1.0, modelparams.nx[2])

    n3 = int(n/3)
    T = ones(Float64, (n, n))
    T[n3:2n3,n3:2n3] = -1.0

    state = ModelState2d(T, -copy(T), modelparams)

    #state.phi = exact_phi(state)
    n3 = n/3
    u1 = linspace(-n3, 2n3, n)
    u2 = linspace(2n3, -n3, n)
    state.phi = state.params.dx[1] * min(u1 .* ones(n)', u2 .* ones(n)',
                                         u1' .* ones(n), u2' .* ones(n))
    return state, physics
end
