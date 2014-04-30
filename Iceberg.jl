module Iceberg
import Base.show
#using Temperature
#using Levelset

# maybe this should be a dictionary?
#PhysicalParams = {:kapl => 1.0,
#                  :kaps => 1.0,
#                  :cpl => 1.0,
#                  :cps => 1.0,
#                  :rhol => 1.0,
#                  :rhos => 1.0,
#                  :tmelt => 0.0}

type PhysicalParams
    Lf::Float64
    kapl::Float64
    kaps::Float64
    cpl::Float64
    cps::Float64
    rhol::Float64
    rhos::Float64
    tmelt::Float64
end

type ModelParams
    dt::Float64
    dx::Tuple
    nt::Integer
    nx::Tuple
end

type ModelState
    temp::Array
    phi::Array
    #Tmat::AbstractArray
end

show(A::ModelState) = @sprintf("Problem of size %s", size(A.phi))

# temporary function to initialize the problem domain
function initialize1d(n=32)

    n2 = int(floor(n/2))
    mpar = ModelParams(1e-5, (1.0/n,), 100, (n,))
    ppar = PhysicalParams(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0)

    # set up a simple temperature field
    T = Array(Float64, n)
    T[1:n2] = 1.0
    T[n2+1:end] = -1.0

    # set up a global phi
    x = linspace(0.5*mpar.dx[1], mpar.dx[1]*n-0.5*mpar.dx[1], n)
    phi = x .- mpar.dx[1]*n2

    return ModelState(T, phi), mpar, ppar

end

# solve the temperature evolution equation
function tsolve!(state::ModelState, mpar::ModelParams, ppar::PhysicalParams)

    # this is a prototype functions that only works for one dimension with a
    # t=1.0 Dirichlet BC on the left
    #
    # simple explicit finite differences are used, so the timestep choice
    # requires care

    n = mpar.nx[1]
    L = spdiagm((ones(n-1), -2ones(n), ones(n-1)), (-1, 0, 1))
    # Apply Dirichlet boundary conditions
    # Warning: elementwise modification of sparse matrices is slow. This could
    # be done during matrix assembly, at loss of clarity
    L[1,2] = 0.0
    L[end,end-1] = 0.0
    L[1,1] = 0.0
    L[end,end] = 0.0

    # Compute a field alpha
    alphas = ppar.kaps / (ppar.cps * ppar.rhos)
    alphal = ppar.kapl / (ppar.cpl * ppar.rhol)
    alpha = [state.phi[i] < 0.0 ? alphal : alphas for i=1:mpar.nx[1]]

    dT = mpar.dt / mpar.dx[1]^2 * alpha .* (L * state.temp[:])
    state.temp += reshape(dT, mpar.nx)
    return state.temp
end

# Pin the temperature of the fluid phase
function mixed_fluid_phase!(state::ModelState, mpar::ModelParams, ppar::PhysicalParams)
    state.temp[state.phi .< 0] = 1.0
end

# computes the normal velocity based on temperature gradient
# 1d only right now
function front_velocity(state::ModelState, mpar::ModelParams, ppar::PhysicalParams)

    # need to look right at the front, and compute the gradient on either side.
    # This is trivial for 1D, but it would be a good idea to think about how I
    # can generalize this to ND.
    Tabs = abs(state.temp .- ppar.tmelt)
    zidx = find(Tabs .== minimum(Tabs))[1]

    if (zidx < 3) | (zidx > mpar.nx[1] - 2)
        warn("Zero level set includes the domain boundary")
        return zeros(Float64, mpar.nx)
    elseif zidx - 1 > ppar.tmelt
        idxsSolid = [zidx - 2, zidx - 1]
        idxsLiquid = [zidx + 1, zidx + 2]
    else
        idxsSolid = [zidx + 1, zidx + 2]
        idxsLiquid = [zidx - 2, zidx - 1]
    end

    gradSolid = 1.0/mpar.dx[1] * dot(state.temp[idxsSolid], [-1, 1])
    gradLiquid = 1.0/mpar.dx[1] * dot(state.temp[idxsLiquid], [-1, 1])

    vel = 1.0/ppar.Lf * (ppar.kaps * gradSolid - ppar.kapl * gradLiquid)
    return vel * ones(Float64, mpar.nx)
end

# reinitialize the level set function
# 1d only right now
#function lsupdate!(state::ModelState, mpar::ModelParams, ppar::PhysicalParams)
#
#    Tabs = abs(state.temp .- ppar.tmelt)
#    zidx = find(Tabs .== minimum(Tabs))[1]
#    
#    for idx=1:mpar.nx[1]
#        x = idx * mpar.dx[1]
#        state.phi[idx] = -sign(state.temp[idx]) * mpar.dx[1] * abs(idx - zidx)
#    end
#
#end

# master function that runs the model
function run()


end

end #module
