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

abstract ModelState

type ModelState1D   <: ModelState
    temp::Array
    phi::Array
    params::ModelParams
end

type ModelState2D   <: ModelState
    temp::Array
    phi::Array
    params::ModelParams
end

print(io::IO, A::ModelState) = @sprintf("Problem of size %s", size(A.phi))

# temporary function to initialize the problem domain
function initialize1d(n=32)

    n2 = int(floor(n/2))
    modelparams = ModelParams(1e-5, (1.0/n,), 100, (n,))
    physics = PhysicalParams(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0)

    # set up a simple temperature field
    T = Array(Float64, n)
    T[1:n2] = 1.0
    T[n2+1:end] = -1.0

    # set up a global phi
    x = linspace(0.5*modelparams.dx[1], modelparams.dx[1]*n-0.5*modelparams.dx[1], n)
    phi = x .- modelparams.dx[1]*n2

    return ModelState1D(T, phi, modelparams), physics

end

# solve the temperature evolution equation
function tsolve!(state::ModelState1D, phys::PhysicalParams)

    # this is a prototype functions that only works for one dimension
    #
    # simple explicit finite differences are used, so the timestep choice
    # requires care

    n = state.params.nx[1]
    L = spdiagm((ones(n-1), -2ones(n), ones(n-1)), (-1, 0, 1))
    # Apply Dirichlet boundary conditions
    # Warning: elementwise modification of sparse matrices is slow. This could
    # be done during matrix assembly, at loss of clarity
    L[1,2] = 0.0
    L[end,end-1] = 0.0
    L[1,1] = 0.0
    L[end,end] = 0.0

    mskSolid = state.phi .> 0.0

    alpha = Array(Float64, state.params.nx)
    alpha[mskSolid] = phys.kaps / (phys.cps * phys.rhos)
    alpha[!mskSolid] = phys.kapl / (phys.cpl * phys.rhol)

    geo = state.params.dt / state.params.dx[1]^2
    state.temp[:] += geo * alpha .* (L * state.temp[:])
    state.temp[~mskSolid] = max(state.temp[~mskSolid], phys.tmelt)
    state.temp[mskSolid] = min(state.temp[mskSolid], phys.tmelt)

    return state.temp
end

# computes the normal velocity based on temperature gradient
# 1d only right now
function front_velocity(state::ModelState1D, phys::PhysicalParams)

    # need to look right at the front, and compute the gradient on either side.
    # This is trivial for 1D, but it would be a good idea to think about how I
    # can generalize this to ND.
    zidx = front_position(state, phys)
    if (zidx < 3) | (zidx > state.params.nx[1] - 2)
        warn("Zero level set includes the domain boundary")
        return zeros(Float64, state.params.nx)
    elseif zidx - 1 > phys.tmelt
        idxsSolid = [zidx - 2, zidx - 1]
        idxsLiquid = [zidx + 1, zidx + 2]
    else
        idxsSolid = [zidx + 1, zidx + 2]
        idxsLiquid = [zidx - 2, zidx - 1]
    end

    gradSolid = 1.0/state.params.dx[1] * dot(state.temp[idxsSolid], [-1, 1])
    gradLiquid = 1.0/state.params.dx[1] * dot(state.temp[idxsLiquid], [-1, 1])

    vel = 1.0/phys.Lf * (phys.kaps * gradSolid - phys.kapl * gradLiquid)
    return vel * ones(Float64, state.params.nx)
end

# return the position of the freezing front
function front_position(state::ModelState1D, phys::PhysicalParams)

    Tabs = abs(state.temp .- phys.tmelt)
    zidx = find(Tabs .== minimum(Tabs))[1]
    return zidx

end

# reinitialize the level set function
# 1d only right now
#function lsupdate!(state::ModelState, phys::PhysicalParams)
#
#    Tabs = abs(state.temp .- phys.tmelt)
#    zidx = find(Tabs .== minimum(Tabs))[1]
#    
#    for idx=1:state.params.nx[1]
#        x = idx * state.params.dx[1]
#        state.phi[idx] = -sign(state.temp[idx]) * state.params.dx[1] * abs(idx - zidx)
#    end
#
#end

# master function that runs the model
function run()


end

end #module
