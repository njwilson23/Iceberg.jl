module Iceberg
using ib_types
import heat: tsolve!

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
    φ = x .- modelparams.dx[1]*n2

    return ModelState1d(T, φ, modelparams), physics
end

# initialize the scenario in Hill (1987) on pp 6-14
function initialize1d_hill(n=32)

    n2 = int(floor(n/2))
    modelparams = ModelParams(1e-5, (1.0/n,), 100, (n,))
    physics = PhysicalParams(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0)

    # set up a simple temperature field
    T = Array(Float64, n)
    T[1] = -1.0
    T[2:end] = 0.0

    # set up a global phi
    #x = linspace(0.5*modelparams.dx[1], modelparams.dx[1]*n-0.5*modelparams.dx[1], n)
    x = linspace(0.0, modelparams.dx[1]*(modelparams.nx[1]-1), modelparams.nx[1])
    φ = x .- modelparams.dx[1]*2

    return ModelState1d(T, -φ, modelparams), physics
end

# computes the normal velocity based on temperature gradient
# 1d only right now
function front_velocity(state::ModelState1d, phys::PhysicalParams)

    # need to look right at the front, and compute the gradient on either side.
    # This is trivial for 1d, but it would be a good idea to think about how I
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

    ∇sol = 1.0/state.params.dx[1] * dot(state.temp[idxsSolid], [-1, 1])
    ∇liq = 1.0/state.params.dx[1] * dot(state.temp[idxsLiquid], [-1, 1])

    vel = 1.0/phys.Lf * (phys.kaps * ∇sol - phys.kapl * ∇liq)
    return vel * ones(Float64, state.params.nx)
end

# return the position of the freezing front
function front_position(state::ModelState1d, phys::PhysicalParams)

    Tabs = abs(state.temp .- phys.tmelt)
    zidx = find(Tabs .== minimum(Tabs))[1]
    return zidx

end

# calculate ∇T
function gradT(state::ModelState)

    T = state.temp
    dx = diff(T, 1)
    dy = diff(T, 2)
    # not finished!

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

end #module
