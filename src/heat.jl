# Functions to evolve the temperature field for the Stefan problem by solving
# the heat equation.

# The function `tsolve!()` uses a stable Crank-Nicolsen scheme that is second
# order in time and space.

# The function `tsolve_explicit()` is less stable and less accurate, and is
# retained for large problems on coarse grids and comparison purposes.

function assemble_mat(state::ModelState1d)
    n = state.params.nx[1]
    L = spdiagm((ones(n-1), -2ones(n), ones(n-1)), (-1, 0, 1))
    # Apply Dirichlet boundary conditions
    L[1,2] = 0.0
    L[end,end-1] = 0.0
    L[1,1] = 0.0
    L[end,end] = 0.0
    return L
end

function assemble_mat(state::ModelState2d)
    m, n = state.params.nx
    L1 = spdiagm((ones(n-1), -2ones(n), ones(n-1)), (-1, 0, 1))
    L2 = spdiagm((ones(m-1), -2ones(m), ones(m-1)), (-1, 0, 1))
    # Apply Dirichlet boundary conditions
    for L_ in (L1, L2)
        L_[1,2] = 0.0
        L_[end,end-1] = 0.0
        L_[1,1] = 0.0
        L_[end,end] = 0.0
    end
    L = kron(L1, speye(m)) + kron(speye(n), L2)
    return L
end

# solve the temperature evolution equation using a 2nd order in time and space
# Crank-Nicolsen scheme
function tsolve!(state::ModelState, phys::PhysicalParams)
    n = prod(state.params.nx)
    L = assemble_mat(state)

    mskSolid = state.phi .> 0.0

    alpha = Array(Float64, n)
    alpha[mskSolid[:]] = phys.kaps / (phys.cps * phys.rhos)
    alpha[!mskSolid[:]] = phys.kapl / (phys.cpl * phys.rhol)

    # Solve (1 - 2/(2h^2)L) * Tnext = (1 + 2/(2h^2)) * T
    k2h2L = state.params.dt / (2*state.params.dx[1]^2) * L * spdiagm(alpha)
    state.temp[:] = (spdiagm(ones(n)) - k2h2L) \ 
                    ((spdiagm(ones(n)) + k2h2L) * state.temp[:])

    state.temp[~mskSolid] = max(state.temp[~mskSolid], phys.tmelt)
    state.temp[mskSolid] = min(state.temp[mskSolid], phys.tmelt)

    return state.temp
end

# solve the temperature evolution equation using a 1st order in time and 2nd
# order in space Euler scheme
function tsolve_explicit!(state::ModelState, phys::PhysicalParams)
    L = assemble_mat(state)

    mskSolid = state.phi .> 0.0

    alpha = Array(Float64, state.params.nx)
    alpha[mskSolid] = phys.kaps / (phys.cps * phys.rhos)
    alpha[!mskSolid] = phys.kapl / (phys.cpl * phys.rhol)

    kh2 = state.params.dt / state.params.dx[1]^2
    state.temp[:] += kh2 * alpha .* (L * state.temp[:])
    state.temp[~mskSolid] = max(state.temp[~mskSolid], phys.tmelt)
    state.temp[mskSolid] = min(state.temp[mskSolid], phys.tmelt)

    return state.temp
end

# computes the normal velocity based on temperature gradient
function front_velocity(state::ModelState1d, phys::PhysicalParams)

    # need to look right at the front, and compute the gradient on either side.
    # This is trivial for 1d, but it would be a good idea to think about how I
    # can generalize this to ND.
    ζ = front_indices(state)
    velocities = Array(Float64, length(ζ))
    dx = state.params.dx[1]

    for (i, z) in enumerate(ζ)

        if (z > 2) & (z < state.params.nx[1] - 1)

            if state.phi[z-1] > 0.0
                idxsSolid = [z - 2, z - 1]
                idxsLiquid = [z + 1, z + 2]
            else
                idxsSolid = [z + 1, z + 2]
                idxsLiquid = [z - 2, z - 1]
            end

        elseif (z <= 2)

            if state.phi[z] > 0.0
                idxsSolid = [1, 1]
                idxsLiquid = [z + 1, z + 2]
            else
                idxsSolid = [z + 1, z + 2]
                idxsLiquid = [1, 1]
            end

        elseif (z >= state.params.nx[1] - 1)

            n = state.params.nx[1]

            if state.phi[z] > 0.0
                idxsSolid = [n, n]
                idxsLiquid = [z - 1, z - 2]
            else
                idxsSolid = [z - 1, z - 2]
                idxsLiquid = [n, n]
            end

        end

        ∇sol = 1.0/dx * dot(state.temp[idxsSolid], [-1, 1])
        ∇liq = 1.0/dx * dot(state.temp[idxsLiquid], [-1, 1])
        velocities[i] = (phys.kaps * ∇sol - phys.kapl * ∇liq) / phys.Lf

    end

    if length(ζ) > 1
        return interp1d([1:state.params.nx[1]], ζ, velocities)
    elseif length(ζ) > 0
        return velocities[1] * ones(Float64, state.params.nx[1])
    else
        return zeros(Float64, state.params.nx[1])
    end
end

# 2d front velocities within a band around the zero level set
#
# requires the level set function φ to be reasonably close to the signed
# distance function in order to accurately determine the bandwidth
#
# Works by computing the band gradient in φ and the global gradient in
# temperature. For each point in the band, computes a potential front velocity
# based on the local temperature (potentially liquid) and the up-φ-gradient
# temperature (potentially solid)
function front_velocity(state::ModelState2d, phys::PhysicalParams)

    # make band global for now
    idxBand = 1:prod(state.params.nx)
    #idxBand = find(-2.0 < φ < 2.0)

    dTdx, dTdy = grad(state.temp, state.params.dx)
    dφdx, dφdy = grad(state.phi, state.params.dx)
    dTdn = dTdx .* dφdx + dTdy .* dφdy

    κl = phys.kapl
    κs = phys.kaps
    stefanvel(Tl_n::Float64, Ts_n::Float64) = (κl*Tl_n + κs*Ts_n) / phys.Lf

    V = zeros(Float64, state.params.nx)

    for i=2:state.params.nx[1]-1
        for j=2:state.params.nx[2]-1
            upstrDir = angle([dφdx[i,j], dφdy[i,j]])
            if upstrDir < 0.25pi || upstrDir >= 1.75pi
                V[i,j] = stefanvel(dTdn[i,j], dTdn[i,j+1])
            elseif 0.25pi <= upstrDir < 0.75pi
                V[i,j] = stefanvel(dTdn[i,j], dTdn[i-1,j])
            elseif 0.75pi <= upstrDir < 1.25pi
                V[i,j] = stefanvel(dTdn[i,j], dTdn[i,j-1])
            elseif 1.25pi <= upstrDir < 1.75pi
                V[i,j] = stefanvel(dTdn[i,j], dTdn[i+1,j])
            end
        end
    end
    return V
end

unitvec(vec::Vector) = vec / norm(vec)
function angle(vec::Vector)
    if vec[1] > 1e-8
        if vec[2] >= 0
            return atan(vec[2] / vec[1])
        elseif vec[2] < 0
            return atan(vec[2] / vec[1]) + 2pi
        end
    elseif vec[1] < -1e-8
        return atan(vec[2] / vec[1]) + pi
    elseif vec[2] > 0
        return 0.5pi
    else
        return pi
    end
end

# 2d front velocities using Chen's (1997) approximation where 
#   cps = cpl = kaps = kapl = 1
# then advects derivatives in temperature in the normal direction
function front_velocity_chen(state::ModelState2d, phys::PhysicalParams)

    #dtdnLiquid = nan*Array(Float64, state.params.nx)
    #dtdnLiquid = nan*Array(Float64, state.params.nx)
    #isSolid = state.phi .>= 0.0

    # Normal directions
    φ = state.phi
    φx, φy = grad(state.phi, state.params.dx)
    # Temperature gradient
    dtdx, dtdy = grad(state.temp, state.params.dx)

    # Extend the derivatives in the normal direction using upwinding
    # dT/dx
    d2tdx2 = zeros(Float64, size(dtdx))
    k = 0.9 * state.params.dx[2]
    for i = 1:state.params.nx[2]
        d2tdx2[:,2:end-1] = 0.5 * (dtdx[:,3:end] - dtdx[:,1:end-2])
        dtdx -= k * sign(φ .* φx) * d2tdx2
    end

    d2tdy2 = zeros(Float64, size(dtdy))
    k = 0.9 * state.params.dx[2]
    for i = 1:state.params.nx[2]
        d2tdy2[2:end-1,:] = 0.5 * (dtdy[3:end,:] - dtdy[1:end-2,:])
        dtdy -= k * sign(φ .* φx) * d2tdy2
    end

    # Temperature gradient normal to level sets
    dtdn = dtdx*φx + dtdy*φy
    return dtdx

    #dtdnSolid[isSolid] = dtdn[isSolid]
    #dtdnLiquid[~isSolid] = dtdn[~isSolid]
    return -dtdn / phys.Lf

    # Alternative algorithm:
    #
    # look a every undefined point in the solid domain.
    # if there is a neighbouring undefined point, it borders an interface
    # associate the solid border cells with a downstream (in φ) liquid value
    # and perform the complimentary operation
    # compute the velocity of the front where there are paired derivatives

end

