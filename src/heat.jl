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

