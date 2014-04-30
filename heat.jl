module heat

using ib_types

function assemble_mat(state::ModelState1d)
    n = state.params.nx[1]
    L = spdiagm((ones(n-1), -2ones(n), ones(n-1)), (-1, 0, 1))
    # Apply Dirichlet boundary conditions
    # Warning: elementwise modification of sparse matrices is slow. This could
    # be done during matrix assembly, at loss of clarity
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
    for L_ in (L1, L2)
        L_[1,2] = 0.0
        L_[end,end-1] = 0.0
        L_[1,1] = 0.0
        L_[end,end] = 0.0
    end
    L = kron(L1, eye(m)) + kron(eye(n), L2)
    return L
end

# solve the temperature evolution equation
function tsolve!(state::ModelState, phys::PhysicalParams)
    # this is a prototype function that only works for one dimension
    #
    # simple explicit finite differences are used, so the timestep choice
    # requires care
    L = assemble_mat(state)

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

end
