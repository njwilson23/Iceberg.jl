# Level set code to handle the computation, reinitialization, and evolution of
# level set functions (φ)

function exact_phi(state::ModelState1d)
    zeta = Iceberg.front_indices(state)
    if length(zeta) > 0
        dx = state.params.dx[1]
        nx = state.params.nx[1]
        dtable = Float64[abs(x-xf) for x=linspace(0.0, (nx-1)*dx, nx),
                                       xf in Iceberg.front_positions(state)]
        sdist = -sign(state.temp) .* minimum(dtable, 2)
    else
        sdist = nan * Array(Float64, state.params.nx[1])
    end
    return sdist[:]
end

function exact_phi(state::ModelState2d)
    error("Not implemented")
end

# reinitialize following the simple relaxation scheme of Rouy and Tourin
# A viscosity solutions approach to shape-from-shading (1992)
function reinitialize!(state::ModelState1d, iters::Int; k=0.01, kd=0.03)
    φ = state.phi
    sφ = sign(φ)
    dx = state.params.dx[1]
    n = state.params.nx[1]

    # Matrix for numerical diffusion
    D = spdiagm((ones(n-1), -2ones(n), ones(n-1)), (-1,0,1))
    D[1,2] = 2.0
    D[end,end-1] = 2.0

    # Neumann boundary conditions for signed distance function
    bc1 = kd*k*(φ[1]   < 0 ? 2/dx : -2/dx)
    bc2 = kd*k*(φ[end] < 0 ? 2/dx : -2/dx)
    for it=1:iters
        dφ = sφ .* (1.0 .- abs(gradient(φ, dx)))
        φ[:] .+= k * dφ
        # Numerical diffusion step
        φ[1] -= bc1
        φ[end] -= bc2
        φ[:] = (spdiagm(ones(n), 0) - D*kd*k/(dx^2)) \ φ
    end
end

# warning: boundary condition code assumes a square array
function reinitialize!(state::ModelState2d, iters::Int; k=0.01, kd=0.005)
    φ = state.phi
    sφ = sign(φ)
    dx = state.params.dx
    m, n = state.params.nx

    # Matrix for numerical diffusion
    D1 = spdiagm((ones(n-1), -2ones(n), ones(n-1)), (-1, 0, 1))
    D2 = spdiagm((ones(m-1), -2ones(m), ones(m-1)), (-1, 0, 1))
    for D_ in (D1, D2)
        D_[1,2] = 2.0
        D_[end,end-1] = 2.0
    end
    D = kron(D1, speye(m)) + kron(speye(n), D2)

    # Boundary indexes in lexigraphic order
    Ω = [1:m, m+1:m:m*(n-2)+1, 2m:m:m*(n-1), m*n-m+1:m*n]
    sort!(Ω)

    # Neumann boundary conditions for signed distance function
    bcs = Array(Float32, length(Ω))
    bcs[φ[Ω] .< 0.0] = 2.0/dx[1] * kd*k
    bcs[φ[Ω] .>= 0.0] = -2.0/dx[1] * kd*k

    for it=1:iters
        dφ = sφ .* (1.0 .- absgrad(φ, dx))
        φ[:,:] .+= k * dφ
        # Numerical diffusion step
        φ[Ω] -= bcs
        # warning: also assumes square array, i.e. dx[1] == dx[2]
        φ[:] = (spdiagm(ones(n^2), 0) - D*kd*k/(dx[1].^2)) \ φ[:]
    end
end

# compute |∇φ| for a field φ
function absgrad(A::Vector, h::Tuple)
    abs(gradient(A, h[1]))
end

# compute |∇φ| for a field φ
# needs to be extended to Nd
function absgrad(A::Array, h::Tuple)
    dx, dy = grad(A, h)
    return sqrt(dx.^2 + dy.^2)
end

# compute ∇φ for a field φ
function grad(A::Array, h::Tuple)
    m, n = size(A)
    dx = Array(Float64, (m,n))
    dy = Array(Float64, (m,n))
    for irow=1:m
        dx[irow,:] = gradient(reshape(A[irow,:], n), h[1])
    end
    for icol = 1:n
        dy[:,icol] = gradient(reshape(A[:,icol], m), h[2])
    end
    return dx, dy
end

# return the index positions of the freezing fronts
function front_indices(state::ModelState1d)

    φ = state.phi
    sφ = sign(φ)
    ζ = Int[]
    for i=2:length(state.temp)
        if sφ[i] != sφ[i-1]
            if sφ[i] == 0           # skip this case because the next
                pass                # condition will catch it next iteration
            elseif sφ[i-1] == 0
                push!(ζ, i-1)
            elseif abs(φ[i]) > abs(φ[i-1])
                push!(ζ, i-1)
            else
                push!(ζ, i)
            end
        end
    end

    return ζ
end

# return the interpolated position of the freezing fronts
function front_positions(state::ModelState1d)
    φ = state.phi
    ζ = front_indices(state)
    if length(ζ) == 0
        return Float64[]
    end
    X = Array(Float64, size(ζ))
    for (i, z) in enumerate(ζ)
        if 1 < z < state.params.nx[1]
            dφ = (φ[z+1] - φ[z-1]) / 2state.params.dx[1]
        elseif z > 1
            dφ = (φ[z] - φ[z-1]) / state.params.dx[1]
        else
            dφ = (φ[z+1] - φ[z]) / state.params.dx[1]
        end

        X[i] = (z-1)*state.params.dx[1] - (φ[z] / dφ)
    end
    return X
end

# Linear interpolation of a 1-d sequence
function interp1d(x::Vector, xp::Vector, yp::Vector)
    dx₁ = (yp[2] - yp[1]) / (xp[2] - xp[1])
    dx₂ = (yp[end] - yp[end-1]) / (xp[end] - xp[end-1])

    lb = 1
    y = Array(Float64, length(x))
    for i=1:length(x)
        if x[i] < xp[1]
            y[i] = yp[1] - (xp[1] - x[i]) * dx₁
        elseif x[i] >= xp[end]
            y[i] = (x[i] - xp[end]) * dx₂ + yp[end]
        else
            ip = lb
            for j = lb:length(xp)
                if xp[j] > x[i]
                    ip = j-1
                    lb = j
                    break
                end
            end

            y[i] = (x[i]-xp[ip]) / (xp[ip+1]-xp[ip]) * (yp[ip+1]-yp[ip]) + yp[ip]
        end
    end
    return y
end

# update the level set function using conservation law methods
function lsupdate!(state::ModelState, dphi::Array)
    error("Not implemented")
end
