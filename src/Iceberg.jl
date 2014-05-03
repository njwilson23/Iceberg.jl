module Iceberg
using ib_types
import heat: tsolve!

export ModelParams, ModelState1d, ModelParams2d,
       PhysicalParams

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

compute_lsfunc(state::ModelState) = compute_lsfunc(state.temp, state.params.dx)
function compute_lsfunc(temp::Vector, dx::Tuple)

    # find zero crossings

end

# reinitialize following the simple relaxation scheme of Rouy and Tourin
# A viscosity solutions approach to shape-from-shading (1992)
function reinitialize!(phi::Vector, iters::Int)
    sphi = sign(phi)
    for it=1:iters
        dphi = sphi .* (1.0 .- abs(gradient(phi)))
        phi[:] .+= dphi
        phi[:] = conv([0.25, 0.5, 0.25], phi)[2:end-1]
    end
end

function reinitialize!(phi::Array, iters::Int)
    sphi = sign(phi)
    convwin = [0.125 0.25 0.125;
                0.25  0.5  0.25;
               0.125 0.25 0.125] / 2
    for it=1:iters
        dphi = sphi .* (1.0 .- absgrad(phi, (1.0, 1.0)))
        phi[:,:] = conv2(convwin, phi+dphi)[2:end-1, 2:end-1]
    end
end

function absgrad(A::Vector, h::Tuple)
    abs(gradient(A, h[1]))
end

# needs to be extended to Nd
function absgrad(A::Array, h::Tuple)

    m, n = size(A)
    dx = Array(Float64, (m,n))
    dy = Array(Float64, (m,n))

    for irow=1:m
        dx[irow,:] = gradient(reshape(A[irow,:], n), h[1])
    end
    for icol = 1:n
        dy[:,icol] = gradient(reshape(A[:,icol], m), h[2])
    end

    return sqrt(dx.^2 + dy.^2)
end


# computes the normal velocity based on temperature gradient
# 1d only right now
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
    else
        return velocities[1] * ones(Float64, state.params.nx[1])
    end
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
function front_positions(state::ModelState1d, phys::PhysicalParams)
    ζ = front_indices(state)
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

# update the level set function by calculating interface velocities and solving
# the hyperbolic evolution equation
function lsupdate!(state::ModelState, phys::PhysicalParams)

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
