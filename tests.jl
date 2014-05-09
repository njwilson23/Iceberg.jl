module tests
using Base.Test
using Iceberg

function test_initialize1d()
    println("test_initialize1d")

    state, physics = Iceberg.initialize1d_step(64)
    @test length(state.params.dx) == 1
    @test length(state.params.nx) == 1

end

function test_interp1d()
    println("test_interp1d")

    xp = [0, 5, 7, 8, 12]
    yp = [0, 2, 6, 3, 7]
    xi = [2.5, 6, 7.5, 11]
    y_ans = [1, 4, 4.5, 6]
    yi = Iceberg.interp1d(xi, xp, yp)
    @test y_ans == yi
end

function test_front_indices()
    println("test_front_indices")

    # single front
    state = ModelState1d(linspace(-1, 1, 50),       # T
                         linspace(25, -25, 50),     # φ
                         ModelParams(1e-5, (0.04,), (50,)))
    ζ = Iceberg.front_indices(state)
    @test ζ == [26]

    # multiple fronts
    x = linspace(0, 4pi, 100)
    state = ModelState1d(sin(x),    # T
                         -sin(x),   # φ
                         ModelParams(1e-5, (0.04,), (50,)))
    ζ = Iceberg.front_indices(state)
    @test ζ == [1, 26, 50, 75]
end

function test_front_velocity()
    println("test_front_velocity")

    # Single front with velocity = (κs dTs - κl dTl) / L
    #                            = (1 * 1 - 0.5 * 1) / 1
    #                            = 0.5
    physics = PhysicalParams(1.0, 0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0)
    x = linspace(0, 1, 101)
    state = ModelState1d(x.-0.5,    # T
                         0.5.-x,    # φ
                         ModelParams(1e-5, (0.01,), (101,)))
    u = Iceberg.front_velocity(state, physics)
    @test_approx_eq u 0.5*ones(101)
end

# ensure that the 1D solution is tolarably close to an analytical solution
function test_hill_onephase()
    println("test_hill_onephase")

    x_numerical = Float64[]
    state, physics = Iceberg.initialize1d_hill(200)
    state.params.dt = 1e-4

    tend = 0.5
    nt = tend / state.params.dt

    for i=1:nt
        Iceberg.tsolve!(state, physics)
        vel = Iceberg.front_velocity(state, physics)
        state.phi += state.params.dt * vel
        
        zeta = Iceberg.front_indices(state)[1]
        push!(x_numerical, zeta[1] * state.params.dx[1])
    end

    t = state.params.dt * [1:nt]
    x_analytical = sqrt(2*0.769 * physics.kaps * t / (physics.rhos*physics.cps))

    residual = x_analytical - x_numerical
    @printf("\tL1 norm: %1.3f\n", sum(abs(residual)))
    @printf("\tL∞ norm: %1.3f (required to be less than 0.01)\n",
            maximum(abs(residual)))
    @test maximum(abs(residual)) < 0.01

end

function test_iceblock()
    println("test_iceblock")

    state, physics = Iceberg.initialize2d_square(64)
    state.params.dt = 5e-3
    physics.cps = 2.1
    physics.cpl = 4.2
    physics.Lf = 335.0

    tend = 0.5
    nt = tend / state.params.dt

    for i = 1:nt
        Iceberg.tsolve!(state, physics)
        vel = Iceberg.front_velocity(state, physics)
        state.phi += state.params.dt * vel
        if i%25 == 0
            Iceberg.reinitialize!(state, 2)   
        end
    end

    # Right now this just verifies that the test runs and provides non-NaN
    # results. Actual validation of results are TODO
    @test countnz(isnan(state.phi)) == 0
    @test countnz(isnan(state.temp)) == 0

end

function test_front_positions()
    println("test_front_positions")
    
    n = 256
    x = linspace(0, 2pi, n)
    dx = 2pi/n
    phi = -cos(2x)
    state = ModelState1d(0*phi, phi, ModelParams(0.01, (dx,),  (n,)))
    X = Iceberg.front_positions(state)
    Xtrue = [pi/4, 3pi/4, 5pi/4, 7pi/4]
    residual = abs(Xtrue - X)
    @test sum(residual) < 0.05

end

# solve the Hill front propagation problem in two dimensions

# Note: until the heat solver can handle Neumann boundary conditions on the
# upper and lower boundaries, this is expected to propagate slightly faster
# than the analytical solution.
function test_hill_onephase_2d(n=32)
    println("test_hill_onephase_2d")

    function find_interface_index(φ::Vector)
        index = find(abs(φ) .== minimum(abs(φ)))
        if length(index) > 0
            return index[1]
        else
            return NaN
        end
    end

    state, physics = Iceberg.initialize2d_frontv(n)
    state.params.dt = 1e-4
    state.params.dx = (1.0/n, 0.5/n)


    tend = 0.5
    nt = int(floor(tend / state.params.dt))
    x_numerical = zeros(Float64, nt)

    for i=1:nt
        Iceberg.tsolve!(state, physics)
        vel = Iceberg.front_velocity(state, physics)
        state.phi += state.params.dt * vel
        
        ζ = find_interface_index(state.phi[int(floor(n/2)),:][:]) - 2
        x_numerical[i] = state.params.dx[1] * ζ

        if i%1 == 0
            Iceberg.reinitialize!(state, 1)
        end
    end

    t = state.params.dt * [1:nt]
    x_analytical = sqrt(2*0.769 * physics.kaps * t / (physics.rhos*physics.cps))
    residual = x_numerical - x_analytical

    @printf("\tL1 norm: %1.3f\n", sum(abs(residual)))
    @printf("\tL∞ norm: %1.3f\n", maximum(abs(residual)))
    @test maximum(abs(residual)) < 0.05

end

test_initialize1d()
test_interp1d()
test_front_indices()
test_front_positions()
test_front_velocity()

test_hill_onephase()
test_hill_onephase_2d()
test_iceblock()

end #module
