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
    @printf("L∞ norm: %1.3f\n", sum(abs(residual)))
    @printf("L1 norm: %1.3f (required to be less than 0.02)\n",
            maximum(abs(residual)))
    @test maximum(abs(residual)) < 0.02

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

test_initialize1d()
test_interp1d()
test_front_indices()
test_front_positions()
test_front_velocity()
test_hill_onephase()

end #module
