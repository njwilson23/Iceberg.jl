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
    # single front
    println("test_front_indices")
    state = ModelState1d(linspace(-1, 1, 50),       # T
                         linspace(25, -25, 50),     # φ
                         ModelParams(1e-5, (0.04,), (50,)))
    ζ = Iceberg.front_indices(state)
    @test ζ == [26]

    # multiple fronts
    x = linspace(0, 4pi, 100)
    state = ModelState1d(sin(x),       # T
                         -sin(x),     # φ
                         ModelParams(1e-5, (0.04,), (50,)))
    ζ = Iceberg.front_indices(state)
    @test ζ == [1, 26, 50, 75]
end


# ensure that the 1D solution is tolarably close to an analytical solution
function test_tsolve1d()


end

test_initialize1d()
test_interp1d()
test_front_indices()

end #module
