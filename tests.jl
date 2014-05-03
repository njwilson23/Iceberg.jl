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
    xp = [0, 5, 7, 8, 12]
    yp = [0, 2, 6, 3, 7]
    xi = [2.5, 6, 7.5, 11]
    y_ans = [1, 4, 4.5, 6]
    yi = Iceberg.interp1d(xi, xp, yp)
    @test y_ans == yi
end


# ensure that the 1D solution is tolarably close to an analytical solution
function test_tsolve1d()


end

test_initialize1d()
test_interp1d()

end #module
