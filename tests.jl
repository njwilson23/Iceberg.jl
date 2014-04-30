module tests
using Base.Test
using Iceberg

function test_initialize1d()
    println("test_initialize1d")

    state, physics = Iceberg.initialize1d(64)
    @test length(state.params.dx) == 1
    @test length(state.params.nx) == 1

end

# ensure that the 1D solution is tolarably close to an analytical solution
function test_tsolve1d()


end

test_initialize1d()

end #module
