module tests
using Base.Test
using Iceberg

function test_initialize1d()
    println("test_initialize1d")

    state, physics = Iceberg.initialize1d(64)
    @test length(state.params.dx) == 1
    @test length(state.params.nx) == 1

end

test_initialize1d()


end #module
