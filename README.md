Iceberg.jl
==========

Performs ice-seawater interface calculations using level set methods. Level set scheme is based in part off that described by Chen, et al. (1997). Model state in passed around using a subtype of the `ModelState` abstract type. Physical constants and parameters are contained in a `PhysicalParameters` type. See `initialize1d()` for an example.

## Model dimension

Models in one, two, and three dimensions will be supported. Dimension is indicated based on the concrete type of the `ModelState` type passed to functions. As of #9eee576048c3d53ee3b9a6ab8767f808a4224335, one-dimensional models are implemented.

## Field evolution

Methods to solve the heat equation are contained in `heat.jl`. The publicly facing method `tsolve!` is imported into the main Iceberg namespace.

Interface fluxes can be calculated using the `velocity` function. In the future, there will be functions to handle level set function evolution and reinitialization.

## Validation

The function `initialize1d_hill()` recreates a scenario easily compared to the analytical solution of Hill (1987), Section 1.3

## To do

- extend to higher dimensions
- supplement explicit heat equation solver with an implicit scheme
- add methods to evolve and reinitialize level set function
- add simple Navier-Stokes field solver for fluid phase

## References

Chen, S., B. Merriman, S. Osher and P. Smereka. _A simple level set method for solving Stefan problems_, 1997. Journal of Computational Physics, 135, 8-29.

Hill, J. M. _One-dimensional Stefan Problems: an introduction_, 1987. Pitman Monographs and Surveys in Pure and Applied Mathematics, 31. John Wiley and Sons, New York.
