Iceberg.jl
==========

[![Build Status](https://travis-ci.org/njwilson23/Iceberg.jl.svg?branch=master)](https://travis-ci.org/njwilson23/Iceberg.jl)

Performs ice-seawater interface calculations using level set methods. The level
set scheme is similar to that described by Chen, et al. (1997). Model state is
passed around using a subtype of the `ModelState` abstract type. Physical
constants and parameters are contained in a `PhysicalParameters` type. See
`initialize1d()` for an example.

Requires a recent build of Julia (i.e. 0.3).

## Model dimension

Models in one, two, and three dimensions will be supported. Dimension is
indicated based on the concrete type of the `ModelState` type passed to
functions. Currently, one-dimensional models are implemented.

## Field evolution

Methods to solve the heat equation are contained in `heat.jl`.

Interface fluxes can be calculated using the `front_velocity()` function. Level
set reinitialization uses a global method described by Chen et al. (1997) and
Osher and Fedkiw (2003). In the future, there will be functions to handle level
set function evolution with a hyperbolic equation solver.

## Validation

The function `initialize1d_hill()` recreates a scenario easily compared to the
analytical solution of Hill (1987), Section 1.3. An instance of this scenario
with `n=200` is a part of the test suite.

## To do

- extend to higher dimensions
- supplement explicit heat equation solver with an implicit scheme
- add methods to evolve the level set function
- add simple Navier-Stokes field solver for fluid phase

## References

Chen, S., B. Merriman, S. Osher and P. Smereka. _A simple level set method for
solving Stefan problems_, 1997. Journal of Computational Physics, 135, 8-29.

Hill, J. M. _One-dimensional Stefan Problems: an introduction_, 1987. Pitman
Monographs and Surveys in Pure and Applied Mathematics, 31. John Wiley and Sons,
New York.

Osher, S. and R. Fedkiw. _Level Set Methods and Dynamic Implicit Surfaces_,
2003. Applied Mathematical Sciences, 153. Springer-Verlag, New York.
