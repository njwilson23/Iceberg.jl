
# maybe this should be a dictionary?
#PhysicalParams = {:kapl => 1.0,
#                  :kaps => 1.0,
#                  :cpl => 1.0,
#                  :cps => 1.0,
#                  :rhol => 1.0,
#                  :rhos => 1.0,
#                  :tmelt => 0.0}

immutable PhysicalParams
    Lf::Float64
    kapl::Float64
    kaps::Float64
    cpl::Float64
    cps::Float64
    rhol::Float64
    rhos::Float64
    tmelt::Float64
end

immutable ModelParams
    dt::Float64
    dx::Tuple
    nx::Tuple
end

abstract ModelState

type ModelState1d   <: ModelState
    temp::Array
    phi::Array
    params::ModelParams
end

type ModelState2d   <: ModelState
    temp::Array
    phi::Array
    params::ModelParams
end

print(io::IO, A::ModelState) = @sprintf("Problem of size %s", size(A.phi))

