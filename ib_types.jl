module ib_types
import Base.show
export ModelState, ModelState1d, ModelState2d
export ModelParams, PhysicalParams

# maybe this should be a dictionary?
#PhysicalParams = {:kapl => 1.0,
#                  :kaps => 1.0,
#                  :cpl => 1.0,
#                  :cps => 1.0,
#                  :rhol => 1.0,
#                  :rhos => 1.0,
#                  :tmelt => 0.0}

type PhysicalParams
    Lf::Float64
    kapl::Float64
    kaps::Float64
    cpl::Float64
    cps::Float64
    rhol::Float64
    rhos::Float64
    tmelt::Float64
end

type ModelParams
    dt::Float64
    dx::Tuple
    nt::Integer
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


end
