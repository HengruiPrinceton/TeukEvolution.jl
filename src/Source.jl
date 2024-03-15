"""
Analytical Source Term for the Teukolsky Equation
"""
module Source

include("Fields.jl")
include("Radial.jl")
include("Sphere.jl")

using .Fields: Field
import .Radial
import .Sphere
export Evaluate_Source!

function Evaluate_Source!(
    Rvals::Vector{Float64}, # Value of R
    t::Float64,   # time to evaluate
    S::Array{ComplexF64,2} # source term to be evaluated
)
    nx, ny = size(S)
    omega = 0.5 # - im*0.2
    for j = 1:ny
        for i = 1:nx
            S[i,j] = Rvals[i]*exp(-im*omega*t)
        end
    end
end


end