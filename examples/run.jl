include("../src/TeukEvolution.jl")

import .TeukEvolution as TE

# write the parameter file

params = Dict(
    "outdir" => "ultra_a0.999",
    "nx" => 1024,   # number of x radial grid points
    "ny" => 40,    # number of y collocation points
    "nt" => 640000, # number of time steps
    "ts" => 400,   # save every ts time steps
    "psi_spin" => -2, # spin-weight of linear evolution scalar
    "id_kind" => "gaussian",
    "runtype" => "linear_field",
    "m_vals" => [2],   # m angular values
    "id_l_ang" => [2],
    "id_ru" => [4.0],
    "id_rl" => [2.5],
    "id_width" => [1.5],
    "ingoing" => true,
    # format: for each m value: [real part, imaginary part]
    "id_amp" => [[0.1,0]],
    "cl" => 1.0, # compactification scale
    "cfl" => 0.5, # CFL number
    "bhs" => 0.5, # black hole spin 
    "bhm" => 1, # black hole mass
    # save with l decomposition
    "save_l" => false,
    "l_to_save" => 2,
    "precision" => Float64, # precision the code is compiled at
)

@time TE.launch(params)
