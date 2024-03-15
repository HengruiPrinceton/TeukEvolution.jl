module TeukEvolution

include("Fields.jl")
include("Io.jl")
include("Radial.jl")
include("Sphere.jl")
include("Id.jl")
include("Evolution.jl")

using .Fields: Field, Initialize_Field
import .Io
import .Radial
import .Sphere
import .Id
using .Evolution: Evo_lin_f, Initialize_Evo_lin_f, Evolve_lin_f!

import TOML

"""
    launch(params::Dict{String,Any})

    Launches time evolution of Teukolsky code, reading in paramter
    dictionary file (see examples).
"""
function launch(params::Dict{String,Any})::Nothing
    println("Launching run")
    println("params=$params")

    nx = params["nx"]
    ny = params["ny"]
    nt = params["nt"]
    ts = params["ts"]

    psi_spin = params["psi_spin"]

    runtype = params["runtype"]

    prec = params["precision"]

    cl = convert(prec, params["cl"])
    cfl = convert(prec, params["cfl"])
    bhs = convert(prec, params["bhs"])
    bhm = convert(prec, params["bhm"])

    outdir = params["outdir"]
    ##===================
    ## Derived parameters
    ##===================
    minr = bhm * (1 + sqrt(1 - (bhs / bhm)^2) ) # horizon (uncompactified)
    maxR = 1 / minr # dt should not depend on cl
    dr = maxR / (nx - 1)
    dt = min(cfl * dr * bhm^2, 6 / ny^2) # make the time step roughly proportional to mass instead of inversely proportional
    println(dt)
    println("Number of threads: $(Threads.nthreads())")

    println("Setting up output directory")
    if !isdir(outdir)
        mkdir(outdir)
    else
        rm(outdir, recursive = true)
        mkdir(outdir)
    end
    println("Writing input parameters")
    cp("./run.jl", outdir*"/run.jl")

    println("Initializing constant fields")
    Rv = Radial.R_vals(nx, dr)
    Yv = Sphere.Y_vals(ny)
    Cv = Sphere.cos_vals(ny)
    Sv = Sphere.sin_vals(ny)
    Mv = params["m_vals"]
    time = 0.0

    ##=================
    ## Dynamical fields 
    ##=================
    println("Initializing linear psi4")
    lin_f = Initialize_Field(
        name = "lin_f",
        spin = psi_spin,
        boost = psi_spin,
        falloff = 1,
        Mvals = Mv,
        nx = nx,
        ny = ny,
    )
    lin_p = Initialize_Field(
        name = "lin_p",
        spin = psi_spin,
        boost = psi_spin,
        falloff = 1,
        Mvals = Mv,
        nx = nx,
        ny = ny,
    )

    ##=======================================
    ## Fixed fields (for evolution equations) 
    ##=======================================
    println("Initializing psi4 evolution operators")
    evo_psi4 = Initialize_Evo_lin_f(
        Rvals = Rv,
        Cvals = Cv,
        Svals = Sv,
        Mvals = Mv,
        bhm = bhm,
        bhs = bhs,
        cl = cl,
        spin = psi_spin,
    )

    ##=============
    ## Initial data
    ##=============
    println("Initial data")

    if params["id_kind"] == "gaussian"
        for (mi, mv) in enumerate(Mv)
            Id.set_gaussian!(
                lin_f[mv],
                lin_p[mv],
                psi_spin,
                mv,
                params["id_l_ang"][mi],
                params["id_ru"][mi],
                params["id_rl"][mi],
                params["id_width"][mi],
                params["id_amp"][mi][1] + params["id_amp"][mi][2] * im,
                cl,
                Rv,
                Yv,
                params["ingoing"],
                dr,
                nx,
                ny,
            )
            Io.save_csv(t = 0.0, mv = mv, outdir = outdir, f = lin_f[mv])
            Io.save_csv(t = 0.0, mv = mv, outdir = outdir, f = lin_p[mv])
        end
    elseif params["id_kind"] == "qnm"
        overtone_n = params["id_overtone_n"]
        overtone_amp = params["id_qnm_amp"]
        for mv in Mv
            for (ni,n) in enumerate(overtone_n)
                Id.set_qnm!(
                    lin_f[mv],
                    lin_p[mv],
                    psi_spin,
                    mv,
                    n,
                    params["id_filename"],
                    overtone_amp[ni],
                    params["id_m"],
                    Rv,
                    Yv,
                )
            end
            Io.save_csv(t = 0.0, mv = mv, outdir = outdir, f = lin_f[mv])
            Io.save_csv(t = 0.0, mv = mv, outdir = outdir, f = lin_p[mv])
        end
    else
        throw(DomainError(params["id_kind"], "Unsupported `id_kind` in parameter file"))
    end

    ##===================
    ## Time evolution 
    ##===================
    println("Beginning evolution")

    for tc = 1:nt
        # time of evolution
        t = (tc - 1) * dt / bhm

        Threads.@threads for mv in Mv
            Evolve_lin_f!(lin_f[mv], lin_p[mv], Rv, evo_psi4[mv], dr, dt, t)

            lin_f_n = lin_f[mv].n
            lin_p_n = lin_p[mv].n
            lin_f_np1 = lin_f[mv].np1
            lin_p_np1 = lin_p[mv].np1

            for j = 1:ny
                for i = 1:nx
                    lin_f_n[i, j] = lin_f_np1[i, j]
                    lin_p_n[i, j] = lin_p_np1[i, j]
                end
            end
        end    
        if tc % ts == 0
            t = tc * dt / bhm
            println("time/bhm ", t)
            Threads.@threads for mv in Mv
                Io.save_csv(t = t, mv = mv, outdir = outdir, f = lin_f[mv])
                Io.save_csv(t = t, mv = mv, outdir = outdir, f = lin_p[mv])
            end
        end
    end
    println("Finished evolution")
    return nothing
end

end
