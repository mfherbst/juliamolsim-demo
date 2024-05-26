using AtomsIO
using DFTK
using LazyArtifacts
using JSON3
using Printf
setup_threading()

function run_calculation(system)
    @assert all(isequal(:Al), atomic_symbol(system))

    DFTK.reset_timer!(DFTK.timer)
    system = attach_psp(system; Al=artifact"pd_nc_sr_pbe_stringent_0.4.1_upf/Al.upf")
    model  = model_PBE(system; smearing=Smearing.Gaussian(), temperature=1e-3)
    kgrid  = kgrid_from_maximal_spacing(model, 0.13)
    basis  = PlaneWaveBasis(model; Ecut=26, kgrid)

    if mpi_master()
        show(stdout, "text/plain", basis)
        println()
    end
    scfres   = self_consistent_field(basis; tol=1e-10)
    forces   = compute_forces_cart(scfres)

    if length(system) > 45  # Performance issue in DFTK
        stresses = "nothing"
    else
        stresses = compute_stresses_cart(scfres)
    end
    if mpi_master()
        println(DFTK.timer)
    end
    (; scfres, forces, stresses)
end

function run_extxyz(file)
    @assert endswith(file, ".extxyz")
    systems = load_trajectory(file)
    bn, _   = splitext(file)

    for (i, system) in enumerate(systems)
        istr = @sprintf "%04i" i
        outfile = joinpath("dftk_output", bn * "-$istr.json")
        isfile(outfile) && continue

        if mpi_master()
            println("#")
            println("# $i / $(length(systems)) ================================")
            println("#")
        end

        res = run_calculation(system)

        extra_data = Dict("forces" => res.forces, "stresses" => res.stresses)
        save_scfres(outfile, res.scfres; extra_data)
    end
end
