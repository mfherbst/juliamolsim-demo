using AtomsIO
using DFTK
using LazyArtifacts
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
    stresses = compute_stresses_cart(scfres)
    if mpi_master()
        println(DFTK.timer)
    end
    (; energy=scfres.energies.total, forces, stresses)
end

function run_extxyz(file)
    @assert endswith(file, ".extxyz")
    systems = load_trajectory(file)
    bn, _   = splitext(file)
    outfile = bn * "_out.extxyz"

    systems_plus_data = empty(systems)
    for (i, system) in enumerate(systems)
        if mpi_master()
            println("#")
            println("# $i / $(length(systems)) ================================")
            println("#")
        end

        res = run_calculation(system)
        push!(systems_plus_data,
              FlexibleSystem(system; res.energy, res.forces, res.stresses))
        if mpi_master()
            save_trajectory(outfile, systems_plus_data)
        end
    end
end
