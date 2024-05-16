#!/bin/sh
#=
BN=$(basename "$0" .jl)
mpiexecjl --project -np 8 julia --project -t1 $BN.jl | tee $BN.log
exit $?
=#
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
    scfres = self_consistent_field(basis; tol=1e-10)
    forces = compute_forces_cart(scfres)
    if mpi_master()
        println(DFTK.timer)
    end
    (; energy=scfres.energies.total, forces)
end

systems = load_trajectory("Al_bulk.extxyz")
for system in systems
    if mpi_master()
        println("#")
        println("# =========================================")
        println("#")
    end
    run_calculation(system)
end
