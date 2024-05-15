#!/bin/sh
#=
BN=$(basename "$0" .jl)
julia --project -t4 $BN.jl | tee $BN.log
exit $?
=#
using AtomsIO
using DFTK
setup_threading()

function run_calculation(system)
    @assert all(isequal(:Al), atomic_symbols(system))

    model  = model_PBE(system; smearing=Smearing.Gaussian(), temperature=1e-3)
    kgrid  = kgrid_from_maximal_spacing(model, 0.13)
    basis  = PlaneWaveBasis(model; Ecut=26, kgrid)

    if mpi_master()
        show(stdout, "text/plain", basis)
        println()
    end

    scfres = self_consistent_field(basis; tol=1e-10)
    forces = compute_forces_cart(scfres)
    (; energy=scfres.energies.total, forces)
end

systems = load_trajectory("Al_bulk.extxyz")
for system in systems
    run_calculation(system)
end
