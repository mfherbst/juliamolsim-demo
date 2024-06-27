#!/bin/bash
#=
BN=$(basename "$0" .jl)
mpiexecjl --project -np 6 julia --project -t1 $BN.jl
exit $?
=#
using AtomsBase
using AtomsIO
using DFTK
using LazyArtifacts
using Unitful

systems = load_trajectory("Al_structures.extxyz")
all_data = map(systems) do system
    system   = attach_psp(system; Al=artifact"pd_nc_sr_pbe_stringent_0.4.1_upf/Al.upf")
    model    = model_PBE(system; smearing=Smearing.Gaussian(), temperature=1e-3)
    basis    = PlaneWaveBasis(model; Ecut=26, kgrid=(12, 12, 12))
    scfres   = self_consistent_field(basis; tol=1e-10)
    forces   = compute_forces_cart(scfres)
    stresses = compute_stresses_cart(scfres)

    properties = (; energy=ustrip(u"eV", scfres.energies.total * u"hartree"),
                    virial=ustrip.(u"eV", stresses * u"eV"))
    atoms = map(system, forces) do at, at_force
        Atom(; pairs(at)..., force=ustrip.(u"eV/Å", at_force * u"hartree/bohr"))
    end
    FlexibleSystem(system; atoms, properties...)
end
save_trajectory("Al_training_data.extxyz", all_data)
