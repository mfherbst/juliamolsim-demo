#!/bin/bash
#=
BN=$(basename "$0" .jl)
julia --project -t1 $BN.jl |& tee $BN.log
exit $?
=#
include("calculations.jl")
using LinearAlgebra

lattice   = 20.0Matrix(I, 3, 3)
atoms     = [ElementPsp(:Al; psp=load_psp(artifact"pd_nc_sr_pbe_stringent_0.4.1_upf/Al.upf"))]
positions = [zeros(3)]
system = periodic_system(lattice, atoms, positions)
model  = model_PBE(lattice, atoms, positions;
                   smearing=Smearing.Gaussian(), temperature=1e-3)
basis  = PlaneWaveBasis(model; Ecut=26, kgrid=[1, 1, 1])

if mpi_master()
    show(stdout, "text/plain", basis)
    println()
end
scfres = self_consistent_field(basis; tol=1e-10, damping=0.5)
@show scfres.energies.total
