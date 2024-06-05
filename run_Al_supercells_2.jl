#!/bin/bash
#=
BN=$(basename "$0" .jl)
mpiexecjl --project -np 6 julia --project -t1 $BN.jl |& tee $BN.log
exit $?
=#
include("calculations.jl")
run_extxyz("Al_supercells_2.extxyz")
