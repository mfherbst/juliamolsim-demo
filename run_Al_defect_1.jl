#!/bin/sh
#=
BN=$(basename "$0" .jl)
mpiexecjl --project -np 8 julia --project -t1 $BN.jl | tee $BN.log
exit $?
=#
include("calculations.jl")
run_extxyz("Al_defect_1.extxyz")
